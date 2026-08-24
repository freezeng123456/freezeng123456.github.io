---
title: "Chapter 5: Conclusions"
description: Paper conclusions, method selection, full experiment inventory, GPU optimization, and reporting standards
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization/chapter-5-unified-view
tags:
  - parallel-in-time
  - methodology
---

> [!note] Content boundary
> Section 5 of the paper (p. 481) contains no numbered subsections or numbered equations. This page first explains the source conclusion paragraph by paragraph. The subsequent algebraic synthesis, method selection, Python experiments, T4 GPU performance, and reporting standards are marked as site supplements and do not use 5.x numbering.

## Section 5: Conclusions

The paper separates the two classes by whether the solution forgets information. A parabolic equation forgets a large amount of information in time and therefore has a solution that is local in time. Many PinT methods can exploit this property, including Parareal, STMG, ParaExp, ParaDiag, and domain-decomposition waveform relaxation.

Hyperbolic equations retain very fine solution features over very long times, which narrows the effective method set. The paper highlights ParaExp, ParaDiag, and SWR, with particular emphasis on the relationship between SWR and tent pitching.

For further study, the authors recommend the monograph by Gander and Lunet (2024, SIAM), which provides historical context, self-contained convergence analyses, and short executable MATLAB programs for individual PinT methods. The code used to produce the results in the paper is public in [wushulin/ActaPinT](https://github.com/wushulin/ActaPinT).

> [!note] Site supplement: from the conclusion to method selection
> The following first filter is this site's extension of the conclusion above, not the paper's wording: determine how quickly the problem forgets high-frequency information, then select the parallel structure; an algorithm's name or historical category cannot replace this dynamical assessment. The phrase "temporal memory" used later on this page is likewise this site's coinage. The paper's own wording is that parabolic problems "tend to forget a lot of information in time and thus have solutions that are local in time" and that hyperbolic problems "retain very fine solution features over very long times".

## Site synthesis: solving one all-at-once system

For a linear all-at-once discretization

$$
A\boldsymbol U=\boldsymbol b,
$$

many PinT iterations can be written as

$$
\boldsymbol U^{k+1}
=\boldsymbol U^k+M^{-1}
(\boldsymbol b-A\boldsymbol U^k).
$$

$M^{-1}$ is a parallel approximation to $A^{-1}$. Each method chooses a different locality, hierarchy, or transform:

| Method    | Main source of $M^{-1}$                           | Concurrent work                              | Carrier of long-range information            |
| --------- | ------------------------------------------------- | -------------------------------------------- | -------------------------------------------- |
| SWR       | space–time subdomain inverses and transmission    | complete subdomain waveforms                 | interface waveforms and characteristic cones |
| PIDC/RIDC | integral residual correction                      | window or correction-level pipeline          | high-order error equation                    |
| ParaExp   | local inhomogeneous solves and exponential action | local forced responses and homogeneous tails | exact $e^{tA}$ propagation                   |
| ParaDiag  | circulant/diagonalizable time operator            | shifted spatial systems after FFT            | global temporal frequencies                  |
| Parareal  | lower-triangular coarse propagator inverse        | fine solves on large intervals               | sequential coarse prediction                 |
| PFASST    | multilevel preconditioner for collocation         | SDC sweeps on time steps                     | FAS coarse collocation correction            |
| MGRiT     | temporal multilevel cycle                         | F relaxation and coarse levels               | C points and overlapping relaxation          |
| STMG      | full space–time multilevel cycle                  | time-block Jacobi                            | coarse space–time grids                      |

This common form poses three questions. Which error modes does $M^{-1}$ reduce? Does it preserve phase, mean value, and shock position? During one application of $M^{-1}$, which operations are genuinely concurrent and which remain sequential?

## Site synthesis: method-selection map

![Map for selecting parallel-in-time methods](assets/diagrams/pint/en/method-selection.svg)

### Strong dissipation and low-to-moderate temporal order

Heat and sufficiently diffusive reaction–diffusion systems are natural candidates for MGRiT and STMG. Parareal provides a quick test of coarse-propagator quality and a nonintrusive baseline. STMG has greater scalability potential but requires access to the all-at-once operator, smoother, and grid transfers.

### High-order collocation

PFASST is attractive when high temporal order is required and a collocation solve is expensive. Nodes, SDC sweeps, the coarse collocation level, and spatial parallel resources need to be designed together. High formal order does not guarantee an efficient temporal pipeline.

### Large linear or linearized systems

ParaExp accurately transports long-range linear information when exponential action scales. ParaDiag removes temporal forward substitution through FFTs when complex shifted spatial systems have capable solvers. ParaDiag-I is limited by conditioning of the temporal eigenvectors. ParaDiag-II must also balance $\alpha$, outer Krylov convergence, and roundoff.

### Transport, waves, and low-viscosity nonlinearity

Prioritize structures that represent characteristics and phase, including SWR/OSWR, tent pitching, ParaExp, $\alpha$-ParaDiag, and phase-aware coarse propagation. Standard Parareal, MGRiT, and STMG remain useful diagnostic baselines. If fine/coarse phase mismatch grows with frequency, increasing the interval count is likely to amplify the difficulty.

### Six questions before choosing a method

1. Do important modes lie near the negative real axis or close to the unit circle/imaginary axis?
2. Do boundaries permit outflow, recirculate a periodic signal, or reflect waves?
3. Does nonlinearity continuously create shocks and high frequencies?
4. Which reusable component is available: a time stepper, shifted spatial solver, exponential action, or all-at-once operator?
5. Is the goal lower iteration count, higher single-node throughput, or multi-node strong/weak scaling?
6. How much intrusive modification, global transformation, and all-time storage is acceptable?

## Site synthesis: parameter reference

| Parameter                            | Location                       | Direct role                              | Quantities to monitor together                            |
| ------------------------------------ | ------------------------------ | ---------------------------------------- | --------------------------------------------------------- |
| SWR overlap and Robin $p$            | subdomain interface            | controls waveform transfer               | window length, viscosity, interface cost                  |
| IDC node count $M$ and corrections   | error equation                 | limit formal order and pipeline depth    | regularity, fill and drain time                           |
| Parareal intervals $N$ and ratio $J$ | fine/coarse propagation        | set concurrency and mismatch             | iteration count, sequential coarse cost, phase error      |
| MGRiT coarsening and CF count        | time hierarchy                 | set overlap contraction and fine work    | total factor at equal fine-solve cost                     |
| ParaDiag $\alpha$                    | cyclic head–tail approximation | smaller values improve the approximation | $\epsilon/\alpha$ roundoff and shifted-solve stability    |
| STMG damping $\eta$                  | time-block Jacobi              | controls high-frequency smoothing        | time integrator, cycle cost, spatial coarsening condition |

Each parameter interacts with the physical spectrum, discrete stability function, and machine cost. Minimizing iteration count in isolation can move the work into a much more expensive iteration.

## Site reproduction: experiment inventory

The Python project now provides four whole-article group entries: `schematics`, `section2`, `section3`, and `section4`. Together they run six computed schematics, 39 numerical figures, and Tables 3.1–3.2, for 47 local computations. The earlier eight baselines and combined `paper_validation` entry remain useful regression targets but no longer describe the project's full coverage.

| Group entry  | Coverage                                     | Local computations |
| ------------ | -------------------------------------------- | -----------------: |
| `schematics` | Figures 1.1, 3.2, 3.4, 3.7, 4.1, and 4.7     |                  6 |
| `section2`   | Figures 2.1–2.4, 24 panels                   |                  4 |
| `section3`   | fifteen numerical figures and Tables 3.1–3.2 |                 17 |
| `section4`   | Figures 4.2–4.6 and 4.8–4.22                 |                 20 |

| Python output               | Page location                  | Machine-readable record                                      |
| --------------------------- | ------------------------------ | ------------------------------------------------------------ |
| `solution_heat_ade`         | Chapter 2, advection–diffusion | [[assets/pint/data/solution_heat_ade.json\|JSON]]            |
| `solution_burgers`          | Chapter 2, Burgers             | [[assets/pint/data/solution_burgers.json\|JSON]]             |
| `solution_wave`             | Chapter 2, wave                | [[assets/pint/data/solution_wave.json\|JSON]]                |
| `parareal_heat_ade`         | Chapter 4, Parareal            | [[assets/pint/data/parareal_heat_ade.json\|JSON]]            |
| `parareal_burgers`          | Chapter 4, Parareal            | [[assets/pint/data/parareal_burgers.json\|JSON]]             |
| `mgrit_heat_ade`            | Chapter 4, MGRiT               | [[assets/pint/data/mgrit_heat_ade.json\|JSON]]               |
| `iterative_paradiag_ade`    | Chapter 3, ParaDiag            | [[assets/pint/data/iterative_paradiag_ade.json\|JSON]]       |
| `stmg_heat_ade`             | Chapter 4, STMG                | [[assets/pint/data/stmg_heat_ade.json\|JSON]]                |
| Figure 3.15 validation      | Chapter 3, ParaDiag            | [[assets/pint/data/figure_3_15_validation.json\|JSON]]       |
| Figure 4.5 validation       | Chapter 4, Parareal            | [[assets/pint/data/figure_4_5_validation.json\|JSON]]        |
| Figures 4.9–4.10 validation | Chapter 4, MGRiT               | [[assets/pint/data/figure_4_10_validation.json\|JSON]]       |
| Figure 4.19 validation      | Chapter 4, STMG                | [[assets/pint/data/figures_4_19_4_20_validation.json\|JSON]] |
| Figure 4.20 validation      | Chapter 4, STMG                | [[assets/pint/data/figures_4_19_4_20_validation.json\|JSON]] |
| T4 GPU validation           | this chapter, GPU acceleration | [[assets/pint/data/gpu_benchmark_t4.json\|JSON]]             |

The cross-experiment summary is [[assets/pint/data/paper_validation_summary.json\|paper_validation_summary.json]].

All 24 upstream MATLAB scripts have now been ported, and experiments described only in the paper text have also been reconstructed. Two gaps must remain explicit in any completeness claim: the Figure 3.14 NKA curves at $\nu=0.02$ are unmatched, and Table 4.1 is a historical supercomputer measurement from a different three-dimensional parallel code and can only be quoted. See the [[en/computational-mathematics/knowledge-notes/time-parallelization/reproduction-audit-2026-08-24|ActaPinT-Python whole-article reproduction audit]].

## Site reproduction: formal run workflow

```bash
python3.11 -m pip install -e ".[test]"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
actapint schematics --output-dir results/schematics
actapint section2 --output-dir results/section2
actapint section3 --output-dir results/section3
actapint section4 --output-dir results/section4
```

The formal run on July 31, 2026 used Python 3.11, NumPy 2.4.6, SciPy 1.17.1, and Matplotlib 3.11.1. The CPU paper suite took about 4 minutes 24 seconds wall time and stayed below 280 MiB peak resident memory. The code path includes CPU sparse factorization, Krylov methods, and FFTs.

On August 24, 2026, reinstalling from the open dependency ranges produced 101 passes, one GPU skip, and two numerical-tolerance failures among 104 tests; the latest public CI run also fails. The historical formal environment and the current dependency environment must therefore be reported separately. An old green run cannot replace a dependency lock; the audit gives the two failure details.

Each experiment writes:

1. editable SVG and high-resolution PNG produced from the same analysis arrays;
2. JSON containing grid, physical parameters, tolerance, stopping convention, and metrics;
3. deterministic seeds where an all-at-once random initial vector is needed;
4. a separate validation summary for paper-matched experiments.

## Site reproduction: GPU acceleration and profiling

Function-level profiling attributes 43.06 of 62.95 seconds in the quick paper suite to Figure 4.5, including 38.11 seconds in Burgers fine propagation. All FFTs take only 0.007 seconds and all GMRES calls total 0.251 seconds. The first CUDA backend therefore batches the 40 independent Burgers fine propagations in each Parareal iteration.

The CuPy backend keeps spatial operators resident on the GPU and assembles and solves 40 independent Newton systems in batches. Causal coarse propagation and the remaining experiments stay on the CPU, yielding a hybrid CPU/GPU implementation.

| T4 double-precision test                 |      CPU |     GPU | Speedup |
| ---------------------------------------- | -------: | ------: | ------: |
| 40 Burgers fine propagators, 32 substeps |  2.893 s | 0.246 s |  11.76× |
| combined `paper_validation` suite        | 263.57 s | 67.92 s |   3.88× |

The maximum absolute CPU/GPU difference for one batch is $2.33\times10^{-15}$. In the earlier `paper_validation` entry, Figure 4.5 retains stopping iterations ADE 14/24/35 and Burgers 14/21/25; its Burgers run follows the upstream-script $\frac12(u^2)_x$. The whole-article `section4` entry uses $(u^2)_x$ to match the paper panel and gives 14/24/36, so the two must not be conflated. The 3.88× value belongs to the earlier combined validation entry, not to all 39 numerical figures in Sections 2–4. Continuing beyond the $10^{-10}$ target toward machine precision produces normal rounding differences between CPU SuperLU and GPU batched LU because their floating-point operation orders differ.

```bash
python3.11 -m pip install -e ".[gpu,test]"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  actapint paper_validation --backend gpu \
  --output-dir results/paper-gpu
```

### Next optimization opportunities

1. Replace dense batched LU with a cyclic tridiagonal CUDA kernel so storage and work scale linearly with $N_x$.
2. Keep the complete Parareal state on the GPU and investigate parallel-prefix/scan coarse propagation to reduce host transfers and the sequential tail.
3. Implement batched complex shifted banded solves on larger ParaDiag grids. FFT and GMRES work on the current $100\times100$ grid is too small for stable GPU benefit.
4. Separate spatial concurrency, temporal concurrency, and communication overlap in multi-GPU strong- and weak-scaling tests.

The machine-readable record is [[assets/pint/data/gpu_benchmark_t4.json\|gpu_benchmark_t4.json]]. Current data demonstrate single-GPU kernel and end-to-end acceleration and do not establish multi-GPU scaling.

## Site reproduction: interpretation boundaries

- Current experiments measure numerical convergence and single-node CPU/GPU performance, without time-dimensional MPI strong or weak scaling.
- Site figures do not report MPI rank count, network volume, initialization cost, or cross-node wall-clock speedup.
- MATLAB and NumPy random generators do not produce identical initial arrays; convergence factors, phases, and final states are better comparison targets.
- Sparse factorization, FFT ordering, and GMRES reductions can differ normally near $10^{-14}$ to $10^{-16}$.
- The invalid expression `nu=0.002max;` in `MGRiT_Heat_ADE.m` is interpreted as $\nu=0.002$ from its branch context and the paper.
- Source Figure 4.9 annotates the maximum Parareal factor at $\nu=0.01$ as $0.9986$, while equation (4.5b) and the upstream-script stability function give $1.0501$ in Python. The site retains both and records this as a reproduction discrepancy.
- STMG validation retains the upstream backward-Euler residual convention and stores a consistent postsmoothing residual separately in JSON.
- Convergence to the sequential fine solution establishes algorithmic consistency and gives no wall-clock speedup guarantee.
- The large-scale STMG data in Table 4.1 belong to the cited three-dimensional parallel implementation and are not measurements from the present Python project.

## Site standard: minimum reporting requirements for future experiments

An algorithmic-convergence experiment should state:

- PDE, boundaries, initial data, and source;
- spatial and temporal discretizations and the fine-grid reference;
- fine/coarse propagators, interval count, and coarsening;
- error norm, residual definition, stopping threshold, and iteration cap;
- parameter sweeps and observed failure points.

A parallel-performance experiment should additionally state:

- CPU/GPU models, precision, and rank/thread/device allocation;
- fine/coarse cost ratio, concurrent tasks per iteration, and load balance;
- initialization, transfer, communication, synchronization, and I/O time;
- total wall time, speedup, efficiency, and baseline implementation;
- strong or weak scaling as the interval count, spatial size, and device count vary.

Error decay by iteration supports a numerical-convergence claim but cannot support a parallel-efficiency claim on its own. Method comparisons must also normalize fine-propagation count or total work so curves with different per-iteration costs are not ranked directly.

## Site-wide coverage table

| Source range                  | Site chapter    | Completeness statement                                                                            |
| ----------------------------- | --------------- | ------------------------------------------------------------------------------------------------- |
| Section 2, pp. 388–396        | Chapter 2       | four models, all boundary settings, every Figure 2.1–2.4 observation group, and PinT implications |
| Sections 3.1–3.2, pp. 396–405 | Chapter 3.1–3.2 | history, WR/SWR, Theorems 3.1–3.2, OSWR, MTP/UTP                                                  |
| Sections 3.3–3.4, pp. 405–415 | Chapter 3.3–3.4 | IDC/PIDC/RIDC derivation and regularity tests, linear and nonlinear ParaExp                       |
| Section 3.5, pp. 415–443      | Chapter 3.5     | ParaDiag-I/II, Theorems 3.5–3.9, BVM, NKA, circulant and $\alpha$-circulant experiments           |
| Sections 4.1–4.4, pp. 443–460 | Chapter 4.1–4.4 | Parareal, PFASST, MGRiT, Theorems 4.1–4.6, Figures 4.1–4.11                                       |
| Sections 4.5–4.6, pp. 460–481 | Chapter 4.5–4.6 | both diagonalized Parareal variants, STMG, Theorems 4.7–4.9, Figures 4.12–4.22, Table 4.1         |
| Section 5, p. 481             | Chapter 5       | hyperbolic/parabolic conclusion, recommended methods, monograph, and public code                  |

Each chapter ends with a more granular source-page audit. Site supplements occupy explicitly labeled sections and are not blended into claims attributed to the paper.

## Summary

The starting point for PinT selection is the temporal memory of the dynamics. Strong diffusion permits coarse temporal levels to represent the remaining slow modes. Transport, waves, and low-viscosity nonlinear systems require preservation of phase, characteristics, and shock position. The unified all-at-once view helps compare algorithms, while a successful implementation must satisfy three conditions together: iterations to tolerance remain controlled, concurrent work dominates runtime, and communication and memory scale. A complete reproduction reports evidence for each layer separately.

## Source

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Section 5, p. 481.
