---
title: Time Parallelization for Hyperbolic and Parabolic Problems
description: Chapter-by-chapter notes and reproducible experiments for parallel-in-time methods
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization
tags:
  - computational-mathematics
  - parallel-in-time
---

These notes follow M. J. Gander, S.-L. Wu, and T. Zhou, _Time Parallelization for Hyperbolic and Parabolic Problems_, Acta Numerica 34 (2025), pp. 385-489. The survey distinguishes methods that remain effective for propagative problems from methods designed primarily for dissipative problems.

> [!info] Source graphics and license
> The close-reading pages include every original Figure 1.1–4.22 and Tables 3.1, 3.2, and 4.1: 48 graphical and tabular assets extracted from the paper PDF as scalable SVGs. Adjacent explanations are newly written. The paper is distributed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); authors, title, and DOI appear under “Primary sources.” Redrawn diagrams and Python reproductions are labeled separately from source graphics.

## Reading sequence

1. [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-1-why-parallelize-in-time|Chapter 1: Why Parallelize in Time?]] follows the abstract and introduction paragraph by paragraph, covering the hardware context, the causal chain, four historical lineages, the two-way classification, and the all-at-once system.
2. [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-2-model-problems|Chapter 2: Model Problems Linking the Parabolic and Hyperbolic World]] follows Sections 2.1–2.4 exactly, reads every source panel across the heat, advection–diffusion, Burgers, and wave equations, and reports the site's numerical supplements separately.
3. [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-3-hyperbolic-methods|Chapter 3: Effective PinT Methods for Hyperbolic Problems]] derives SWR, PIDC/RIDC, ParaExp, and ParaDiag and includes ParaDiag-II experiments for heat, ADE, and wave problems.
4. [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-4-parabolic-methods|Chapter 4: PinT Methods Designed for Parabolic Problems]] analyzes Parareal, PFASST, MGRiT, two diagonalized variants, and STMG, with recomputed convergence studies.
5. [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-5-unified-view|Chapter 5: Conclusions]] explains the source conclusion paragraph by paragraph and marks method selection, experiment inventory, GPU performance, and reporting rules as site supplements.

## Paper-coverage progress

“Paragraph-level complete” means that the claims, equations, figures, historical links, qualifications, and section relationships have all been checked against the source. Existing experiments remain in place after the corresponding source discussion.

| Source range           | Website chapter | Current status                                                                                                                                               |
| ---------------------- | --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Abstract and Section 1 | Chapter 1       | **Paragraph-level complete**: every argument on pp. 385–388, equation (1.1), and the original Figure 1.1 are covered                                         |
| Sections 2.1–2.4       | Chapter 2       | **Paragraph-level complete**: (2.1)–(2.7), four models, every boundary setting, original Figures 2.1–2.4, and three supplemental experiments                 |
| Sections 3.1–3.5.2     | Chapter 3       | **Paragraph-level complete**: five close-reading pages cover (3.1)–(3.68), Theorems 3.1–3.9, Figures 3.1–3.18, and both tables                               |
| Sections 4.1–4.6       | Chapter 4       | **Paragraph-level complete**: four close-reading pages cover (4.1)–(4.44), Theorems 4.1–4.9, Figures 4.1–4.22, and Table 4.1                                 |
| Section 5              | Chapter 5       | **Paragraph-level complete**: paper conclusions are covered; unified analysis, GPU work, and the experiment ledger are explicitly marked as site supplements |

> [!warning] Implementation coverage is not exception-free reproduction
> The [[en/computational-mathematics/knowledge-notes/time-parallelization/reproduction-audit-2026-08-24|whole-article experiment audit]] of August 24, 2026 confirms generators or computed redraws for all 45 numbered figures and local computations for two of the three tables. The two Figure 3.14 NKA curves remain unmatched, Table 4.1 is an external supercomputer result that is only quoted, and a current clean installation exposes two numerical-tolerance failures. “Paragraph-level complete” describes the site's coverage of the source text; it must not be read as exception-free numerical and engineering reproduction.

## Method map

| Method    | Parallel unit                             | Mechanism                              | Natural regime                 |
| --------- | ----------------------------------------- | -------------------------------------- | ------------------------------ |
| SWR       | overlapping space-time subdomains         | waveform transmission                  | transport and wave problems    |
| PIDC/RIDC | correction levels and time nodes          | deferred correction and pipelining     | initial-value problems         |
| ParaExp   | inhomogeneous and homogeneous subproblems | exponential propagation                | linear systems                 |
| ParaDiag  | all-at-once temporal matrix               | circulant approximation and FFT        | linear or linearized systems   |
| Parareal  | coarse time intervals                     | coarse prediction plus fine correction | moderate to strong dissipation |
| PFASST    | collocation nodes across time steps       | SDC with multilevel correction         | high-order integration         |
| MGRiT     | a hierarchy of temporal grids             | relaxation and coarse-grid correction  | long time intervals            |
| STMG      | the full space-time grid                  | space-time smoothing and coarsening    | parabolic systems              |

## Three organizing principles

### How causal constraints enter a parallel computation

For a one-step method, $u_{n+1}=\Phi_{\Delta t}(u_n)$ forms a sequential recurrence. A parallel-in-time method exposes the coupling among all temporal unknowns and constructs a parallel approximation to the inverse. Concurrency then occurs inside an iteration, decomposition, or transform.

### Dissipation determines whether a coarse representation is informative

Parabolic dynamics damp high-frequency error. A coarse propagator can therefore reproduce the slowly varying components that remain relevant over long intervals. Hyperbolic dynamics preserve phase information. Small phase errors in the coarse problem may then accumulate across time intervals.

### Iteration count is not parallel efficiency

A simplified cost model is

$$
T_{\mathrm{parallel}}
\approx K(C_G+C_F/P)+C_{\mathrm{comm}}+C_{\mathrm{setup}},
$$

where $K$ is the iteration count, $C_G$ and $C_F$ are the coarse and fine propagation costs, and $P$ is the temporal concurrency. The numerical experiments in these notes measure convergence, not end-to-end parallel speedup.

> [!note] Numerical provenance
> Original paper figures on these pages are extracted from the source PDF with their coordinates and panels intact. Experiments labeled “site reproduction” were regenerated on 31 July 2026 from the Python project on the new experiment server. The initial formal results used the SciPy CPU path. A subsequent CuPy/T4 hybrid backend batches the independent Burgers fine propagators on the GPU, preserves the Figure 4.5 stopping iterations within the earlier `paper_validation` entry, and reduces that combined suite from 263.57 to 67.92 seconds. Its Burgers flux differs from the current whole-article `section4` entry, and the timing is not the total wall time of all four whole-article groups.
>
> These T4 values are the archived single-GPU record from July 31, 2026. Its scope and the current dependency drift are documented in the [[en/computational-mathematics/knowledge-notes/time-parallelization/reproduction-audit-2026-08-24|reproduction audit]].

## Primary sources

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), Acta Numerica 34 (2025), pp. 385-489.
- Original MATLAB examples: [wushulin/ActaPinT](https://github.com/wushulin/ActaPinT).
- Python conversion, extensions, and formal results: [freezeng123456/ActaPinT-Python](https://github.com/freezeng123456/ActaPinT-Python).
