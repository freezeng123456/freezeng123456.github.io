---
title: ActaPinT-Python Whole-Article Reproduction Audit
description: Chapter-by-chapter audit of implementation coverage, paper agreement, test health, T4 evidence, and remaining gaps
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization/reproduction-audit-2026-08-24
tags:
  - parallel-in-time
  - reproducible-computing
  - numerical-experiments
---

This page audits commit [`9bf1293`](https://github.com/freezeng123456/ActaPinT-Python/commit/9bf12930a5fbd0403abf8a436c511a90de4025fd) of [ActaPinT-Python](https://github.com/freezeng123456/ActaPinT-Python), against M. J. Gander, S.-L. Wu, and T. Zhou, _Time Parallelization for Hyperbolic and Parabolic Problems_, and the 24 MATLAB scripts pinned at upstream commit `68431b7`. It distinguishes the existence of an entry point, successful artifact generation, quantitative agreement with the paper, and parallel-performance evidence.

> [!summary] Verdict
> The project supplies a numerical generator or computed redraw for all 45 numbered figures and local computations for two of the three tables. All 24 upstream MATLAB scripts have been ported. Its **implementation coverage is therefore close to whole-article complete**, but “every experiment has been reproduced” is too strong: the nearest-Kronecker curves of Figure 3.14 at $\nu=0.02$ remain unmatched, Table 4.1 is a 262,144-core result from a different three-dimensional parallel implementation and is only quoted, and a clean installation with the current unlocked dependency range fails two numerical-tolerance tests.

The machine-readable audit is [[assets/pint/data/full_reproduction_audit_20260824.json|full_reproduction_audit_20260824.json]].

## 1. Evidence levels

The audit uses four levels:

1. **entry-point coverage**: the figure or table has an independent CLI entry;
2. **artifact coverage**: the entry writes parameter/metric JSON and an SVG, PNG, or table;
3. **paper agreement**: the output is compared with printed values, digitized markers, or pixel statistics;
4. **performance reproduction**: hardware, precision, wall time, and a baseline are reported in addition to numerical correctness.

A visually similar curve reaches only the second level. Copying Table 4.1 into the site does not mean that its fourth-level scaling experiment was rerun.

## 2. Whole-article coverage

| Chapter   |                                Paper content |                                      Local entries | Audit result                                                                             |
| --------- | -------------------------------------------: | -------------------------------------------------: | ---------------------------------------------------------------------------------------- |
| Chapter 1 |                                  1 schematic |                                                  1 | Figure 1.1 has no numerical experiment; it is redrawn from computed marker heights       |
| Chapter 2 |    4 numerical figures, 24 space–time panels |                                                  4 | all implemented and compared panel by panel with the paper rasters                       |
| Chapter 3 | 15 numerical figures, 3 schematics, 2 tables |                                                 20 | all have entries; the $\nu=0.02$ NKA curves of Figure 3.14 are only partially reproduced |
| Chapter 4 |  20 numerical figures, 2 schematics, 1 table |        22 figure entries; no local table generator | Figures 4.1–4.22 are covered; Table 4.1 is quoted only                                   |
| Chapter 5 |           no figures, tables, or experiments |                                     not applicable | prose conclusion, with no experimental gap                                               |
| **Total** |         **45 figures + 3 tables = 48 items** | **47 items have a local generator or computation** | **implementation coverage 47/48, with one local curve gap and one external-table gap**   |

The 47/48 count measures the existence of a local computation, not perfect pointwise agreement. Figure 3.14 has an entry and is therefore counted among the 47, while its unmatched curves remain explicit exceptions.

## 3. Numerical evidence by chapter

### Chapter 2: model problems

All 24 panels of Figures 2.1–2.4 are generated for the heat, advection–diffusion, Burgers, and wave equations. The comparison is calibrated against paper panels rather than based on visual inspection: the overall median smoothed correlation is $0.996$ and the overall median decile mismatch is $0.0026$. Matching the shock location in Figure 2.3 requires $(u^2)_x$ instead of the $\frac12(u^2)_x$ printed in both the article and the upstream script. This is a diagnosed source inconsistency, not a parameter change that should be hidden.

### Chapter 3: methods effective for hyperbolic problems

SWR, UTP, PIDC/RIDC, ParaExp, ParaDiag I/BVM/NKA, and ParaDiag II all have independent entries. Representative quantitative results are:

| Result                                                           |                                                                                                Reproduced value |
| ---------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------: |
| Figure 3.15, heat GMRES                                          |                                                                                                    2 iterations |
| Figure 3.15, ADE at $\nu=10^{-3}$                                |                                                                                                    3 iterations |
| Figure 3.15, ADE at $\nu=10^{-6}$                                |                                                                                                   13 iterations |
| Whole-article `section3` entry, Figure 3.15 wave GMRES callbacks |                                                               115 total; first below $10^{-11}$ at callback 103 |
| Earlier `paper_validation` entry, same wave case                 | below preconditioned residual $10^{-11}$ at iteration 89; reaches its SciPy true-residual stop at iteration 103 |

All Figure 3.1 OSWR counts, the twelve condition numbers in Table 3.1, and the main Figure 3.15 spectral extrema are also tested. The sequential Jacobian-solve column of Table 3.2 is reproduced exactly.

### Chapter 4: methods designed for parabolic problems

Every numbered figure for Parareal, PFASST, MGRiT, both diagonalized-Parareal variants, and STMG has an entry. Headline results are:

| Result                                                          | Reproduced value                                                      |
| --------------------------------------------------------------- | --------------------------------------------------------------------- |
| Whole-article `section4`, Figure 4.5 ADE reaches $10^{-10}$     | 14 / 24 / 35 iterations                                               |
| Whole-article `section4`, Figure 4.5 Burgers reaches $10^{-10}$ | 14 / 24 / 36 iterations with $(u^2)_x$, which matches the paper panel |
| Earlier GPU `paper_validation`, Figure 4.5 Burgers              | 14 / 21 / 25 iterations with the upstream-script $\frac12(u^2)_x$     |
| Figure 4.9, long-time factors at $\nu=0.002$                    | Parareal 1.421082; MGRiT 1.281200                                     |
| Figure 4.19, best sampled damping                               | Heat 0.500; ADE ($\nu=0.01$) 0.372                                    |

The 124 digitized error markers of Figure 4.10 are reproduced within 6% overall. Figure 4.9 prints a maximum Parareal factor of 0.9986 at $\nu=0.01$, while the formula and the upstream stability function give about 1.0501; the project retains both and treats this as a source conflict.

This also reveals that the project's two validation paths are not interchangeable. The earlier `paper_validation` command is a representative subset that predates the whole-article groups; `section3` and `section4` are the later full entries. They differ in the wave-GMRES callback/stopping convention and in the Burgers flux used for Figure 4.5. The README and older site text placed some of these counts side by side without naming the entry. This audit always records the entry and flux.

## 4. Items that prevent an unconditional “complete reproduction” claim

| Item                                   | Current status                                                 | Why it is not complete                                                                                                                                    |
| -------------------------------------- | -------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Figure 3.14 NKA at $\nu=0.02$          | average-Jacobian curves match pointwise; NKA speed-up does not | the weights from (3.48) are nearly one, and direct optimization of all 26 diagonal weights still reaches only 0.389 rather than the published 0.280/0.161 |
| Table 4.1 STMG strong/weak scaling     | the site quotes the source table                               | it belongs to the Gander–Neumüller (2016) three-dimensional implementation and historical machine, not to this Python code                                |
| Rough-source panels of Figures 3.5/3.6 | trend and parameters match                                     | empirical scale factors 2.447 and 1.469 remain necessary; the printed error definition does not produce them                                              |
| Low-step end of Figure 3.11(a)         | main trend matches                                             | a bounded 0–3.5-decade roundoff difference remains                                                                                                        |
| Final markers of Figure 4.10           | almost all markers match                                       | the last two or three heat and $\nu=0.1$ points lie at the roundoff floor and flatten in the paper while Python continues to fall                         |

These bounded and documented differences do not invalidate the corresponding algorithm entries, but they rule out the phrase “pointwise, exception-free reproduction of the whole article.”

## 5. Current engineering reproducibility

On August 24, 2026, a clean installation from the open ranges in `pyproject.toml` selected Python 3.12.4, NumPy 2.5.2, SciPy 1.18.1, and Matplotlib 3.11.1. Pytest collected 104 cases: 101 passed, one CUDA case was skipped on a non-GPU host, and two failed.

| Failed test                                  |                                                                     Observation | Assessment                                                                                                                                               |
| -------------------------------------------- | ------------------------------------------------------------------------------: | -------------------------------------------------------------------------------------------------------------------------------------------------------- |
| fixed-ratio endpoint of Figure 3.10          | 0.217770 versus a digitized paper marker near 0.18 and a 20% relative tolerance | exceeds the upper edge by about 0.00177; a brittle digitization/version tolerance                                                                        |
| purely imaginary spectrum of the wave system |      maximum spurious real part $3.24\times10^{-7}$ at spectral radius about 40 | relative leakage about $8\times10^{-9}$; the analytic structure is unchanged, but the absolute $10^{-8}$ threshold is too strict for the new eigensolver |

The latest public [GitHub Actions run](https://github.com/freezeng123456/ActaPinT-Python/actions/runs/30839901657) also fails in its Python 3.12 job. The project declares lower bounds such as `numpy>=1.26` and `scipy>=1.11` but no lock or upper bounds. A previously passing environment therefore does not guarantee that a future clean installation passes. A verified lock file and scale-aware numerical tests are needed.

Rechecking on the same machine with the historically recorded NumPy 2.4.6 and SciPy 1.17.1 makes the Figure 3.10 tolerance test pass, while the absolute pure-imaginary-spectrum threshold still fails in the macOS arm64 dense eigensolver. This separates the two causes: the first is primarily dependency drift, while the second also needs a platform-independent, scale-aware assertion.

All four formal groups nevertheless complete in the same open-dependency environment:

| Group        | Local computations |                      Wall time |
| ------------ | -----------------: | -----------------------------: |
| `schematics` |                  6 |                        16.70 s |
| `section2`   |                  4 |                        38.13 s |
| `section3`   |                 17 |                       784.58 s |
| `section4`   |                 20 |                       629.08 s |
| **Total**    |             **47** | **1468.49 s (24 min 28.49 s)** |

The run writes 47 SVGs, 47 PNGs, and 48 JSON files—the extra JSON is the Chapter 2 summary—for about 30 MiB. The largest absolute numerical difference from the committed Chapter 2 metrics is $4.17\times10^{-11}$; Figure 3.1 counts and the Table 3.2 sequential column are identical; Figure 4.9 factors differ by at most about $8.2\times10^{-13}$. Divergent and roundoff-amplified Chapter 3 tails cannot be summarized responsibly with one maximum over every JSON leaf, and the project reports already exclude those quantities from stable headline metrics.

## 6. Correct scope of the T4 evidence

The existing July 31, 2026 Tesla T4 double-precision record validates a hybrid CPU/GPU backend. Forty independent Burgers fine propagators are batched on the GPU; causal coarse propagation and the rest of the paper suite remain on the CPU.

| T4 double-precision test                 |      CPU |     GPU | Speedup |
| ---------------------------------------- | -------: | ------: | ------: |
| 40 Burgers fine propagators, 32 substeps |  2.893 s | 0.246 s |  11.76× |
| combined paper-validation suite          | 263.57 s | 67.92 s |   3.88× |

The maximum absolute CPU/GPU difference for one batch is $2.33\times10^{-15}$, and Figure 4.5 stopping iterations within the earlier `paper_validation` entry are unchanged. That GPU record uses the upstream-script $\frac12(u^2)_x$ and must not be conflated with the 14/24/36 counts from the whole-article `section4` entry, which uses $(u^2)_x$ to match the paper panel. The record demonstrates numerical agreement and end-to-end benefit for one single-T4 hotspot. It does not establish multi-GPU, MPI, or Table 4.1 scaling.

## 7. Final assessment

- If the question is whether every chapter has a code path, **yes**: Chapters 1 and 5 contain no numerical experiments, and every numbered figure in Chapters 2–4 has a local entry or computed redraw.
- If the question is whether every published experimental value is reproduced without exception, **no**: the two Figure 3.14 NKA curves and Table 4.1 are explicit gaps, with three additional bounded discrepancies.
- If the question is whether any current clean environment is guaranteed to pass in one command, **no**: the latest CI and this unlocked-dependency installation both expose two numerical-tolerance failures.
- If the question is whether the GPU evidence validates whole-article parallel scaling, **no**: it covers a hybrid single-GPU Burgers batching hotspot only.

The most accurate description is therefore: **whole-article implementation coverage is essentially complete and most figures and tables have quantitative matches; one local curve family, one external supercomputer table, and present-day dependency drift remain unresolved.**

## Primary sources

- [Paper DOI](https://doi.org/10.1017/S0962492924000072)
- [Original MATLAB repository](https://github.com/wushulin/ActaPinT)
- [Python reproduction repository](https://github.com/freezeng123456/ActaPinT-Python)
- [Whole-article reproduction index](https://github.com/freezeng123456/ActaPinT-Python/blob/main/docs/PAPER_REPRODUCTION.md)
- [Per-item results](https://github.com/freezeng123456/ActaPinT-Python/blob/main/docs/RESULTS.md)
