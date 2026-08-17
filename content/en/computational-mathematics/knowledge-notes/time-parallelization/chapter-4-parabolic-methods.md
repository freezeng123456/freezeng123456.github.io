---
title: "Chapter 4: PinT Methods Designed for Parabolic Problems"
description: Complete analysis of Parareal, PFASST, MGRiT, diagonalized Parareal, and space–time multigrid
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization/chapter-4-parabolic-methods
tags:
  - parallel-in-time
  - parabolic-PDE
---

> [!note] Reading scope
> This chapter covers Section 4 of the paper (pp. 443–481) and retains its exact hierarchy: the Section 4 introduction, 4.1 historical development, 4.2 Parareal, 4.3 PFASST, 4.4 MGRiT, the two diagonalization-based Parareal variants in 4.5, and 4.6 STMG. Python reproductions, parameter comparisons, and coverage audits are marked as site supplements and do not take paper section numbers. Theoretical contraction factors, measured iteration curves, and wall-clock performance carry distinct labels.

## Source-to-page map

This page retains the chapter-level synthesis and site reproductions. The equation-by-equation, theorem-by-theorem, and figure-by-figure arguments are split into the following pages:

| Paper section      | Source pages | Close-reading page                                                                                                                      | Coverage                                                 |
| ------------------ | ------------ | --------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------- |
| Section 4, 4.1–4.2 | pp. 443–452  | [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-4-1-parareal\|Historical development and Parareal]]         | (4.1)–(4.9), Theorems 4.1–4.4, Figures 4.1–4.5           |
| 4.3–4.4            | pp. 452–461  | [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-4-2-pfasst-mgrit\|PFASST and MGRiT]]                        | (4.10)–(4.13), Theorems 4.5–4.6, Figures 4.6–4.11        |
| 4.5                | pp. 460–471  | [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-4-3-diagonalized-parareal\|Diagonalization-based Parareal]] | (4.14)–(4.29), Theorems 4.7–4.8, Figures 4.12–4.17       |
| 4.6                | pp. 472–481  | [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-4-4-stmg\|Space–time multigrid]]                            | (4.30)–(4.44), Theorem 4.9, Figures 4.18–4.22, Table 4.1 |

## Section 4 introduction and 4.1 Historical development

### What parabolic methods exploit

The methods in Chapter 3 can address long-range propagation directly, but each has a constraint: SWR needs transmission operators, ParaExp is principally linear, and nonlinear ParaDiag solves Newton systems over long windows. Parabolic equations offer another route. Diffusion rapidly suppresses high-frequency error, so a late state becomes weakly dependent on fine details from the distant past. A coarse temporal model may still retain the dynamically important slow modes.

This chapter studies Parareal, PFASST, MGRiT, two forms of diagonalized Parareal, and space–time multigrid (STMG). The first three construct most of their hierarchy in time and normally retain one spatial grid. STMG coarsens space and time together and employs a block Jacobi smoother that is parallel across time.

![Coarse/fine hierarchy and parallel work in parabolic PinT methods](assets/diagrams/pint/en/parabolic-multilevel-map.svg)

Parareal has roots in multiple shooting, Nievergelt-style precomputation, and Saha's method. It subsequently inspired variants such as PITA, PFASST, and MGRiT. These methods communicate long-range information through a low-cost global approximation and repair fine details with concurrent local work. As a problem approaches the hyperbolic limit, phase and propagation errors on the coarse level cease to decay automatically, and convergence degrades.

## 4.2 Parareal

### Algorithm and concurrency

Partition $[0,T]$ into $N$ large intervals. Let $\mathcal F$ be an accurate, expensive fine propagator and $\mathcal G$ a cheap coarse propagator. Parareal updates

$$
U_{n+1}^{k+1}
=\mathcal G(U_n^{k+1})
+\mathcal F(U_n^k)-\mathcal G(U_n^k). \tag{4.1}
$$

All $\mathcal F(U_n^k)$ evaluations are concurrent within iteration $k$. The new values $\mathcal G(U_n^{k+1})$ still advance sequentially in $n$. The coarse term establishes the latest causal prediction, while $\mathcal F(U_n^k)-\mathcal G(U_n^k)$ repairs the previous iteration's fine/coarse discrepancy.

Parareal can be read as an approximate Newton method for multiple shooting or as coarse-propagator preconditioning of the lower-triangular all-at-once system. It is nonintrusive: existing fine and coarse integrators can participate if they expose an interface that advances one large interval from a supplied initial state.

### Linear modal analysis

For spatial mode $\lambda_\ell$, let one fine step have stability function $R_f(\lambda_\ell\Delta t)$, with $J$ fine steps per large interval, and let the coarse stability function be $R_g(\lambda_\ell\Delta T)$. Theorem 4.1 writes the error iteration as a strictly lower-triangular Toeplitz matrix. Strict triangularity gives finite-step exactness: exact arithmetic reaches the sequential fine solution in at most $N$ iterations, and the first $k$ large intervals are exact after iteration $k$.

Finite termination alone does not create speedup. If the tolerance requires nearly $N$ iterations, each containing a sequential coarse sweep, little temporal advantage remains.

Theorem 4.2 supplies two complementary estimates. For a small number of intervals, a combinatorial bound captures superlinear acceleration with $k$. Over a long horizon, the more informative modal linear factor is

$$
\rho_\ell
=\frac{|R_g-R_f^J|}{1-|R_g|}. \tag{4.5b}
$$

The numerator measures fine/coarse mismatch; the denominator measures the coarse propagator's damping margin. Figure 4.2 shows that a short-horizon superlinear phase and a long-horizon nearly linear phase can coexist.

### Nonlinear finite-step behavior

Under Lipschitz stability of the coarse propagator and a local error of order $p$, Theorem 4.3 establishes finite-step and superlinear bounds for nonlinear Parareal. The mechanism remains the interval-by-interval growth of the exact region. Nonlinearity makes the constants depend on solution regularity and propagator Lipschitz factors, so the observed iteration count can be highly sensitive to the window and physical parameters.

### The time integrator sets the limiting factor

Theorem 4.4 analyzes parabolic modes on the negative real axis. With an L-stable fine method, backward Euler coarse propagation, and enough fine steps per interval, the worst long-time factor is about $0.3$. If the fine method is only A-stable, high-frequency modes are insufficiently damped and the required $J$ grows with the most dangerous frequency. Figure 4.3 compares the trapezoidal rule with two SDIRK fine propagators. It demonstrates that the fine method and fine/coarse step ratio materially change the convergence region even after the coarse method has been fixed. Switching to Radau IIA as the **coarse** propagator at the continuous level (with $\mathcal F$ the exact propagator) reduces the theoretical worst factor to roughly $0.068$.

### From the heat equation toward the hyperbolic limit

Figures 4.4–4.5 use periodic ADE and Burgers problems with $T=4$, $\Delta T=0.1$, $\Delta x=1/128$, $J=32$, backward Euler coarse propagation, and a second-order L-stable SDIRK fine method. As viscosity falls, ADE develops increasing fine/coarse phase mismatch. Burgers' equation also changes local transport speed and shock position. Both slow markedly, and Burgers exhibits approximate divergence around $\nu\le10^{-3}$.

The wave equation is more severe. Parareal struggles to control phase unless coarse propagation is almost as accurate as fine propagation; such a coarse solver often loses its cost advantage. Semi-Lagrangian or phase-optimized coarse solvers help linear advection. A broadly effective coarse model for nonlinear hyperbolic dynamics remains open.

## 4.3 PFASST

### From SDC to a temporal hierarchy

PFASST combines Parareal-style concurrency across large intervals with high-order SDC collocation. On each time step, choose $M_f$ fine collocation nodes and write

$$
\boldsymbol U_f
=\boldsymbol U_{0,f}
+\Delta t\,Q_f\boldsymbol f(\boldsymbol U_f). \tag{4.10}
$$

Direct solution of the dense collocation system is expensive. SDC preconditions it with an easily solved lower-triangular approximation based on implicit Euler; each sweep removes part of the collocation residual. PFASST adds a coarse level with $M_c$ nodes and pipelines fine and coarse sweeps across different time steps.

Lagrange interpolation defines restriction and prolongation between fine and coarse nodes. The coarse level uses the full approximation scheme (FAS), which transfers the fine residual as a $\tau$ correction into the coarse collocation equation. A PFASST iteration contains concurrent fine sweeps, fine-to-coarse transfer, a sequential or pipelined coarse sweep, and coarse-to-fine correction. It can be understood through either Parareal or multigrid on collocation equations.

### Numerical observation and limitations

Figure 4.6 uses periodic heat and ADE problems with $T=3$, $\Delta x=1/128$, $\Delta t=1/64$, source parameter $\sigma=1000$, three Radau IIA nodes on the fine level, and two on the coarse level. The heat equation converges rapidly. As ADE viscosity decreases, high-frequency propagation across time steps is represented less faithfully on the coarse collocation level, and convergence slows.

PFASST is attractive when high temporal order is needed and each collocation solve is expensive. Its performance depends on nodes, SDC sweeps, coarse cost, node transfer, pipeline fill, and the allocation of spatial parallel resources. Formal collocation order alone does not establish parallel efficiency.

## 4.4 MGRiT

### F points, C points, and FCF relaxation

MGRiT marks every $J$th point of the time grid as a C point and labels the intervening points F. F relaxation performs fine propagation concurrently between neighboring C points. C relaxation updates coarse points, and coarse-grid correction carries information over long distances. A two-level FCF iteration performs F relaxation followed by C and another F relaxation, yielding an overlapping update related to Parareal.

One FCF iteration uses roughly two sets of fine propagation and therefore costs more than one Parareal iteration. The extra CF segment supplies overlap, so each iteration may make two large intervals exact; at most about $\lceil N/2\rceil$ iterations reach the sequential fine solution. General $F(CF)^\nu$ trades additional fine work for stronger contraction.

### Convergence factors and cost-fair comparison

Theorem 4.5 gives the long-time modal factor for two-level FCF:

$$
\rho_{\mathrm{MGRiT},\ell}
=|R_f^J|
\frac{|R_g-R_f^J|}{1-|R_g|}.
$$

It contains one extra $|R_f^J|$ relative to Parareal, so dissipative modes contract further. Each additional CF segment adds a similar fine-propagation factor and another set of fine solves.

A fair comparison places one FCF iteration beside two Parareal iterations with comparable fine work. Theorem 4.6 gives representative worst factors for L-stable fine methods. With backward Euler coarse propagation, Parareal and one FCF have factors about $0.2984$ and $0.1115$. A second-order Lobatto IIIC combination gives about $0.0817$ and $0.0197$. One FCF has the better per-iteration factor, yet it can be slightly worse than the square of the Parareal factor at equal fine-solve cost. Figure 4.8 maps this comparison in the complex plane.

### ADE and Burgers experiments

Figures 4.9–4.10 use $T=5$, $J=20$, $\Delta T=1/8$, $\Delta x=1/160$, backward Euler coarse propagation, and SDIRK22 fine propagation. Every heat mode lies well inside the convergent region. As ADE viscosity decreases from $0.1$ to $0.01$ and $0.002$, dangerous modes approach the high-phase, low-damping region. At $\nu=0.002$, Parareal and FCF have measured linear-stage factors around $1.4211$ and $1.2812$, so both exhibit transient growth.

MGRiT retains finite-step decay. That property states that the strictly triangular error eventually vanishes; it does not remove amplification during the preceding iterations. A method that spends the tolerance-relevant phase in transient growth has not become a scalable solver.

Figure 4.11 shows the same pattern for nonlinear Burgers. With adequate diffusion, FCF improves the factor per iteration. At matched fine-propagation cost, one FCF is often comparable to two Parareal iterations. Both degrade when coarse propagation fails to represent nonlinear shock motion.

## 4.5 Diagonalization-based Parareal

The paper diagonalizes two different parts of Parareal. The first replaces sequential coarse-grid correction across the $N$ large intervals by an FFT-diagonalizable head–tail system. The second constructs a cheap coarse propagator from an $\alpha$-circulant system over the $J$ fine points inside each large interval.

### 4.5.1 Diagonalization-based CGC

Standard Parareal applies its coarse correction sequentially in time. The diagonalized variant introduces

$$
U_0^{k+1}=u_0+\alpha U_N^{k+1}
$$

and modifies the correction right-hand side consistently so the original initial-value problem remains the fixed point. With linear backward Euler coarse propagation, this produces an $\alpha$-circulant all-at-once matrix that FFTs reduce to $N$ independent shifted coarse spatial systems.

As $\alpha\to0$, the standard sequential coarse correction is recovered, while transform roundoff grows. Theorem 4.7 states that if ordinary Parareal has factor $\rho$, then

$$
\alpha\le \frac{\rho}{1+\rho}
$$

retains the same convergence scale. The result is proved for negative real spectra; for complex spectra there is only numerical evidence. In Figure 4.12, the heat equation has $\rho\approx0.22$ and threshold about $0.18$; ADE at $\nu=0.1$ (complex spectrum) has $\rho\approx0.39$ and threshold about $0.28$, in agreement with the prediction.

The nonlinear version applies quasi-Newton to the all-at-once coarse equation and builds a block $\alpha$-circulant system from an average Jacobian. The Burgers results in Figure 4.13 exhibit the same threshold behavior.

Direct insertion of the head–tail condition into MGRiT changes its fixed point and can diverge. The paper supplies a consistent condition in which head–tail coupling uses the difference between the latest two iterates. At small $\alpha$, it retains the original MGRiT convergence rate and can also be used in Parareal.

### 4.5.2 Diagonalization-based coarse solver

The second construction uses the same integrator and small step for fine and coarse propagation. Fine propagation performs $J$ sequential steps. Coarse propagation puts the $J$ unknown time points into one $\alpha$-circulant all-at-once system and solves them concurrently. At $\alpha=0$, coarse and fine propagators coincide, so Parareal converges in one iteration but offers no cost advantage. At $\alpha>0$, the diagonalizable approximation can be roughly $J$ times cheaper in parallel.

Theorem 4.8 exposes a sharp equation-class distinction. A parabolic problem with negative real spectrum has

$$
\rho=\alpha,
$$

independent of the number $N$ of large intervals. A wave problem with imaginary spectrum satisfies

$$
\rho\le \frac{2\alpha N}{1+\alpha}.
$$

Figure 4.14 verifies the $\alpha$ slope for the heat equation. Figure 4.15 shows wave convergence slowing with $N$ at $\alpha=0.01$. A superlinear phase remains: increasing $N$ from 24 to 960 costs only about two additional iterations to reach discretization accuracy. The linear bound is sharpest when $\alpha N$ is small.

Figures 4.16–4.17 apply the method to Burgers. The factor is approximately linear in $\alpha$, and $\alpha=10^{-3}$ is robust with respect to the interval count. Smaller $\alpha$ also reduces viscosity sensitivity. This second diagonalized coarse propagator works across parabolic and hyperbolic cases and has broader scope than the first variant.

## 4.6 Space–time multigrid

### All-at-once system and block Jacobi smoothing

For $\boldsymbol u'=A\boldsymbol u+\boldsymbol g$, a general one-step integrator reads

$$
r_1\boldsymbol u_{n+1}=r_2\boldsymbol u_n+\widetilde{\boldsymbol f}_n,
\qquad n=0,\ldots,N_t-1, \tag{4.30}
$$

where $r_1,r_2$ are matrix polynomials in $\Delta tA$. Stacking all time points gives the all-at-once system

$$
K\boldsymbol U
=\begin{bmatrix}
r_1\\-r_2&r_1\\&\ddots&\ddots\\&&-r_2&r_1
\end{bmatrix}
\boldsymbol U
=\boldsymbol b. \tag{4.31}
$$

STMG applies damped block Jacobi:

$$
\boldsymbol U^{j+1}
=\boldsymbol U^j
+\eta\left(I_t\otimes r_1\right)^{-1}
(\boldsymbol b-K\boldsymbol U^j),
\qquad r_1=I_x-\theta\Delta t A. \tag{4.32}
$$

$I_t\otimes r_1$ is block diagonal in time, so every spatial block can be solved concurrently. The damping parameter $\eta$ controls removal of high-frequency error.

A two-level cycle performs $s_1$ presmoothing steps, restricts the residual in space and time, solves on a grid with $\Delta X=2\Delta x$ and $\Delta T=2\Delta t$, prolongs the correction, and applies $s_2$ postsmoothing steps. Recursive use yields full STMG and assigns long-wave error in both dimensions to coarse levels.

Earlier parabolic multigrid often used point Gauss–Seidel in the time direction. That smoother is sequential forward substitution: it coarsens space effectively but does not scale in time. Time-parallel block Jacobi is a defining ingredient of modern STMG.

### Local Fourier analysis

For the one-dimensional heat equation, centered spatial differences, and backward Euler, local Fourier analysis decomposes error into space–time frequencies. Theorem 4.9 gives the optimal damping

$$
\eta=\frac12.
$$

High temporal frequencies are damped by at most $1/\sqrt2$, permitting temporal coarsening. For the normalized heat equation in the paper, $\Delta t/\Delta x^2\ge1/\sqrt2$ places high spatial frequencies under the same bound and permits spatial coarsening as well. Restoring a diffusion coefficient gives the nondimensional ratio $\nu\Delta t/\Delta x^2$. The ADE symbol contains an imaginary part; the paper still identifies $\eta=1/2$ as a sound starting point for backward Euler.

### Damping, smoothing count, and the time integrator

Figure 4.19 scans $\eta$ for two-level backward-Euler STMG. Both the heat equation and ADE with $\nu=0.01$ perform well near $1/2$. Figure 4.20 fixes $\eta=1/2$ and compares one versus three block Jacobi smoothing steps. Additional smoothing raises the cost of a cycle and materially reduces the cycle count. ADE remains slower than heat, although its viscosity sensitivity weakens with more smoothing and a superlinear phase appears.

Figure 4.21 changes the time integrator to the trapezoidal rule. The heat equation has convergence difficulty across the damping scan even with ten smoothing steps. ADE converges and improves with more smoothing; the better sampled damping is around $0.8$. STMG therefore depends on the stability function of the time integrator. The backward-Euler damping result cannot simply be transferred to the trapezoidal rule.

### Large-scale scaling results

Table 4.1 reports modern STMG on a three-dimensional heat equation. In weak scaling, cores increase from 1 to 262,144, time steps from 2 to 524,288, and degrees of freedom from 59,768 to 15,667,822,592. The iteration count remains seven and total time changes only from about 28.8 to 30.0 seconds. The reference forward solve, parallel only in space, grows to roughly 4,988,060 seconds.

For strong scaling, a problem with 512 time steps and 15,300,608 unknowns falls from 7,635.2 seconds on one core to 30.0 seconds on 256 cores. A larger case with 524,288 time steps falls from 15,205.9 seconds on 512 cores to 30.0 seconds on 262,144 cores. These data come from the three-dimensional implementation cited by the paper and demonstrate the weak- and strong-scaling potential of STMG for parabolic systems.

### Nonlinear FAS-STMG

For $\boldsymbol u'=\boldsymbol f(\boldsymbol u)$, the $\theta$ method gives

$$
K(\boldsymbol U)
=(B\otimes I_x)\boldsymbol U
-\Delta t(\widetilde B\otimes I_x)\boldsymbol f(\boldsymbol U)
=\boldsymbol b. \tag{4.42}
$$

A nonlinear block Jacobi smoother solves independent temporal-block corrections, for example by an inner Newton iteration. Local Fourier analysis no longer provides a direct optimized damping value, so the parameter is problem dependent. The coarse level uses FAS: restrict the current approximation and residual, solve a nonlinear coarse equation with consistency correction, prolong the coarse-solution difference, and postsmooth.

Figure 4.22 uses two block Jacobi smoothing steps for Burgers' equation. Nonlinear STMG converges rapidly with adequate diffusion and deteriorates as viscosity falls. The best experimental damping is $\eta=1/4$, illustrating the difference from linear backward-Euler theory.

Overall, STMG is among the most effective PinT solvers currently available for parabolic problems and has demonstrated large-scale scalability. It is more intrusive than Parareal because it needs access to the all-at-once discretization, smoother, and grid transfers. Robustness for hyperbolic problems and across time integrators remains an active issue.

## Site numerical supplement: Python reproduction results

### Parareal baselines and Figure 4.5

| Problem             | Formal parameters                     | Iterations |   Final maximum error |
| ------------------- | ------------------------------------- | ---------: | --------------------: |
| advection–diffusion | $N_x=128$, $N=40$, $J=32$, $\nu=0.02$ |         39 | $6.106\times10^{-15}$ |
| Burgers             | $N_x=128$, $N=40$, $J=32$, $\nu=1$    |         16 | $3.544\times10^{-12}$ |

![Parareal convergence baseline for advection–diffusion](assets/pint/parareal-ade-baseline.svg)

![Parareal convergence baseline for viscous Burgers](assets/pint/parareal-burgers-baseline.svg)

With $T=4$, $\Delta T=0.1$, $\Delta x=1/128$, and $J=32$, the iteration counts required to reduce maximum error below $10^{-10}$ are:

| Equation            | $\nu=1$ | $\nu=0.1$ | $\nu=0.02$ |
| ------------------- | ------: | --------: | ---------: |
| advection–diffusion |      14 |        24 |         35 |
| Burgers             |      14 |        21 |         25 |

![Parareal convergence for ADE and Burgers as diffusion weakens](assets/pint/parareal-figure-4-5.svg)

These counts measure convergence toward the sequential fine solution. GPU performance appears in [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-5-unified-view#site-reproduction-gpu-acceleration-and-profiling|Chapter 5]].

### MGRiT baseline and Figures 4.9–4.10

The near-hyperbolic baseline uses $N_x=160$, $N=40$, $J=20$, $T=5$, and $\nu=0.002$. At the end of the reported window, Parareal has maximum error $2.895\times10^2$, while two-level MGRiT eventually reaches $5.551\times10^{-16}$ through finite-step exactness.

![Baseline comparison of Parareal and MGRiT on near-hyperbolic ADE](assets/pint/mgrit-baseline.svg)

For cost comparison between one FCF and two Parareal iterations, the modal long-time factors are:

| Problem          | Parareal factor | FCF-MGRiT factor |
| ---------------- | --------------: | ---------------: |
| heat             |          0.2824 |           0.0835 |
| ADE, $\nu=0.1$   |          0.4453 |           0.2719 |
| ADE, $\nu=0.01$  |          1.0501 |           0.9021 |
| ADE, $\nu=0.002$ |          1.4211 |           1.2812 |

![Modal long-time factors for Heat and ADE at several viscosities](assets/pint/mgrit-figure-4-9.svg)

![Equal-fine-solve-cost convergence of Parareal and FCF-MGRiT](assets/pint/mgrit-figure-4-10.svg)

> [!note] One numerical difference from Figure 4.9
> The source figure annotates the maximum Parareal factor at $\nu=0.01$ as $0.9986$. Direct evaluation in the Python conversion, using equation (4.5b), the printed parameters, and the stability function in upstream `MGRiT_Heat_ADE.m`, gives $1.0501$. The MGRiT value $0.9021$ for the same case and the other two pairs agree closely with the source figure. This page retains the computed reproduction value and records the discrepancy explicitly instead of overwriting it with the source annotation. Both values lie near one and support the same qualitative conclusion: long-time convergence is already close to its critical regime at this viscosity.

### STMG damping and smoothing validation

The trapezoidal-rule baseline uses $N_x=N_t=255$, $\nu=10^{-3}$, and three pre- and postsmoothing steps. In the sampled scan, the lowest error after 15 cycles occurs near $\eta=0.98$.

![STMG damping scan with the trapezoidal rule](assets/pint/stmg-baseline.svg)

The paper-grid backward-Euler validation gives:

| Problem         | Best sampled $\eta$ after 15 iterations |
| --------------- | --------------------------------------: |
| heat            |                                   0.500 |
| ADE, $\nu=0.01$ |                                   0.372 |

![Damping scan for backward-Euler STMG](assets/pint/stmg-figure-4-19.svg)

The finite-grid minimum for ADE is $0.372$, while the paper's analysis and plot use $1/2$ as a robust choice. The two statements concern different quantities: error after 15 cycles on one finite grid and a high-frequency smoothing bound.

At fixed $\eta=0.5$, three pre- and postsmoothing steps require fewer cycles than one:

![STMG convergence with one versus three pre- and postsmoothing steps](assets/pint/stmg-figure-4-20.svg)

A performance decision should compare

$$
\text{work to tolerance}
=(\text{cost per cycle})\times(\text{cycle count})
+\text{communication and memory cost}.
$$

## Site method comparison: parameters and failure modes

| Method                | Key parameter                                   | Tradeoff controlled                            | Typical failure mode                                                      |
| --------------------- | ----------------------------------------------- | ---------------------------------------------- | ------------------------------------------------------------------------- |
| Parareal              | interval count $N$, step ratio, $\mathcal G$    | concurrency, coarse cost, propagation accuracy | iteration count approaches $N$; phase error grows                         |
| PFASST                | collocation nodes, SDC sweeps, coarse nodes     | high order, pipeline depth, coarse cost        | coarse collocation misrepresents weakly damped propagation                |
| MGRiT                 | coarsening, F/FCF/$F(CF)^\nu$                   | overlapping contraction versus fine work       | no equal-cost gain; transient amplification                               |
| diagonalized Parareal | $\alpha$, location of diagonalization           | sequential coarse fraction versus roundoff     | small $\alpha$ is ill-conditioned; large $\alpha$ weakens convergence     |
| STMG                  | $\eta$, smoothing count, coarsening, integrator | high-frequency smoothing versus cycle cost     | integrator mismatch; degradation for low viscosity or hyperbolic dynamics |

## Source coverage audit

| Source location                 | This page | Material covered                                                                                                |
| ------------------------------- | --------- | --------------------------------------------------------------------------------------------------------------- |
| Sections 4 and 4.1, pp. 443–444 | 4.1       | parabolic temporal locality, method scope, history, and the distinction of STMG                                 |
| Section 4.2, pp. 444–452        | 4.2       | Parareal update, Theorems 4.1–4.4, Figures 4.1–4.5, nonlinear behavior and hyperbolic degradation               |
| Section 4.3, pp. 452–455        | 4.3       | collocation equations, SDC, FAS transfer, PFASST iteration, Figure 4.6                                          |
| Section 4.4, pp. 455–460        | 4.4       | FCF structure, Theorems 4.5–4.6, Figures 4.7–4.11, cost-fair comparison, Burgers                                |
| Section 4.5, pp. 460–471        | 4.5       | two diagonalization locations, Theorems 4.7–4.8, Figures 4.12–4.17, nonlinearity and consistent MGRiT condition |
| Section 4.6, pp. 472–481        | 4.6       | all-at-once STMG, block Jacobi, Theorem 4.9, Figures 4.18–4.22, Table 4.1, nonlinear FAS                        |

## Summary

These methods turn diffusion-induced temporal locality into an algorithmic advantage. Parareal corrects a cheap global prediction with concurrent fine solves. PFASST combines SDC collocation with FAS. MGRiT strengthens coarse-grid correction through temporal hierarchy and overlapping relaxation. Diagonalized variants reduce the sequential coarse component. STMG treats spatial and temporal scales together. Their parabolic performance can be excellent, while still depending on spectral matching, time integrators, cost-fair comparisons, and complete hardware overhead. As viscosity falls, phase and long-range memory become dominant and standard coarse-grid mechanisms progressively lose their advantage.

## Source

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Section 4, pp. 443–481.
