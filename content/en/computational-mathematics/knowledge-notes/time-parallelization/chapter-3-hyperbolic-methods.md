---
title: "Chapter 3: Effective PinT Methods for Hyperbolic Problems"
description: Complete derivations and numerical interpretation of SWR, parallel deferred correction, ParaExp, and ParaDiag
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization/chapter-3-hyperbolic-methods
tags:
  - parallel-in-time
  - hyperbolic-PDE
---

> [!note] Reading scope
> This chapter covers Section 3 of the paper (pp. 396–443) and retains its exact hierarchy: the Section 3 introduction, 3.1 historical development, 3.2 SWR, 3.3 IDC, 3.4 ParaExp, and the two ParaDiag branches in 3.5.1–3.5.2. Equations, theorems, and paper experiments follow the source argument. Python results, parameter comparisons, and coverage audits are marked as site supplements and do not take paper section numbers.

## Source-to-page map

This page retains the chapter-level synthesis and the site's reproduction. The complete equation-by-equation, theorem-by-theorem, and figure-by-figure reading is split into the following pages:

| Paper section      | Source pages | Close-reading page                                                                                                                                        | Coverage                                                          |
| ------------------ | ------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------- |
| Section 3, 3.1–3.2 | pp. 396–405  | [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-3-1-history-and-swr\|Historical development and Schwarz waveform relaxation]] | (3.1)–(3.4), Theorems 3.1–3.2, Figures 3.1–3.3                    |
| 3.3                | pp. 405–411  | [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-3-2-idc\|Parallel integral deferred correction]]                              | (3.5)–(3.12), Theorem 3.3, Figures 3.4–3.6                        |
| 3.4                | pp. 412–415  | [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-3-3-paraexp\|ParaExp]]                                                        | (3.13)–(3.21), Theorem 3.4, Figures 3.7–3.8                       |
| 3.5, 3.5.1         | pp. 415–430  | [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-3-4-paradiag-i\|Direct ParaDiag]]                                             | (3.22)–(3.48), Theorems 3.5–3.7, Figures 3.9–3.14, Tables 3.1–3.2 |
| 3.5.2              | pp. 431–442  | [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-3-5-paradiag-ii\|Iterative ParaDiag]]                                         | (3.49)–(3.68), Theorems 3.8–3.9, Figures 3.15–3.18                |

## Section 3 introduction and 3.1 Historical development

### Why these methods can handle long-range propagation

Hyperbolic equations carry fine structures along characteristics over long distances in time. An effective PinT method must preserve path and phase information or solve the coupling across the full time domain directly. The paper places four families in this category: Schwarz waveform relaxation (SWR), parallel integral deferred correction (IDC), ParaExp, and ParaDiag.

Their origins differ. SWR combines domain decomposition with waveform relaxation. IDC grows out of defect correction and high-order integration. ParaExp exploits exponential propagation in linear systems. ParaDiag uses a diagonalizable all-at-once time matrix. Most of these techniques also work on parabolic problems. Their common feature is a parallel mechanism that does not require a strongly dissipative coarse propagator.

Historically, SWR connects closely to mapped and unmapped tent pitching. Parallel IDC includes window-level PIDC and correction-level RIDC pipelines. ParaExp began as a linear direct construction and later acquired nonlinear iterative variants. ParaDiag developed direct, waveform-relaxation, stationary-iteration, and Krylov-preconditioned forms.

## 3.2 Schwarz waveform relaxation

### From classical waveform relaxation to space–time decomposition

For $\boldsymbol u'=A\boldsymbol u+\boldsymbol g$, classical waveform relaxation splits $A=M+N$ and iterates

$$
\boldsymbol u^{k+1\prime}-M\boldsymbol u^{k+1}
=N\boldsymbol u^k+\boldsymbol g.
$$

A block-diagonal or colored $M$ permits concurrent component solves. The performance depends directly on the algebraic splitting, and a poor choice may converge slowly or diverge. Classical spatial domain decomposition also solves a new elliptic problem at each time step and normally ties every subdomain to the same temporal discretization.

SWR decomposes the continuous spatial domain first and solves an entire time window on each subdomain. Neighboring subdomains exchange a time-dependent interface waveform over overlaps or artificial boundaries. Each subdomain may use a locally suitable discretization, and the transmission condition can be designed from the PDE's propagation mechanism.

![Schwarz waveform relaxation exchanges interface waveforms over a complete time window](assets/diagrams/pint/en/schwarz-waveform-relaxation.svg)

The iterated unknown is an interface function. Dirichlet transmission communicates only solution values. Robin, Ventcel, and convolution conditions also approximate normal fluxes or a fuller Dirichlet-to-Neumann map. Versions with optimized transmission conditions are usually called OSWR.

### 3.2.1 First-order parabolic problems: overlap and Robin optimization

The paper illustrates the mechanism on a one-dimensional advection–diffusion equation. Two subdomains overlap by a width $l$. At each iteration, both complete time-window problems are solved in parallel with Robin transmission

$$
(\partial_x+p)u_1^{k+1}=(\partial_x+p)u_2^k,
\qquad
(\partial_x-p)u_2^{k+1}=(\partial_x-p)u_1^k. \tag{3.1}
$$

The parameter $p>0$ controls the interface operator; $p\to\infty$ gives Dirichlet exchange. The analysis proceeds in Laplace/Fourier space. Theorem 3.1 supplies an optimized Robin parameter and a worst-frequency estimate for the continuous problem. The Dirichlet case satisfies a bound of the form

$$
\rho_{\mathrm D}\leq
\exp\!\left(-\frac{l\pi}{\nu T}\right).
$$

The bound displays three trends. More overlap accelerates convergence. A longer time window increases the interface-coupling burden. Smaller $\nu$ makes directed transport carry data through the overlap faster, so SWR can improve for this particular model as diffusion decreases.

Figure 3.1 uses $L=8.2$, $T=5$, $\Delta t=0.01$, $\Delta x=0.02$, $l=2\Delta x$, and a Gaussian initial condition. The curves confirm faster convergence at lower viscosity and a strong advantage for optimized Robin transmission. In the four-subdomain test at $\nu=0.1$, Dirichlet and optimized Robin transmission require approximately 92 and 28 iterations. The corresponding continuous two-subdomain predictions are about 32 and 4. Continuous versus discrete operators, two versus four subdomains, and the boundary configuration account for the difference; the theoretical estimate should not be treated as an exact iteration forecast for that discrete experiment.

Higher-order Ventcel conditions can improve asymptotic convergence further. A temporal convolution can approximate an exact transparent boundary and deliver a mesh-independent theoretical factor, at the cost of a more complicated and temporally nonlocal interface operator.

### 3.2.2 Second-order hyperbolic problems: finite-step propagation

For the wave equation, SWR has a direct geometric interpretation. Each iteration advances correct interface information across one finite propagation distance. Theorem 3.2 states that Dirichlet transmission on two overlapping subdomains becomes exact over the full window once

$$
k>\frac{Tc}{(\beta-\alpha)L},
$$

where the physical overlap is $(\beta-\alpha)L$ and $c$ is the wave
speed. The journal omits $L$ from the denominator, which agrees with
the corrected form only for $L=1$. Finite propagation speed drives the
result: one iteration expands the exact region by one characteristic
cone, and enough cones cover the full space–time domain.

Figure 3.2 visualizes the cones. A portion of each subdomain solution already agrees with the exact solution, and the next interface exchange enlarges that portion. The same geometry motivates red–black SWR, in which colored space–time blocks are processed concurrently and some redundant work buys additional parallelism.

### Tent pitching, MTP, and UTP

Tent pitching builds inclined space–time elements according to finite propagation speed. Once data on a tent's lower boundary is available, its interior can be computed independently. Mapped tent pitching (MTP) maps each inclined tent to a regular cylinder so a standard solver can be applied. Mapping adds implementation work and can reduce the observed order.

Unmapped tent pitching (UTP) retains the original geometry and can be interpreted as red–black SWR or restricted additive Schwarz on the all-at-once system. A residual determines how far a tent can advance, eliminating the explicit mapping and its associated order loss. Strictly independent tents are unavailable for parabolic equations because their propagation speed is infinite. With weak diffusion, SWR/UTP may still work by iteratively correcting the cross-tent influence.

## 3.3 Time-parallel IDC

### The residual error equation

Consider

$$
u'(t)=f(u(t),t),\qquad u(0)=u_0.
$$

For a current approximation $u^k(t)$, define the integral residual

$$
r^k(t)=u_0+\int_0^t f(u^k(s),s)\,ds-u^k(t). \tag{3.6}
$$

Let $e^k=u-u^k$. Substitution into the integral equation gives an integral error relation; differentiation gives

$$
e^{k\prime}(t)
=f(u^k(t)+e^k(t),t)-f(u^k(t),t)+r^{k\prime}(t). \tag{3.8}
$$

IDC obtains a low-order predictor, discretely solves this error equation, and updates $u^{k+1}=u^k+e^k$. The paper derives the correction for a $\theta$ method, then introduces quadrature weights from interpolatory integration to obtain recurrence (3.11). Integrating the residual avoids direct numerical differentiation of a noisy defect.

![Correction and pipeline structure of IDC, PIDC, and RIDC](assets/diagrams/pint/en/idc-pipeline.svg)

Theorem 3.3 states that a base method of order $p$ on $M$ uniform nodes reaches

$$
\mathcal O\!\left(\Delta t^{\min\{M,(k+1)p\}}\right)
$$

after $k$ corrections. Node count eventually limits the order. Spectral deferred correction (SDC) with $J$ Gauss–Lobatto nodes can reach order $2J-1$. PFASST in Chapter 4 places SDC inside a multilevel parallel structure.

### Why standard IDC remains sequential

A long horizon is usually divided into windows. The predictor and corrections finish on one window before its endpoint initializes the next. Standard IDC therefore processes windows sequentially, while each correction also contains a recurrence across nodes. High-order accuracy alone creates no temporal concurrency.

### 3.3.1 Pipeline IDC (PIDC)

PIDC executes different correction sweeps on different windows at the same time. Once the pipeline is full, the predictor works on a later window, the first correction works on the preceding window, and higher corrections trail behind. The correction count roughly determines concurrency; fill and drain phases reduce efficiency for short runs.

There is an accuracy risk. A later window begins from an endpoint that keeps changing as the previous window receives higher corrections. A high correction level may therefore start from rough, unsettled initial data. IDC order theory requires adequate temporal regularity. Irregular initial data, weak diffusion, and persistent high-frequency structure can remove the expected order gain.

Figure 3.5 uses periodic advection–diffusion with $\Delta x=1/64$, $T=3$, window length $\Delta T=0.1$, $M=5$ nodes per window, and backward Euler as the base method. A narrow source with $\sigma=1000$ has low regularity, while $\sigma=5$ is smooth. The results separate into three regimes: repeated IDC/PIDC corrections do not reliably raise order for low-regularity data; smooth data with strong diffusion benefits clearly; smooth data with weak diffusion still suffers because long-lived high frequencies undermine the ideal correction assumptions.

### 3.3.2 Revisionist IDC (RIDC)

RIDC maintains a sliding $M$-node window for every correction level. The levels advance across successive time steps like an assembly line, without waiting for a whole time window to finish. This design reduces global synchronization and keeps different correction levels active on separate cores.

RIDC improves scheduling while retaining the regularity requirement. Figure 3.6 again shows limited order improvement under low regularity and weak diffusion. For hyperbolic problems, preservation of high frequencies is physically valuable and numerically hazardous for deferred-correction order recovery.

## 3.4 ParaExp

### Exact decomposition for linear systems

ParaExp constructs a direct parallel solution for

$$
\boldsymbol u'(t)=A\boldsymbol u(t)+\boldsymbol g(t).
$$

Partition $[0,T]$ into $N$ intervals $[T_{n-1},T_n]$. First solve the zero-initial-value inhomogeneous problems concurrently:

$$
\boldsymbol v_n'(t)=A\boldsymbol v_n(t)+\boldsymbol g(t),
\qquad \boldsymbol v_n(T_{n-1})=0. \tag{3.13}
$$

Then propagate every interface contribution through a homogeneous tail:

$$
\boldsymbol w_n'(t)=A\boldsymbol w_n(t),
\quad t\in(T_{n-1},T],
\qquad \boldsymbol w_n(T_{n-1})=\boldsymbol v_{n-1}(T_{n-1}), \tag{3.14}
$$

with $\boldsymbol v_0(T_0)=\boldsymbol u_0$. Linearity gives the exact reconstruction

$$
\boldsymbol u(t)=\boldsymbol v_n(t)+
\sum_{j=1}^{n}\boldsymbol w_j(t),
\qquad t\in[T_{n-1},T_n]. \tag{3.15}
$$

![ParaExp separates local forced responses from global homogeneous propagation](assets/diagrams/pint/en/paraexp-decomposition.svg)

The central operation is

$$
\boldsymbol w_n(t)=e^{(t-T_{n-1})A}\boldsymbol v_{n-1}(T_{n-1}).
$$

An exponential action can jump directly to any later time. Its cost is controlled by the matrix, tolerance, and approximation method and need not grow linearly with the number of intermediate time steps. Rational Krylov or polynomial/Chebyshev methods are common for large sparse systems. Scaling-and-squaring with Padé approximation is appropriate for smaller dense matrices. A cited wave-equation study reported parallel efficiency up to roughly 80%; the value is specific to its exponential implementation, partition, and machine.

### Nonlinear extension and its limits

For

$$
\boldsymbol u'=A\boldsymbol u+B(\boldsymbol u)+\boldsymbol g, \tag{3.17}
$$

the linear homogeneous part continues to use exponential propagation, while $B(\boldsymbol u)$ enters iterative inhomogeneous subproblems. Theorem 3.4 has two consequences. After iteration $k$, the first $k$ time intervals agree with the sequential fine solution. At coarse time points, the iteration is identical to a simplified Parareal method whose coarse propagator is the linear homogeneous evolution and whose fine propagator resolves the full nonlinearity.

Figure 3.8 compares ParaExp with standard Parareal for Burgers' equation using spatial step $0.01$, $T=2$, fine step $0.01/20$, and standard-Parareal coarse step $0.01$. ParaExp is faster at large viscosity because linear diffusion dominates. Standard Parareal becomes faster as nonlinear transport gains importance. ParaExp diverges at $\nu=0.02$, and standard Parareal also fails deeper in the hyperbolic regime. The linear ParaExp construction is powerful; the nonlinear splitting quality sets the range of its extension.

## 3.5 ParaDiag: diagonalization in time

### Two ParaDiag routes

ParaDiag converts an all-at-once system into independent spatial systems. The paper distinguishes:

1. **ParaDiag-I**, which exactly diagonalizes a specially designed time discretization and acts as a direct solver;
2. **ParaDiag-II**, which approximates the original time matrix by a circulant or $\alpha$-circulant one and uses the result as a stationary iteration or Krylov preconditioner.

Both need a well-conditioned temporal eigenvector matrix and efficient solvers for complex shifted spatial systems. Both use three stages: transform in time, solve all shifted spatial problems concurrently, and apply the inverse transform.

![ParaDiag time transform, independent spatial solves, and inverse transform](assets/diagrams/pint/en/paradiag-three-stage.svg)

### 3.5.1 Direct ParaDiag methods (ParaDiag-I)

#### Backward Euler on a geometric time mesh

For the linear system (2.1), variable-step backward Euler gives

$$
K=B\otimes I_x-I_t\otimes A. \tag{3.23}
$$

When $\Delta t_n=\mu^{n-1}\Delta t_1$, the matrix $B$ admits $B=VDV^{-1}$. Hence

$$
K^{-1}
=(V\otimes I_x)
(D\otimes I_x-I_t\otimes A)^{-1}
(V^{-1}\otimes I_x). \tag{3.25}
$$

The middle block diagonal contains $N_t$ independent shifted spatial problems.

Writing $\mu=1+\rho$ exposes the main conflict. A larger $\rho$ produces substantial step variation and truncation error of order $\mathcal O(\rho^2)$. A very small $\rho$ makes $B$ approach a nondiagonalizable Jordan structure, amplifying roundoff like $\epsilon\rho^{-(N_t-1)}$. Theorem 3.5 balances these contributions and gives the scale of an optimal $\rho$. Figures 3.9–3.10 confirm the U-shaped error curve and the severe limit on useful time steps in double precision. Increasing $N_t$ eventually lets conditioning and roundoff dominate.

#### Wave equation and the trapezoidal rule

The wave equation is converted to a first-order system and integrated by a variable-step trapezoidal rule to retain its energy behavior. The all-at-once system remains diagonalizable. Theorem 3.6 yields the same competition between truncation and $\epsilon\rho^{-(N_t-1)}$ roundoff. Figure 3.11 and Table 3.1 show rapid growth of the eigenvector condition number; the error starts increasing beyond roughly $N_t=32$.

#### Boundary-value methods improve conditioning

A boundary-value method (BVM) applies centered differences at the first $N_t-1$ nodes and backward Euler to close the final node:

$$
\frac{\boldsymbol u_{n+1}-\boldsymbol u_{n-1}}{2\Delta t}
=A\boldsymbol u_n+\boldsymbol g_n.
$$

The complete diagonalizable system remains second order even though the terminal formula is first order. Theorem 3.7 gives $\operatorname{Cond}(V)=\mathcal O(N_t^2)$, a substantial improvement over the geometric-step construction. Figure 3.12 shows second-order decay on a uniform mesh without the earlier rapid deterioration. A direct second-order wave formulation can use $B^2\otimes I_x-I_t\otimes A$ and avoid the doubled storage of a first-order conversion.

#### Nonlinear ParaDiag-I

Newton linearization of the nonlinear all-at-once equations produces a different Jacobian block at every time point. Replacing these blocks by an average Jacobian recovers a Kronecker approximation and gives a quasi-Newton method whose shifted Jacobian systems remain concurrent after the time transform. Multiple sequential windows are preferable when the Jacobian changes too much over a long horizon.

Figure 3.13 and Table 3.2 show the behavior for Burgers' equation. At $\nu=0.1$, the count of parallel Jacobian solves is far below the sequential Newton count. Lower viscosity makes convergence more sensitive to $T$, and the method fails for $\nu=0.002$ on longer windows.

A richer approximation solves the nearest-Kronecker-product problem

$$
\min_{\Phi_k\ \mathrm{diagonal}}
\left\|
\nabla F(\boldsymbol U^k)-\Phi_k\otimes A_k
\right\|. \tag{3.47}
$$

Here $\Phi_k$ is diagonal in time. The NKA strategy can compute it
offline on a coarse spatial model and retain the amplitude of the
Jacobian's time variation. Figure 3.14 shows a clear improvement,
especially for longer windows such as $T=1.3$.

### 3.5.2 Iterative ParaDiag methods (ParaDiag-II)

#### Strang circulant preconditioner

The all-at-once matrix for a linear multistep method is

$$
K=B_1\otimes I_x-B_2\otimes\Delta t A.
$$

ParaDiag-II replaces the Toeplitz-like matrices $B_1,B_2$ by Strang circulants $C_1,C_2$:

$$
P=C_1\otimes I_x-C_2\otimes\Delta t A.
$$

The discrete Fourier matrix simultaneously diagonalizes the circulants, so applying $P^{-1}$ again consists of an FFT, independent shifted spatial solves, and an inverse FFT. A close approximation supports stationary iteration. When its contraction is weak, $P$ is usually more effective inside GMRES. Krylov convergence can remain rapid with clustered eigenvalues even when the stationary spectral radius exceeds one.

Theorem 3.8 gives a structural finite-step result for symmetric negative $A$: only a limited collection of eigenvalues differs from one, so exact-arithmetic GMRES has a finite termination bound. The bound grows with spatial dimension and alone does not imply rapid convergence for large or nonsymmetric systems.

Figure 3.15 uses $T=2$, $\Delta t=1/50$, and $\Delta x=1/100$. The circulant preconditioner is excellent for the heat equation and advection–diffusion at $\nu=10^{-3}$. Spectral clustering deteriorates as viscosity falls. For the wave equation, the nonunit eigenvalues spread in the complex plane and require many more outer iterations. This behavior directly reflects Chapter 2: diffusion localizes the head–tail mismatch introduced by cyclic closure, while persistent propagation carries that mismatch throughout the horizon.

#### $\alpha$-circulants and a waveform-relaxation interpretation

Continuous head–tail waveform relaxation derives the same
$\alpha$-circulant matrix. The coupling joins the initial value at
iteration $k$ to the terminal values from iterations $k$ and $k-1$ on
the **same time window**; it does not open a new window. After
discretization, the eigenvector matrix has the form
$\Lambda_\alpha F^*$ and can still be applied with FFTs.

$\alpha=1$ recovers the standard circulant approximation. Purely periodic propagation can make this choice singular. Values $0<\alpha<1$ break exact cyclic closure and can greatly improve advection- and wave-dominated cases, as the stationary results in Figure 3.16 demonstrate.

For stable one-step methods and symmetric two-step methods, Theorem 3.9 gives

$$
\frac{1}{1+\alpha}
\le |\lambda(P_\alpha^{-1}K)|
\le \frac{1}{1-\alpha},
\qquad
\rho(I-P_\alpha^{-1}K)
\le \frac{\alpha}{1-\alpha}.
$$

These bounds are independent of the spatial grid and time-step count under the stated stability assumptions. (The printed paper reverses the endpoints as $1/(1-\alpha)\le|\lambda|\le1/(1+\alpha)$, which is an empty interval for $0<\alpha<1$; the close-reading page gives the check and the derivation.) A Numerov threshold experiment confirms their dependence on stability: $\gamma=1/120$ satisfies the condition, whereas the slightly unstable $1/120.01$ invalidates the predicted behavior.

Smaller $\alpha$ improves the iteration factor while increasing the
transform condition number to
$\alpha^{-(N_t-1)/N_t}\le1/\alpha$; thus
$O(\epsilon/\alpha)$ is a convenient conservative roundoff bound.
Figure 3.18 displays this tradeoff. Updating through an error equation
prevents repeated transformation of a large solution vector and reduces
roundoff contamination. A more general multistep Volterra analysis
likewise places the preconditioned eigenvalues at
$1+\mathcal O(\alpha)$.

#### Nonlinear ParaDiag-II

Nonlinear problems use Newton–Krylov. Each Newton step produces an all-at-once Jacobian, and GMRES is preconditioned by $P_\alpha$ built from an average Jacobian. Stationary iteration often fails on long nonlinear windows, while Krylov methods can still exploit spectral clustering. Shorter windows make Jacobian blocks more similar; NKA retains additional temporal variation. Both improve preconditioner quality.

## Site numerical supplement: Python validation of Figure 3.15

### Baseline experiment

The recomputed baseline uses $N_x=N_t=100$, $T=2$, $\nu=10^{-6}$, $\alpha=1$, and GMRES tolerance $10^{-12}$. It converges in 13 iterations with true relative residual $1.152\times10^{-14}$.

![GMRES convergence baseline for ParaDiag-II on advection–diffusion](assets/pint/paradiag-baseline.svg)

### Paper-grid validation

| Problem                            |                                         Python result | Interpretation                                                  |
| ---------------------------------- | ----------------------------------------------------: | --------------------------------------------------------------- |
| heat                               |                                      2 Krylov updates | eigenvalues tightly clustered near one                          |
| advection–diffusion, $\nu=10^{-3}$ |                                                     3 | circulant and original all-at-once systems are close            |
| advection–diffusion, $\nu=10^{-6}$ |                                                    13 | weak diffusion transports the cyclic-closure mismatch           |
| wave                               | preconditioned residual below $10^{-11}$ at update 89 | nonunit eigenvalues spread along $\operatorname{Re}\lambda=0.5$ |

![ParaDiag-II spectra and GMRES convergence for heat, advection–diffusion, and wave equations](assets/pint/paradiag-figure-3-15.svg)

The 3- and 13-update ADE results match Figures 3.15(c,d). The wave run uses the paper's $\gamma=1/100$ and $\alpha=1$. Single-threaded SciPy crosses a preconditioned residual of $10^{-11}$ at update 89, close to the paper curve ending around 88. SciPy requires 103 updates under its true-relative-residual threshold of $10^{-12}$. MATLAB and SciPy differ in residual normalization, restart, and stopping conventions, so the spectral geometry and convergence phase are more meaningful than a single raw count.

## Site method comparison: parameters, applicability, and implementation cost

| Method      | Key parameter or choice                                | Property controlled                         | Main risk                                                               |
| ----------- | ------------------------------------------------------ | ------------------------------------------- | ----------------------------------------------------------------------- |
| SWR         | overlap, window length, transmission operator          | rate of characteristic information transfer | poor tuning or expensive interface operator                             |
| PIDC/RIDC   | nodes, correction levels, window and pipeline depth    | formal order and concurrency                | loss of order under low regularity                                      |
| ParaExp     | exponential-action method, linear/nonlinear split      | cost of homogeneous tails                   | expensive action or inaccurate nonlinear split                          |
| ParaDiag-I  | time discretization, $N_t$, eigenvector conditioning   | direct concurrency                          | conflict between truncation and roundoff                                |
| ParaDiag-II | $\alpha$, outer Krylov method, shifted-solve tolerance | clustering and stability                    | small $\alpha$ amplifies roundoff; large $\alpha$ weakens approximation |

The latest public [ActaPinT-Python](https://github.com/freezeng123456/ActaPinT-Python) commit now gives independent entries, SVG/PNG figures, and JSON metrics for all fifteen numerical figures, three computed schematics, and Tables 3.1–3.2 in Section 3. SWR, PIDC/RIDC, ParaExp, direct ParaDiag, BVM/NKA, ParaDiag II, and wave-domain decomposition are no longer migration-list placeholders. One exception remains: the average-Jacobian curves of Figure 3.14 match, but the two published nearest-Kronecker speed-up curves at $\nu=0.02$ cannot be recovered from the documented weights. See the [[en/computational-mathematics/knowledge-notes/time-parallelization/reproduction-audit-2026-08-24|whole-article reproduction audit]] for itemized evidence and current test health.

## Source coverage audit

| Source location                 | This page | Material covered                                                                                                                          |
| ------------------------------- | --------- | ----------------------------------------------------------------------------------------------------------------------------------------- |
| Sections 3 and 3.1, pp. 396–398 | 3.1       | scope and origins of hyperbolic-effective methods, tent-pitching context                                                                  |
| Section 3.2, pp. 398–405        | 3.2       | WR versus SWR, Dirichlet/Robin/Ventcel/convolution transmission, Theorems 3.1–3.2, Figures 3.1–3.3, MTP/UTP                               |
| Section 3.3, pp. 405–411        | 3.3       | IDC residual and recurrence, Theorem 3.3, SDC, PIDC/RIDC schedules, Figures 3.4–3.6, regularity limitation                                |
| Section 3.4, pp. 412–415        | 3.4       | linear ParaExp decomposition and exponential action, nonlinear iteration, Theorem 3.4, Figures 3.7–3.8                                    |
| Sections 3.5–3.5.1, pp. 415–430 | 3.5.1     | ParaDiag-I stages, geometric steps, Theorems 3.5–3.7, wave/BVM variants, nonlinear quasi-Newton and NKA, Figures 3.9–3.14, Tables 3.1–3.2 |
| Section 3.5.2, pp. 431–442      | 3.5.2     | Strang and $\alpha$-circulant preconditioners, Theorems 3.8–3.9, Figures 3.15–3.18, stability/roundoff tradeoff, nonlinear Newton–Krylov  |

## Summary

The four families preserve long-range information in different ways. SWR transports interface waveforms along characteristic geometry. IDC pipelines high-order residual corrections. ParaExp applies the matrix exponential to linear homogeneous responses. ParaDiag transforms the global temporal coupling into concurrent spatial solves. Each avoids exclusive reliance on a dissipative coarse propagator and introduces its own constraint: SWR needs effective transmission, IDC needs temporal regularity, ParaExp needs an efficient exponential action and a meaningful split, and ParaDiag needs stable diagonalization plus scalable shifted solvers.

## Source

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Section 3, pp. 396–443.
