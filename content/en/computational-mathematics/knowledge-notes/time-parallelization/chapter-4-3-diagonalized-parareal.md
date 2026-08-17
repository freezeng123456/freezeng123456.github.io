---
title: "4.5: Diagonalization-Based Parareal"
description: The two complete routes of parallel coarse-grid correction and an interval-local diagonalized coarse propagator
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization/chapter-4-3-diagonalized-parareal
tags:
  - parallel-in-time
  - Parareal
  - ParaDiag
---

> [!note] Reading scope
> This page follows Section 4.5 (pp. 460–471). It covers equations (4.14)–(4.29), Theorems 4.7–4.8, Remark 4.2, and Figures 4.12–4.17. Both variants use diagonalization, but at different locations: the first parallelizes the global coarse-grid correction; the second defines a special coarse propagator inside each coarse interval.

## 4.5 Diagonalization-based Parareal

### Distinguishing the two constructions

- **Diagonalized CGC (Section 4.5.1; Wu 2018, Wu and Zhou 2019):** modifies the serial correction across the $N_t$ coarse points. Its concurrency lies across coarse points, and its convergence mechanism remains close to standard Parareal, so its main range is parabolic.
- **Diagonalized coarse propagator (Section 4.5.2; Gander and Wu 2020):** retains the outer Parareal form and uses ParaDiag on the $J$ fine points inside every $[T_n,T_{n+1}]$. Coarse and fine propagation use the same integrator and step size. The paper notes that this coarse propagator transports **all frequency components** over a very long time, so it can also handle hyperbolic problems.

## 4.5.1 Diagonalization-based CGC

### From serial CGC to head–tail coupling

Standard coarse-grid correction is

$$
\boldsymbol u_{n+1}^{k+1}
=\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})
+\boldsymbol b_{n+1}^k,
\qquad n=0,\ldots,N_t-1, \tag{4.14}
$$

$$
\boldsymbol b_{n+1}^k
=\mathcal F(T_n,T_{n+1},\boldsymbol u_n^k)
-\mathcal G(T_n,T_{n+1},\boldsymbol u_n^k),
$$

starting from $\boldsymbol u_0^{k+1}=\boldsymbol u_0$. Wu (2018) replaced this initial condition by

$$
\boldsymbol u_0^{k+1}=\alpha\boldsymbol u_{N_t}^{k+1}+\boldsymbol u_0.
$$

To preserve the original fixed point, define

$$
\widetilde{\boldsymbol u}_n^k=
\begin{cases}
\boldsymbol u_0,&n=0,\\
\boldsymbol u_n^k,&n\ge1,
\end{cases}
$$

and iterate

$$
\left\{
\begin{aligned}
\boldsymbol u_{n+1}^{k+1}
={}&\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})
+\mathcal F(T_n,T_{n+1},\widetilde{\boldsymbol u}_n^k)
-\mathcal G(T_n,T_{n+1},\boldsymbol u_n^k),\\
\boldsymbol u_0^{k+1}&=\alpha\boldsymbol u_{N_t}^{k+1}+\boldsymbol u_0.
\end{aligned}
\right. \tag{4.15}
$$

This predates the more natural ParaDiag-II condition (3.55); the modified value $\widetilde{\boldsymbol u}_0^k$ maintains consistency at convergence.

### Linear all-at-once system and three-stage solve

For $\boldsymbol u'=A\boldsymbol u$ and backward Euler coarse propagation, substitution of the head–tail condition yields

$$
(C_\alpha\otimes I_x-I_t\otimes\Delta TA)\boldsymbol U^{k+1}
=\boldsymbol g^k, \tag{4.16}
$$

where

$$
C_\alpha=
\begin{bmatrix}
1&&&-\alpha\\-1&1\\&\ddots&\ddots\\&&-1&1
\end{bmatrix},
$$

$$
\boldsymbol g^k=
\begin{bmatrix}
\boldsymbol u_0+(I_x-\Delta TA)\boldsymbol b_1^k\\
(I_x-\Delta TA)\boldsymbol b_2^k\\\vdots\\
(I_x-\Delta TA)\boldsymbol b_{N_t}^k
\end{bmatrix}.
$$

The paper writes the diagonalized solve as

$$
\left\{
\begin{aligned}
\boldsymbol U^{a,k+1}&=(F\otimes I_x)\boldsymbol g^k,\\
(\lambda_nI_x-\Delta TA)\boldsymbol u_n^{b,k+1}
&=\boldsymbol u_n^{a,k+1},
&&n=1,\ldots,N_t,\\
\boldsymbol U^{k+1}&=(F^*\otimes I_x)\boldsymbol U^{b,k+1}.
\end{aligned}
\right. \tag{4.17}
$$

Here $\{\lambda_n\}$ are the eigenvalues of $C_\alpha$ and $F$ is the discrete Fourier matrix, as in (3.51) and (3.50).

> [!note] Site supplement: the diagonal scaling for $\alpha\ne1$
> As printed, (4.17) uses only $F$ and $F^*$, which diagonalize $C_\alpha$ only for $\alpha=1$. By (3.59) of the paper's Section 3.5.2, the general case is $C_\alpha=V_\alpha D_\alpha V_\alpha^{-1}$ with $V_\alpha=\Lambda_\alpha F^*$ and $\Lambda_\alpha=\operatorname{diag}(1,\alpha^{-1/N_t},\ldots,\alpha^{-(N_t-1)/N_t})$, so an implementation needs the corresponding diagonal scaling.

As in Section 3.5.2, the essential stages are an FFT-like transform, independent shifted spatial solves, and the inverse transform.

### Theorem 4.7: threshold for matching standard Parareal

As $\alpha\to0$, equation (4.15) approaches standard CGC, while alpha-circulant roundoff increases, particularly in single or half precision. Theorem 4.7 is due to Wu (2018): let $\rho$ be the convergence factor of standard Parareal (4.14) and $\rho_{\mathrm{new}}$ that of the new variant (4.15), with $\mathcal G$ a stable integrator. Then

$$
\rho_{\mathrm{new}}=\rho,
\qquad\text{if}\qquad
\alpha\le\frac{\rho}{1+\rho}.
$$

The paper notes that this was proved for linear problems $\boldsymbol u'=A\boldsymbol u+\boldsymbol g$ where $A$ has negative real eigenvalues; for other cases, such as complex eigenvalues, numerical results suggest it holds as well but there is no proof.

The practical choice is the threshold itself: reducing alpha further does not improve the asymptotic rate and increases roundoff exposure. Typically both $\rho$ and alpha are $O(10^{-1})$, a regime in which the roundoff incurred by the diagonalization is negligible.

![Original Figure 4.12: standard and diagonalized CGC for heat and ADE](assets/papers/time-parallelization/source-figures/figure-4-12.svg)

The test uses periodic data, $u_0(x)=\sin(2\pi x)$, backward Euler coarse and SDIRK22 fine propagation, $T=4$, $J=10$, $\Delta T=0.1$, and $\Delta x=1/128$. Panel (a) is heat, with $\rho\approx0.22$ and threshold $0.18$: $\alpha=0.25,0.4$ are slower than standard CGC, while $\alpha=0.1$ tracks it. Panel (b) is ADE at $\nu=0.1$, with $\rho\approx0.39$ and threshold $0.28$: $\alpha=0.1,0.25$ track standard CGC and $\alpha=0.4$ is clearly slower. Note that the ADE semi-discretization has complex eigenvalues, exactly the case Theorem 4.7 does not cover and for which only numerical evidence is available.

### Nonlinear all-at-once quasi-Newton solve

With backward Euler coarse propagation, define

$$
\boldsymbol b_{n+1}^k
=\mathcal F(T_n,T_{n+1},\widetilde{\boldsymbol u}_n^k)
-\mathcal G(T_n,T_{n+1},\boldsymbol u_n^k).
$$

The correction is

$$
(C_\alpha\otimes I_x)\boldsymbol U^{k+1}
-\Delta TF(\boldsymbol U^{k+1})=\boldsymbol g^k. \tag{4.18}
$$

The $n$th block of $F$ is $f(\boldsymbol u_n^{k+1}-\boldsymbol b_n^k)$. The inner quasi-Newton iteration is

$$
P_\alpha^{k+1,l}\Delta\boldsymbol U^{k+1,l}
=\boldsymbol g^k-(C_\alpha\otimes I_x)\boldsymbol U^{k+1,l}
+\Delta TF(\boldsymbol U^{k+1,l}),
$$

$$
\boldsymbol U^{k+1,l+1}
=\boldsymbol U^{k+1,l}+\Delta\boldsymbol U^{k+1,l}, \tag{4.19a}
$$

$$
P_\alpha^{k+1,l}
=C_\alpha\otimes I_x-I_t\otimes\Delta TA^{k+1,l}, \tag{4.19b}
$$

The exact block Jacobian and its separable approximation are

$$
\nabla F(\boldsymbol U^{k+1,l})
=\operatorname{blkdiag}_{n=1}^{N_t}
\nabla f(\boldsymbol u_n^{k+1,l}-\boldsymbol b_n^k),
$$

$$
A^{k+1,l}
=\frac1{N_t}\sum_{n=1}^{N_t}
\nabla f(\boldsymbol u_n^{k+1,l}-\boldsymbol b_n^k),
\qquad
\nabla F\approx I_t\otimes A^{k+1,l}.
$$

> [!warning] Source check: the average-Jacobian index
> The journal prints
> $J^{-1}\sum_{j=1}^J\nabla f(\boldsymbol u_n-\boldsymbol b_n)$:
> the summand contains no $j$, and the Section 4.5.1 all-at-once system
> has $N_t$ blocks rather than $J$. The display above follows the block
> Jacobian definition.

The preconditioner has the structure of (4.16), so every increment uses
(4.17). A nearest Kronecker approximation could replace the average
Jacobian but is not pursued here.

![Original Figure 4.13: the two CGCs for Burgers' equation at two viscosities](assets/papers/time-parallelization/source-figures/figure-4-13.svg)

The left and right panels of Figure 4.13 use Burgers' equation at $\nu=1$ and $0.01$, respectively, with the same problem set-up and discretization parameters as the heat and ADE runs. Each compares $\alpha=0.4,0.25,0.1$ with standard CGC. Alpha $0.4$ is slowest in both viscosity regimes; at $\nu=1$ the $\alpha=0.1$ curve is closest to the standard one, while at $\nu=0.01$ it is $\alpha=0.25$ that essentially coincides with the standard curve and $\alpha=0.1$ that runs slightly below it. The conclusion the paper draws is that the influence of $\alpha$ on the convergence rate remains as in the linear case, and Wu (2018, Section 4) shows the rate mirrors Parareal with standard CGC when $\alpha$ is chosen appropriately small. The paper states no nonlinear analogue of the $\rho/(1+\rho)$ threshold.

### Remark 4.2: MGRiT needs a consistent head–tail condition

Directly transplanting (4.15) into MGRiT diverges for every alpha. The consistent condition is

$$
\boldsymbol u_1^{k+1}
=\alpha(\boldsymbol u_{N_t}^{k+1}-\boldsymbol u_{N_t}^k)+\boldsymbol u_1,
$$

which leads to

$$
\left\{
\begin{aligned}
\boldsymbol u_0^{k+1}&=\boldsymbol u_0,\\
\boldsymbol u_1^{k+1}
&=\alpha(\boldsymbol u_{N_t}^{k+1}-\boldsymbol u_{N_t}^k)+\boldsymbol u_1,\\
\boldsymbol u_{n+1}^{k+1}
&=\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})
+\widetilde{\boldsymbol b}_{n+1}^k.
\end{aligned}
\right. \tag{4.20}
$$

Here $\widetilde{\boldsymbol b}_{n+1}^k=\mathcal F(T_n,T_{n+1},\widetilde{\boldsymbol s}_n^k)-\mathcal G(T_n,T_{n+1},\widetilde{\boldsymbol s}_n^k)$ and $\widetilde{\boldsymbol s}_n^k=\mathcal F(T_{n-1},T_n,\widetilde{\boldsymbol u}_{n-1}^k)$. Note that $\widetilde{\boldsymbol u}_n^k$ differs from the Section 4.5.1 definition: for the MGRiT variant the paper redefines it as

$$
\widetilde{\boldsymbol u}_n^k=
\begin{cases}
\boldsymbol u_n,&n=0,1,\\
\boldsymbol u_n^k,&n\ge2,
\end{cases}
$$

so at $n=1$ it uses the converged value $\boldsymbol u_1$ rather than the current iterate $\boldsymbol u_1^k$, consistent with the $\boldsymbol u_1$ appearing in the head–tail condition above. For small alpha this variant matches original MGRiT, with the same threshold mechanism as Theorem 4.7.

## 4.5.2 Diagonalization-based coarse solver

### Fine propagation and a head–tail coarse propagator

Inside each coarse interval, both propagators use the same linear-theta method and $\Delta t=\Delta T/J$. Fine propagation advances sequentially:

$$
\boldsymbol v_{j+1}-\boldsymbol v_j
=\Delta t[\theta f(\boldsymbol v_{j+1})
+(1-\theta)f(\boldsymbol v_j)],
\quad j=0,\ldots,J-1,
\quad \boldsymbol v_0=\boldsymbol u_n. \tag{4.21}
$$

$\theta=1$ is backward Euler and $\theta=1/2$ is trapezoidal. The
linear-theta method exposes the structure; an $s$-stage Runge–Kutta
generalization is given in the appendix of Gander and Wu (2020). The
special coarse propagator $\mathcal F_\alpha^*$ changes only the
condition to

$$
\boldsymbol v_0=\alpha\boldsymbol v_J+(1-\alpha)\boldsymbol u_n. \tag{4.22}
$$

The $J$ fine steps are now a head–tail system solvable in parallel by ParaDiag.

### Nonlinear all-at-once system and quasi-Newton iteration

With $\boldsymbol V=(\boldsymbol v_1^\top,\ldots,\boldsymbol v_J^\top)^\top$,

$$
\underbrace{(C_\alpha\otimes I_x)\boldsymbol V
-\Delta tF(\boldsymbol V)}_{K(\boldsymbol V)}
=\boldsymbol b(\boldsymbol u_n), \tag{4.23}
$$

$$
\boldsymbol b(\boldsymbol u_n)
=((1-\alpha)\boldsymbol u_n^\top,0,\ldots,0)^\top.
$$

The paper defines $\boldsymbol V$ and $\boldsymbol b(\boldsymbol u_n)$ in an unnumbered inline display; the tag (4.24) belongs to the following group, which defines $C_\alpha$ and $F(\boldsymbol V)$:

$$
C_\alpha=
\begin{bmatrix}
1&&&-\alpha\\
-1&1\\
&\ddots&\ddots\\
&&-1&1
\end{bmatrix},
\qquad
F(\boldsymbol V)=
\begin{bmatrix}
\theta f(\boldsymbol v_1)+(1-\theta)f(\alpha\boldsymbol v_J+(1-\alpha)\boldsymbol u_n)\\
\theta f(\boldsymbol v_2)+(1-\theta)f(\boldsymbol v_1)\\
\vdots\\
\theta f(\boldsymbol v_J)+(1-\theta)f(\boldsymbol v_{J-1})
\end{bmatrix}. \tag{4.24}
$$

The $-\alpha$ entry in the corner of $C_\alpha$ is exactly the head–tail condition (4.22): replacing it by $0$ makes $C_\alpha$ strictly lower triangular and turns the all-at-once system back into sequential fine propagation. The first block of $F$ differs in shape from the others for the same reason, since $\boldsymbol v_0$ has been replaced by $\alpha\boldsymbol v_J+(1-\alpha)\boldsymbol u_n$. Note that the weights $\theta$ and $1-\theta$ are already built into $F$, which matters for the Jacobian derived below.

The quasi-Newton update is

$$
P_\alpha(\boldsymbol V^l)\Delta\boldsymbol V^l
=\boldsymbol b(\boldsymbol u_n)-K(\boldsymbol V^l),
\qquad
\boldsymbol V^{l+1}=\boldsymbol V^l+\Delta\boldsymbol V^l, \tag{4.25a}
$$

$$
P_\alpha(\boldsymbol V^l)
=C_\alpha\otimes I_x
-\Delta t\widetilde C_{\alpha,\theta}\otimes\overline{\nabla f}(\boldsymbol V^l), \tag{4.25b}
$$

where

$$
\widetilde C_{\alpha,\theta}=
\begin{bmatrix}
\theta&&&(1-\theta)\alpha\\
1-\theta&\theta\\&\ddots&\ddots\\&&1-\theta&\theta
\end{bmatrix},
$$

Let
$J_0=\nabla f(\alpha\boldsymbol v_J^l+(1-\alpha)\boldsymbol u_n)$ and
$J_j=\nabla f(\boldsymbol v_j^l)$. Differentiating the already weighted
$F(\boldsymbol V)$ in (4.24) gives

$$
\mathcal J_F(\boldsymbol V^l)=
\begin{bmatrix}
\theta J_1&&&\alpha(1-\theta)J_0\\
(1-\theta)J_1&\theta J_2\\
&\ddots&\ddots\\
&&(1-\theta)J_{J-1}&\theta J_J
\end{bmatrix}.
$$

Therefore

$$
\nabla K=C_\alpha\otimes I_x-\Delta t\,\mathcal J_F.
$$

> [!warning] Source check: the overloaded $\nabla F$
> The source writes
> $(\widetilde C_{\alpha,\theta}\otimes I_x)\nabla F$, although (4.24)
> has already placed the theta and head–tail weights inside $F$.
> Multiplying by $\widetilde C_{\alpha,\theta}$ again would double
> those weights if $\nabla F$ retains that definition. The display
> above differentiates (4.24) directly and removes the ambiguity.

and the separable approximation uses the special average

$$
\overline{\nabla f}(\boldsymbol V^l)
=\frac1J\left[
\sum_{j=1}^{J-1}\nabla f(\boldsymbol v_j^l)
+\nabla f\!\left(
\alpha\boldsymbol v_J^l+(1-\alpha)\boldsymbol u_n
\right)
\right].
$$

The final block is evaluated at the head–tail state rather than at
$\boldsymbol v_J^l$. The two temporal matrices are simultaneously
diagonalizable. The outer iteration is

$$
\boldsymbol u_{n+1}^{k+1}
=\mathcal F_\alpha^*(T_n,T_{n+1},\boldsymbol u_n^{k+1})
+\mathcal F(T_n,T_{n+1},\boldsymbol u_n^k)
-\mathcal F_\alpha^*(T_n,T_{n+1},\boldsymbol u_n^k). \tag{4.26}
$$

### Linear system, concurrency, and limiting cases

For $f(\boldsymbol u)=A\boldsymbol u$,

$$
(C_\alpha\otimes I_x
-\widetilde C_{\alpha,\theta}\otimes\Delta tA)\boldsymbol V
=\boldsymbol b(\boldsymbol u_n), \tag{4.27}
$$

$$
\boldsymbol b(\boldsymbol u_n)
=([(I_x+\Delta t(1-\theta)A)(1-\alpha)\boldsymbol u_n]^\top,0,\ldots,0)^\top.
$$

We use $\widetilde C_{\alpha,\theta}$ consistently; the source
alternates it with $\widetilde C_{\theta,\alpha}$ for the same matrix.
Let $H_J=(0,\ldots,0,1)\in\mathbb R^{1\times J}$. The coarse output is
$\mathcal F_\alpha^*=(H_J\otimes I_x)\boldsymbol V=\boldsymbol v_J$,
equivalently

$$
\left\{
\begin{aligned}
\boldsymbol v_{j+1}-\boldsymbol v_j
&=\Delta tA[\theta\boldsymbol v_{j+1}+(1-\theta)\boldsymbol v_j],\\
\boldsymbol v_0&=\alpha\boldsymbol v_J+(1-\alpha)\boldsymbol u_n^k.
\end{aligned}
\right. \tag{4.28}
$$

At $\alpha=0$, coarse propagation equals sequential fine propagation and outer Parareal converges in one iteration with no speedup. For $0<\alpha<1$, all $J$ points are solved at once; with enough spatial-solve resources, wall time is approximately $1/J$ of sequential fine propagation.

### Theorem 4.8: parabolic and hyperbolic spectra

Theorem 4.8 is due to Gander and Wu (2020). For
$\boldsymbol u'=A\boldsymbol u+\boldsymbol g$,
$\boldsymbol u(0)=\boldsymbol u_0$, with
$A\in\mathbb C^{N_x\times N_x}$, suppose both $\mathcal F$ and
$\mathcal F_\alpha^*$ use a stable one-step Runge–Kutta method and
$\{\boldsymbol u_n^k\}$ is iteration $k$ of (4.26). Let

$$
e^k=\max_{1\le n\le N_t}\|\boldsymbol u_n-\boldsymbol u_n^k\|_\infty,
$$

where $\{\boldsymbol u_n\}$ is the converged solution rather than the exact PDE solution.

Then

$$
e^k\le\rho^ke^0,
\qquad
\rho=
\begin{cases}
\alpha,&\sigma(A)\subset\mathbb R_-,\\[4pt]
\dfrac{2\alpha N_t}{1+\alpha},&\sigma(A)\subset i\mathbb R.
\end{cases} \tag{4.29}
$$

The heat-equation factor is independent of the number of coarse intervals. The imaginary-spectrum bound grows linearly with $N_t$ and may be loose when $\alpha N_t$ is large.

![Original Figure 4.14: sharp rho=alpha prediction for the heat equation](assets/papers/time-parallelization/source-figures/figure-4-14.svg)

The heat test uses homogeneous Dirichlet data, $u_0=\sin^2(2\pi x)$, trapezoidal integration, $\Delta T=1/12$, $J=10$, and $\Delta x=1/100$. The left and right panels use $N_t=36$ and $72$; each compares $\alpha=10^{-1},10^{-2},10^{-3}$. The measured dashed curves are almost parallel to the theoretical dotted curves, and doubling $N_t$ does not change the slope set by $\rho=\alpha$.

![Original Figure 4.15: joint influence of alpha and the coarse-interval count on the wave equation](assets/papers/time-parallelization/source-figures/figure-4-15.svg)

The wave equation is first reduced to the first-order system $\boldsymbol w'=\boldsymbol{Aw}$ with $\boldsymbol w=(\boldsymbol u^\top,(\boldsymbol u')^\top)^\top$ and $\boldsymbol A=\begin{bmatrix}0&I_x\\A&0\end{bmatrix}$, which gives $\sigma(\boldsymbol A)\subset i\mathbb R$, the second branch of Theorem 4.8. The test uses periodic data and $\boldsymbol w(0)=(\sin^2(2\pi\boldsymbol x_h)^\top,\boldsymbol 0^\top)^\top$. Panel (a) fixes $\alpha=0.01$ and compares $N_t=24,48,96$, so increasing the interval count visibly slows convergence. Panel (b) fixes $\alpha=10^{-4}$ and compares $N_t=24,48,96,960$; increasing the count from 24 to 960 costs only about two extra iterations to reach $\max\{\Delta t^2,\Delta x^2\}$. The panels isolate the joint effect of $\alpha N_t$.

![Original Figure 4.16: the linear bound is sharp for small alpha Nt and conservative when the product is large](assets/papers/time-parallelization/source-figures/figure-4-16.svg)

Figure 4.16(a) fixes $N_t=24$ and compares $\alpha=10^{-2}$ with $10^{-4}$; panel (b) fixes $\alpha=10^{-4}$ and compares $N_t=24$ with $960$. Only the small-product case $\alpha=10^{-4},N_t=24$ closely follows the dotted bound from (4.29). The other two measured curves exhibit superlinear decay, making the linear bound visibly conservative.

![Original Figure 4.17: iterations needed to reach 1e-8 on Burgers' equation](assets/papers/time-parallelization/source-figures/figure-4-17.svg)

The Burgers test uses periodic data, $u_0=\sin^2(2\pi x)$, $\Delta T=0.1$, $J=10$, and $\Delta x=1/100$; the three curves use $\nu=1,0.01,10^{-4}$. Panel (a) fixes $N_t=40$, showing that small alpha accelerates convergence and reduces viscosity sensitivity. Panel (b) fixes $\alpha=10^{-3}$; from $N_t=10$ to $160$, the count stays between two and five iterations. Nonlinear theory gives $\rho=O(\alpha)$ under exact solution of (4.23) and suitable Lipschitz assumptions.

## Final comparison

| Property                  | Diagonalized CGC (4.15)    | Diagonalized coarse propagator (4.26)                                |
| ------------------------- | -------------------------- | -------------------------------------------------------------------- |
| diagonalization direction | global $N_t$ coarse points | $J$ fine points within each coarse interval                          |
| modified component        | CGC                        | coarse propagator                                                    |
| integrators               | coarse and fine may differ | identical integrator and step size                                   |
| main range                | parabolic problems         | parabolic and hyperbolic problems                                    |
| alpha rule                | $\alpha\le\rho/(1+\rho)$   | parabolic factor $\alpha$; hyperbolic bound $2\alpha N_t/(1+\alpha)$ |

## Equation, theorem, and figure audit

| Source item                            | Paper section | Coverage                                                                 |
| -------------------------------------- | ------------- | ------------------------------------------------------------------------ |
| (4.14)–(4.17)                          | 4.5.1         | standard/head–tail CGC, linear all-at-once matrix, three stages          |
| Theorem 4.7, Figure 4.12               | 4.5.1         | alpha threshold, roundoff tradeoff, heat and ADE                         |
| (4.18)–(4.19), Figure 4.13             | 4.5.1         | nonlinear system, average-Jacobian quasi-Newton, Burgers                 |
| Remark 4.2, (4.20)                     | 4.5.1         | consistent MGRiT head–tail condition and convergent variant              |
| (4.21)–(4.26)                          | 4.5.2         | equal-integrator fine/coarse propagation, nonlinear system, outer update |
| (4.27)–(4.28)                          | 4.5.2         | linear form, terminal extraction, alpha-zero limit, J-way concurrency    |
| Theorem 4.8, (4.29), Figures 4.14–4.17 | 4.5.2         | negative-real/imaginary bounds and all heat, wave, Burgers figures       |

## Source

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Section 4.5, pp. 460–471.
