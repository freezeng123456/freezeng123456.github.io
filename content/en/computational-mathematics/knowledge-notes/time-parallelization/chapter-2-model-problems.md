---
title: "Chapter 2: Model Problems Linking the Parabolic and Hyperbolic World"
description: Starting from "can we compute the late solution first?" to read the heat, advection–diffusion, Burgers, and wave equations
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization/chapter-2-model-problems
tags:
  - parallel-in-time
  - PDE
---

> [!note] Scope of this page
> This page follows Section 2 of the paper (printed pp. 388–396) and covers equations (2.1)–(2.7) together with every panel of Figures 2.1–2.4. The source figures are extracted directly from the paper PDF with axes and panel labels unchanged. The paper is distributed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); full attribution appears at the end. Any interpretation, diagram, or experiment produced for this site is labeled explicitly.

## What this chapter is for

Section 2 presents no algorithms. It does one thing: it builds a measuring stick. The reason Chapters 3 and 4 split the methods into two groups comes entirely from here.

Causality says the late solution is determined by the early one, which makes the time direction look inherently sequential. Section 2 immediately counters with a question: **can we compute $t\in(1.7,2.2)$ before knowing the solution for $t<1.7$?** For the heat equation with Dirichlet boundaries the answer is yes, and the paper punctuates it with an exclamation mark.

That is the pivot of the whole chapter. Causality is a logical dependence, but the **numerically effective dependence** can be far shorter. If dissipation has erased old information, its actual influence on the late solution falls below the discretization error, and "we must know the early solution first" stops being true in any computational sense.

So the chapter measures one quantity throughout: **how long the solution's memory is, and which frequencies that memory keeps**. Those two answers decide whether PinT can succeed.

## Understand the experimental design before reading the figures

This is the easiest thing to miss on a first reading and the most rewarding to notice. The four figures are not four unrelated computations but one controlled experiment: the paper reuses exactly the same source term and exactly the same oscillatory initial condition across all four equations, so panels can be compared directly.

Only three knobs actually change: the equation itself (heat, then advection–diffusion, then Burgers, then wave, adding transport and nonlinearity in turn), the boundary condition (Dirichlet, Neumann, or periodic, deciding whether information can leave the domain), and the diffusion strength $\nu$ (taking $1$, $10^{-2}$, and $5\times10^{-4}$, deciding how fast fine scales decay).

There is also one trap: **the last column of each figure is a different experiment**. Figures 2.1(d), 2.2(d) and (h), 2.3(d) and (h), and 2.4(d) all switch the source off and use the oscillatory initial condition instead. The earlier columns ask how long the response to an external forcing survives; the last column asks how long the initial information itself survives. Do not compare across the two groups.

## Two notations, and one important exception

Most time-parallel methods are described and analyzed on the ODE system obtained after spatial semidiscretization. The linear form is

$$
\begin{aligned}
\boldsymbol u'(t)&=A\boldsymbol u(t)+\boldsymbol g(t),
&&t\in(0,T],\\
\boldsymbol u(0)&=\boldsymbol u_0,
\end{aligned}
\tag{2.1}
$$

where $A\in\mathbb R^{N_x\times N_x}$ comes from semidiscretizing the PDE in space. The nonlinear form is

$$
\begin{aligned}
\boldsymbol u'(t)&=\boldsymbol f(\boldsymbol u(t),t),
&&t\in(0,T],\\
\boldsymbol u(0)&=\boldsymbol u_0,
\end{aligned}
\tag{2.2}
$$

with $\boldsymbol f:\mathbb R^{N_x}\times\mathbb R\to\mathbb R^{N_x}$ nonlinear in its first argument. The Burgers-type semidiscretization can be written $\boldsymbol f(\boldsymbol u(t),t)=A\boldsymbol u(t)+B\boldsymbol u^2(t)+\boldsymbol g(t)$.

**The exception is worth remembering**: domain-decomposition methods do not use this ODE notation and instead solve and analyze continuous space–time subproblems. The difference becomes visible immediately when reading SWR in Chapter 3.

Every model is posed on the one-dimensional unit interval $\Omega=(0,1)$. The paper states explicitly that this is not a real restriction, since the applicability and convergence properties of PinT methods generally do not depend on the spatial dimension. Dimension changes the cost of a single spatial solve, but the three mechanisms studied here are unchanged: decay of past information, survival of fine scales, and the escape or return of information through a boundary.

## 2.1 Heat equation: the boundary condition decides how fast the problem forgets

The parabolic reference model is

$$
\partial_tu(x,t)=\partial_{xx}u(x,t)+g(x,t),
\qquad (x,t)\in\Omega\times(0,T], \tag{2.3}
$$

with initial data $u(x,0)=u_0(x)$. The text first introduces homogeneous Dirichlet or Neumann boundaries, and the figure adds periodic boundaries to test whether heat can leave the domain.

The first three panels use zero initial data and a source that is localized in both space and time:

$$
g(x,t)=10\sum_{j=1}^{4}
\exp\!\left(
-\sigma\left[(t-t_j)^2+(x-0.5)^2\right]
\right), \tag{2.4}
$$

with $(t_1,t_2,t_3,t_4)=(0.1,0.6,1.35,1.85)$ and $\sigma=200$. The source sits at $x=0.5$ and its four time centres are well separated. Taking $\sigma=200$ gives characteristic widths of about $1/\sqrt{200}\approx0.07$ in both space and time, so it appears as four bright bands in the space–time plot.

![Source Figure 2.1: heat-equation solutions for three boundary conditions and oscillatory initial data](assets/papers/time-parallelization/source-figures/figure-2-1.svg)

The equation and the source are identical in all four panels; **only the boundary condition changes**. The vertical axis is time $t\in[0,3]$ and the horizontal axis is space $x\in[0,1]$.

**Panel (a) uses homogeneous Dirichlet boundaries, and the thing to look at is that the four bright blobs never touch.** Each heating event forms an isolated blob that decays back to the dark blue background before the next one arrives, and the background stays dark blue throughout, which means the heat genuinely leaves. Computing the response to the fourth source on $t\in(1.7,2.2)$ therefore does not require the complete earlier solution. The paper offers a good analogy: predicting your living-room temperature a week or a month ahead mainly requires knowing whether the heater will be on and the windows closed, while the detailed earlier temperature field matters little.

**Panel (b) switches to homogeneous Neumann boundaries, and the thing to look at is the staircase in the background colour.** This is probably the most immediately readable picture in the chapter: the background steps from dark blue to blue, cyan, yellow-green, and yellow, moving up one level at each heating event. That staircase is the accumulated heat. Diffusion still smooths the local structure of each pulse, but the spatial mean cannot escape through the boundary, so the level only rises. The solution on $t\in(1.7,2.2)$ therefore still contains contributions from the first three sources, and **causality is genuinely binding here**. The paper's analogy is a perfectly insulated room, where you must know how often and how long the heater ran because the heat stays forever. It adds a practical remark: a real room always leaks heat slowly, which a Robin boundary condition would model better.

**Panel (c) uses periodic boundaries and is almost indistinguishable from (b).** The same staircase appears: local spatial variation is smoothed rapidly while the constant component keeps accumulating, and $t\in(1.7,2.2)$ is still influenced by the first three sources. Periodicity returns heat from one end to the other, so it cannot leave either.

**Panel (d) is a different experiment and also a frequency-selection test.** Here the source is switched off ($g=0$), the boundary stays periodic, and the initial condition is

$$
u_0(x)=\sin^2\!\left(8\pi(1-x)^2\right).
$$

Only a very thin band of fine spatial stripes survives near $t\approx0$; above it the panel becomes an almost uniform block of colour. That block is **not an absence of information**. It is the spatial mean of the initial condition, the one component that survives after the periodic heat equation rapidly damps every nonconstant Fourier mode. The paper adds a quantitative comparison: the surviving constant is roughly the same as the constant accumulated from the first two sources in panels (b) and (c). The apparently featureless upper part of the panel is precisely the record of low-frequency information persisting.

Taken together, the Dirichlet heat problem is highly local in time, a result discussed in Gander, Ohlberger and Rave (2024); this temporal locality can also be compared with spatial locality in solvation models from computational chemistry, see Ciaramella and Gander (2017, 2018a, 2018b). Neumann and periodic problems remain viable for PinT provided the algorithm can carry low-frequency components such as the constant mode over long times, and a coarse grid is the standard mechanism.

> **One sentence to remember**: high frequencies disappear on their own through dissipation and need no communication, while low frequencies never disappear and must be transported over long times by a coarse grid. Everything Parareal, PFASST, MGRiT, and STMG rely on in Chapter 4 sits in that sentence.

## 2.2 Advection–diffusion: $\nu$ turns the dial toward hyperbolic, but the boundary decides whether you can see it

Adding directed transport on the same interval gives

$$
\partial_tu(x,t)+\partial_xu(x,t)
-\nu\partial_{xx}u(x,t)
=g(x,t),
\qquad (x,t)\in\Omega\times(0,T], \tag{2.5}
$$

with $u(x,0)=u_0(x)$ and $\nu>0$. The transport speed is $1$ and $\nu$ controls the diffusion strength. Homogeneous Dirichlet and periodic boundaries form the main contrast; a footnote in the paper states that Neumann conditions would add no new qualitative behaviour, so no separate set of panels is shown.

![Source Figure 2.2: eight advection–diffusion solutions with Dirichlet and periodic boundaries](assets/papers/time-parallelization/source-figures/figure-2-2.svg)

The top row (a–d) uses zero Dirichlet boundaries and the bottom row (e–h) uses periodic boundaries. Panels (a–c) and (e–g) take $u_0=0$ with source (2.4) and diffusion $\nu=1$, $10^{-2}$, and $5\times10^{-4}$ from left to right. Panels (d) and (h) switch the source off, reuse the oscillatory initial condition from Figure 2.1(d), and take $\nu=5\times10^{-4}$.

**Reading the top row from left to right shows diffusion giving way to transport.** In panel (a) with $\nu=1$ diffusion dominates and the bands do not tilt, closely resembling Figure 2.1(a). In panel (b) the bands lean clearly to the right with slope set by the transport speed $1$, showing that transport now controls where the signal is. In panel (c) the trajectories are narrower and sharper, close to pure translation along characteristics, so fine scales survive longer. All three share one feature: **each band disappears once it reaches $x=1$**. The bottom row shows the other fate. Panel (e) retains only low frequencies such as the constant; panel (f) lets the signal return from $x=1$ to $x=0$, so the slanted trajectories of early pulses persist much longer; panel (g) recirculates equally sharp trajectories many times, and information created long before the current time still determines the fine structure of the present solution.

**Panels (d) and (h) are the pair to stare at.** Same equation, same $\nu=5\times10^{-4}$, same oscillatory initial condition, and **only the boundary condition differs**. In (d) the fine stripes occupy only the lower-left triangle, travel along the diagonal, leave through $x=1$, and by roughly $t=1$ the panel is a clean dark blue with nothing left. In (h) the same stripes wrap around repeatedly and fill the entire panel through $t=3$. This pair cleanly separates the effect of the outflow boundary from the effect of weak diffusion.

That yields the most important methodological warning in the chapter:

> [!warning] A Dirichlet outflow boundary hides long-range propagation
> For the Dirichlet problem every component eventually diffuses or leaves the domain regardless of $\nu$, so even at $\nu=5\times10^{-4}$ the interval $t\in(1.25,2.5)$ can be computed before the exact earlier solution is finished. The difficulty only appears when **periodic boundaries and small diffusion hold at the same time**. The practical consequence is direct: benchmarking a PinT method on an advection-dominated problem with Dirichlet outflow produces falsely optimistic results. The paper explicitly flags this as a point to watch when testing.

Under periodic boundaries the conclusion depends on $\nu$. Panel (e) retains only constants and other low frequencies, so a coarse grid may still work. Panels (f) and (g) retain progressively finer information as diffusion weakens, and panel (h) shows that high-frequency initial data can travel very far. A later interval can then no longer be determined independently of the earlier state, and the difficulty is most pronounced in the hyperbolic limit $\nu\to0$.

> [!note] Site supplement: phase and damping
> For a continuous Fourier mode $e^{ikx}$, propagation carries both the phase factor $e^{-ikt}$ and the diffusive damping $e^{-\nu k^2t}$, so mode $k$ has a survival timescale of about $1/(\nu k^2)$ and smaller $\nu$ keeps higher $k$ alive longer. The key point is that damping and phase are independent: a coarse propagator can be perfectly stable with the right amplitude and still place its correction at the wrong spatial location if its phase speed is off. This frequency-domain reading is provided to interpret the figure; the formula does not appear in Section 2.2 of the paper.

## 2.3 Burgers: nonlinearity adds the regeneration of high frequencies

Adding nonlinearity to the same comparison gives

$$
\begin{aligned}
\partial_tu(x,t)-\nu\partial_{xx}u(x,t)
+\frac12\partial_x\!\left(u^2(x,t)\right)
&=g(x,t),
&& (x,t)\in\Omega\times(0,T],\\
u(x,0)&=u_0(x),
&&x\in\Omega,
\end{aligned}
\tag{2.6}
$$

with $\nu>0$. The boundary conditions, source, initial data, and three diffusion values match Figure 2.2 panel for panel, so the two figures can be overlaid directly.

![Source Figure 2.3: eight Burgers solutions with Dirichlet and periodic boundaries](assets/papers/time-parallelization/source-figures/figure-2-3.svg)

**First compare (a) and (e) with the heat equation: at $\nu=1$ diffusion still dominates.** The four inputs in panel (a) stay separated and decay, while panel (e) retains an accumulated low-frequency background. Nonlinear transport has not yet produced visible sharp fronts.

**Next compare (b) and (c) with Figures 2.2(b) and (c): the shape of the bands has changed.** In the linear case the leading and trailing edges stay roughly parallel; here each band becomes an asymmetric wedge whose upper edge is markedly steeper. The reason is that regions of larger amplitude travel faster, so the waveform itself deforms. Panel (f) in the bottom row recirculates these deformed structures through the periodic boundary, so early inputs leave a persistent background and slanted trajectories. As diffusion weakens further, **even a smooth source** produces very steep edges, that is, viscous shocks carrying high spatial frequencies. In panel (c) these fronts eventually leave the Dirichlet domain, while panel (g) recirculates them so that fine scales travel far in space and time. Compared with linear transport, the new phenomena are shape change and shock generation.

**Finally compare (h) with Figure 2.2(h), the most telling pair.** The linear 2.2(h) shows many parallel stripes of uniform width, because modes simply translate. The Burgers 2.3(h) also begins with fine stripes near the bottom, but they merge rapidly into a few extremely sharp slanted fronts separated by broad smooth regions, which is shock coalescence. Panel (d) instead empties gradually under Dirichlet outflow, leaving only a weak tail at late times.

The distinction is **qualitative** and worth thinking through. In a linear problem high frequencies only decay, and smooth data can never grow fine structure by itself. Nonlinear transport actively steepens gradients and moves energy toward high frequencies. A coarse grid has a resolution ceiling, and Burgers dynamics keeps manufacturing structure above it. Even if the coarse grid could represent the initial data, it cannot represent what the dynamics creates, and it gets the shock **location** wrong as well.

The Dirichlet case still keeps its precomputability, since every component eventually diffuses or leaves, so $t\in(1.25,2.5)$ can be computed first. As $\nu\to0$ the problem approaches a hyperbolic limit with natural shocks, later solutions depend on the full frequency content of earlier ones, and the precomputation opportunity disappears. Periodic boundaries with small diffusion remain the essential stress-test combination — the wave equation in the next section propagates fine information over long times without needing either of those two extra conditions.

## 2.4 Second-order wave equation: the difficulty becomes unconditional

The hyperbolic reference model is

$$
\begin{aligned}
\partial_{tt}u(x,t)&=c^2\partial_{xx}u(x,t)+g(x,t),
&& (x,t)\in(0,1)\times(0,T],\\
u(x,0)&=u_0(x),
&&x\in(0,1),\\
\partial_tu(x,0)&=0,
&&x\in(0,1),
\end{aligned}
\tag{2.7}
$$

with $c>0$; the figure uses $c^2=0.2$. (Figure 2.4 floats ahead of the Section 2.4 heading in the paper layout, but its content and discussion belong to Section 2.4.)

![Source Figure 2.4: wave-equation solutions for three boundary conditions and oscillatory initial data](assets/papers/time-parallelization/source-figures/figure-2-4.svg)

**In panels (a), (b), and (c) the thing to look at is the V-shaped cone opening from $(x,t)=(0.5,0)$.** All three use zero initial displacement, zero initial velocity, and source (2.4). Each localized excitation sends a wave in **both** directions, and the edges of the cone are the characteristics. Waves reflect at the boundaries, so the region $t\gtrsim1.5$ is a superposition of four excitations after several reflections — the bright yellow patch near $t\approx2$ in (a) is exactly that. Neumann reflection in (b) changes how waves combine at the boundary without removing long-range propagation, and scalloped ripples from reflection are visible near the top. Periodicity in (c) merely reconnects the paths; the preservation of detail is unchanged.

**Panel (d) is the picture of having nowhere to escape.** It switches off the source and uses periodic boundaries, zero initial velocity, and the same oscillatory initial condition $u_0(x)=\sin^2(8\pi(1-x)^2)$. Multiple frequencies travel in both directions and superpose repeatedly, filling all of $0<t<3$ with a dense interference texture formed by crossing left- and right-going characteristics, with no sign of decay. The paper notes that this detailed dependence on every initial frequency would occur equally under Dirichlet or Neumann boundaries.

The key argument of this section is a causal explanation, and it is the sharpest point in the chapter:

> Advection–diffusion and Burgers **need** periodic boundaries to expose the difficulty because the transport term is first order and propagation has a direction — information travels one way, and a Dirichlet boundary is an exit. The wave equation propagates in several directions, and Dirichlet and Neumann boundaries **reflect rather than absorb**, so there is no exit at all.

Add one more fact: the wave equation has no dissipative term, so no mode decays and phase error only accumulates rather than being smoothed away. A coarse temporal model's phase or wave-speed error therefore cannot be rescued by the dynamics, and an effective method must transport fine detail deliberately. Chapter 3 turns to SWR, parallel IDC, ParaExp, and ParaDiag for exactly this reason.

## Site synthesis: the four models as one memory spectrum

Placing the four source-free panels with oscillatory data side by side — Figures 2.1(d), 2.2(h), 2.3(h), and 2.4(d) — gives the clearest statement of the chapter's conclusion. From the same initial condition $u_0=\sin^2(8\pi(1-x)^2)$ come four fates: **a uniform block of colour, dense parallel stripes, a few sharp shock fronts, and a full-domain interference grid.**

![The four model equations span temporal locality to long-range memory](assets/diagrams/pint/en/model-memory-spectrum.svg)

| Model               | Dominant mechanism                         | Information retained over long times                                                                                | PinT implication in the paper                                |
| ------------------- | ------------------------------------------ | ------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------ |
| heat                | diffusion                                  | rapid forgetting with Dirichlet boundaries; constants and other low frequencies with Neumann or periodic boundaries | communicate slow modes accurately on coarse levels           |
| advection–diffusion | directed transport and diffusion           | phase and fine scales under periodic boundaries and small $\nu$                                                     | move progressively finer information across long times       |
| Burgers             | nonlinear transport, shocks, and diffusion | high-frequency fronts that persist and are regenerated under periodic boundaries and small $\nu$                    | track shock location and shape in a coarse representation    |
| wave                | bidirectional propagation and reflection   | amplitude, phase, and propagation paths under every boundary type                                                   | use a temporal mechanism designed for long-range propagation |

The three criteria compress into three sentences: **boundary conditions decide whether information can leave; the diffusion parameter decides how fast fine scales decay; nonlinearity decides whether high frequencies are regenerated.** When analyzing any later algorithm, these three questions test whether its experimental setup actually reaches the difficulty.

> [!tip] Site synthesis: the one-sentence version
> A coarse propagator is essentially a low-pass filter. Parabolic dynamics happens to be a low-pass filter too, since it discards high frequencies on its own, so the coarse solver and the true dynamics discard the same information and agree on what matters. Hyperbolic dynamics is not a low-pass filter; it carries every frequency unchanged far into the future, so the coarse solver discards precisely the information that determines the answer. That is why Chapter 4's methods fail on Chapter 3's problems, and why the paper abandons the classification by algorithmic lineage in favour of a classification by dynamics.

## Three questions to test yourself

1. Why can the Dirichlet heat equation be solved late-first while the Neumann heat equation cannot? The answer lies in whether information can leave the domain, not in the equation itself.
2. Why must periodic boundaries be used when testing PinT methods on transport-dominated problems? Because Dirichlet outflow hides the difficulty and the method's weakness together.
3. What exactly does Burgers add over advection–diffusion? Not vague "more nonlinearity" but the regeneration of high frequencies, which a coarse grid's resolution ceiling can never catch up with.

## Site numerical supplement: three recomputed solutions

The following Python experiments contrast diffusion plus outflow with persistent propagation. They use independent initial data, grids, and source-free settings. They do not reproduce the source-driven experiments in Figures 2.2–2.4. Their numerical values are reported separately from the paper figures.

### Advection–diffusion solution

The experiment uses homogeneous Dirichlet boundaries, $\Delta t=\Delta x=10^{-3}$, $T=3$, $\nu=5\times10^{-4}$, unit advection speed, and initial data $\sin^2(8\pi(1-x)^2)$. The final $L^\infty$ norm is $7.022\times10^{-99}$.

![Recomputed space–time solution of the advection–diffusion equation](assets/pint/model-advection-diffusion.svg)

The near-zero final state results from the present combination of diffusion and Dirichlet outflow. A periodic problem retains recirculating signals.

### Viscous Burgers solution

The experiment uses $\Delta t=\Delta x=1/400$, $T=3$, $\nu=5\times10^{-4}$, and homogeneous Dirichlet boundaries. The maximum over space and time is $1.045940$, and the final $L^\infty$ norm is $0.325871$.

![Recomputed space–time solution of the viscous Burgers equation](assets/pint/model-burgers.svg)

### Wave solution

The source-free experiment uses the trapezoidal rule with $\Delta t=\Delta x=1/400$, $T=3$, and $c^2=0.2$. It uses the same oscillatory initial displacement as the previous experiments and zero initial velocity. The final displacement has $L^\infty$ norm $0.948217$.

![Recomputed space–time solution of the wave equation](assets/pint/model-wave.svg)

### Numerical summary

| Supplemental experiment | Grid and horizon                   |                      Reported metric |
| ----------------------- | ---------------------------------- | -----------------------------------: |
| advection–diffusion     | $\Delta x=\Delta t=10^{-3}$, $T=3$ | final $L^\infty=7.022\times10^{-99}$ |
| viscous Burgers         | $\Delta x=\Delta t=1/400$, $T=3$   |            final $L^\infty=0.325871$ |
| wave equation           | $\Delta x=\Delta t=1/400$, $T=3$   |            final $L^\infty=0.948217$ |

These quantities characterize the discrete dynamics. They do not measure PinT convergence or hardware speedup.

## Quick cross-reference

Use the table below when reading section by section alongside the paper.

| Source location                     | Page section                     | Equations and figures        | Main question                                                   |
| ----------------------------------- | -------------------------------- | ---------------------------- | --------------------------------------------------------------- |
| Section 2 introduction, pp. 388–389 | Two notations and one exception  | (2.1)–(2.2)                  | Why four PDEs link parabolic and hyperbolic dynamics            |
| Section 2.1, pp. 389–391            | 2.1 Heat equation                | (2.3)–(2.4), Figure 2.1(a–d) | How diffusion and boundaries create temporal locality           |
| Section 2.2, pp. 391–393            | 2.2 Advection–diffusion equation | (2.5), Figure 2.2(a–h)       | How far fine scales travel as diffusion weakens                 |
| Section 2.3, pp. 393–395            | 2.3 Burgers' equation            | (2.6), Figure 2.3(a–h)       | How nonlinearity generates and transports shocks                |
| Section 2.4, pp. 395–396            | 2.4 Second-order wave equation   | (2.7), Figure 2.4(a–d)       | Why multidirectional propagation and reflection preserve memory |

## Source

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Section 2, pp. 388–396.
