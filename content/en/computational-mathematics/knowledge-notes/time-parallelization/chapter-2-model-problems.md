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

## Starting from a counter-intuitive question

Section 2 contains no algorithms. It settles a more basic question, namely which evolution problems are suited to time parallelization at all, and the reason Chapters 3 and 4 group methods by dynamics rather than by algorithmic lineage comes entirely from here.

The starting point looks like it contradicts common sense. Causality tells us that the late solution is determined by the early one, which makes the time direction appear inherently sequential, yet the paper opens by asking whether we can compute $t\in(1.7,2.2)$ before knowing the solution for $t<1.7$. For the heat equation with Dirichlet boundaries the answer is yes, and the paper marks it with an exclamation point.

There is no contradiction with causality. Causality describes a logical dependence, whereas what matters in computation is the numerically effective dependence: once dissipation has erased old information, its influence on the late solution falls below the discretization error, and the requirement to know the early solution first no longer holds in any computational sense. The rest of the chapter measures one thing throughout, which is how long the solution's memory lasts and which frequencies that memory keeps. Those two answers decide directly whether time parallelization can work.

## The four figures form a single controlled experiment

The four figures look like four unrelated computations, but they are in fact one carefully controlled experiment. The paper reuses the same source term and the same oscillatory initial condition across all four equations and varies only three things: the equation itself, moving from the heat equation through advection–diffusion and Burgers to the wave equation and adding transport and nonlinearity along the way; the boundary condition, switching among Dirichlet, Neumann, and periodic to decide whether information can leave the domain; and the diffusion strength $\nu$, taking $1$, $10^{-2}$, and $5\times10^{-4}$ to decide how fast fine scales decay. Because everything else is held fixed, the panels can be compared directly against one another.

One thing is easy to misread. The last column of each figure does not belong to the same class of experiment: Figures 2.1(d), 2.2(d) and (h), 2.3(d) and (h), and 2.4(d) all switch the source off and use the oscillatory initial condition instead. The earlier columns ask how long the response to an external forcing survives, while the last column asks how long the initial information itself survives, and the two should not be compared with each other.

Two pieces of notation are needed before the models themselves. Most time-parallel methods are described and analyzed on the system of ordinary differential equations obtained after spatial semidiscretization. The linear case is written

$$
\begin{aligned}
\boldsymbol u'(t)&=A\boldsymbol u(t)+\boldsymbol g(t),
&&t\in(0,T],\\
\boldsymbol u(0)&=\boldsymbol u_0,
\end{aligned}
\tag{2.1}
$$

where $A\in\mathbb R^{N_x\times N_x}$ comes from semidiscretizing the PDE in space, and the nonlinear case is written

$$
\begin{aligned}
\boldsymbol u'(t)&=\boldsymbol f(\boldsymbol u(t),t),
&&t\in(0,T],\\
\boldsymbol u(0)&=\boldsymbol u_0,
\end{aligned}
\tag{2.2}
$$

where $\boldsymbol f:\mathbb R^{N_x}\times\mathbb R\to\mathbb R^{N_x}$ is nonlinear in its first argument and the Burgers-type semidiscretization can be written $\boldsymbol f(\boldsymbol u(t),t)=A\boldsymbol u(t)+B\boldsymbol u^2(t)+\boldsymbol g(t)$. Domain-decomposition methods are the exception to this notation. They are not reduced to a system of ODEs but are solved and analyzed directly on continuous space–time subdomains, and the difference becomes visible as soon as one reaches SWR in Chapter 3.

All four models are posed on the one-dimensional unit interval $\Omega=(0,1)$. The paper explains that this is not a real restriction, since the applicability and convergence properties of time-parallel methods generally do not depend on the spatial dimension. Dimension changes the cost of a single spatial solve considerably, but the three mechanisms this chapter cares about are unaffected: how old information decays, whether fine scales survive, and how information leaves or returns through a boundary.

## 2.1 Heat equation

The parabolic reference model is

$$
\partial_tu(x,t)=\partial_{xx}u(x,t)+g(x,t),
\qquad (x,t)\in\Omega\times(0,T], \tag{2.3}
$$

with initial data $u(x,0)=u_0(x)$. The text itself mentions only homogeneous Dirichlet and Neumann boundaries; the figure adds periodic boundaries in order to compare whether heat can leave the domain.

The first three panels use zero initial data, and all the heat comes from a source that is localized in both space and time:

$$
g(x,t)=10\sum_{j=1}^{4}
\exp\!\left(
-\sigma\left[(t-t_j)^2+(x-0.5)^2\right]
\right), \tag{2.4}
$$

with $(t_1,t_2,t_3,t_4)=(0.1,0.6,1.35,1.85)$ and $\sigma=200$. The source sits fixed at $x=0.5$ and its four time centres are well separated. Taking $\sigma=200$ puts the characteristic width at roughly $1/\sqrt{200}\approx0.07$ in both space and time, so it appears in the space–time plot as four non-overlapping bright bands.

![Source Figure 2.1: heat-equation solutions for three boundary conditions and oscillatory initial data](assets/papers/time-parallelization/source-figures/figure-2-1.svg)

The equation and the source are identical in all four panels and only the boundary condition changes. The vertical axis is time $t\in[0,3]$ and the horizontal axis is space $x\in[0,1]$.

Panel (a) uses homogeneous Dirichlet boundaries. The four bright blobs never touch: each temperature peak has already decayed back into the dark blue background before the next heating event arrives, and the background stays dark blue from beginning to end, which means the heat really does flow out through the boundary. Given that, computing the response to the fourth source on $t\in(1.7,2.2)$ does not require an accurate solution at earlier times. The paper makes the point with an everyday example: predicting your living-room temperature a week or a month into winter depends on whether the heater will be running and the windows closed, while the detailed temperature field in the room right now is almost irrelevant.

Panel (b) switches to homogeneous Neumann boundaries, and what changes in the picture is the background colour. It rises from dark blue through blue, cyan, yellow-green, and yellow, stepping up one level at each heating event. Diffusion still smooths away the local structure left by each pulse, but the spatial mean cannot escape through the boundary, so heat only enters and never leaves, and the background level becomes a reading of the accumulated heat. The solution on $t\in(1.7,2.2)$ therefore still carries contributions from the first three sources, and causality is genuinely binding here. The paper's analogy is a perfectly insulated room: to know how warm it is now you must know how often and for how long the heater ran, because the heat that entered stays. The paper adds that a real room always leaks heat slowly, which a Robin boundary condition models more realistically.

Panel (c) uses periodic boundaries, and its long-time behaviour is almost indistinguishable from the Neumann case, showing the same stepwise rise in level. Local variation is smoothed quickly, the constant component keeps accumulating, and $t\in(1.7,2.2)$ is still influenced by the first three sources. Only the mechanism differs: periodic boundaries let heat travel around from one end to the other rather than reflecting it back, but either way it cannot escape.

Panel (d) belongs to the other class of experiment. Here the source is switched off, the boundary stays periodic, and the initial condition is

$$
u_0(x)=\sin^2\!\left(8\pi(1-x)^2\right).
$$

Only a very thin band of stripes survives near $t\approx0$, and above it the panel becomes an almost uniform block of colour. That block does not mean the information has vanished. It is the one component that survives after the periodic heat equation rapidly damps every nonconstant Fourier mode, namely the spatial mean of the initial condition. The paper adds a quantitative comparison: this surviving constant is roughly the same as the constant accumulated from the first two sources in panels (b) and (c). The upper part of the picture looks featureless, but it is precisely the record of low-frequency information persisting.

Putting the four panels together settles the heat equation. The Dirichlet case is highly local in time, as analyzed in Gander, Ohlberger and Rave (2024), and this temporal locality can be compared with the spatial locality of solvation models in computational chemistry, see Ciaramella and Gander (2017, 2018a, 2018b). The Neumann and periodic cases still leave room for time parallelization, provided the algorithm can carry low-frequency components such as the constant mode accurately over long times, and a coarse grid is the mechanism the paper points to. Put differently, high-frequency error disappears on its own through dissipation and needs no communication, while low-frequency components never decay by themselves and must be transported across time. Everything Parareal, PFASST, MGRiT, and STMG do in Chapter 4 rests on that premise.

## 2.2 Advection–diffusion equation

Adding directed transport on the same interval gives

$$
\partial_tu(x,t)+\partial_xu(x,t)
-\nu\partial_{xx}u(x,t)
=g(x,t),
\qquad (x,t)\in\Omega\times(0,T], \tag{2.5}
$$

with $u(x,0)=u_0(x)$ and $\nu>0$. The transport speed is $1$ and $\nu$ controls the diffusion strength. Homogeneous Dirichlet and periodic boundaries form the main contrast, and a footnote in the paper explains that Neumann conditions would produce no new qualitative behaviour, so no separate set of panels is shown.

![Source Figure 2.2: eight advection–diffusion solutions with Dirichlet and periodic boundaries](assets/papers/time-parallelization/source-figures/figure-2-2.svg)

The top row (a–d) uses zero Dirichlet boundaries and the bottom row (e–h) uses periodic boundaries. Panels (a–c) and (e–g) take $u_0=0$ with source (2.4) and diffusion $\nu=1$, $10^{-2}$, and $5\times10^{-4}$ from left to right, while (d) and (h) switch the source off, reuse the oscillatory initial condition from Figure 2.1(d), and take $\nu=5\times10^{-4}$.

Reading the top row from left to right shows diffusion gradually giving way to transport. In panel (a) with $\nu=1$ diffusion still dominates, the bands do not tilt, and the picture resembles the Dirichlet heat case. By panel (b) the bands lean clearly to the right with a slope equal to the transport speed $1$, and the position of the signal is now set by transport. In panel (c) the trajectories are narrower and cleaner, close to pure translation along characteristics, and fine scales survive longer. All three share one feature: once a band reaches $x=1$ it leaves the domain and does not come back.

The bottom row shows the other outcome. Panel (e) resembles the periodic heat case, retaining only low frequencies such as the constant. Panel (f) returns the signal from $x=1$ to $x=0$, so the slanted trajectories left by early pulses run much longer, and diffusion remains visible as the stripes broaden during propagation. Panel (g) sends equally sharp trajectories around the periodic boundary again and again, so information created long before the present moment still determines the fine structure of the current solution.

Panels (d) and (h) are the pair worth comparing most closely in this section. The equation, the value $\nu=5\times10^{-4}$, and the oscillatory initial condition are identical, and the only difference is the boundary condition. In (d) the fine stripes occupy only the triangular region in the lower left, travel along the diagonal, leave through $x=1$, and above roughly $t=1$ the panel is a clean dark blue with nothing left. In (h) the same stripes wrap around repeatedly and fill the entire panel through $t=3$. This pair cleanly separates the role of the outflow boundary from the role of weak diffusion: the fine scales vanish in (d) not because diffusion smoothed them away but because the boundary carried them off.

> [!warning] Dirichlet outflow hides long-range propagation
> For the Dirichlet problem every component eventually diffuses away or leaves the domain regardless of $\nu$, so even at $\nu=5\times10^{-4}$ the interval $t\in(1.25,2.5)$ can be computed before the exact earlier solution is finished. The real difficulty appears only when periodic boundaries and small diffusion hold at the same time. This has a direct consequence for benchmarking: testing a transport-dominated problem with Dirichlet outflow yields an over-optimistic conclusion, and the paper explicitly lists it as something to watch for when testing.

Under periodic boundaries the conclusion changes with $\nu$. Panel (e) retains only constants and other low frequencies over long times, so a coarse grid may still suffice. As diffusion weakens, the information retained in (f) and (g) becomes progressively finer, and by (h) high-frequency initial data travels very far. A later time interval can then no longer be determined independently of the earlier state, and the difficulty is sharpest in the hyperbolic limit $\nu\to0$.

> [!note] Site supplement: phase and damping
> For a continuous Fourier mode $e^{ikx}$, propagation carries both the phase factor $e^{-ikt}$ and the diffusive damping $e^{-\nu k^2t}$, so mode $k$ has a survival timescale of about $1/(\nu k^2)$ and smaller $\nu$ keeps higher $k$ alive longer. The point is that damping and phase are independent of each other: a coarse propagator can be perfectly stable and get the amplitude right, yet if its phase speed is off, the correction it returns lands at the wrong spatial location. This frequency-domain reading is provided to interpret the figure; the formula does not appear in Section 2.2 of the paper.

## 2.3 Burgers' equation

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

with $\nu>0$. The boundary conditions, source, initial data, and three diffusion values match Figure 2.2 panel for panel, so the two figures can be read on top of each other.

![Source Figure 2.3: eight Burgers solutions with Dirichlet and periodic boundaries](assets/papers/time-parallelization/source-figures/figure-2-3.svg)

At $\nu=1$ diffusion still dominates and nonlinear transport has not yet shown its effect: the four inputs in panel (a) stay separated and decay individually, while panel (e) retains an accumulated low-frequency background, so the two resemble the Dirichlet and periodic heat cases respectively.

Placing (b) and (c) next to the corresponding panels of Figure 2.2 reveals the first change that nonlinearity brings, which is that the shape of the bands is different. In the linear case the leading and trailing edges of a band stay roughly parallel, whereas here each band becomes an asymmetric wedge with a markedly steeper upper edge. The reason is that regions of larger amplitude travel faster, so the waveform keeps deforming as it propagates. Panel (f) in the bottom row recirculates these already deformed structures through the periodic boundary, so early inputs leave a persistent background and slanted trajectories behind them. Once diffusion weakens further, even a smooth source produces very steep edges in the solution, that is, viscous shocks carrying high spatial frequencies. In panel (c) these fronts eventually leave the Dirichlet domain, while panel (g) brings them back around so that fine scales travel far in space and time.

The second change requires comparing (h) with Figure 2.2(h). The linear 2.2(h) consists of many parallel stripes of uniform width, because the individual modes simply translate without interacting. The Burgers panel 2.3(h) also begins with fine stripes near the bottom, but they merge quickly into a few extremely sharp slanted fronts separated by broad smooth regions, which is shocks coalescing. Panel (d), by contrast, empties gradually under Dirichlet outflow and leaves only a weak tail at later times.

This difference is qualitative and worth thinking through. In a linear problem high frequencies only decay, and smooth data can never grow fine structure by itself; nonlinear transport instead actively steepens gradients and moves energy toward high frequencies. A coarse grid has a fixed resolution ceiling, and Burgers dynamics keeps manufacturing structure above it, so even if the coarse grid could represent the initial data it cannot represent what the dynamics creates, and it will get the shock location wrong as well.

The Dirichlet case still keeps the precomputability that outflow provides, since every component eventually diffuses away or leaves the domain, and $t\in(1.25,2.5)$ can still be computed first. As $\nu\to0$ the problem approaches a hyperbolic limit in which shocks are natural, later solutions depend on the full frequency content of earlier ones, and the opportunity to precompute disappears. Exposing this difficulty still requires periodic boundaries and small diffusion together, whereas the wave equation in the next section needs neither of those extra conditions.

## 2.4 Second-order wave equation

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

with $c>0$, and the figure uses $c^2=0.2$. Figure 2.4 floats ahead of the Section 2.4 heading in the paper layout, but its content and discussion belong to Section 2.4.

![Source Figure 2.4: wave-equation solutions for three boundary conditions and oscillatory initial data](assets/papers/time-parallelization/source-figures/figure-2-4.svg)

Panels (a), (b), and (c) all use zero initial displacement, zero initial velocity, and source (2.4), and what they share visually is a V-shaped cone opening from $(x,t)=(0.5,0)$. Each localized excitation sends one wave in each direction, and the two edges of the cone are the characteristics. The waves reflect once they reach a boundary, so the region $t\gtrsim1.5$ is a superposition of four excitations after several reflections, which is how the bright yellow patch near $t\approx2$ in panel (a) forms. Neumann reflection in panel (b) changes the way waves combine at the boundary without weakening long-range propagation, and scalloped ripples from reflection remain visible near the top. The periodic boundary in panel (c) merely reconnects the propagation paths, leaving the long-term preservation of detail unchanged.

Panel (d) switches the source off and uses periodic boundaries, zero initial velocity, and the same oscillatory initial condition $u_0(x)=\sin^2(8\pi(1-x)^2)$. The many frequencies in the initial data travel in both directions and superpose repeatedly, filling the whole interval $0<t<3$ with a dense interference texture woven from left- and right-going characteristics, with no sign of decay anywhere. The paper notes that this detailed dependence on every initial frequency is equally present under Dirichlet and Neumann boundaries.

The paper's explanation for why the wave equation exposes the difficulty without needing periodic boundaries is direct. The transport term in advection–diffusion and Burgers is first order and propagation has a definite direction, so information travels one way and a Dirichlet boundary becomes an exit, which is why the exit has to be closed off with periodicity before the difficulty shows. The wave equation propagates in several directions, and Dirichlet and Neumann boundaries reflect rather than absorb, so no exit exists at all. On top of that the equation has no dissipative term, no mode decays, and phase error simply accumulates instead of being smoothed away by the dynamics. Once a coarse temporal model is off in phase or wave speed, the dynamics will not correct it, so an effective method has to take responsibility for transporting fine detail over long times. Chapter 3 turns to SWR, parallel IDC, ParaExp, and ParaDiag precisely to address this.

## Site synthesis: the four models as one memory spectrum

Placing the four source-free panels with oscillatory data side by side, namely Figures 2.1(d), 2.2(h), 2.3(h), and 2.4(d), presents the chapter's conclusion in its most direct form. The same initial condition $u_0=\sin^2(8\pi(1-x)^2)$ meets four fates across the four equations: first it collapses into a uniform block of colour, then it becomes dense parallel stripes, then it condenses into a few sharp shock fronts, and finally it fills the whole space–time domain with an interference grid.

![The four model equations span temporal locality to long-range memory](assets/diagrams/pint/en/model-memory-spectrum.svg)

| Model               | Dominant mechanism                         | Information retained over long times                                                                                | PinT implication in the paper                                |
| ------------------- | ------------------------------------------ | ------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------ |
| heat                | diffusion                                  | rapid forgetting with Dirichlet boundaries; constants and other low frequencies with Neumann or periodic boundaries | communicate slow modes accurately on coarse levels           |
| advection–diffusion | directed transport and diffusion           | phase and fine scales under periodic boundaries and small $\nu$                                                     | move progressively finer information across long times       |
| Burgers             | nonlinear transport, shocks, and diffusion | high-frequency fronts that persist and are regenerated under periodic boundaries and small $\nu$                    | track shock location and shape in a coarse representation    |
| wave                | bidirectional propagation and reflection   | amplitude, phase, and propagation paths under every boundary type                                                   | use a temporal mechanism designed for long-range propagation |

Judging whether a test configuration actually reaches the difficulty takes only three questions: does the boundary condition let information leave, how fast does the diffusion parameter make fine scales decay, and will nonlinearity regenerate high frequencies. All three can be applied directly when analyzing any of the algorithms that follow.

The differences among the four models ultimately come down to the relationship between the coarse and fine representations, and the following summary is this site's own framing rather than the paper's. A coarse propagator is essentially a low-pass filter. Parabolic dynamics happens to be a low-pass filter as well, discarding high frequencies of its own accord, so the coarse solver and the true dynamics discard the same information and agree on what matters. Hyperbolic dynamics has no such property and carries every frequency unchanged far into the future, so the coarse solver discards exactly the information that determines the answer. That is why the methods of Chapter 4 fail on the problems of Chapter 3, and why the paper abandons a classification by algorithmic lineage in favour of one by dynamics.

## Site numerical supplement: three recomputed solutions

The following Python experiments contrast diffusion combined with outflow against persistent propagation. They use their own initial data, grids, and source-free settings rather than reproducing the source-driven experiments of Figures 2.2–2.4, and their numbers are reported separately from the paper figures.

### Advection–diffusion solution

The experiment uses homogeneous Dirichlet boundaries, $\Delta t=\Delta x=10^{-3}$, $T=3$, $\nu=5\times10^{-4}$, unit advection speed, and initial data $\sin^2(8\pi(1-x)^2)$. The final $L^\infty$ norm is $7.022\times10^{-99}$.

![Recomputed space–time solution of the advection–diffusion equation](assets/pint/model-advection-diffusion.svg)

The near-zero final state follows from this particular combination of diffusion and Dirichlet outflow, whereas a periodic problem would retain recirculating signals.

### Viscous Burgers solution

The experiment uses $\Delta t=\Delta x=1/400$, $T=3$, $\nu=5\times10^{-4}$, and homogeneous Dirichlet boundaries. The maximum over space and time is $1.045940$, and the final $L^\infty$ norm is $0.325871$.

![Recomputed space–time solution of the viscous Burgers equation](assets/pint/model-burgers.svg)

### Wave solution

The source-free experiment uses the trapezoidal rule with $\Delta t=\Delta x=1/400$, $T=3$, and $c^2=0.2$. It reuses the oscillatory initial displacement of the previous experiments together with zero initial velocity. The final displacement has $L^\infty$ norm $0.948217$.

![Recomputed space–time solution of the wave equation](assets/pint/model-wave.svg)

### Numerical summary

| Supplemental experiment | Grid and horizon                   |                      Reported metric |
| ----------------------- | ---------------------------------- | -----------------------------------: |
| advection–diffusion     | $\Delta x=\Delta t=10^{-3}$, $T=3$ | final $L^\infty=7.022\times10^{-99}$ |
| viscous Burgers         | $\Delta x=\Delta t=1/400$, $T=3$   |            final $L^\infty=0.325871$ |
| wave equation           | $\Delta x=\Delta t=1/400$, $T=3$   |            final $L^\infty=0.948217$ |

These quantities describe the dynamics of the discrete solutions. They do not measure time-parallel iteration speed or hardware speedup.

## Quick cross-reference

Use the table below when reading section by section alongside the paper.

| Source location                     | Page section                                | Equations and figures        | Main question                                                   |
| ----------------------------------- | ------------------------------------------- | ---------------------------- | --------------------------------------------------------------- |
| Section 2 introduction, pp. 388–389 | The four figures as a controlled experiment | (2.1)–(2.2)                  | Why four PDEs link parabolic and hyperbolic dynamics            |
| Section 2.1, pp. 389–391            | 2.1 Heat equation                           | (2.3)–(2.4), Figure 2.1(a–d) | How diffusion and boundaries create temporal locality           |
| Section 2.2, pp. 391–393            | 2.2 Advection–diffusion equation            | (2.5), Figure 2.2(a–h)       | How far fine scales travel as diffusion weakens                 |
| Section 2.3, pp. 393–395            | 2.3 Burgers' equation                       | (2.6), Figure 2.3(a–h)       | How nonlinearity generates and transports shocks                |
| Section 2.4, pp. 395–396            | 2.4 Second-order wave equation              | (2.7), Figure 2.4(a–d)       | Why multidirectional propagation and reflection preserve memory |

## Source

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Section 2, pp. 388–396.
