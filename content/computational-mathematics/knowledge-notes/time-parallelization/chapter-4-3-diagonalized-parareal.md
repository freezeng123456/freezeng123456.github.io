---
title: 4.5：基于对角化的 Parareal
description: 并行粗网格校正与区间内对角化粗传播子的两条完整路线
lang: zh
translation: en/computational-mathematics/knowledge-notes/time-parallelization/chapter-4-3-diagonalized-parareal
tags:
  - 时间并行
  - Parareal
  - ParaDiag
---

> [!note] 阅读范围
> 本页对应论文 Section 4.5（pp. 460–471），覆盖公式 (4.14)–(4.29)、Theorems 4.7–4.8、Remark 4.2 和 Figures 4.12–4.17。两种方法都使用对角化，作用位置和适用范围不同：第一种并行化跨粗点的 CGC，第二种在每个粗区间内构造可并行的特殊粗传播子。

## 4.5 基于对角化的 Parareal

### 两条路线先分清

- **对角化 CGC（Section 4.5.1，Wu 2018；Wu 与 Zhou 2019）**：修改 Parareal 跨 $N_t$ 个粗点的顺序粗校正。并行宽度来自粗时间点；收敛机制仍接近标准 Parareal，主要适合抛物问题。
- **对角化粗传播子（Section 4.5.2，Gander 与 Wu 2020）**：保留标准 Parareal 粗校正的外形，在每个 $[T_n,T_{n+1}]$ 内，用 ParaDiag 同时处理 $J$ 个细步。粗、细传播使用同一个积分器和步长；这种粗传播能长时间输运**全部频率分量**，因此也能处理双曲问题。

## 4.5.1 基于对角化的 CGC

### 从顺序 CGC 到首尾耦合

标准 Parareal 的粗网格校正为

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

并从 $\boldsymbol u_0^{k+1}=\boldsymbol u_0$ 顺序推进。Wu (2018) 改用

$$
\boldsymbol u_0^{k+1}=\alpha\boldsymbol u_{N_t}^{k+1}+\boldsymbol u_0.
$$

为了让收敛极限仍满足原初值问题，需定义

$$
\widetilde{\boldsymbol u}_n^k=
\begin{cases}
\boldsymbol u_0,&n=0,\\
\boldsymbol u_n^k,&n\ge1,
\end{cases}
$$

并使用

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

这里的首尾条件早于 ParaDiag-II 的自然条件 (3.55)，形式也略有不同；修正 $\widetilde{\boldsymbol u}_0^k$ 保证了极限一致性。

### 线性全时间系统与三步解法

对 $\boldsymbol u'=A\boldsymbol u$，粗层取后向 Euler。把首尾条件代入第一条粗方程，得到

$$
(C_\alpha\otimes I_x-I_t\otimes\Delta TA)\boldsymbol U^{k+1}
=\boldsymbol g^k, \tag{4.16}
$$

其中

$$
C_\alpha=
\begin{bmatrix}
1&&&-\alpha\\
-1&1\\
&\ddots&\ddots\\
&&-1&1
\end{bmatrix},
$$

$$
\boldsymbol g^k=
\begin{bmatrix}
\boldsymbol u_0+(I_x-\Delta TA)\boldsymbol b_1^k\\
(I_x-\Delta TA)\boldsymbol b_2^k\\
\vdots\\
(I_x-\Delta TA)\boldsymbol b_{N_t}^k
\end{bmatrix}.
$$

对角化时间矩阵后得到

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

这里 $\{\lambda_n\}$ 是 $C_\alpha$ 的特征值，$F$ 是离散 Fourier 矩阵，分别对应 (3.51) 和 (3.50)。

> [!note] 本站补充：$\alpha\ne1$ 时的对角缩放
> (4.17) 按原文排印只写了 $F$ 与 $F^*$，而这只在 $\alpha=1$ 时对角化 $C_\alpha$。按论文 Section 3.5.2 的 (3.59)，一般情形是 $C_\alpha=V_\alpha D_\alpha V_\alpha^{-1}$，其中 $V_\alpha=\Lambda_\alpha F^*$、$\Lambda_\alpha=\operatorname{diag}(1,\alpha^{-1/N_t},\ldots,\alpha^{-(N_t-1)/N_t})$，即实现中还需要相应的对角缩放。

与 Section 3.5.2 相同，核心是 FFT 变换、独立移位空间求解和逆变换。CGC 因而能在全部粗点上同时完成。

### Theorem 4.7：保持标准 Parareal 速度的阈值

当 $\alpha\to0$，(4.15) 回到标准 CGC；过小 $\alpha$ 又会放大 $\alpha$-循环对角化的舍入误差，在单精度或半精度下尤其明显。Theorem 4.7 引自 Wu（2018）：设标准 Parareal (4.14) 的收敛因子为 $\rho$，新变体 (4.15) 的收敛因子为 $\rho_{\mathrm{new}}$，粗求解器 $\mathcal G$ 为稳定积分器，则

$$
\rho_{\mathrm{new}}=\rho,
\qquad\text{只要}\qquad
\alpha\le\frac{\rho}{1+\rho}.
$$

Theorem 4.7 只对 $A$ 具有负实特征值的线性问题
$\boldsymbol u'=A\boldsymbol u+\boldsymbol g$ 得到证明；复谱目前只有
数值证据。

因此实用选择为阈值本身 $\alpha=\rho/(1+\rho)$：再减小不会改善渐近速度，只会增加舍入风险。通常 $\rho=O(10^{-1})$，所以 $\alpha$ 也在 $10^{-1}$ 量级，此时对角化引入的舍入误差可以忽略。

![原论文 Figure 4.12：标准与对角化 CGC 在热方程和 ADE 上的误差](assets/papers/time-parallelization/source-figures/figure-4-12.svg)

实验使用周期边界、$u_0(x)=\sin(2\pi x)$、后向 Euler 粗层、SDIRK22 细层、$T=4$、$J=10$、$\Delta T=0.1$、$\Delta x=1/128$。(a) 是热方程，测得 $\rho\approx0.22$，阈值约 $0.18$；因此 $\alpha=0.25,0.4$ 慢于标准 CGC，$\alpha=0.1$ 与其重合。(b) 是 $\nu=0.1$ 的 ADE，$\rho\approx0.39$、阈值约 $0.28$；这里 $\alpha=0.1,0.25$ 都能跟上标准 CGC，$\alpha=0.4$ 明显变慢。注意 ADE 的半离散矩阵具有复特征值，正是 Theorem 4.7 未覆盖、只有数值证据支持的情形。

### 非线性全时间准 Newton

粗层继续用后向 Euler，定义

$$
\boldsymbol b_{n+1}^k
=\mathcal F(T_n,T_{n+1},\widetilde{\boldsymbol u}_n^k)
-\mathcal G(T_n,T_{n+1},\boldsymbol u_n^k).
$$

粗校正与首尾条件形成

$$
(C_\alpha\otimes I_x)\boldsymbol U^{k+1}
-\Delta T F(\boldsymbol U^{k+1})=\boldsymbol g^k. \tag{4.18}
$$

$F$ 的第 $n$ 个块为 $f(\boldsymbol u_n^{k+1}-\boldsymbol b_n^k)$，$\boldsymbol g^k$ 由 $\boldsymbol b_1^k+\boldsymbol u_0,\boldsymbol b_2^k,\ldots$ 组成。内层准 Newton 为

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

精确 Jacobian 与可分离近似分别是

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

> [!warning] 原文公式核对：平均 Jacobian 的求和指标
> 期刊版把这里印成
> $J^{-1}\sum_{j=1}^J\nabla f(\boldsymbol u_n-\boldsymbol b_n)$：
> 被加项不含 $j$，且 4.5.1 的全时间系统有 $N_t$ 个块而非 $J$ 个。
> 上式按完整块对角 Jacobian 的定义修正。

预条件矩阵与 (4.16) 同构，增量仍由 (4.17) 并行求解。NKA 可以
替换平均 Jacobian，但此处不再展开。

![原论文 Figure 4.13：两种黏性 Burgers 方程上的两类 CGC](assets/papers/time-parallelization/source-figures/figure-4-13.svg)

Figure 4.13 把同一问题推进到非线性 Burgers：左右面板取
$\nu=1$ 与 $0.01$，并比较 $\alpha=0.4,0.25,0.1$ 和标准 CGC。
两种黏性下 $\alpha=0.4$ 都最慢；$\nu=1$ 时 $\alpha=0.1$
最接近标准曲线，$\nu=0.01$ 时则是 $\alpha=0.25$ 几乎重合。
这说明 $\alpha$ 的影响延续到非线性问题；Wu（2018, Section 4）
证明了“足够小的 $\alpha$ 可保持标准 CGC 的速度”，但没有给出
$\rho/(1+\rho)$ 的非线性阈值。

### Remark 4.2：MGRiT 需要一致的首尾条件

把 (4.15) 的条件直接搬到 MGRiT 会对任意 $\alpha$ 发散。收敛的一致条件应为

$$
\boldsymbol u_1^{k+1}
=\alpha(\boldsymbol u_{N_t}^{k+1}-\boldsymbol u_{N_t}^k)+\boldsymbol u_1.
$$

由此得到

$$
\left\{
\begin{aligned}
\boldsymbol u_0^{k+1}&=\boldsymbol u_0,\\
\boldsymbol u_1^{k+1}
&=\alpha(\boldsymbol u_{N_t}^{k+1}-\boldsymbol u_{N_t}^k)+\boldsymbol u_1,\\
\boldsymbol u_{n+1}^{k+1}
&=\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})
+\widetilde{\boldsymbol b}_{n+1}^k,\quad n=1,\ldots,N_t-1.
\end{aligned}
\right. \tag{4.20}
$$

这里 $\widetilde{\boldsymbol b}_{n+1}^k=\mathcal F(T_n,T_{n+1},\widetilde{\boldsymbol s}_n^k)-\mathcal G(T_n,T_{n+1},\widetilde{\boldsymbol s}_n^k)$，$\widetilde{\boldsymbol s}_n^k=\mathcal F(T_{n-1},T_n,\widetilde{\boldsymbol u}_{n-1}^k)$。注意这里的 $\widetilde{\boldsymbol u}_n^k$ 与 Section 4.5.1 不同；MGRiT 变体重新定义为

$$
\widetilde{\boldsymbol u}_n^k=
\begin{cases}
\boldsymbol u_n,&n=0,1,\\
\boldsymbol u_n^k,&n\ge2,
\end{cases}
$$

即 $n=1$ 处取的是收敛值 $\boldsymbol u_1$ 而非当前迭代 $\boldsymbol u_1^k$，这与上面首尾条件里出现的 $\boldsymbol u_1$ 一致。小 $\alpha$ 时，该变体与原 MGRiT 同速，Theorem 4.7 的阈值机制仍适用。Parareal 本身也可使用同样一致的差分首尾条件。

## 4.5.2 基于对角化的粗传播子

### 细传播与首尾耦合粗传播

第二条路线在每个大区间内使用相同的线性-$\theta$ 积分器和相同的 $\Delta t=\Delta T/J$。细传播顺序执行

$$
\boldsymbol v_{j+1}-\boldsymbol v_j
=\Delta t[\theta f(\boldsymbol v_{j+1})
+(1-\theta)f(\boldsymbol v_j)],
\quad j=0,\ldots,J-1,
\quad \boldsymbol v_0=\boldsymbol u_n. \tag{4.21}
$$

$\theta=1$ 为后向 Euler，$\theta=1/2$ 为梯形规则。这里以线性
$\theta$ 方法展示结构；$s$ 级 Runge–Kutta 的推广见 Gander 与
Wu（2020）的附录。特殊粗传播 $\mathcal F_\alpha^*$ 使用同一差分式，
改用

$$
\boldsymbol v_0=\alpha\boldsymbol v_J+(1-\alpha)\boldsymbol u_n. \tag{4.22}
$$

$J$ 个细步因此构成首尾耦合系统，可由 ParaDiag 同时求解。

### 非线性全时间系统与准 Newton

令 $\boldsymbol V=(\boldsymbol v_1^\top,\ldots,\boldsymbol v_J^\top)^\top$，则

$$
\underbrace{(C_\alpha\otimes I_x)\boldsymbol V
-\Delta tF(\boldsymbol V)}_{K(\boldsymbol V)}
=\boldsymbol b(\boldsymbol u_n), \tag{4.23}
$$

$$
\boldsymbol b(\boldsymbol u_n)
=((1-\alpha)\boldsymbol u_n^\top,0,\ldots,0)^\top.
$$

原文把 $\boldsymbol V$ 与 $\boldsymbol b(\boldsymbol u_n)$ 放在一段未编号的行内公式中，编号 (4.24) 属于随后定义 $C_\alpha$ 与 $F(\boldsymbol V)$ 的那组公式：

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

$C_\alpha$ 的右上角 $-\alpha$ 就是 (4.22) 的首尾耦合条件：若把它换成 $0$，$C_\alpha$ 退回严格下三角，全时间系统重新变成顺序细传播。$F$ 的首块之所以与其余块形状不同，也是同一原因——$\boldsymbol v_0$ 已被 $\alpha\boldsymbol v_J+(1-\alpha)\boldsymbol u_n$ 替换。注意 $\theta$ 与 $1-\theta$ 的权重已经写进 $F$ 内部，这一点在下面推导 Jacobian 时是关键。

准 Newton 更新为

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

其中

$$
\widetilde C_{\alpha,\theta}=
\begin{bmatrix}
\theta&&&(1-\theta)\alpha\\
1-\theta&\theta\\
&\ddots&\ddots\\
&&1-\theta&\theta
\end{bmatrix},
$$

令
$J_0=\nabla f(\alpha\boldsymbol v_J^l+(1-\alpha)\boldsymbol u_n)$、
$J_j=\nabla f(\boldsymbol v_j^l)$。由 (4.24) 中已经带权的
$F(\boldsymbol V)$ 直接求导，得到

$$
\mathcal J_F(\boldsymbol V^l)=
\begin{bmatrix}
\theta J_1&&&\alpha(1-\theta)J_0\\
(1-\theta)J_1&\theta J_2\\
&\ddots&\ddots\\
&&(1-\theta)J_{J-1}&\theta J_J
\end{bmatrix}.
$$

因此

$$
\nabla K=C_\alpha\otimes I_x-\Delta t\,\mathcal J_F.
$$

> [!warning] 原文公式核对：$\nabla F$ 的记号重载
> 正式版写成
> $(\widetilde C_{\alpha,\theta}\otimes I_x)\nabla F$，但 (4.24)
> 已把 $\theta$、$1-\theta$ 和 head–tail 权重放进 $F$。
> 若 $\nabla F$ 仍指这个带权向量的 Jacobian，再乘一次
> $\widetilde C_{\alpha,\theta}$ 就会重复计权。上式直接对 (4.24)
> 求导，避免该歧义。

而可分离近似使用特殊平均

$$
\overline{\nabla f}(\boldsymbol V^l)
=\frac1J\left[
\sum_{j=1}^{J-1}\nabla f(\boldsymbol v_j^l)
+\nabla f\!\left(
\alpha\boldsymbol v_J^l+(1-\alpha)\boldsymbol u_n
\right)
\right].
$$

最后一块取 head–tail 状态的 Jacobian，并不是
$\nabla f(\boldsymbol v_J^l)$。$C_\alpha$ 与
$\widetilde C_{\alpha,\theta}$ 同时对角化，内层空间系统因而可以并行。

外层 Parareal 仍写成

$$
\boldsymbol u_{n+1}^{k+1}
=\mathcal F_\alpha^*(T_n,T_{n+1},\boldsymbol u_n^{k+1})
+\mathcal F(T_n,T_{n+1},\boldsymbol u_n^k)
-\mathcal F_\alpha^*(T_n,T_{n+1},\boldsymbol u_n^k). \tag{4.26}
$$

### 线性系统、并行性与极限

$f(\boldsymbol u)=A\boldsymbol u$ 时，

$$
(C_\alpha\otimes I_x
-\widetilde C_{\alpha,\theta}\otimes\Delta tA)\boldsymbol V
=\boldsymbol b(\boldsymbol u_n), \tag{4.27}
$$

$$
\boldsymbol b(\boldsymbol u_n)
=([(I_x+\Delta t(1-\theta)A)(1-\alpha)\boldsymbol u_n]^\top,0,\ldots,0)^\top.
$$

这里统一使用 $\widetilde C_{\alpha,\theta}$；正式版在
$\widetilde C_{\theta,\alpha}$ 与 $\widetilde C_{\alpha,\theta}$
之间交替，但二者指同一个矩阵。令
$H_J=(0,\ldots,0,1)\in\mathbb R^{1\times J}$，粗传播取
$\mathcal F_\alpha^*=(H_J\otimes I_x)\boldsymbol V=\boldsymbol v_J$。
它等价于

$$
\left\{
\begin{aligned}
\boldsymbol v_{j+1}-\boldsymbol v_j
&=\Delta tA[\theta\boldsymbol v_{j+1}+(1-\theta)\boldsymbol v_j],\\
\boldsymbol v_0&=\alpha\boldsymbol v_J+(1-\alpha)\boldsymbol u_n^k.
\end{aligned}
\right. \tag{4.28}
$$

$\alpha=0$ 时粗传播等于顺序细传播，外层一轮收敛，同时完全失去加速。$0<\alpha<1$ 时，$J$ 个点一次对角化求解；若空间系统有足够资源，墙钟成本约为顺序细传播的 $1/J$。

### Theorem 4.8：抛物谱与双曲谱

Theorem 4.8 引自 Gander 与 Wu（2020）。对线性初值问题
$\boldsymbol u'=A\boldsymbol u+\boldsymbol g$、
$\boldsymbol u(0)=\boldsymbol u_0$，
$A\in\mathbb C^{N_x\times N_x}$，若 $\mathcal F$ 与
$\mathcal F_\alpha^*$ 都使用稳定单步 Runge–Kutta 方法，并令
$\{\boldsymbol u_n^k\}$ 为 (4.26) 的第 $k$ 轮迭代，

$$
e^k=\max_{1\le n\le N_t}\|\boldsymbol u_n-\boldsymbol u_n^k\|_\infty,
$$

其中 $\{\boldsymbol u_n\}$ 是收敛解（不是 PDE 精确解），

则

$$
e^k\le\rho^ke^0,
\qquad
\rho=
\begin{cases}
\alpha,&\sigma(A)\subset\mathbb R_-,\\[4pt]
\dfrac{2\alpha N_t}{1+\alpha},&\sigma(A)\subset i\mathbb R.
\end{cases} \tag{4.29}
$$

热方程因子与粗区间数无关；纯虚谱上界随 $N_t$ 线性增长。后者是上界，$\alpha N_t$ 较大时可能很松。

![原论文 Figure 4.14：热方程上 rho=alpha 的锐利预测](assets/papers/time-parallelization/source-figures/figure-4-14.svg)

热方程采用齐次 Dirichlet 边界、$u_0=\sin^2(2\pi x)$、梯形规则、$\Delta T=1/12$、$J=10$、$\Delta x=1/100$。左、右面板分别取 $N_t=36$ 与 $72$，每幅都比较 $\alpha=10^{-1},10^{-2},10^{-3}$。实测虚线与理论点线几乎平行，$N_t$ 加倍没有改变由 $\rho=\alpha$ 决定的斜率。

![原论文 Figure 4.15：波动方程上 alpha 与粗区间数的共同影响](assets/papers/time-parallelization/source-figures/figure-4-15.svg)

波动方程先化为一阶系统 $\boldsymbol w'=\boldsymbol{Aw}$，其中 $\boldsymbol w=(\boldsymbol u^\top,(\boldsymbol u')^\top)^\top$、$\boldsymbol A=\begin{bmatrix}0&I_x\\A&0\end{bmatrix}$，由此得到 $\sigma(\boldsymbol A)\subset i\mathbb R$，正是 Theorem 4.8 的第二个分支。实验取周期边界、$\boldsymbol w(0)=(\sin^2(2\pi\boldsymbol x_h)^\top,\boldsymbol 0^\top)^\top$。(a) 固定 $\alpha=0.01$，比较 $N_t=24,48,96$，区间数增加会明显减慢；(b) 固定 $\alpha=10^{-4}$，比较 $N_t=24,48,96,960$，从 24 增到 960 只多约两轮达到 $\max\{\Delta t^2,\Delta x^2\}$。两个面板把 $\alpha N_t$ 的联合作用直接分离出来。

![原论文 Figure 4.16：小 alpha Nt 时理论因子较锐利，大乘积时出现超线性](assets/papers/time-parallelization/source-figures/figure-4-16.svg)

Figure 4.16(a) 固定 $N_t=24$，比较 $\alpha=10^{-2}$ 与 $10^{-4}$；(b) 固定 $\alpha=10^{-4}$，比较 $N_t=24$ 与 $960$。只有 $\alpha=10^{-4},N_t=24$ 的小乘积组合紧贴 (4.29) 的点线上界；另两组实测曲线出现超线性下降，线性上界明显偏保守。

![原论文 Figure 4.17：Burgers 方程达到 1e-8 所需迭代数](assets/papers/time-parallelization/source-figures/figure-4-17.svg)

Burgers 实验取周期边界、$u_0=\sin^2(2\pi x)$、$\Delta T=0.1$、$J=10$、$\Delta x=1/100$，三条曲线对应 $\nu=1,0.01,10^{-4}$。(a) 固定 $N_t=40$，小 $\alpha$ 加快收敛并削弱黏性的影响；(b) 固定 $\alpha=10^{-3}$，$N_t=10$ 到 $160$ 时迭代数只在 2–5 轮之间变化。非线性理论在精确求解 (4.23) 和 Lipschitz 条件下给出 $\rho=O(\alpha)$。

## 两种路线的最终对照

| 问题         | 对角化 CGC (4.15)        | 对角化粗传播 (4.26)                                  |
| ------------ | ------------------------ | ---------------------------------------------------- |
| 对角化方向   | 全局 $N_t$ 个粗点        | 每个粗区间的 $J$ 个细点                              |
| 改动位置     | CGC                      | 粗传播子                                             |
| 粗细积分器   | 可不同                   | 相同积分器与步长                                     |
| 主要适用范围 | 抛物问题                 | 抛物与双曲问题                                       |
| 关键参数     | $\alpha\le\rho/(1+\rho)$ | 抛物因子 $\alpha$，双曲上界 $2\alpha N_t/(1+\alpha)$ |

## 公式、定理与图表覆盖核对

| 原文项目                               | 论文小节 | 覆盖状态                                                 |
| -------------------------------------- | -------- | -------------------------------------------------------- |
| (4.14)–(4.17)                          | 4.5.1    | 标准/首尾 CGC、线性全时间矩阵、三步并行解                |
| Theorem 4.7, Figure 4.12               | 4.5.1    | $\alpha$ 阈值、舍入折中、热与 ADE 实验                   |
| (4.18)–(4.19), Figure 4.13             | 4.5.1    | 非线性系统、平均 Jacobian 准 Newton、Burgers 实验        |
| Remark 4.2, (4.20)                     | 4.5.1    | MGRiT 的一致首尾条件及收敛变体                           |
| (4.21)–(4.26)                          | 4.5.2    | 同积分器细/粗传播、非线性全时间系统、准 Newton、外层更新 |
| (4.27)–(4.28)                          | 4.5.2    | 线性化、终点提取、$\alpha=0$ 极限与 $J$ 路并行           |
| Theorem 4.8, (4.29), Figures 4.14–4.17 | 4.5.2    | 负实/纯虚谱界、热、波动与 Burgers 全部原图               |

## 本页原文

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Section 4.5, pp. 460–471.
