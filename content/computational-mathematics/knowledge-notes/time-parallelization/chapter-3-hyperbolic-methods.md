---
title: 第三章：对双曲问题有效的 PinT 方法
description: Schwarz 波形松弛、并行延迟校正、ParaExp 与 ParaDiag 的完整推导和数值解释
lang: zh
translation: en/computational-mathematics/knowledge-notes/time-parallelization/chapter-3-hyperbolic-methods
tags:
  - 时间并行
  - 双曲-PDE
---

> [!note] 阅读范围
> 本章对应论文 Section 3（pp. 396–443）。正文严格保留原论文层级：Section 3 导论、3.1 历史发展、3.2 SWR、3.3 IDC、3.4 ParaExp，以及 3.5.1/3.5.2 两类 ParaDiag。公式、定理和论文实验按论证顺序解释；Python 结果、参数比较与覆盖审计另标为本站补充，不占用论文编号。

## 论文与精读页的对应关系

本页保留章节全景与本站复现实验。逐公式、逐定理和逐图表的完整推导拆在以下页面，便于按论文顺序阅读：

| 论文小节           | 原文页码    | 精读页                                                                                                                      | 覆盖范围                                                          |
| ------------------ | ----------- | --------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------- |
| Section 3、3.1–3.2 | pp. 396–405 | [[computational-mathematics/knowledge-notes/time-parallelization/chapter-3-1-history-and-swr\|历史发展与 Schwarz 波形松弛]] | (3.1)–(3.4), Theorems 3.1–3.2, Figures 3.1–3.3                    |
| 3.3                | pp. 405–411 | [[computational-mathematics/knowledge-notes/time-parallelization/chapter-3-2-idc\|并行积分延迟校正]]                        | (3.5)–(3.12), Theorem 3.3, Figures 3.4–3.6                        |
| 3.4                | pp. 412–415 | [[computational-mathematics/knowledge-notes/time-parallelization/chapter-3-3-paraexp\|ParaExp]]                             | (3.13)–(3.21), Theorem 3.4, Figures 3.7–3.8                       |
| 3.5、3.5.1         | pp. 415–430 | [[computational-mathematics/knowledge-notes/time-parallelization/chapter-3-4-paradiag-i\|直接 ParaDiag]]                    | (3.22)–(3.48), Theorems 3.5–3.7, Figures 3.9–3.14, Tables 3.1–3.2 |
| 3.5.2              | pp. 431–442 | [[computational-mathematics/knowledge-notes/time-parallelization/chapter-3-5-paradiag-ii\|迭代 ParaDiag]]                   | (3.49)–(3.68), Theorems 3.8–3.9, Figures 3.15–3.18                |

## Section 3 导论与 3.1 历史发展

### 这一组方法为何能处理长程传播

双曲方程把精细结构沿特征线带到很远的时间。有效算法必须保留传播
路径与相位，或直接求解跨越整个时间区间的耦合。下面四条路线分别从
界面传播、校正流水线、指数传播和全时间代数入手：SWR、并行 IDC、
ParaExp 与 ParaDiag。

这些方法的来源不同。SWR 继承区域分解与波形松弛；IDC 来自缺陷校正和高阶积分；ParaExp 利用线性系统的指数传播；ParaDiag 利用全时间矩阵的可对角化结构。它们大多也能处理抛物问题。共同点在于并行机制不依赖强耗散型粗传播子。

历史上，SWR 与 mapped tent pitching、unmapped tent pitching 有紧密联系。IDC 的并行版本包括按时间窗排布的 PIDC 和按校正层流水化的 RIDC。ParaExp 属于线性直接方法，随后扩展出非线性迭代形式。ParaDiag 则同时发展出直接对角化、波形松弛、定常迭代和 Krylov 预条件等形态。

## 3.2 Schwarz 波形松弛

### 从经典波形松弛到时空区域分解

对线性系统 $\boldsymbol u'=A\boldsymbol u+\boldsymbol g$，经典波形松弛先作矩阵分裂 $A=M+N$，再迭代求解

$$
\boldsymbol u^{k+1\prime}-M\boldsymbol u^{k+1}
=N\boldsymbol u^k+\boldsymbol g.
$$

若 $M$ 由块对角或着色结构组成，多个分量可以并行计算。这个办法直接依赖代数分裂；不合适的 $M$ 会导致很慢的收敛，甚至发散。经典空间区域分解通常在每个时间步求解一个椭圆问题，也会把所有子域绑定到相同的时间离散。

SWR 先在连续空间区域上分解，再让每个子域一次求解完整时间窗。相邻子域在重叠区或界面交换整条时间波形。各子域可以采用适合局部方程和网格的离散，传输条件也能根据 PDE 的传播机制设计。

![Schwarz 波形松弛在完整时间窗上交换界面波形](assets/diagrams/pint/zh/schwarz-waveform-relaxation.svg)

这种方法的迭代对象是界面函数。Dirichlet 条件只交换解值；Robin、Ventcel 或卷积条件还近似法向通量和更完整的 Dirichlet-to-Neumann 映射。带优化传输条件的版本通常称为 OSWR。

### 3.2.1 一阶抛物问题：重叠宽度与 Robin 参数

以一维对流扩散为例，把区域分成两个重叠子域，重叠宽度记为 $l$。
每次迭代在两个子域上并行求解完整时间窗，并在人工边界施加 Robin 条件

$$
(\partial_x+p)u_1^{k+1}=(\partial_x+p)u_2^k,
\qquad
(\partial_x-p)u_2^{k+1}=(\partial_x-p)u_1^k. \tag{3.1}
$$

$p>0$ 控制传输条件；$p\to\infty$ 对应 Dirichlet 交换。误差分析在 Laplace/Fourier 频域内进行。Theorem 3.1 给出优化 Robin 参数及连续问题的最坏频率收敛估计。Dirichlet 情形满足形如

$$
\rho_{\mathrm D}\leq
\exp\!\left(-\frac{l\pi}{\nu T}\right)
$$

的界。它清楚展示三个趋势：增加重叠会加快收敛；增大时间窗会让跨窗耦合更难；$\nu$ 减小时，定向传播使信息更快穿过重叠区，当前模型中的 SWR 反而可以加速。

Figure 3.1 取 $L=8.2$、$T=5$、$\Delta t=0.01$、$\Delta x=0.02$、$l=2\Delta x$ 和 Gaussian 初值。数值曲线验证了黏性下降时收敛改善，也显示 Robin 条件明显优于 Dirichlet。四子域实验中，$\nu=0.1$ 时 Dirichlet 与优化 Robin 分别约需 92 和 28 次迭代；二子域连续理论预测约为 32 和 4 次。差距来自连续/离散、二子域/四子域以及边界设置的不同，不能把理论界直接当作该离散实验的精确迭代数。

更高阶的 Ventcel 条件可以进一步改善渐近收敛。若使用时间卷积近似精确透明边界，理论上可以获得与网格无关的收敛因子；代价是界面算子更复杂，并带有时间非局部性。

### 3.2.2 二阶双曲问题：有限步传播

对波动方程，SWR 的性质更加直接。每轮迭代都会把正确的界面信息向相邻子域推进一个有限传播距离。Theorem 3.2 表明，在两个重叠子域上使用 Dirichlet 传输时，只要

$$
k>\frac{Tc}{(\beta-\alpha)L},
$$

第 $k$ 轮之后整个时间窗内的界面误差就会消失。实际重叠宽度是
$(\beta-\alpha)L$，$c$ 是波速；期刊版漏掉了分母中的 $L$，在
$L=1$ 时才与上式相同。结论来自有限传播速度：一轮迭代只能把正确性
扩展到一个特征锥覆盖的区域，足够多轮后特征锥覆盖完整时空域。

Figure 3.2 用特征锥展示这一过程。每个子域中已有一部分解与精确解一致，下一轮通过界面数据把这片正确区域继续扩大。这个几何解释也引出红黑 SWR：相邻时空块按颜色并行计算，允许一定冗余工作，以换取更大的并发度。

### Tent pitching、MTP 与 UTP

Tent pitching 按有限传播速度构造倾斜的时空单元。每个 tent 的底边已有数据后，内部计算可以独立进行。Mapped tent pitching（MTP）把倾斜 tent 映射到规则柱体，便于使用标准求解器；映射会增加实现成本，也可能造成阶数下降。

Unmapped tent pitching（UTP）保留原始时空几何，可以理解为红黑 SWR 或全时间系统上的限制加性 Schwarz。残差决定某个 tent 还能向上推进多远，省去显式映射与相关阶数损失。抛物方程具有无限传播速度，无法形成严格独立的 tent；当扩散很小，SWR/UTP 仍可通过少量迭代校正跨 tent 影响。

## 3.3 时间并行 IDC

### IDC 的残差误差方程

考虑初值问题

$$
u'(t)=f(u(t),t),\qquad u(0)=u_0.
$$

对当前近似 $u^k(t)$，积分形式的缺陷或残差可写为

$$
r^k(t)=u_0+\int_0^t f(u^k(s),s)\,ds-u^k(t). \tag{3.6}
$$

令真实误差 $e^k=u-u^k$。把 $u=u^k+e^k$ 代回积分方程，可以得到误差的积分形式；对时间求导后得到

$$
e^{k\prime}(t)
=f(u^k(t)+e^k(t),t)-f(u^k(t),t)+r^{k\prime}(t). \tag{3.8}
$$

IDC 先用低阶方法得到预测解，再离散误差方程并更新
$u^{k+1}=u^k+e^k$。用 $\theta$ 方法离散局部动力学差，再以插值求积
注入积分残差，就得到节点递推 (3.11)，无需直接数值微分噪声残差。

![IDC、PIDC 与 RIDC 的校正和流水线关系](assets/diagrams/pint/zh/idc-pipeline.svg)

Theorem 3.3 说明：基础积分器若为 $p$ 阶，在 $M$ 个均匀节点上完成 $k$ 次校正后，整体阶数达到

$$
\mathcal O\!\left(\Delta t^{\min\{M,(k+1)p\}}\right).
$$

校正阶数最终受节点数限制。使用 Gauss–Lobatto 节点的谱延迟校正（SDC）在 $J$ 个节点上可达到 $2J-1$ 阶，第四章的 PFASST 会把 SDC 放入多层时间并行结构。

### 为什么普通 IDC 仍然串行

长时间区间通常被切成多个 window。一个 window 内完成预测和若干校正后，末值成为下一个 window 的初值。普通 IDC 仍需顺序处理 window，校正层内部也有节点递推，因而高阶精度本身没有自动带来时间并行。

### 3.3.1 流水线 IDC（PIDC）

PIDC 让不同 window 同时执行不同编号的校正 sweep。流水线填满后，预测层处理较晚 window，第一校正层处理前一个 window，更高校正层继续落后。并发度大致由校正层数决定，启动和排空阶段会降低短作业的效率。

这个调度存在一个精度风险。较晚 window 的初值会随着前一 window 的高阶校正不断改变。较高校正层可能从一个尚不光滑、尚未稳定的初值启动。IDC 的升阶理论需要足够的时间正则性；初值不规则、扩散太弱或解含高频结构时，预期高阶会丢失。

Figure 3.5 使用周期对流扩散方程，$\Delta x=1/64$、$T=3$、window 长度 $\Delta T=0.1$、每窗 $M=5$ 个节点和后向 Euler 基础积分器。$\sigma=1000$ 的窄源产生低正则性，$\sigma=5$ 给出较平滑输入。结果可分成三层：低正则性时多次 IDC/PIDC 校正无法稳定实现高阶；平滑且扩散较强时校正明显有效；平滑但扩散很小时，长寿命高频仍会破坏理想升阶。

### 3.3.2 修正型 IDC（RIDC）

RIDC 为每个校正层保留一个滑动的 $M$ 节点窗口。各层在连续时间步上像装配线一样推进，不必等待完整 window 结束。它减少全局同步，适合在多个核上持续运行不同校正层。

RIDC 改善调度，没有取消正则性条件。Figure 3.6 再次表明，低正则性和弱扩散会限制校正层带来的精度提升。对双曲问题，保留高频既是物理优势，也是 IDC 升阶的数值风险。

## 3.4 ParaExp

### 线性问题的精确分解

ParaExp 针对

$$
\boldsymbol u'(t)=A\boldsymbol u(t)+\boldsymbol g(t)
$$

构造直接时间并行。把 $[0,T]$ 分成 $N$ 个子区间 $[T_{n-1},T_n]$。第一步在所有子区间并行求解零初值非齐次问题

$$
\boldsymbol v_n'(t)=A\boldsymbol v_n(t)+\boldsymbol g(t),
\qquad \boldsymbol v_n(T_{n-1})=0. \tag{3.13}
$$

第二步把每个界面上的贡献通过齐次方程向后续时间传播：

$$
\boldsymbol w_n'(t)=A\boldsymbol w_n(t),
\quad t\in(T_{n-1},T],
\qquad \boldsymbol w_n(T_{n-1})=\boldsymbol v_{n-1}(T_{n-1}), \tag{3.14}
$$

其中 $\boldsymbol v_0(T_0)=\boldsymbol u_0$。线性叠加给出精确重构

$$
\boldsymbol u(t)=\boldsymbol v_n(t)+
\sum_{j=1}^{n}\boldsymbol w_j(t),
\qquad t\in[T_{n-1},T_n]. \tag{3.15}
$$

![ParaExp 将局部受迫响应与全局齐次传播分开](assets/diagrams/pint/zh/paraexp-decomposition.svg)

关键点在齐次尾部：

$$
\boldsymbol w_n(t)=e^{(t-T_{n-1})A}\boldsymbol v_{n-1}(T_{n-1}).
$$

矩阵指数作用可以直接跳到任意后续时刻，成本取决于矩阵与目标精度，
通常不随中间步数线性增长。大型稀疏系统常用有理 Krylov 或
Chebyshev 近似，小型稠密矩阵可用 scaling-and-squaring 与 Padé。
Gander 与 Güttel（2013）的波动算例达到约 80% 并行效率，但这一数字
依赖具体指数算法、分区和硬件。

### 非线性扩展及其边界

对

$$
\boldsymbol u'=A\boldsymbol u+B(\boldsymbol u)+\boldsymbol g, \tag{3.17}
$$

非线性版本仍把线性齐次部分交给指数传播，把
$B(\boldsymbol u)$ 留在迭代的局部问题中。Theorem 3.4 给出两点：
第 $k$ 轮后 $[0,T_k]$ 已精确；在粗节点上，该迭代等价于简化
Parareal，粗传播只含线性齐次演化，细传播处理完整非线性。

Figure 3.8 比较 Burgers 方程上的 ParaExp 与标准 Parareal，空间步长 $0.01$、$T=2$，细步长为 $0.01/20$，标准 Parareal 粗步长为 $0.01$。黏性较大时，线性扩散是主要动力学，ParaExp 的粗模型很有效；黏性下降后，非线性输运占比增大，标准 Parareal 一度更快；$\nu=0.02$ 时 ParaExp 发散。进一步靠近双曲极限时，标准 Parareal 也会失效。ParaExp 的线性构造很强，非线性拆分质量决定了扩展版本的适用范围。

## 3.5 ParaDiag：沿时间对角化

### 两条 ParaDiag 路线

ParaDiag 把全时间耦合系统变成多个独立空间系统，分为两类：

1. **ParaDiag-I** 精确对角化经过特殊设计的时间离散，形成直接求解器；
2. **ParaDiag-II** 用循环或 $\alpha$-循环时间矩阵近似原矩阵，再作定常迭代或 Krylov 预条件。

二者都需要时间特征向量矩阵条件良好，并且每个复移位空间系统能够高效求解。它们的三步结构一致：时间方向变换、并行移位求解、逆变换。

![ParaDiag 的时间变换、独立空间求解与逆变换](assets/diagrams/pint/zh/paradiag-three-stage.svg)

### 3.5.1 直接 ParaDiag 方法（ParaDiag-I）

#### 几何变步长的后向 Euler

对线性系统 (2.1)，变步长后向 Euler 的全时间矩阵具有

$$
K=B\otimes I_x-I_t\otimes A. \tag{3.23}
$$

若时间步按 $\Delta t_n=\mu^{n-1}\Delta t_1$ 几何增长，$B$ 可以写成 $B=VDV^{-1}$。于是

$$
K^{-1}
=(V\otimes I_x)
(D\otimes I_x-I_t\otimes A)^{-1}
(V^{-1}\otimes I_x), \tag{3.25}
$$

中间块对角系统包含 $N_t$ 个互不依赖的移位空间问题。

参数 $\mu=1+\rho$ 暴露出直接法的核心矛盾。较大的 $\rho$ 让时间步变化明显，截断误差约为 $\mathcal O(\rho^2)$；$\rho$ 太小时，$B$ 接近不可对角化的 Jordan 结构，舍入误差按 $\epsilon\rho^{-(N_t-1)}$ 放大。Theorem 3.5 给出两者的平衡以及最优 $\rho$ 的尺度。Figures 3.9–3.10 验证了 U 形误差曲线，也显示双精度下可用时间步数十分有限。直接增加 $N_t$ 会让条件数和舍入误差迅速占据主导。

#### 波动方程与梯形规则

波动方程先转成一阶系统，再用变步长梯形规则保持能量性质。对应全时间系统同样可对角化。Theorem 3.6 仍得到截断误差与 $\epsilon\rho^{-(N_t-1)}$ 舍入放大的竞争。Figure 3.11 和 Table 3.1 表明，随着 $N_t$ 增大，特征向量条件数迅速增长，误差在大约 $N_t>32$ 后开始上升。

#### 边值方法缓解病态性

边值方法（BVM）在前 $N_t-1$ 个节点使用中心差分，最后一个节点用后向 Euler 封闭全时间系统：

$$
\frac{\boldsymbol u_{n+1}-\boldsymbol u_{n-1}}{2\Delta t}
=A\boldsymbol u_n+\boldsymbol g_n,
$$

配合终端离散形成可对角化矩阵。尽管末端公式为一阶，整体仍可达到二阶。Theorem 3.7 给出 $\operatorname{Cond}(V)=\mathcal O(N_t^2)$，比几何变步长方案稳定得多。Figure 3.12 展示了随均匀时间步缩小的二阶误差，未出现前述快速恶化。对二阶波动方程还可直接构造含 $B^2\otimes I_x-I_t\otimes A$ 的系统，避免一阶化造成存储量翻倍。

#### 非线性 ParaDiag-I

非线性全时间系统先做 Newton 线性化。各时间点的 Jacobian 不同，
单一 Kronecker 结构随之消失；用平均 Jacobian 近似全部时间块后，
时间变换仍能把移位 Jacobian 问题并行化。若长窗口内的轨道变化过大，
应缩短窗口并顺序处理多个 window。

Figure 3.13 与 Table 3.2 的 Burgers 实验显示：$\nu=0.1$ 时并行 Jacobian 求解数量远少于顺序 Newton；黏性下降后，收敛对总时间 $T$ 越来越敏感；$\nu=0.002$ 且时间窗较长时方法失败。

近似 Jacobian 还可通过最近 Kronecker 乘积问题选择：

$$
\min_{\Phi_k\ \mathrm{diagonal}}
\left\|
\nabla F(\boldsymbol U^k)-\Phi_k\otimes A_k
\right\|. \tag{3.47}
$$

其中 $\Phi_k$ 是对角时间缩放。NKA 可在粗空间模型上离线计算
$\Phi_k$，保留 Jacobian 随时间变化的幅值。Figure 3.14 显示它在
$T=1.3$ 等较长窗口上明显改善准 Newton 收敛。

### 3.5.2 迭代 ParaDiag 方法（ParaDiag-II）

#### Strang 循环预条件器

线性多步法的全时间系统写成

$$
K=B_1\otimes I_x-B_2\otimes\Delta t A.
$$

ParaDiag-II 用 Strang 循环矩阵 $C_1,C_2$ 替代 Toeplitz 型 $B_1,B_2$，构造

$$
P=C_1\otimes I_x-C_2\otimes\Delta t A.
$$

循环矩阵由离散 Fourier 矩阵同时对角化，因此 $P^{-1}$ 仍通过“FFT、独立移位空间求解、逆 FFT”实现。若 $P$ 很接近 $K$，可直接使用定常迭代；收敛因子较差时，$P$ 更适合作 GMRES 预条件器。即使定常迭代的谱半径超过 1，Krylov 方法仍可能利用聚集谱快速收敛。

Theorem 3.8 对对称负定 $A$ 给出一个结构性结果：预条件矩阵只有有限组非单位特征值，GMRES 在精确算术下有有限步上界。该上界随空间维数增长，不能单独保证大问题上的快速收敛。

Figure 3.15 采用 $T=2$、$\Delta t=1/50$、$\Delta x=1/100$。循环预条件器对热方程和 $\nu=10^{-3}$ 的对流扩散方程非常有效；$\nu$ 继续下降时谱聚集变差；波动方程的非单位特征值沿复平面展开，外层迭代明显增加。这组实验直接连接第二章的结论：耗散让循环闭合造成的首尾差异局部化，持续传播则把这项差异带遍整个时间域。

#### $\alpha$-循环与波形松弛解释

连续时间的 head-tail 波形松弛也会导出同一个 $\alpha$-循环矩阵。
这里耦合的是**同一时间窗第 $k$ 轮迭代**的初值、本轮末值和上一轮末值，
不是另开一个新时间窗。离散后的特征向量矩阵具有
$\Lambda_\alpha F^*$ 形式，仍可使用 FFT。

$\alpha=1$ 对应标准循环近似。周期传播问题在这个取值下可能出现奇异性；$0<\alpha<1$ 打破首尾完全周期闭合，并常常显著改善对流和波动问题。Figure 3.16 展示了无需 Krylov 加速的强收敛。

Theorem 3.9 针对稳定单步法和对称两步法给出谱界：

$$
\frac{1}{1+\alpha}
\le |\lambda(P_\alpha^{-1}K)|
\le \frac{1}{1-\alpha},
\qquad
\rho(I-P_\alpha^{-1}K)
\le \frac{\alpha}{1-\alpha}.
$$

这组界与空间网格和时间步数无关，但要求时间积分器满足相应稳定性。（原文印刷版把两端写反成 $1/(1-\alpha)\le|\lambda|\le1/(1+\alpha)$，对 $0<\alpha<1$ 是空区间；精读页给出了核对与推导。）Numerov 参数的临界实验说明，$\gamma=1/120$ 满足条件，略微越过稳定界的 $1/120.01$ 就会破坏预测。

减小 $\alpha$ 会改善迭代因子，同时使时间变换的条件数增至
$\alpha^{-(N_t-1)/N_t}\le1/\alpha$；因此
$O(\epsilon/\alpha)$ 是方便的保守舍入上界。Figure 3.18 展示了这条
折中。采用误差方程更新可以避免反复把大解向量带入病态变换，从而降低
舍入污染。更一般的多步 Volterra 分析也给出特征值偏离 $1$ 为
$\mathcal O(\alpha)$ 的结论。

#### 非线性 ParaDiag-II

非线性问题使用 Newton–Krylov。每次 Newton 线性化产生全时间 Jacobian，GMRES 再用平均 Jacobian 构造的 $P_\alpha$ 预处理。定常迭代在非线性长窗中经常失效，Krylov 方法仍能利用特征值聚集。缩短时间窗会让 Jacobian 更接近，NKA 也能保留更多时间变化，两者都能改善预条件质量。

## 本站数值补充：Figure 3.15 的 Python 验证

### 基线实验

重新计算的基线使用 $N_x=N_t=100$、$T=2$、$\nu=10^{-6}$、$\alpha=1$ 和 GMRES 容差 $10^{-12}$。算法在 13 次迭代后收敛，真实相对残差为 $1.152\times10^{-14}$。

![对流扩散方程 ParaDiag-II 的 GMRES 收敛基线](assets/pint/paradiag-baseline.svg)

### 论文网格验证

| 问题                    |                         Python 结果 | 解释                                               |
| ----------------------- | ----------------------------------: | -------------------------------------------------- |
| 热方程                  |                    2 次 Krylov 更新 | 特征值紧密聚集在 1 附近                            |
| 对流扩散，$\nu=10^{-3}$ |                                3 次 | 循环预条件器与原全时间系统接近                     |
| 对流扩散，$\nu=10^{-6}$ |                               13 次 | 弱扩散使首尾闭合误差长期传播                       |
| 波动方程                | 第 89 次时预条件残差低于 $10^{-11}$ | 非单位特征值沿 $\operatorname{Re}\lambda=0.5$ 展开 |

![热方程、对流扩散方程与波动方程的 ParaDiag-II 谱和 GMRES 收敛](assets/pint/paradiag-figure-3-15.svg)

ADE 的 3 次和 13 次与 Figure 3.15(c,d) 一致。波动实验使用论文的 $\gamma=1/100$、$\alpha=1$；单线程 SciPy 在第 89 次达到预条件残差 $10^{-11}$，与论文曲线约 88 次的终点接近。若按 SciPy 的真实相对残差 $10^{-12}$ 停止，则需 103 次。MATLAB 与 SciPy 的残差归一化、restart 和停止规则不同，因此应比较谱形态和收敛阶段，不宜只比较一个迭代计数。

## 本站方法比较：参数、适用性与实现代价

| 方法        | 关键参数或选择                      | 决定的性质           | 主要风险                            |
| ----------- | ----------------------------------- | -------------------- | ----------------------------------- |
| SWR         | 重叠宽度、时间窗、传输条件          | 特征信息跨子域的速度 | 参数不合适、界面算子过贵            |
| PIDC/RIDC   | 节点数、校正层数、窗口与流水线深度  | 形式阶数和并发度     | 低正则性导致升阶失败                |
| ParaExp     | 指数作用算法、线性/非线性拆分       | 齐次尾部传播成本     | 矩阵指数昂贵、非线性拆分失真        |
| ParaDiag-I  | 时间离散、$N_t$、特征向量条件数     | 直接并发规模         | 截断与舍入误差冲突                  |
| ParaDiag-II | $\alpha$、外层 Krylov、移位求解容差 | 预条件聚集与稳定性   | $\alpha$ 过小放大舍入，过大削弱近似 |

公开的 [ActaPinT-Python](https://github.com/freezeng123456/ActaPinT-Python) 最新提交已经把 Section 3 的 15 张数值图、3 张计算式示意图和 Tables 3.1–3.2 全部接入独立入口，并保存 SVG、PNG 与 JSON。SWR、PIDC/RIDC、ParaExp、直接 ParaDiag、BVM/NKA、ParaDiag II 与波动区域分解不再只是迁移清单。需要保留的例外是 Figure 3.14：平均 Jacobian 曲线可匹配，$\nu=0.02$ 时论文给出的两条 nearest-Kronecker 加速曲线仍无法由所述权重恢复。逐项证据和当前测试状态见 [[computational-mathematics/knowledge-notes/time-parallelization/reproduction-audit-2026-08-24|全文实验复现审计]]。

## 原文覆盖核对

| 原文位置                        | 本页对应 | 已覆盖内容                                                                                                        |
| ------------------------------- | -------- | ----------------------------------------------------------------------------------------------------------------- |
| Section 3 与 3.1，pp. 396–398   | 3.1      | 双曲有效方法的范围、历史来源及 tent pitching 线索                                                                 |
| Section 3.2，pp. 398–405        | 3.2      | WR 与 SWR 差异、Robin/Dirichlet/Ventcel/卷积条件、Theorems 3.1–3.2、Figures 3.1–3.3、MTP/UTP                      |
| Section 3.3，pp. 405–411        | 3.3      | IDC 积分残差与递推、Theorem 3.3、SDC、PIDC/RIDC 调度、Figures 3.4–3.6、正则性限制                                 |
| Section 3.4，pp. 412–415        | 3.4      | ParaExp 线性分解与指数作用、非线性迭代、Theorem 3.4、Figures 3.7–3.8                                              |
| Sections 3.5–3.5.1，pp. 415–430 | 3.5.1    | ParaDiag-I 三步法、几何步长、Theorems 3.5–3.7、波动/BVM、非线性准 Newton 与 NKA、Figures 3.9–3.14、Tables 3.1–3.2 |
| Section 3.5.2，pp. 431–442      | 3.5.2    | Strang 与 $\alpha$-循环预条件、Theorems 3.8–3.9、Figures 3.15–3.18、稳定性/舍入折中、非线性 Newton–Krylov         |

## 小结

本章的四类方法以不同方式保留长程信息。SWR 沿界面和特征传播波形；IDC 把高阶误差校正排成流水线；ParaExp 用矩阵指数直接传播线性齐次响应；ParaDiag 把全时间耦合变换成并行空间求解。它们都能避开强耗散粗传播子的单一依赖，同时各有约束：SWR 需要合适传输条件，IDC 需要时间正则性，ParaExp 依赖指数作用与线性拆分，ParaDiag 需要稳定的时间对角化和可扩展的移位求解器。

## 本章原文

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Section 3, pp. 396–443。
