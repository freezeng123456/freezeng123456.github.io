---
title: 第四章：为抛物问题设计的 PinT 方法
description: Parareal、PFASST、MGRiT、对角化 Parareal 与时空多重网格的完整分析
lang: zh
translation: en/computational-mathematics/knowledge-notes/time-parallelization/chapter-4-parabolic-methods
tags:
  - 时间并行
  - 抛物-PDE
---

> [!note] 阅读范围
> 本章对应论文 Section 4（pp. 443–481）。正文严格保留 Section 4 导论、4.1 历史发展、4.2 Parareal、4.3 PFASST、4.4 MGRiT、4.5 两种对角化 Parareal 和 4.6 STMG 的原始层级。Python 复现、参数比较和覆盖审计统一标为本站补充，不占用论文编号。理论收敛因子、实测迭代曲线和墙钟性能采用不同标签。

## 论文与精读页的对应关系

本页保留章节全景与本站复现实验。逐公式、逐定理和逐图表的论证拆在以下页面：

| 论文小节           | 原文页码    | 精读页                                                                                                                      | 覆盖范围                                                 |
| ------------------ | ----------- | --------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------- |
| Section 4、4.1–4.2 | pp. 443–452 | [[computational-mathematics/knowledge-notes/time-parallelization/chapter-4-1-parareal\|历史发展与 Parareal]]                | (4.1)–(4.9), Theorems 4.1–4.4, Figures 4.1–4.5           |
| 4.3–4.4            | pp. 452–461 | [[computational-mathematics/knowledge-notes/time-parallelization/chapter-4-2-pfasst-mgrit\|PFASST 与 MGRiT]]                | (4.10)–(4.13), Theorems 4.5–4.6, Figures 4.6–4.11        |
| 4.5                | pp. 460–471 | [[computational-mathematics/knowledge-notes/time-parallelization/chapter-4-3-diagonalized-parareal\|基于对角化的 Parareal]] | (4.14)–(4.29), Theorems 4.7–4.8, Figures 4.12–4.17       |
| 4.6                | pp. 472–481 | [[computational-mathematics/knowledge-notes/time-parallelization/chapter-4-4-stmg\|时空多重网格]]                           | (4.30)–(4.44), Theorem 4.9, Figures 4.18–4.22, Table 4.1 |

## Section 4 导论与 4.1 历史发展

### 抛物方法利用了什么

第三章的方法可以直接处理长程传播，但各自带有约束：SWR 需要传输条件，ParaExp 主要适合线性问题，ParaDiag 的非线性版本要在长时间窗上求解 Newton 系统。抛物方程提供了另一条路。扩散会迅速削弱高频误差，让较晚状态对很久以前的精细信息不再敏感。粗时间模型由此可能保留真正重要的慢模态。

本章讨论 Parareal、PFASST、MGRiT、两种对角化 Parareal 和时空多重网格（STMG）。前三类主要在时间上构造粗细层级，通常保持同一个空间网格；STMG 同时粗化空间和时间，并使用可沿时间并行的块 Jacobi 平滑器。

![抛物型时间并行方法的粗细层级与并行位置](assets/diagrams/pint/zh/parabolic-multilevel-map.svg)

历史上，Parareal 与多重打靶、Nievergelt 的预计算思想和 Saha 的方法有关。后续发展出 PITA、PFASST 与 MGRiT 等变体。它们都通过低成本全局近似传递长程信息，再用并行的局部细计算修复细节。进入双曲极限后，粗层的相位和传播误差不再自动消失，收敛便会明显退化。

## 4.2 Parareal

### 算法与并行结构

把 $[0,T]$ 分为 $N$ 个大区间。$\mathcal F$ 是准确且昂贵的细传播子，$\mathcal G$ 是便宜的粗传播子。Parareal 更新为

$$
U_{n+1}^{k+1}
=\mathcal G(U_n^{k+1})
+\mathcal F(U_n^k)-\mathcal G(U_n^k). \tag{4.1}
$$

在第 $k$ 轮中，所有 $\mathcal F(U_n^k)$ 可以同时计算；新的 $\mathcal G(U_n^{k+1})$ 仍按 $n$ 顺序传播。粗项建立最新的因果预测，差值 $\mathcal F(U_n^k)-\mathcal G(U_n^k)$ 修复上轮各区间的粗细失配。

Parareal 可视为多重打靶的近似 Newton 法，也可视为全时间下三角系统的粗传播预条件迭代。它是非侵入式算法：已有的细、粗时间推进器只要暴露“给定初值，推进一个大区间”的接口即可接入。

### 线性模态分析

对空间特征模态 $\lambda_\ell$，设细传播一个小步的稳定函数为 $R_f(\lambda_\ell\Delta t)$，每个大区间有 $J$ 个小步；粗传播稳定函数为 $R_g(\lambda_\ell\Delta T)$。Theorem 4.1 把误差迭代写成严格下三角 Toeplitz 矩阵。严格下三角性带来有限步性质：精确算术下最多 $N$ 轮得到串行细解，且第 $k$ 轮后前 $k$ 个大区间已经精确。

有限步终止不能自动带来加速。若达到容差需要接近 $N$ 轮，每轮又包含一次顺序粗传播，时间并行的收益会很小。

Theorem 4.2 给出两类互补估计。时间区间数较小时，组合型上界体现随 $k$ 加快的超线性收敛；长时间区间中，更有解释力的是逐模态线性因子

$$
\rho_\ell
=\frac{|R_g-R_f^J|}{1-|R_g|}. \tag{4.5b}
$$

分子衡量粗、细传播差异，分母衡量粗传播的耗散裕量。Figure 4.2 表明，短区间的超线性阶段和长区间的近线性阶段可以同时存在。

### 非线性有限步性质

Theorem 4.3 在粗传播具有 Lipschitz 稳定性和 $p$ 阶局部误差时，给出非线性 Parareal 的有限步与超线性误差界。核心机制仍是逐区间建立精确性。非线性会让常数依赖解的正则性和传播子的 Lipschitz 常数，因而实际迭代数可能对时间窗和动力学参数非常敏感。

### 时间积分器怎样决定极限因子

Theorem 4.4 研究负实轴上的抛物模态。若细传播使用 L-稳定方法，粗传播使用后向 Euler，并且每个大区间包含足够多的细步，最坏长时间因子约为 $0.3$。若细方法仅 A-稳定，高频模态不会被充分衰减，所需 $J$ 会随最危险频率增长。Figure 4.3 对比梯形规则与两种 SDIRK 细传播，说明选定粗方法以后，细方法与粗细步数比仍然会显著改变收敛区域。若在连续层面（$\mathcal F$ 取精确传播子）改用 Radau IIA 作为**粗**传播子，理论最坏因子可降到约 $0.068$。

### 从热方程走向双曲极限

Figures 4.4–4.5 使用周期 ADE 与 Burgers 方程，$T=4$、$\Delta T=0.1$、$\Delta x=1/128$、$J=32$，粗传播为后向 Euler，细传播为二阶 L-稳定 SDIRK。黏性下降时，ADE 的粗细相位失配增大；Burgers 方程还会改变局部传播速度和激波位置。两类方程都明显变慢，Burgers 在 $\nu\le10^{-3}$ 附近出现近似发散。

波动方程更加严苛。除非粗传播与细传播几乎一样准确，Parareal 很难控制相位误差；达到这种精度的粗传播往往也失去成本优势。线性平流可以用半 Lagrange 或相位优化粗传播改善，通用非线性双曲粗模型仍然是开放问题。

## 4.3 PFASST

### 从 SDC 到时间多层

PFASST 把 Parareal 的大区间并发与 SDC 的高阶配置求解结合起来。每个时间步上选择 $M_f$ 个细配置节点，配置方程写成

$$
\boldsymbol U_f
=\boldsymbol U_{0,f}
+\Delta t\,Q_f\boldsymbol f(\boldsymbol U_f). \tag{4.10}
$$

直接求解稠密配置系统很昂贵。SDC 使用基于隐式 Euler 的易解下三角近似做预条件 sweep，每次 sweep 消除一部分配置残差。PFASST 再引入含 $M_c$ 个节点的粗层，让不同时间步上的细、粗 sweep 形成流水线。

细、粗节点之间通过 Lagrange 插值限制与延拓。粗层采用全近似格式（FAS），把细层残差以 $\tau$ 校正带入粗配置方程。一个 PFASST 迭代包含并行细 sweep、细到粗传递、顺序或流水化粗 sweep、粗到细校正。这个结构既可从 Parareal 理解，也可看作配置方程上的多重网格。

### 数值观察与限制

Figure 4.6 使用周期热方程与 ADE，$T=3$、$\Delta x=1/128$、$\Delta t=1/64$，热源参数 $\sigma=1000$，细层为三个 Radau IIA 节点，粗层为两个节点。热方程收敛较快；ADE 的黏性下降后，跨时间步的高频传播越来越难被粗配置层表示，收敛随之变慢。

PFASST 适合需要高阶时间精度且单步配置求解昂贵的场景。实际性能取决于配置节点、SDC sweep 数、粗层成本、节点传递、流水线填充以及空间并行资源分配。只报告配置阶数不足以判断并行效率。

## 4.4 MGRiT

### F 点、C 点与 FCF 松弛

MGRiT 在时间网格上每隔 $J$ 个细点选一个 C 点，其余为 F 点。F 松弛在各相邻 C 点之间并行细推进；C 松弛更新粗点；粗网格校正负责长距离传播。两层 FCF 迭代可写成与 Parareal 相似的重叠更新：先做 F 松弛，再做 C 与第二次 F 松弛。

一次 FCF 大约使用两组细传播，因此每轮通常比 Parareal 贵。额外的 CF 段提供重叠，理论上每轮可以让两个大区间进入精确区；最多约 $\lceil N/2\rceil$ 轮达到串行细解。更一般的 $F(CF)^\nu$ 用更多重叠换取更强收缩。

### 收敛因子和公平成本比较

Theorem 4.5 给出两层 FCF 的长时间模态因子

$$
\rho_{\mathrm{MGRiT},\ell}
=|R_f^J|
\frac{|R_g-R_f^J|}{1-|R_g|}.
$$

它比 Parareal 因子多一个 $|R_f^J|$，因此耗散模态会额外收缩。每增加一个 CF 段又会乘一个类似细传播因子，同时增加一组细求解成本。

公平比较应把一次 FCF 与两次 Parareal 放在接近的细求解工作量上。Theorem 4.6 对 L-稳定细方法给出代表性最坏因子。后向 Euler 粗传播下，Parareal 与一次 FCF 约为 $0.2984$ 和 $0.1115$；二阶 Lobatto IIIC 组合约为 $0.0817$ 和 $0.0197$。一次 FCF 的单轮因子更小，但仍可能略差于两次 Parareal 因子的平方。Figure 4.8 以复平面收敛区域进一步展示这种等成本关系。

### ADE 与 Burgers 实验

Figures 4.9–4.10 使用 $T=5$、$J=20$、$\Delta T=1/8$、$\Delta x=1/160$，粗传播为后向 Euler，细传播为 SDIRK22。热方程的逐模态因子远小于 1。随着 ADE 黏性从 $0.1$ 降到 $0.01$ 和 $0.002$，危险模态逐渐靠近高相位、低耗散区域。在 $\nu=0.002$ 时，Parareal 与 FCF 的线性阶段因子分别约为 $1.4211$ 和 $1.2812$，两者都会暂态增长。

MGRiT 的有限步下降仍然存在；它只说明严格下三角误差最终清零，无法消除前几十轮的放大。若实际容差要求在暂态增长阶段就停止，方法没有形成可扩展求解器。

Figure 4.11 对非线性 Burgers 方程得到相同主线。扩散充分时 FCF 每轮优于 Parareal；按细传播次数配平后，一次 FCF 常与两次 Parareal 接近。粗传播无法准确表示非线性激波时，两者都退化。

## 4.5 基于对角化的 Parareal

MGRiT 仍受顺序且可能失真的粗校正限制。Section 4.5 从两个方向消除
这段瓶颈：第一种跨 $N$ 个粗点对角化 CGC，第二种在每个粗区间内部
对角化 $J$ 个细点，构造低成本粗传播。

### 4.5.1 基于对角化的 CGC

标准 Parareal 的粗校正按时间顺序执行。对角化版本加入首尾条件

$$
U_0^{k+1}=u_0+\alpha U_N^{k+1}
$$

并同步修改校正右端，使原始初值问题仍是迭代的不动点。线性后向 Euler 情形产生 $\alpha$-循环全时间矩阵，可以经过 FFT 化为 $N$ 个独立的移位粗空间问题。

$\alpha\to0$ 时恢复标准顺序粗校正，但变换的舍入误差会增强。Theorem 4.7 表明，若标准 Parareal 的收敛因子为 $\rho$，选择

$$
\alpha\le \frac{\rho}{1+\rho}
$$

可保持同量级收敛。该结论是对负实谱证明的；复谱只有数值证据。Figure 4.12 对热方程测得 $\rho\approx0.22$、阈值约 $0.18$；对 $\nu=0.1$ 的 ADE（复谱）测得 $\rho\approx0.39$、阈值约 $0.28$，与理论预测吻合。

非线性版本把全时间粗方程做准 Newton 线性化，再用平均 Jacobian 构造块 $\alpha$-循环系统。Figure 4.13 的 Burgers 实验显示了相同的 $\alpha$ 阈值现象。

直接把这条 head-tail 条件放进 MGRiT 会改变不动点并导致发散。
一致的替代条件改用最近两轮末值之差；小 $\alpha$ 时，它与原 MGRiT
保持相同收敛率，也可反过来用于 Parareal。

### 4.5.2 基于对角化的粗传播子

第二种方法让细、粗传播使用同一个时间积分器和同一个小步长。细传播顺序执行 $J$ 个小步；粗传播把这 $J$ 个未知时间点放进一个 $\alpha$-循环 all-at-once 系统并并行求解。$\alpha=0$ 时粗传播等于细传播，Parareal 一轮收敛，却没有并行成本优势；$\alpha>0$ 产生可对角化近似，粗传播有望约便宜 $J$ 倍。

Theorem 4.8 给出清晰的方程类型差异。负实谱的抛物问题具有因子

$$
\rho=\alpha,
$$

并且与大区间数 $N$ 无关。纯虚谱的波动问题满足上界

$$
\rho\le \frac{2\alpha N}{1+\alpha}.
$$

Figure 4.14 对热方程验证了 $\alpha$ 斜率。Figure 4.15 的波动实验显示，$\alpha=0.01$ 时增加 $N$ 会变慢；实际曲线仍有超线性阶段，从 $N=24$ 增到 $960$ 只多约两轮达到离散误差。小 $\alpha N$ 区域的线性上界最接近实际。

Figures 4.16–4.17 把方法用于 Burgers 方程。误差因子随 $\alpha$ 近似线性，$\alpha=10^{-3}$ 时对时间区间数较鲁棒；减小 $\alpha$ 也能缓解小黏性带来的退化。这种粗传播同时适用于抛物与双曲问题，是两种对角化 Parareal 中适用面更广的一种。

## 4.6 时空多重网格

### 全时间系统与块 Jacobi 平滑

对线性系统 $\boldsymbol u'=A\boldsymbol u+\boldsymbol g$，一般单步积分器写成

$$
r_1\boldsymbol u_{n+1}=r_2\boldsymbol u_n+\widetilde{\boldsymbol f}_n,
\qquad n=0,\ldots,N_t-1, \tag{4.30}
$$

其中 $r_1,r_2$ 是 $\Delta tA$ 的矩阵多项式。叠起全部时间点得到 all-at-once 系统

$$
K\boldsymbol U
=\begin{bmatrix}
r_1\\-r_2&r_1\\&\ddots&\ddots\\&&-r_2&r_1
\end{bmatrix}
\boldsymbol U
=\boldsymbol b. \tag{4.31}
$$

STMG 采用阻尼块 Jacobi：

$$
\boldsymbol U^{j+1}
=\boldsymbol U^j
+\eta\left(I_t\otimes r_1\right)^{-1}
(\boldsymbol b-K\boldsymbol U^j),
\qquad r_1=I_x-\theta\Delta t A. \tag{4.32}
$$

$I_t\otimes r_1$ 在时间上是块对角矩阵，所以所有时间点的空间块可以并行求解。阻尼 $\eta$ 决定高频误差的削弱速度。

两层 cycle 先做 $s_1$ 次预平滑，把残差同时限制到粗空间、粗时间网格，求解粗网格方程，延拓校正，再做 $s_2$ 次后平滑。递归应用便得到完整 STMG。它同时粗化 $\Delta X=2\Delta x$ 与 $\Delta T=2\Delta t$，把空间和时间的长波误差交给粗层。

早期抛物多重网格常用逐点 Gauss–Seidel 沿时间推进。该平滑器本质上是串行前代，虽然空间粗化有效，时间方向难以扩展。块 Jacobi 改变了并行位置，是现代 STMG 的关键。

### 局部 Fourier 分析

对一维热方程、中心空间差分和后向 Euler，局部 Fourier 分析把误差分成时空频率并求出平滑符号。Theorem 4.9 给出最优阻尼

$$
\eta=\frac12.
$$

在高时间频率集合上最大化 $|\rho|$ 后，平滑因子至多为
$1/\sqrt2$，因此可以安全做时间粗化。在归一化热方程中，若
$\Delta t/\Delta x^2\ge1/\sqrt2$，空间高频也被压到同一上界；
恢复扩散系数后，无量纲比变为 $\nu\Delta t/\Delta x^2$。ADE
符号虽含虚部，$\eta=1/2$ 仍是后向 Euler 的稳健起点。

### 阻尼、平滑次数与时间积分器

Figure 4.19 扫描后向 Euler 两层 STMG 的 $\eta$。热方程和 $\nu=0.01$ 的 ADE 都在 $1/2$ 附近表现良好。Figure 4.20 固定 $\eta=1/2$，比较一次与三次块 Jacobi 平滑。更多平滑使每个 cycle 更贵，同时显著减少 cycle 数；ADE 比热方程慢，但平滑次数增加后对黏性的敏感性减弱，并出现超线性阶段。

Figure 4.21 改用梯形规则后，热方程即使增加到十次平滑也很难稳定收敛；ADE 可以收敛，较多平滑会改善速度，采样中的较优阻尼约为 $0.8$。这说明 STMG 的效果依赖时间积分器的稳定函数。后向 Euler 的理论阻尼不能直接移植到梯形规则。

### 大规模扩展结果

Table 4.1 汇总三维热方程上的现代 STMG 数据。弱扩展中，核心数从 1 增到 262,144，时间步从 2 增到 524,288，自由度从 59,768 增到 15,667,822,592；迭代数保持 7，总时间约从 28.8 秒变到 30.0 秒。只沿空间并行的顺序时间推进列随规模增长到约 4,988,060 秒。

强扩展中，一个含 512 个时间步、15,300,608 个自由度的问题从 1 核
的 7,635.2 秒降到 256 核的 30.0 秒；更大的 524,288 时间步问题从
512 核的 15,205.9 秒降到 262,144 核的 30.0 秒。这些数据取自
Gander 与 Neumüller（2016）的三维实现，展示了 STMG 在抛物问题上的
弱、强扩展潜力。

### 非线性 FAS-STMG

非线性半离散系统 $\boldsymbol u'=\boldsymbol f(\boldsymbol u)$ 经过 $\theta$ 方法后形成

$$
K(\boldsymbol U)
=(B\otimes I_x)\boldsymbol U
-\Delta t(\widetilde B\otimes I_x)\boldsymbol f(\boldsymbol U)
=\boldsymbol b. \tag{4.42}
$$

非线性块 Jacobi 在每个时间块并行求解局部校正，内部可使用 Newton。局部 Fourier 分析不再直接适用，阻尼需要根据问题测试。粗层采用 FAS：限制当前解与残差，求解带 $\tau$ 一致性校正的非线性粗问题，延拓粗解差，再做后平滑。

Figure 4.22 对 Burgers 方程使用两次块 Jacobi 平滑。扩散充分时非线性 STMG 收敛很快；黏性降低后明显恶化。实验中较优阻尼为 $\eta=1/4$，体现了非线性与线性后向 Euler 理论的差别。

综合 Section 4.6，STMG 是当前抛物问题上最有力的 PinT 方法之一，
并已展示大规模扩展。它需要访问全时间离散、平滑器和网格传递，侵入性
高于 Parareal。现有图直接展示的是低扩散、近双曲极限下的退化；真正
双曲问题以及不同时间积分器下的鲁棒性仍需进一步研究。

## 本站数值补充：Python 复现结果

### Parareal 基线与 Figure 4.5

| 问题     | 正式参数                              | 迭代次数 |          最终最大误差 |
| -------- | ------------------------------------- | -------: | --------------------: |
| 对流扩散 | $N_x=128$，$N=40$，$J=32$，$\nu=0.02$ |       39 | $6.106\times10^{-15}$ |
| Burgers  | $N_x=128$，$N=40$，$J=32$，$\nu=1$    |       16 | $3.544\times10^{-12}$ |

![对流扩散方程上的 Parareal 收敛基线](assets/pint/parareal-ade-baseline.svg)

![黏性 Burgers 方程上的 Parareal 收敛基线](assets/pint/parareal-burgers-baseline.svg)

固定 $T=4$、$\Delta T=0.1$、$\Delta x=1/128$ 与 $J=32$，把最大误差降到 $10^{-10}$ 以下所需迭代为：

| 方程     | $\nu=1$ | $\nu=0.1$ | $\nu=0.02$ |
| -------- | ------: | --------: | ---------: |
| 对流扩散 |      14 |        24 |         35 |
| Burgers  |      14 |        21 |         25 |

![扩散减弱时 ADE 与 Burgers 方程的 Parareal 收敛](assets/pint/parareal-figure-4-5.svg)

这些迭代数衡量向串行细解的收敛。GPU 性能数据见[[computational-mathematics/knowledge-notes/time-parallelization/chapter-5-unified-view#本站复现gpu-加速与性能剖析|第五章]]。

### MGRiT 基线与 Figures 4.9–4.10

近双曲基线使用 $N_x=160$、$N=40$、$J=20$、$T=5$、$\nu=0.002$。报告窗口末，Parareal 最大误差为 $2.895\times10^2$；两层 MGRiT 最终通过有限步性质降到 $5.551\times10^{-16}$。

![近双曲 ADE 上 Parareal 与 MGRiT 的基线比较](assets/pint/mgrit-baseline.svg)

按“一次 FCF 对两次 Parareal”的细求解成本比较，逐模态长时间因子为：

| 问题             | Parareal 因子 | FCF-MGRiT 因子 |
| ---------------- | ------------: | -------------: |
| 热方程           |        0.2824 |         0.0835 |
| ADE，$\nu=0.1$   |        0.4453 |         0.2719 |
| ADE，$\nu=0.01$  |        1.0501 |         0.9021 |
| ADE，$\nu=0.002$ |        1.4211 |         1.2812 |

![Heat 与不同黏性 ADE 的逐模态长时间收敛因子](assets/pint/mgrit-figure-4-9.svg)

![Parareal 与 FCF-MGRiT 的等细求解成本曲线](assets/pint/mgrit-figure-4-10.svg)

> [!note] Figure 4.9 的一个数值差异
> 原论文图中把 $\nu=0.01$ 的 Parareal 最大因子标为 $0.9986$；Python 转换项目按论文公式 (4.5b)、打印参数和上游 `MGRiT_Heat_ADE.m` 的稳定函数直接计算，得到 $1.0501$。同一算例的 MGRiT 因子 $0.9021$ 以及其余两组因子均与原图基本一致。这里保留复现值并明确记录差异，不用原图标注覆盖计算结果。两个数都位于 $1$ 附近，支持的定性结论相同：这一黏性下长时间收敛已经接近临界状态。

### STMG 阻尼与平滑验证

梯形规则基线使用 $N_x=N_t=255$、$\nu=10^{-3}$ 及三次前、后平滑。采样中，15 个 cycle 后最小误差位于 $\eta\approx0.98$。

![梯形规则下 STMG 阻尼参数的基线扫描](assets/pint/stmg-baseline.svg)

在 Figure 4.19 对应的后向 Euler 网格上复算得到：

| 问题            | 15 次迭代后采样到的最佳 $\eta$ |
| --------------- | -----------------------------: |
| 热方程          |                          0.500 |
| ADE，$\nu=0.01$ |                          0.372 |

![后向 Euler STMG 的阻尼参数扫描](assets/pint/stmg-figure-4-19.svg)

ADE 的离散采样最小点为 $0.372$，而高频理论和 Figure 4.18 把
$1/2$ 作为稳健选择。两者衡量的是特定有限网格上的 15-cycle 误差
和高频平滑上界，结论层次不同。

固定 $\eta=0.5$ 后，三次前、后平滑比一次平滑使用更少 cycle 降低误差：

![一次与三次前后平滑的 STMG 收敛比较](assets/pint/stmg-figure-4-20.svg)

真正的性能决策应比较

$$
\text{达到容差的总工作量}
=(\text{单 cycle 成本})\times(\text{cycle 数})
+\text{通信与内存开销}.
$$

## 本站方法比较：参数与失效方式

| 方法            | 关键参数                             | 参数控制的折中                 | 典型失效方式                    |
| --------------- | ------------------------------------ | ------------------------------ | ------------------------------- |
| Parareal        | 区间数 $N$、粗细步长比、$\mathcal G$ | 并发度、粗成本与传播准确度     | 迭代数接近 $N$，相位误差增长    |
| PFASST          | 配置节点、SDC sweep、粗节点          | 高阶精度、流水线深度与粗层成本 | 弱扩散下粗配置层失真            |
| MGRiT           | 粗化因子、F/FCF/$F(CF)^\nu$          | 重叠收缩与细求解工作量         | 等成本优势消失、暂态放大        |
| 对角化 Parareal | $\alpha$、对角化位置                 | 串行粗校正比例与舍入误差       | $\alpha$ 过小病态，过大改变收敛 |
| STMG            | $\eta$、平滑次数、粗化、积分器       | 高频平滑与 cycle 成本          | 积分器不匹配，双曲/小黏性退化   |

## 原文覆盖核对

| 原文位置                       | 本页对应 | 已覆盖内容                                                                         |
| ------------------------------ | -------- | ---------------------------------------------------------------------------------- |
| Sections 4 与 4.1，pp. 443–444 | 4.1      | 抛物时间局部性、方法范围、历史与 STMG 的区别                                       |
| Section 4.2，pp. 444–452       | 4.2      | Parareal 更新、Theorems 4.1–4.4、Figures 4.1–4.5、非线性与双曲退化                 |
| Section 4.3，pp. 452–455       | 4.3      | 配置方程、SDC、FAS 传递、PFASST 迭代与 Figure 4.6                                  |
| Section 4.4，pp. 455–460       | 4.4      | FCF 结构、Theorems 4.5–4.6、Figures 4.7–4.11、等成本比较与 Burgers                 |
| Section 4.5，pp. 460–471       | 4.5      | 两种对角化位置、Theorems 4.7–4.8、Figures 4.12–4.17、非线性和 MGRiT 一致性条件     |
| Section 4.6，pp. 472–481       | 4.6      | all-at-once STMG、块 Jacobi、Theorem 4.9、Figures 4.18–4.22、Table 4.1、非线性 FAS |

## 小结

这一组方法把扩散造成的时间局部性转化为算法优势。Parareal 用廉价粗传播修正并行细传播；PFASST 在配置节点上叠加 SDC 与 FAS；MGRiT 用时间多层和重叠松弛增强粗校正；对角化变体减少顺序粗传播；STMG 同时处理空间与时间尺度。它们在抛物问题上可以非常高效，性能仍取决于粗细传播的谱匹配、时间积分器、等成本比较和完整硬件开销。黏性下降后，相位与长程记忆重新成为主导，标准粗网格机制会逐步失去优势。

## 本章原文

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Section 4, pp. 443–481。
