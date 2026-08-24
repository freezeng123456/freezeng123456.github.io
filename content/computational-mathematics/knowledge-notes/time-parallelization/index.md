---
title: 双曲与抛物问题的时间并行方法
description: 时间并行方法的分章笔记与可复现实验
lang: zh
translation: en/computational-mathematics/knowledge-notes/time-parallelization
tags:
  - 计算数学
  - 时间并行
---

本专题依据 M. J. Gander、S.-L. Wu 和 T. Zhou 的综述 _Time Parallelization for Hyperbolic and Parabolic Problems_（Acta Numerica 34, 2025, pp. 385-489）整理。原文区分了对传播型问题仍然有效的方法，以及主要为耗散问题设计的方法。

> [!info] 原论文图表与许可
> 精读页收录 Figure 1.1–4.22 和 Tables 3.1、3.2、4.1 的完整原图，共 48 个图表资产，均从论文 PDF 提取为可缩放 SVG；相邻中文说明重新撰写。论文以 [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) 发布，作者、题名和 DOI 见本页“第一手资料”。本站重绘图与 Python 复现实验均另行标注，避免和原图混淆。

## 阅读顺序

1. [[computational-mathematics/knowledge-notes/time-parallelization/chapter-1-why-parallelize-in-time|第一章：为什么要做时间并行？]]从硬件瓶颈和因果链出发，建立四条方法谱系、动力学分类与全时间系统视角。
2. [[computational-mathematics/knowledge-notes/time-parallelization/chapter-2-model-problems|第二章：连接抛物与双曲世界的模型问题]]用热方程、对流扩散、Burgers 和波动方程比较边界条件、耗散与时间记忆。
3. [[computational-mathematics/knowledge-notes/time-parallelization/chapter-3-hyperbolic-methods|第三章：对双曲问题有效的 PinT 方法]]沿“界面传播—校正流水线—指数传播—全时间代数”推导 SWR、PIDC/RIDC、ParaExp 和 ParaDiag。
4. [[computational-mathematics/knowledge-notes/time-parallelization/chapter-4-parabolic-methods|第四章：为抛物问题设计的 PinT 方法]]从粗细传播失配出发，依次分析 Parareal、PFASST、MGRiT、两种对角化路线和 STMG。
5. [[computational-mathematics/knowledge-notes/time-parallelization/chapter-5-unified-view|第五章：结论]]收束为方法选型、证据分层与复现规范，并单列 GPU 性能记录。

## 论文内容覆盖进度

“段落级完成”表示正文论点、公式、图示、历史线索、限定条件和章节关系都已与原文逐项核对。已有实验不会因此删除，它们会放在对应原文解释之后。

| 原文范围           | 网页对应 | 当前状态                                                                                     |
| ------------------ | -------- | -------------------------------------------------------------------------------------------- |
| 摘要与 Section 1   | 第一章   | **段落级完成**：覆盖 pp. 385–388 的全部论点、公式 (1.1) 和 Figure 1.1 原图                   |
| Sections 2.1–2.4   | 第二章   | **段落级完成**：覆盖 (2.1)–(2.7)、四类模型、全部边界设置、Figures 2.1–2.4 原图与三个补充实验 |
| Sections 3.1–3.5.2 | 第三章   | **段落级完成**：五组精读页覆盖 (3.1)–(3.68)、Theorems 3.1–3.9、Figures 3.1–3.18 和两张表     |
| Sections 4.1–4.6   | 第四章   | **段落级完成**：四组精读页覆盖 (4.1)–(4.44)、Theorems 4.1–4.9、Figures 4.1–4.22 和 Table 4.1 |
| Section 5          | 第五章   | **段落级完成**：覆盖论文结论，并将统一视角、GPU 分析和实验清单明确标为本站补充               |

> [!warning] 实现覆盖不等于无例外复现
> 2026-08-24 的 [[computational-mathematics/knowledge-notes/time-parallelization/reproduction-audit-2026-08-24|全文实验复现审计]]确认：45 张编号图全部已有生成器或计算式重绘，3 张表中 2 张可本地计算；但 Figure 3.14 的两条 NKA 曲线仍未匹配，Table 4.1 只是引用的外部超算结果，而且当前开放依赖安装有两个数值容差测试失败。专题的“段落级完成”描述网页对论文内容的覆盖，不应误读为每项数值和工程环境都已无差异复现。

## 方法图谱

| 方法      | 并行单元           | 机制             | 自然适用区域     |
| --------- | ------------------ | ---------------- | ---------------- |
| SWR       | 重叠时空子域       | 波形传输         | 输运与波动问题   |
| PIDC/RIDC | 校正层与时间节点   | 延迟校正与流水线 | 初值问题         |
| ParaExp   | 非齐次与齐次子问题 | 指数传播         | 线性系统         |
| ParaDiag  | 全时间耦合矩阵     | 循环近似与 FFT   | 线性或线性化系统 |
| Parareal  | 粗时间区间         | 粗预测加细校正   | 中等到强耗散     |
| PFASST    | 跨时间步的配置节点 | SDC 与多层校正   | 高阶时间积分     |
| MGRiT     | 时间网格层次       | 松弛与粗网格校正 | 长时间区间       |
| STMG      | 完整时空网格       | 时空平滑与粗化   | 抛物系统         |

## 三条组织原则

### 因果约束如何进入并行计算

单步方法 $u_{n+1}=\Phi_{\Delta t}(u_n)$ 形成串行递推。时间并行方法把所有时间点的耦合显式写出，再构造其逆算子的并行近似。并发工作由此进入迭代、分解或变换过程。

### 耗散决定粗表示是否包含有效信息

抛物动力学会衰减高频误差，因此粗传播子仍能表示长时间尺度上保留下来的慢变分量。双曲动力学保留相位信息，粗问题中的微小相位误差可能沿多个时间区间累积。

### 迭代次数不等于并行效率

一个简化成本模型是

$$
T_{\mathrm{parallel}}
\approx K(C_G+C_F/P)+C_{\mathrm{comm}}+C_{\mathrm{setup}},
$$

其中 $K$ 是迭代次数，$C_G$ 和 $C_F$ 分别是粗、细传播成本，$P$ 是时间并发度。本专题的数值实验测量算法收敛性，不测量端到端并行加速比。

> [!note] 数值来源
> 页面中的原论文图表直接取自论文 PDF，并保留原坐标与面板。标为“本站复现”的实验于 2026 年 7 月 31 日在新的实验服务器上由 Python 项目重新生成。初始正式结果使用 SciPy CPU 路径；随后加入 CuPy/T4 混合后端，将独立 Burgers 细传播批量放到 GPU。早期 `paper_validation` 内 Figure 4.5 的停止迭代保持一致，组合验证套件由 263.57 秒降至 67.92 秒；它使用的 Burgers 通量和当前全文 `section4` 入口不同，也不是四个全文分组的总墙钟。
>
> 上述 T4 数字是 2026-07-31 的已归档单 GPU 记录；其范围和当前依赖漂移见 [[computational-mathematics/knowledge-notes/time-parallelization/reproduction-audit-2026-08-24|复现审计]]。

## 第一手资料

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), Acta Numerica 34 (2025), pp. 385-489.
- 原始 MATLAB 算例：[wushulin/ActaPinT](https://github.com/wushulin/ActaPinT)。
- Python 转换、扩充与正式结果：[freezeng123456/ActaPinT-Python](https://github.com/freezeng123456/ActaPinT-Python)。
