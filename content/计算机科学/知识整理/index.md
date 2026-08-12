---
title: 机器学习知识整理
description: 从算法抽象到训练、推理与硬件执行的机器学习学习笔记
lang: zh
translation: en/computer-science/knowledge-notes
tags:
  - 机器学习
draft: false
---

这一部分按“问题—机制—数据流—性能边界”组织机器学习知识，不堆叠零散链接。

## AI 系统

### [[计算机科学/知识整理/efficient-llm-inference|高效大模型推理：从 Token 调度到 GPU Kernel]]

以 vLLM 为主线，用四张流程图和张量维度表解释：

- continuous batching 如何把调度粒度从请求下沉到 token；
- PagedAttention 如何用 block table 管理动态 KV Cache；
- Llama 3 8B 的数据如何流经 Attention、FFN、LM head 与 sampler；
- prefill、decode、抢占分别受什么硬件与资源约束。

### [[计算机科学/知识整理/how-gpus-work|GPU 如何工作：从 CUDA 生命周期到 SIMT]]

把 GPU 从一个抽象的“加速器”展开为可追踪的执行链路：

- CPU 如何通过运行时、驱动、命令队列和 doorbell 向 GPU 提交工作；
- CUDA 的 Grid、Block 与 Thread 如何映射到 SM、Warp 和执行单元；
- 寄存器、共享内存与线程上限如何共同限制 occupancy；
- 为什么高 occupancy 不保证高性能，以及 SIMD、SIMT 和 warp divergence 的关系。

### [[计算机科学/知识整理/nvidia-gemm-ping-pong|数学博士眼中的 NVIDIA GEMM 中的 ping-pong 流水]]

微信公众号外部来源全文归档，保留原文的公式、代码、图注和 7 张正文图片，聚焦：

- shared memory 双缓冲与多级 software pipeline 如何隐藏操作数搬运延迟；
- Hopper 的 TMA、WGMMA、warp specialization 与 two-consumer ping-pong 如何错峰隐藏 epilogue；
- cooperative 与 ping-pong 的资源权衡，以及 Blackwell 的 TMEM、tcgen05 和 CLC 调度迁移。

## 内容组织

- **基础理论**：优化、泛化、表示学习、在线学习；
- **系统实现**：训练系统、推理服务、数据管线、可观测性；
- **实验方法**：对照、消融、统计不确定性、可复现协议。
