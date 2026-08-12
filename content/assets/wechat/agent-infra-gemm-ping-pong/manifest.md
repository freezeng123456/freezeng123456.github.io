# 《数学博士眼中的NVIDIA GEMM 中的 ping-pong 流水》图片归档

本目录保存文章正文中出现的 7 张图片，按源页面出现顺序编号。图片来自外部数据源，版权归原作者及相关权利人所有；公开发布到本知识库基于用户已确认具备正文图片的转载/公开发布权限。

- 原文：[微信公众号原文](https://mp.weixin.qq.com/s/pWC6dq2MJZgRBCnRttKCFg)
- 作者：Agent&Infra
- 抓取日期：2026-08-12
- 下载格式：PNG；使用源页面提供的图片地址保存到本地，未对图片内容进行编辑

| 本地文件        | 原文图注                                                                                       | 外部图片源                                                                                                                                                                                                      |
| --------------- | ---------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `figure-01.png` | Roofline 模型：从 bandwidth-bound 到 compute-bound。来源：JAX Scaling Book / Google DeepMind   | [mmbiz.qpic.cn 图片地址](https://mmbiz.qpic.cn/mmbiz_png/QPNjE9VGkMU2ha32LpO3biav4st3btZ2SZXO4S652Wiazj41t8vHMdbVDhOX8niaJn7etuuxK8XXGREphJlicZ1yoZwgprPLAF83KvnY8Rtwyiak/640?wx_fmt=png&from=appmsg)           |
| `figure-02.png` | LOAD 与 MMA 的流水重叠：把等待时间压进计算时间。来源：Colfax Research                          | [mmbiz.qpic.cn 图片地址](https://mmbiz.qpic.cn/sz_mmbiz_png/QPNjE9VGkMXMv0uMyprhwjJx5ibWas2JZOYDRVyupgTniahvwcRDQiaW5lJ4niamOfiaibexPm6eADPibibR0vicRmhqJUuWXxick1LZnwgqibD4Lrudick/640?wx_fmt=png&from=appmsg) |
| `figure-03.png` | 两个 alternating SMEM stages：S_0 与 S_1。来源：Colfax Research                                | [mmbiz.qpic.cn 图片地址](https://mmbiz.qpic.cn/mmbiz_png/QPNjE9VGkMW4rbWsWyJVLtxxcDUl6gKNRqibnjDala2VWAt7czv7eq2jMicGLRjwKW7aibwYcwNT9XsLSm8FDicjB1ndqcaZ9WRhT3Pv38qXg24/640?wx_fmt=png&from=appmsg)            |
| `figure-04.png` | CUTLASS GEMM mainloop 的 software pipeline。来源：NVIDIA CUTLASS Documentation                 | [mmbiz.qpic.cn 图片地址](https://mmbiz.qpic.cn/mmbiz_png/QPNjE9VGkMW9dFZuicHVW7XGcvOsUgG3qjPx9A90bQEZsp2QYmeiaiaf3305ajLxiaACtu8LNKYMyJFRO2iaNcuibiaHrWNpibY32Mqnib658lYic968s/640?wx_fmt=png&from=appmsg)      |
| `figure-05.png` | Hopper Ping-Pong 的 full async pipeline：TMA、SMEM、barrier 与 Tensor Core。来源：PyTorch Blog | [mmbiz.qpic.cn 图片地址](https://mmbiz.qpic.cn/mmbiz_png/QPNjE9VGkMUktD5L4cyaav2cX5vURjJkbH25y19F755BB9AP2icm4q6ibkJ3Z0ib6vpe2H5mXiaxRpo8CSzibPOibChhazKqcxSMtgHONWCfGGcbA/640?wx_fmt=png&from=appmsg)          |
| `figure-06.png` | CUTLASS Ping-Pong Architecture：1 个 producer，2 个 consumer。来源：PyTorch Blog               | [mmbiz.qpic.cn 图片地址](https://mmbiz.qpic.cn/sz_mmbiz_png/QPNjE9VGkMVMkRMR0osIniazHBG9DI9fgeeO2O2clUgZiabia92EQfDOe8IQuKAWceBJ7OCKW33PaZ88MjXwlBRJxPiakgfnU2icPiacmBENwgfys/640?wx_fmt=png&from=appmsg)       |
| `figure-07.png` | Blackwell Tensor Memory layout and addressing。来源：Colfax Research                           | [mmbiz.qpic.cn 图片地址](https://mmbiz.qpic.cn/mmbiz_png/QPNjE9VGkMXDoial6ImEvU6VZMUVTiaMtjxgiaJZ8qYXbd7DhehMDLeSjvcWy6h7a6DwL9HsNibmTMktiaRw8RmHGgbb0kyzuofmn3n0XBrTwYEE/640?wx_fmt=png&from=appmsg)           |
