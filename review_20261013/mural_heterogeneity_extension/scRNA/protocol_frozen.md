# GSE183320 壁细胞异质性：结果生成前冻结的事后扩展协议

固定日期：2026-10-13（上海）。已有壁细胞总体 PDE5A 分布已知；下述亚群/状态/轨迹结果尚未生成。

1. 使用作者公开的 normalized expression，保留其数值，**不再 log、不反推原始 UMI**；仅作者注释 SMC 与 pericyte。分别建图，避免把两类细胞的差异写成去分化。
2. 预先列定四个描述性状态 panel：contractile、matrix_fibrosis、SASP_associated、inflammatory。它们是研究者选择的先验 panel，不是完整通路活性指标；SASP panel 不能单独鉴定衰老细胞。PDE5A 不进入任何 panel。
3. 每类细胞基因至少 10 个细胞检出。HVG 使用 Scanpy seurat flavor（适于 log-normalized data），最多 2,000，按 donor 分层选取；排除 PDE5A、全部 panel 基因、线粒体/核糖体基因。PCA 30 维、kNN 15、Leiden 主分辨率 1.5，敏感性 0.8/1.2/2.0，seed=20261013。缩放上限 10；不以 PDE5A、疾病标签或显著性选择分辨率。
4. 亚群只用中性编号；每群至少 50 个细胞、至少 5 位各贡献 >=10 个细胞的供体、最大单供体占比 <=40%，才标记 donor-supported。各群仍全部保留、输出供体构成。
5. marker 是 **cluster-vs-rest 的供体配对 log-normalized均值差**，不是 count pseudobulk。每位供体至少 cluster 10、rest 20 个细胞，至少 5 位有效供体；用双侧 Wilcoxon 配对秩检验，BH 校正全部 cluster×gene 探索。输出全部 marker，不能将同一数据上定义后检验的亚群当成独立验证。
6. 状态分数是 Scanpy expression-matched control score（ctrl_size=50, n_bins=25）；输出真实 panel 覆盖。全基因非负 normalized 表达保留；不将缺失基因补零计为已测。
7. DPT 仅针对每种细胞最大连通分量；若不足该类 80% 则暂停该类轨迹。使用 diffusion map 15 维；root 为 contractile 得分最高 10% 细胞中离其 PCA 质心最近者。另用该 top 10% 中第二近的 root 作预定根敏感性。PDE5A 及全部状态panel已从邻接图排除。方向表示从较高收缩状态出发的转录差异，**不是时间、细胞年龄、衰老或去分化的证明**。
8. 统计单位为 donor。每 donor/type >=30 个有 DPT 细胞且 >=5 PDE5A检出，计算 PDE5A 与 DPT/四状态分数的 Spearman 以及控制 log1p(nCount_RNA) 的偏秩相关。供体等权平均 rho，用 9,999 次双侧 donor sign-flip 作探索检验（要求供体效应围绕零对称），BH family 为 2 types×5 predictors；报告每 donor 值、有效 n、95% t 区间、健康限制、留一供体区间。主分析为深度调整 rho；未调整/替代根为敏感性。不把相关 p 值当作轨迹方向验证。
9. 无 raw counts：不运行 count-based differential expression；无 spliced/unspliced：不运行 RNA velocity；无纵向数据：不声称衰老末端。健康限制复用冻结的图，作为表型限制，不宣称独立重建。
10. 不改变既有供体重叠审计。GSE135922 的非共享 3 位供体只能支持壁细胞背景，不作为本次高分辨率亚群的自动独立重复。
