# PDE5A 壁细胞配体–受体分析：事后扩展固定方案

冻结日期：2026-10-13；冻结时尚未计算新 ligand–receptor 结果。本方案是已知 PDE5A 壁细胞定位后的事后扩展，不是预注册或验证性研究。

## 输入与边界

- GSE183320 作者 RNA log-normalized 表达矩阵及原供体、原细胞注释。不使用 integrated assay、scale.data 或推造的 raw counts。
- SMC 与 pericyte 分开，主要接收群为 endothelial（artery/vein/choriocapillaris/ec_unknown 合并）及 macrophage（macrophage-inf/macrophage-res 合并）。供体始终为独立推断单位。
- 使用 LIANA 1.10.0 的 consensus 全资源；导出原始资源并记录 SHA256。PDE5A 不作为配体–受体候选基因或打分输入。不能把酶 PDE5A 画成分泌配体。
- 作者 log-normalized 数据不再次归一化。不使用有限候选基因代替全测量基因背景。

## 固定分组与最低数量

- 主分析：在每 donor × 原壁细胞类型内，PDE5A>0 的细胞以其非零中位数划分；high 为严格大于中位数，low_positive 为 0<表达≤中位数。ties 全部放入 low_positive，不能随机打散或为了凑 n 改阈值。
- 敏感性：同一 high 与 PDE5A=0（undetected）比较。未检出不称为生物学阴性或真正低表达。
- 每 donor × 类型的 high、比较组各≥10 细胞，每同供体接收细胞群≥20。否则不可评价，不记为无通讯或零。
- 独立供体≥5 才进行供体层面的差值检验。每组细胞数、非零中位数、供体身份、疾病状态和 nCount/nFeature 差异全部保留。

## LIANA 与供体层面推断

- 用 LIANA geometric_mean 获取每个供体每 sender–receiver 的表达支持评分；固定 expr_prop=0.10、min_cells=10、use_raw=False、return_all_lrs=True。全部复合体亚单位必须可检测；记录未通过比例门槛者，不能当作零。
- 配对比较同一供体、同一原壁细胞类型、同一 receiver 的同一 ligand–receptor。只有两组都通过表达支持门槛才进入该供体配对。禁止把低组缺失值填零制造差异。
- 对每个可评价配对差值 geometric_mean(high)−geometric_mean(comparator)做双侧供体符号翻转检验；n≤15 枚举全部 2^n 符号，n>15 用 99,999 次固定种子 Monte Carlo。均值效应、样本SD、供体 n、正方向供体比例同时报告；95%不确定区间通过供体重采样 9,999 次（固定种子）。
- BH 校正覆盖主分析所有可评价 SMC/pericyte × endothelial/macrophage × 全 consensus LR；敏感性单独成家族。绝不把排名或 cell-level permutation P 当作供体重复。
- 健康供体限制、剔除明确共享供体1/2/11/12、留一供体方向稳定性为固定敏感性。小于5供体报告不可评价。不能用敏感性较小P替代主分析。
- 深度差异先报告供体内 high−comparator 的 log1p(nCount_RNA)/log1p(nFeature_RNA)标准化差异；结果须承认捕获深度混杂。若高低分组数量不足以完成有效匹配，不能写成已排除深度因素。
- 预先固定深度敏感性：同 donor×壁细胞类型×比较内，以两组 pooled nFeature_RNA 的四分位箱分层，各箱两组抽取相同数量 min(n_high,n_comparator)，固定种子20261013抽样一次；匹配后各组仍≥10。该单次固定抽样仅为技术敏感性，结果不得冒充全面批次校正。对匹配组重新调用 LIANA；两项比较及各自深度匹配分别校正，保留全部结果。
- LIANA 的细胞层置换不运行（n_perms=None）；由供体配对检验提供不确定性。两个高低sender共享同一receiver，因此差异主要来自sender配体表达，不能写成接收细胞发生了真实反应。

## 解释规则

- 候选通讯需 BH q<0.05、均值差值>0、至少70%供体方向一致；仍仅称 transcript-supported candidate interaction。
- TGF-beta、PDGF、Notch、血管活性等只能在全资源结果之后作生物学归纳；禁止先筛热门通路再只对其校正。所有可评价结果完整导出。
- PDE5A 高低表达不是他达拉非暴露或抑制实验；候选 LR 共表达不能建立传递方向、蛋白结合或功能激活。
- 已完成 Xenium 主要邻近分析无 BH 显著正富集，应在叙事中保留。独立 RNA 结果不能把该空间阴性覆盖成空间验证。

## 实现校正记录（2026-10-13，首次独立评分QA后）

原计划采用 return_all_lrs=True 以保留未通过门槛者。独立从 RNA 矩阵复算发现 LIANA 1.10.0 的该分支在复合体最小亚单位汇总前去重，导致部分多亚单位评分不符（例如 PGF→FLT1_KDR；单体吻合）。最终实现改为 return_all_lrs=False，正确组装复合体后手工与完整资源左连接，缺失标为不可评价、评分NaN。所有分组、阈值、资源、效应和检验规则不变；初次评分及统计结果作废，最终结果以修复后独立评分QA通过的版本为准。
