# GSE135922：PDE5A壁细胞定位的供体与区域复核

## 本次结论

在按实验纠正cluster编号后，975个mural-like细胞来自7位供体，342个细胞的原作者归一化PDE5A表达大于0。合并每位供体黄斑与周边的阳性细胞数/总细胞数后，供体未加权均值±样本SD为 **31.47±10.33%（n=7）**。无记录眼病限制为 **29.30±6.03%（n=5）**。这里的“无记录眼病”沿用GEO文字，不等于无组织学异常的健康眼。

7位供体的壁细胞PDE5A检测比例为25.00%、24.53%、22.54%、51.25%、35.20%、36.55%、25.23%；每位有53–240个壁细胞。阈值>=1/5/10/20/30/50均保留7位供体，故供体合并区域的结果完全相同。留一供体均值为28.18–32.96%；去掉最高的donor4后仍为28.18%。这支持该定位不完全由一个供体贡献，不构成疾病、药物响应或患者易感性验证。

## 原有汇总中的实质性问题

GSE135922含两个实验，其cluster标签不是同一个命名空间：

- 非富集donor1–3：11个cluster，mural-like是cluster5，180 cells。
- CD31富集donor4–7：13个cluster，mural-like是cluster1，795 cells。

原`R/08_singlecell_GSE135922.R`按相同cluster整数直接合并两个实验。例如旧cluster1的828 cells为富集实验795个mural-like加上非富集实验33个其他细胞；旧cluster5则混入大量富集实验内皮。不能将旧表中的同号cluster作为整个队列统一细胞类型。原文件未被覆盖。

本次先保存`fixed_protocol.json`，再仅查看RGS5、ACTA2、PDGFRB、NOTCH3及辅助标志物，按实验独立确定mural-like，并在`annotation.json`冻结后才统计新的PDE5A结果。两实验所选cluster的四个mural marker各自均排第一。没有按PDE5A值挑cluster，也没有把mural-like强拆为SMC和pericyte。PECAM1在全部14个源文件的标题中缺失，在初始表头审核后从辅助panel删除；未填补为0。

这是此前已看过项目汇总结果后的回顾性复核。规则先于本次新供体分析冻结，不是前瞻性预注册。

## 分层与区域结果

| 分层 | 供体n | 壁细胞n | 供体均值±SD阳性率 |
|---|---:|---:|---:|
| 全部 | 7 | 975 | 31.47±10.33% |
| 无记录眼病 | 5 | 664 | 29.30±6.03% |
| 非富集实验 | 3 | 180 | 24.02±1.31% |
| CD31富集实验 | 4 | 795 | 37.06±10.72% |

两实验处理流程、取样及供体组成不同；富集组较高不能作为富集处理的效应或AMD差异。GEO记录donor3/4为AMD；同行复核取得的原论文Table1进一步说明为neovascular AMD。其余1/2/5/6/7符合“无记录眼病”限制，但原文亦记载部分组织学沉积，不能统称完全健康。

按每个donor-region>=20 cells：黄斑6供体（538 cells）为34.38±12.02%，周边5供体（387 cells）为30.56±9.34%。这些两个均值来自不完全相同的供体集合，不能直接用来断言区域差异。两区域都符合阈值的只有donor4–7（4对、全部为CD31富集实验），黄斑减周边的差值为5.48±6.48个百分点，仅作配对描述。完整阈值敏感性与配对差值另见CSV。区域不是额外独立供体。

## 跨数据集边界

GEO样本级`Reanalysis of`已确认GSE135922原donor4–7的黄斑样本被GSE183320复用，因此两个数据集不能被称为完全独立复制。具体对应与来源证据见同级`../metadata_audit/documented_sample_reuse.csv`。在主分析完成后，依据新增元数据排除全部明确复用donor4–7，剩余donor1–3仍为 **24.02±1.31%（n=3，180 cells）**，三位均有PDE5A阳性壁细胞。这三位属于同一非富集实验；结果不能外推为大型独立验证队列。

跨论文的donor5/6年龄、性别及后来AMD标签存在冲突。原GEO“无记录眼病”分析仍保留原始口径n=5，但新增保守敏感性同时排除donor5/6后，保留donor1/2/7，结果为 **24.92±0.36%（n=3，323 cells）**。这项附加分析由元数据问题触发，不是依据PDE5A值选择性排除；原冻结方案保持不动，新增规则、输入SHA256和结果见`metadata_addendum.md`及`audit_metadata_sensitivity.json`。

这些非CSCR、未建立他达拉非暴露的组织只提供PDE5A的细胞定位背景。检测比例受捕获效率、测序深度、细胞回收和原作者处理影响；表达>0不是蛋白、药理抑制、局部空间接触或通透性的测量。区域差异不等于基于物理坐标的空间转录组分析。

## 文件与复现

- `PDE5A_mural_donor_pooled_regions.csv`：画主供体点图；每位供体的阳性细胞分子与总细胞分母、归一化表达均值均可追溯。
- `PDE5A_mural_donor_summary_thresholds.csv`：主结果、限制与>=1/5/10/20/30/50阈值。
- `PDE5A_mural_donor_region.csv`：每位供体每区域的数据点，含不满足主阈值的点并明确标注，不删除。
- `PDE5A_mural_region_summary_thresholds.csv`、`PDE5A_mural_paired_region_differences.csv`：区域与配对结果。
- `PDE5A_mural_leave_one_donor_out.csv`：留一供体敏感性。
- `marker_cluster_means_by_experiment.csv`、`annotation.json`：不以PDE5A决定细胞注释的证据。
- `audit_markers.json`、`audit_quantify.json`：输入SHA256、规则/脚本SHA256、软件与运行时间。
- `metadata_sensitivity_summary.csv`、`metadata_sensitivity_donor_points.csv`、`metadata_sensitivity_leave_one_out.csv`：排除明确复用供体及排除临床标签冲突供体的附加敏感性。

在项目根目录使用已有Python环境运行：

```powershell
& 'C:/Users/18328/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' 'review_20261013/mural_PDE5A_validation/GSE135922/run_donor_region_validation.py' --stage markers
& 'C:/Users/18328/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' 'review_20261013/mural_PDE5A_validation/GSE135922/run_donor_region_validation.py' --stage quantify
& 'C:/Users/18328/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' 'review_20261013/mural_PDE5A_validation/GSE135922/run_metadata_sensitivity.py'
```

首次运行前应保留并审核`annotation.json`；脚本不会自行根据PDE5A修改注释。此次marker阶段约27.66秒；最终quantify阶段约27.57秒，Windows峰值工作集约403.64 MiB。无额外包安装，无原始文件或旧分析覆盖。原文件扩展名为`.tsv.gz`，实际为含引号的空格分隔R表并带额外行号；脚本显式补足行号列，避免`usecols`引起列错位。

## 可用英文表述

We reassessed PDE5A localization in GSE135922 while retaining donor, anatomical region and experimental processing. Because the unenriched and CD31-enriched experiments used different cluster labels, mural-like populations were identified separately using RGS5, ACTA2, PDGFRB and NOTCH3 before calculating the new PDE5A summaries. The selected populations corresponded to cluster 5 in the unenriched experiment and cluster 1 in the enriched experiment. Macular and peripheral counts were pooled within each donor by summing positive cells and total cells. With a minimum of 20 cells per donor, all seven donors were eligible. PDE5A positivity averaged 31.47% (sample SD, 10.33%; 975 cells), and the corresponding estimate among five donors with no noted ophthalmologic disease was 29.30% (SD, 6.03%; 664 cells). Results were unchanged across donor-level minimum-cell thresholds of 1, 5, 10, 30 and 50; leave-one-donor-out means ranged from 28.18% to 32.96%. These findings support reproducible mural-cell localization within this contextual dataset. They do not establish tadalafil responsiveness, and reuse of samples in GSE183320 prevents treating the two accessions as wholly independent replication.

Metadata review identified reuse of macular samples from four GSE135922 donors in GSE183320. Excluding all four shared donors retained three donors with a mean PDE5A positivity of 24.02% (SD, 1.31%; 180 cells). A separate conservative sensitivity analysis excluded both donor labels affected by cross-publication clinical-label uncertainty from the original no-noted-disease subset; positivity was 24.92% (SD, 0.36%; three donors and 323 cells). These metadata-driven analyses were added after the primary summary and are reported as sensitivity analyses rather than prespecified confirmatory tests.
