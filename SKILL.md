---
name: metab-analysis
description: 代谢组学差异分析全套流程。输入丰度矩阵数据目录，自动运行 PCA → OPLS-DA → 火山图 → 聚类热图 → KEGG 通路富集。触发词：分析代谢组、代谢组分析、代谢组学分析、代谢组差异、代谢物差异、代谢组数据处理、代谢组重新分析、跑代谢组、做代谢组。
---

# 代谢组学差异分析

一键运行 PCA、OPLS-DA、火山图、聚类热图、KEGG 通路富集全套分析，输出表格 + PDF。

## 触发条件

用户说"分析代谢组"、"代谢组分析"、"代谢组学分析"、"代谢组差异分析"、"代谢物差异"、"代谢组数据处理"、"代谢组重新分析"、"跑代谢组"、"做代谢组"等，且指明了数据目录时自动触发。

> 如果你的数据之前已跑过代谢组分析流程（如 MetaboAnalyst、XCMS、Biotree 等），想用 R 重新分析，同样触发本技能。

## 前置：找到 R

R 4.4.1 安装在 `D:/RRR/R/R-4.4.1/`。调用方式：

```bash
"D:/RRR/R/R-4.4.1/bin/Rscript.exe" -e "source('<script>.R')"
```

或在 R 脚本中首行写入 `setwd("<数据目录>")` 后直接执行：

```bash
"D:/RRR/R/R-4.4.1/bin/Rscript.exe" "<脚本路径>"
```

如果你的 R 路径不同，用 `where Rscript` 或搜索 `R-4*` 目录找到它，然后以此为准。

## 分析流程

### Step 0 — 确认输入

向用户确认：
- **数据目录**：包含 `metab_abund_named.txt`（或名称为 `metab_abund.txt` 等 tab 分隔丰度矩阵）的目录路径
- **对比组**：默认自动提取列名前缀（如 `WT` vs `KO_2`、`KO_6`），若用户指定则使用用户指定
- **输出目录**：默认同数据目录

### Step 1 — 复制脚本到数据目录

将 `scripts/` 下的 4 个 R 脚本复制到目标数据目录（或用 `setwd()` 直接指向数据目录，读脚本时加绝对路径）。

### Step 2 — 依次运行

按顺序运行以下 4 个脚本，每一步完成后检查是否有错误再继续：

| 顺序 | 脚本 | 功能 | 输出 | 依赖 |
|---|---|---|---|---|
| 1 | `02_OPLSDA_analysis.R` | OPLS-DA 建模 + t检验 + 差异筛选 | `OPLSDA_Diff_Results_*.csv` | — |
| 2 | `01_PCA_analysis.R` | PCA 降维（分组 + 三组合并） | `01_PCA_*.pdf` | — |
| 3 | `03_Volcano_plot.R` | 火山图 | `02_Volcano_*.pdf` | Step1 |
| 4 | `04_Heatmap_plot.R` | 差异代谢物聚类热图 | `03_Heatmap_*.pdf` | Step1 |
| 5 | `05_KEGG_targets.R` | 提取显著差异代谢物 KEGG ID | `04_KEGG_Target_*.csv` | Step1 |
| 6 | `06_KEGG_bubble.R` | KEGG 通路富集 + Top10 气泡图 | `05_KEGG_*.csv`, `05_KEGG_Bubble_*.pdf` | Step5 |

> Step 1 必须先跑。Step 2-4 可并行。Step 5-6 串行（5→6）。Step 6 需联网获取 KEGG 数据库。

### Step 3 — 报告结果

输出关键数字：每个对比组的差异代谢物数量（上调/下调）、模型 Q²，并列出生成的文件。

## 脚本适配规则

当数据格式与原脚本假设不一致时，按以下规则适配：

### 列名识别
- 对照组：`WT`、`CK`、`Control`、`NC`
- 处理组：`KO_`、`T`（如 `T1`、`T2`）、`Treat`
- QC：`QC`
- 列名模式：`<组名>_<编号>`（如 `WT_1`、`KO_2_3`）

### 数据格式
- **tab 分隔** → `read.table(..., sep="\t")`
- **csv 逗号分隔** → `read.csv(...)`
- 第 1 列：代谢物 ID（行名）
- 其余列：样本丰度值

### 关键阈值
- FC cutoff：1.5 倍（即 `log10(1.5) ≈ 0.176`）
- P-value cutoff：0.05
- VIP cutoff：1.0
- OPLS-DA：predI=1, orthoI=1, permI=200

### 常见问题处理
1. **`comment.char` 问题**：SMILES 中 `#` 被误作注释 → 读文件时加 `comment.char = ""`
2. **VIP 向量长度不匹配**：`getVipVn()` 排除近零方差变量 → 用 `NA` 填充缺失项再替换为 0
3. **PCA 零方差**：`prcomp` 报错 → `apply(x, 1, var) > 1e-10` 过滤
4. **椭圆点数不足**：减少 `stat_ellipse` 或忽略 warning
5. **R 未安装** → 提示用户安装或自行安装（优先用已有的 `D:/RRR/R/R-4.4.1/`）

## KEGG 通路富集（Step 5-6）

Step 6 需要联网获取 KEGG 数据库。**物种代码**可通过命令行参数传入：

```bash
Rscript 06_KEGG_bubble.R <数据目录> [organism_code]
```

默认物种为 **dosa**（水稻 Oryza sativa japonica）。常见物种代码：

| 代码 | 物种 |
|---|---|
| `dosa` | 水稻 (Oryza sativa japonica) |
| `hsa` | 人类 (Homo sapiens) |
| `mmu` | 小鼠 (Mus musculus) |
| `ath` | 拟南芥 (Arabidopsis thaliana) |

气泡图取 p.adjust 最小的 Top 10 通路，全量结果保存为 CSV。不设 p.adjust 硬阈值，保留所有富集通路供查看。
