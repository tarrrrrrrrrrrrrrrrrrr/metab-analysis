---
name: metab-analysis
description: Use when the user asks for metabolomics analysis, differential metabolite analysis, KEGG enrichment, or transcriptome-metabolome joint analysis, especially Chinese prompts such as 代谢组分析、代谢物差异、联合分析、转录组代谢组相关网络、基因-代谢物网络、Cytoscape节点边表、R≥0.9网络图、论文组学图.
---

# Metabolomics And Multi-Omics Analysis

This skill supports two related workflows:

1. **Metabolomics differential analysis** from an abundance matrix: PCA, OPLS-DA, differential metabolites, volcano plots, heatmaps, KEGG target extraction and KEGG bubble plots.
2. **Transcriptome-metabolome joint analysis** from existing DEG, metabolite matrix and RNA matrix files: Pearson correlation network, evidence tables, Cytoscape node/edge tables and a publication-style focused network figure.

Use it for paper-oriented omics workflows where the user expects reusable R code, auditable CSV outputs and figure files.

## Runtime

Prefer the user's local R if present:

```powershell
D:\RRR\R\R-4.4.1\bin\x64\Rscript.exe
```

If that path is absent, resolve R from a shortcut or `Get-Command Rscript`. Do not silently switch to another language for figure generation when the user requested R.

## Workflow A: Metabolomics Differential Analysis

Use this when the user has a metabolite abundance matrix such as `metab_abund_named.txt`.

Run the existing scripts in order:

| Step | Script | Purpose | Main outputs |
|---|---|---|---|
| 1 | `02_OPLSDA_analysis.R` | OPLS-DA, VIP, t-test and differential metabolite screening | `OPLSDA_Diff_Results_*.csv` |
| 2 | `01_PCA_analysis.R` | PCA and QC visualization | `01_PCA_*.pdf` |
| 3 | `03_Volcano_plot.R` | Volcano plots | `02_Volcano_*.pdf` |
| 4 | `04_Heatmap_plot.R` | Differential metabolite heatmaps | `03_Heatmap_*.pdf` |
| 5 | `05_KEGG_targets.R` | Extract significant KEGG compounds | `04_KEGG_Target_*.csv` |
| 6 | `06_KEGG_bubble.R` | KEGG enrichment and bubble plots | `05_KEGG_*.csv`, `05_KEGG_Bubble_*.pdf` |

One-command entry:

```powershell
& "D:\RRR\R\R-4.4.1\bin\x64\Rscript.exe" "scripts\run_all.R" "C:\path\to\data" dosa
```

Default organism is `dosa` for rice. Change it when the manuscript species is not rice.

## Workflow B: Transcriptome-Metabolome Joint Network

Use this when the user has already produced DEG and metabolite results and asks to redo a joint analysis, replace old paths/files, use `R=0.9`, create a core gene-metabolite network or prepare a paper figure.

### Expected Input Files

The input directory should contain these files, or equivalent files with the same columns:

| File | Required columns |
|---|---|
| `Final_Scaled_RNA_Matrix_Seedling_DXH.csv` | first column gene id, then sample columns such as `CK.1`, `CK.2`, `CK.3`, `T1.1`, `T1.2`, `T1.3` |
| `DESeq2_T1_vs_CK_Significant_DEGs.csv` | `Gene_ID`, `log2FoldChange`, `pvalue`, `padj`, `baseMean` |
| `DEMs_Matrix_Seedling_DXH_Final.csv` | `metab_id`, `Metabolite`, sample columns, `Log2FC_T1_vs_CK`, `Pval_T1_vs_CK`, `FDR_T1_vs_CK`, `VIP` |
| `DX-H_KEGG_Target_Unique_T1_vs_CK.csv` | `Metabolite`, `KEGG_ID`, `log2FC`, `Pvalue`, `VIP`, `Status`, `Mode` |

If file names differ, copy or adapt the script arguments rather than editing the source data.

### Reanalysis Command

```powershell
& "D:\RRR\R\R-4.4.1\bin\x64\Rscript.exe" `
  "scripts\07_multiomics_network.R" `
  "C:\path\to\input" `
  "C:\path\to\output" `
  0.90
```

This script:

- aligns shared RNA and metabolite sample columns;
- uses Pearson correlation across paired samples;
- filters edges with `abs(R) >= threshold` and `P < 0.05`;
- joins DEG statistics and metabolite statistics;
- prefers KEGG target-table `log2FC/Pvalue/VIP` for manuscript evidence when available;
- exports all correlations, filtered edges, evidence edges, Cytoscape nodes, focused target-axis edges and input matrices.

### Network Figure Command

```powershell
& "D:\RRR\R\R-4.4.1\bin\x64\Rscript.exe" `
  "scripts\08_plot_focused_network_old_layout.R" `
  "C:\path\to\output" `
  "SmartLeader"
```

The plotting script follows the old focused-network style:

- `ggraph + tidygraph`;
- force-directed layout;
- red/blue edges for positive/negative correlations;
- gene circles filled by gene log2FC with black outlines;
- metabolite diamonds in green with black outlines;
- label leader lines only where `ggrepel` judges they are helpful;
- PDF, SVG, TIFF and PNG preview exports.

Use `SmartLeader` for the recommended version. `NoLeader` is useful when labels are very close and the user wants no black guide lines.

## Evidence And Writing Rules

- Treat correlation networks as candidate co-variation evidence, not proof of regulation.
- State the exact sample columns used for correlation. If RNA has 3+3 samples but metabolomics has more replicates, only paired shared sample columns enter the network.
- Report both threshold and edge count, e.g. `abs(R) >= 0.90, P < 0.05`.
- Distinguish metabolite `Pvalue + VIP` evidence from strict global FDR evidence.
- Keep generated CSV edge/node tables alongside figures so the manuscript figure is auditable.

## Common Problems

| Problem | Fix |
|---|---|
| Old scripts still contain another project name such as `Spike` or `YT1` | Replace with command-line arguments and current input/output paths |
| `Rscript` not in PATH | Resolve from `D:\RRR\R\R-4.4.1\bin\x64\Rscript.exe` or the user's R shortcut |
| PowerShell expands `$Weight` or `$Target` inside inline R | Put verification code in a `.R` file or quote carefully |
| PDF export fails with Cairo stream error | The existing PDF may be open; export with a new suffix |
| Too many label leader lines | Use `SmartLeader` or `NoLeader` mode in the plotting script |

