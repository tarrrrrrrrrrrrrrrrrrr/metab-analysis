---
name: metab-analysis
description: Use when the user asks for metabolomics analysis, differential metabolite analysis, KEGG enrichment, or transcriptome-metabolome joint analysis, including prompts about daixiezufenxi, daixiewu chayi, lianhe fenxi, gene-metabolite networks, Cytoscape node/edge tables, R>=0.9 correlation networks, or paper-ready omics figures.
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

Use this when the user has already produced DEG and metabolite results, or has a narrative analysis text (联合分析思路/*.txt) with named metabolites and axes, and asks to build a gene-metabolite correlation network.

### Complete Step-by-Step SOP

**Do NOT run 07/08 scripts directly.** This workflow requires per-project adaptation. Follow every step:

#### B.0 — Read the narrative file first

If the user provides a narrative text (联合分析思路/*.txt), parse:
- Each axis name and its storyline
- All metabolite names listed under each axis
- Any hub gene IDs mentioned

This becomes the `axis_map` in Step B.3.

#### B.1 — Survey the data directory

```bash
ls "path/to/data"
```

Identify these files (names vary by project):

| Role | Look for | Key columns check |
|---|---|---|
| RNA expression matrix | `Final_Scaled_RNA_Matrix_*.csv` | Gene ID col + sample cols (CK.1, CK.2, CK.3, T2.1, T2.2, T2.3) |
| DEG statistics | `Treat2_Significant_Genes.csv` or `*_Significant_DEGs.csv` | `Gene_ID`, `log2FoldChange`, `pvalue`, `padj` |
| Metabolite abundance | ⚠️ **ALWAYS use POS_Diff + NEG_Diff full data** (see B.2) | Sample cols + `Log2FC_T2_vs_CK`, `Pval_T2_vs_CK`, `FDR_T2_vs_CK` |
| KEGG target table | `04_KEGG_Target_Unique_*.csv` | `Metabolite`, `KEGG_ID`, `log2FC`, `Pvalue`, `VIP` |

**CRITICAL: Check sample column naming.** RNA matrices often use hyphens (`CK-1`, `CK-2`), while metabolomics files use dots (`CK.1`, `CK.2`). Normalize to dots before analysis: `names(rna) <- gsub("-", ".", names(rna), fixed = TRUE)`

Extract exact sample column names with:
```bash
head -1 "file.csv" | tr ',' '\n' | grep -E "CK|T[0-9]"
```

Shared samples are typically 3+3 (e.g. `CK.1, CK.2, CK.3, T2.1, T2.2, T2.3`). Do NOT include extra replicates (CK.4-CK.6, T2.4-T2.6) in the correlation.

#### B.2 — Build the FULL DEM matrix (MANDATORY)

⚠️ **Never use a pre-filtered "Ready" or "Optimized" DEM file** (e.g. `DEMs_Matrix_*_Ready.csv`, `DEMs_Matrix_*_Optimized.csv`). These are curated subsets (often 35-83 metabolites) and will exclude narrative metabolites that exist in the raw data, causing the network to miss key nodes.

**Always merge the full POS_Diff and NEG_Diff files:**

```r
# Identify the full abundance files — typically named like:
#   POS_Diff_CK_vs_T2_6.csv, NEG_Diff_CK_vs_T2_6.csv
# or with project-specific suffixes like _Leaf_SCI.csv

col_shared <- c("metab_id", "Metabolite", samples, 
                "Log2FC_T2_vs_CK", "Pval_T2_vs_CK", "FDR_T2_vs_CK")

load_diff <- function(fname) {
  df <- read.csv(file.path(input_dir, fname), check.names = FALSE, stringsAsFactors = FALSE)
  for (mc in setdiff(col_shared, names(df))) df[[mc]] <- NA_real_
  df %>% select(any_of(col_shared))
}

pos_dem <- load_diff("POS_Diff_CK_vs_T2_6.csv")
neg_dem <- load_diff("NEG_Diff_CK_vs_T2_6.csv")

dem <- bind_rows(pos_dem, neg_dem) %>%
  distinct(Metabolite, .keep_all = TRUE)
```

This yields ~1,200 metabolites instead of 35-83.

#### B.3 — Build the axis_map and VERIFY

```r
axis_map <- tibble::tribble(
  ~Target,                      ~Narrative_Axis,
  "Plantamajoside",             "Axis 1: Phenylpropanoid & Lignin",
  "Eugenyl Acetate",            "Axis 1: Phenylpropanoid & Lignin",
  "Uric Acid",                  "Axis 2: Glutathione & Redox",
  # ... every metabolite from the narrative file
)
```

**Then verify every metabolite exists in DEM:**

```r
in_dem <- axis_map$Target %in% dem$Metabolite
for (i in seq_len(nrow(axis_map))) {
  cat(sprintf("  %-35s %s  %s\n", axis_map$Target[i], 
              ifelse(in_dem[i], "✅ IN DEM", "❌ MISSING"), 
              axis_map$Narrative_Axis[i]))
}
```

Report coverage to the user. Metabolites marked ❌ exist in KEGG target table but lack abundance values — they cannot enter the network. Flag them explicitly so the user can adjust the narrative or supplement data.

#### B.4 — Write the adapted 07 script

Write a project-specific `07_multiomics_network_<PROJ>.R` into the output directory. Template:

```r
# Key parameters to adapt:
# - input_dir, output_dir
# - File names (RNA, DEG, POS_Diff, NEG_Diff, KEGG)
# - Sample vector (check exact names, normalize dots)
# - axis_map content
# - prefix for output file names
# - DEG column mapping (Gene ID vs Gene_ID, etc.)

# After axis_map, log coverage:
cat("\n=== Axis metabolite coverage ===\n")
for (i in seq_len(nrow(axis_map))) {
  cat(sprintf("  %-40s %s  %s\n", axis_map$Target[i], 
              ifelse(in_dem[i], "✅", "❌"), axis_map$Narrative_Axis[i]))
}

# After target_edges, log per-axis counts:
cat(sprintf("\nRNA: %d | DEM: %d | Pairs: %d | |R|>=%.2f: %d | Target: %d\n", ...))
for (ax in unique(target_edges$Narrative_Axis)) {
  ae <- filter(target_edges, Narrative_Axis == ax)
  cat(sprintf("  %s: %d edges (%s)\n", ax, nrow(ae), paste(unique(ae$Target), collapse=", ")))
}
```

#### B.5 — Run the 07 script

```bash
"D:/RRR/R/R-4.4.1/bin/Rscript.exe" "path/to/07_multiomics_network_<PROJ>.R" \
  "path/to/input" "path/to/output" 0.90
```

Default threshold: `|R| >= 0.90, P < 0.05`. If user asks for a different threshold, change the last argument (e.g. `0.85`).

#### B.6 — Write and run the 08 network plot script

Write a project-specific `08_plot_focused_network_<PROJ>.R` with:

```r
# Adapt:
# - output_dir, prefix
# - short_label() function to alias long metabolite names
# - plot title (include project context: T2 vs CK, Root vs Leaf, etc.)
```

```bash
"D:/RRR/R/R-4.4.1/bin/Rscript.exe" "path/to/08_plot_focused_network_<PROJ>.R" \
  "path/to/output" "Compact"
```

If PDF is locked (cairo error), use a suffix like `_v2` and rename after.

#### B.7 — Clean up

- Delete old/intermediate network outputs (older suffixes, v2, v3)
- Rename final versions to the standard name (no suffix)
- Keep the project-specific 07/08 scripts for reproducibility

### Network Figure Defaults

- `ggraph + tidygraph`, `layout = "fr"`, niter=5000, seed=123
- Edge: `scale_edge_width_continuous(range = c(0.12, 0.55))`, alpha 0.42
- Gene node: circle `shape = 21`, black stroke 0.24, fill = Gene Log2FC (blue-white-red gradient)
- Metabolite node: dark-green square `shape = 22`, black stroke 0.32, fill = `#276231`
- Labels: `geom_node_text(repel = TRUE, size = 2.35, fontface = "bold.italic")`, segment.size = 0.11, segment.alpha = 0.48
- Output: PDF (cairo_pdf), SVG (svglite), TIFF (ragg 600dpi LZW), PNG preview (ragg 180dpi)
- Plot size: Compact mode 10×9 inches
- `Top20 per metabolite` by |R| — keeps the network readable

### Gene Uniqueness

Never duplicate a gene. A shared gene appears once and connects to all its retained metabolites. When node tables have real-name labels, enforce uniqueness by displayed label, not raw ID. Use `GENE::` prefix and `group_by(Source, Target)` to merge duplicates. This prevents cases like two `5_8S_rRNA` labels appearing as separate nodes.

### Evidence And Writing Rules

- Correlation networks are candidate co-variation evidence, not proof of regulation.
- State exact sample columns used. Use only paired shared columns (3+3, not 3+6).
- Report both threshold and edge count: `Pearson |R| >= 0.90, P < 0.05, n = 6 samples (df = 4)`.
- Distinguish metabolite `Pvalue + VIP` evidence from strict global FDR evidence.
- Keep CSV edge/node tables alongside figures for auditability.

### Common Problems

| Problem | Fix |
|---|---|
| Narrative metabolites missing from network | Check DEM source — use POS_Diff + NEG_Diff full, not "Ready" subset |
| RNA and DEM sample names don't match | Normalize: `names(rna) <- gsub("-", ".", names(rna), fixed = TRUE)` |
| "Ready" DEM has only 35-83 metabolites | It's a pre-filtered subset; always merge POS_Diff + NEG_Diff |
| KEGG table has metabolites absent from abundance | Flag as ❌ in axis_map coverage report; can't enter correlation |
| DEG column "Gene ID" (space) vs "Gene_ID" (underscore) | Rename: `if ("Gene ID" %in% names(deg)) names(deg)[names(deg) == "Gene ID"] <- "Gene_ID"` |
| PDF export fails: cairo error "writing to output stream" | File is open/locked; use a suffix `_v2` then rename after |
| `Rscript` not in PATH | Use `D:\RRR\R\R-4.4.1\bin\x64\Rscript.exe` or resolve from shortcut |
| Only 2-3 axes appear when narrative has 4 | Metabolites for missing axes are in KEGG table only (not DEM); report to user |
