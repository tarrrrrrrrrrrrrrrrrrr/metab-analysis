suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(reshape2)
})

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1) args[1] else getwd()
output_dir <- if (length(args) >= 2) args[2] else input_dir
cor_threshold <- if (length(args) >= 3) as.numeric(args[3]) else 0.90
sample_arg <- if (length(args) >= 4) args[4] else "CK.1,CK.2,CK.3,T1.1,T1.2,T1.3"
samples <- trimws(strsplit(sample_arg, ",")[[1]])

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_required <- function(file_name, ...) {
  path <- file.path(input_dir, file_name)
  if (!file.exists(path)) stop("Missing required file: ", path)
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, ...)
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

rna <- read_required("Final_Scaled_RNA_Matrix_Seedling_DXH.csv")
if (names(rna)[1] == "" || grepl("^Unnamed", names(rna)[1])) names(rna)[1] <- "Gene_ID"
deg <- read_required("DESeq2_T1_vs_CK_Significant_DEGs.csv")
dem <- read_required("DEMs_Matrix_Seedling_DXH_Final.csv")
kegg <- read_required("DX-H_KEGG_Target_Unique_T1_vs_CK.csv")

missing_rna <- setdiff(samples, names(rna))
missing_dem <- setdiff(samples, names(dem))
if (length(missing_rna) > 0) stop("RNA matrix missing samples: ", paste(missing_rna, collapse = ", "))
if (length(missing_dem) > 0) stop("DEM matrix missing samples: ", paste(missing_dem, collapse = ", "))

rna_mat <- rna %>%
  select(Gene_ID, all_of(samples)) %>%
  distinct(Gene_ID, .keep_all = TRUE)
rna_values <- as.matrix(rna_mat[, samples])
mode(rna_values) <- "numeric"
rownames(rna_values) <- rna_mat$Gene_ID

dem_mat <- dem %>%
  select(metab_id, Metabolite, all_of(samples), Log2FC_T1_vs_CK, Pval_T1_vs_CK, FDR_T1_vs_CK, VIP) %>%
  distinct(Metabolite, .keep_all = TRUE)
dem_values_raw <- as.matrix(dem_mat[, samples])
mode(dem_values_raw) <- "numeric"
rownames(dem_values_raw) <- dem_mat$Metabolite
dem_values_scaled <- t(scale(t(dem_values_raw)))
dem_values_scaled[is.na(dem_values_scaled)] <- 0

r_mat <- cor(t(rna_values), t(dem_values_scaled), method = "pearson", use = "pairwise.complete.obs")
n_mat <- matrix(length(samples), nrow = nrow(r_mat), ncol = ncol(r_mat), dimnames = dimnames(r_mat))
t_stat <- r_mat * sqrt((n_mat - 2) / pmax(1 - r_mat^2, 1e-15))
p_mat <- 2 * pt(-abs(t_stat), df = n_mat - 2)

all_edges <- melt(r_mat, varnames = c("Source", "Target"), value.name = "Weight") %>%
  mutate(
    Pvalue = as.vector(p_mat),
    Abs_Weight = abs(Weight),
    Direction = ifelse(Weight >= 0, "Positive", "Negative")
  ) %>%
  arrange(desc(Abs_Weight), Pvalue)

edge_filtered <- all_edges %>%
  filter(Abs_Weight >= cor_threshold, Pvalue < 0.05) %>%
  arrange(desc(Abs_Weight), Pvalue)

deg_info <- deg %>%
  transmute(
    Source = Gene_ID,
    Gene_Log2FC = safe_num(log2FoldChange),
    Gene_pvalue = safe_num(pvalue),
    Gene_padj = safe_num(padj),
    Gene_baseMean = safe_num(baseMean)
  )

dem_info <- dem_mat %>%
  transmute(
    Target = Metabolite,
    metab_id = metab_id,
    Metab_Log2FC_DEM = safe_num(Log2FC_T1_vs_CK),
    Metab_Pvalue_DEM = safe_num(Pval_T1_vs_CK),
    Metab_FDR = safe_num(FDR_T1_vs_CK),
    Metab_VIP = safe_num(VIP)
  )

kegg_info <- kegg %>%
  transmute(
    Target = Metabolite,
    KEGG_ID = KEGG_ID,
    Metab_Log2FC_KEGG = safe_num(log2FC),
    Metab_Pvalue_KEGG = safe_num(Pvalue),
    Metab_VIP_KEGG = safe_num(VIP),
    Metab_Status_KEGG = Status,
    Metab_Mode_KEGG = Mode
  )

axis_map <- tibble::tribble(
  ~Target, ~Narrative_Axis,
  "1-O-Feruloyl-Beta-D-Glucose", "Axis 1: Phenolic glycosylation defense reserve",
  "Caffeic Acid 3-Glucoside", "Axis 1: Phenolic glycosylation defense reserve",
  "Adenosine Monophosphate", "Axis 2: AMP energy-cost and membrane protection",
  "4-Hydroxybenzoate", "Axis 3: Substrate extraction and phenolic precursor flux",
  "D-Sorbitol", "Axis 5: Carbon fixation and sugar reserve adjustment"
)

edge_evidence <- edge_filtered %>%
  left_join(deg_info, by = "Source") %>%
  left_join(dem_info, by = "Target") %>%
  left_join(kegg_info, by = "Target") %>%
  left_join(axis_map, by = "Target") %>%
  mutate(
    Metab_Log2FC = ifelse(!is.na(Metab_Log2FC_KEGG), Metab_Log2FC_KEGG, Metab_Log2FC_DEM),
    Metab_Pvalue = ifelse(!is.na(Metab_Pvalue_KEGG), Metab_Pvalue_KEGG, Metab_Pvalue_DEM),
    Metab_VIP_Final = ifelse(!is.na(Metab_VIP_KEGG), Metab_VIP_KEGG, Metab_VIP),
    Narrative_Axis = ifelse(is.na(Narrative_Axis), "Other candidate metabolite", Narrative_Axis)
  )

target_edges <- edge_evidence %>%
  filter(Narrative_Axis != "Other candidate metabolite") %>%
  arrange(Narrative_Axis, desc(Abs_Weight), Pvalue)

node_gene <- edge_evidence %>%
  distinct(Node_ID = Source) %>%
  left_join(deg_info, by = c("Node_ID" = "Source")) %>%
  mutate(Label = Node_ID, Log2FC = Gene_Log2FC, Node_Type = "Hub_Gene", Shape = "Ellipse") %>%
  select(Node_ID, Label, Log2FC, Node_Type, Shape, Gene_padj)

node_metab <- edge_evidence %>%
  distinct(Node_ID = Target) %>%
  left_join(dem_info, by = c("Node_ID" = "Target")) %>%
  left_join(kegg_info, by = c("Node_ID" = "Target")) %>%
  left_join(axis_map, by = c("Node_ID" = "Target")) %>%
  mutate(
    Label = Node_ID,
    Log2FC = ifelse(!is.na(Metab_Log2FC_KEGG), Metab_Log2FC_KEGG, Metab_Log2FC_DEM),
    Node_Type = "Metabolite",
    Shape = "Diamond"
  ) %>%
  select(Node_ID, Label, Log2FC, Node_Type, Shape, Metab_Pvalue_DEM, Metab_Pvalue_KEGG, Metab_FDR, Metab_VIP, Metab_VIP_KEGG, Narrative_Axis)

node_table <- bind_rows(
  node_gene %>% mutate(Metab_Pvalue_DEM = NA_real_, Metab_Pvalue_KEGG = NA_real_, Metab_FDR = NA_real_, Metab_VIP = NA_real_, Metab_VIP_KEGG = NA_real_, Narrative_Axis = NA_character_),
  node_metab %>% mutate(Gene_padj = NA_real_)
) %>%
  select(Node_ID, Label, Log2FC, Node_Type, Shape, Narrative_Axis, Gene_padj, Metab_Pvalue_DEM, Metab_Pvalue_KEGG, Metab_FDR, Metab_VIP, Metab_VIP_KEGG)

target_node_ids <- unique(c(target_edges$Source, target_edges$Target))
target_nodes <- node_table %>% filter(Node_ID %in% target_node_ids)

write.csv(all_edges, file.path(output_dir, "GH_All_Gene_Metabolite_Correlations.csv"), row.names = FALSE)
write.csv(edge_filtered, file.path(output_dir, "GH_Cytoscape_Edges_R0.90.csv"), row.names = FALSE)
write.csv(edge_evidence, file.path(output_dir, "GH_Evidence_Edges_R0.90.csv"), row.names = FALSE)
write.csv(target_edges, file.path(output_dir, "GH_Target_Axis_Evidence_Edges_R0.90.csv"), row.names = FALSE)
write.csv(node_table, file.path(output_dir, "GH_Cytoscape_Nodes_R0.90.csv"), row.names = FALSE)
write.csv(target_nodes, file.path(output_dir, "GH_Target_Axis_Nodes_R0.90.csv"), row.names = FALSE)
write.csv(rna_mat, file.path(output_dir, "GH_RNA_Input_Matrix_6samples.csv"), row.names = FALSE)
write.csv(dem_mat, file.path(output_dir, "GH_DEM_Input_Matrix_6samples.csv"), row.names = FALSE)

summary_table <- tibble::tibble(
  Metric = c("RNA genes", "Metabolites", "Samples used", "Correlation threshold", "All possible gene-metabolite pairs", "Edges with |R| >= threshold and P < 0.05", "Target-axis evidence edges"),
  Value = c(nrow(rna_values), nrow(dem_values_scaled), length(samples), cor_threshold, nrow(all_edges), nrow(edge_filtered), nrow(target_edges))
)
write.csv(summary_table, file.path(output_dir, "GH_Reanalysis_Summary.csv"), row.names = FALSE)
print(summary_table)

