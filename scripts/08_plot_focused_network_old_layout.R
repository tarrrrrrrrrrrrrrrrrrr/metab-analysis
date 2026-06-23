suppressPackageStartupMessages({
  library(dplyr)
  library(ggraph)
  library(tidygraph)
  library(stringr)
  library(showtext)
  library(svglite)
  library(ragg)
})

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[[1]] else getwd()
mode <- if (length(args) >= 2) args[[2]] else "Compact"
prefix <- if (length(args) >= 3) args[[3]] else NA_character_
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

showtext_auto()
while (!is.null(dev.list())) dev.off()

find_one <- function(pattern) {
  files <- list.files(output_dir, pattern = pattern, full.names = TRUE)
  if (!is.na(prefix)) files <- files[grepl(paste0("^", prefix, "_"), basename(files))]
  if (length(files) == 0) stop("No file matching: ", pattern)
  files[[1]]
}

edge_file <- find_one("_Target_Axis_Evidence_Edges_R0\\.90\\.csv$")
node_file <- find_one("_Target_Axis_Nodes_R0\\.90\\.csv$")
if (is.na(prefix)) prefix <- sub("_Target_Axis_Evidence_Edges_R0\\.90\\.csv$", "", basename(edge_file))

edges_raw <- read.csv(edge_file, stringsAsFactors = FALSE, check.names = FALSE)
nodes_raw <- read.csv(node_file, stringsAsFactors = FALSE, check.names = FALSE)

short_label <- function(x) {
  x <- as.character(x)
  x <- gsub("Cyanidin 3-Galactoside", "Cyanidin-3-Gal", x, fixed = TRUE)
  x <- gsub("Isovitexin 2''-O-Glucoside", "Isovitexin-2G", x, fixed = TRUE)
  x <- gsub("Phenylacetylglutamine", "Phenylacetyl-Gln", x, fixed = TRUE)
  x <- gsub("Eugenyl Acetate", "Eugenyl-Ac", x, fixed = TRUE)
  x <- gsub("Vitamin K3", "Vit K3", x, fixed = TRUE)
  x <- ifelse(nchar(x) > 20, paste0(substr(x, 1, 17), "..."), x)
  x
}

edges_focused <- edges_raw %>%
  mutate(
    Source = str_trim(as.character(Source)),
    Target = str_trim(as.character(Target)),
    Weight = as.numeric(Weight),
    Pvalue = if ("Pvalue" %in% names(.)) as.numeric(Pvalue) else NA_real_,
    Abs_Weight = abs(Weight),
    Direction = ifelse(Weight > 0, "Positive", "Negative")
  ) %>%
  filter(Abs_Weight >= 0.90) %>%
  { if (any(!is.na(.$Pvalue))) filter(., Pvalue < 0.05) else . } %>%
  group_by(Target) %>%
  slice_max(order_by = Abs_Weight, n = 20, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(Target, desc(Abs_Weight), Pvalue)

active_ids <- unique(c(edges_focused$Source, edges_focused$Target))
gene_degree <- edges_focused %>% distinct(Source, Target) %>% count(Source, name = "Gene_Target_Count")
metab_degree <- edges_focused %>% count(Target, name = "Metab_Gene_Count")

gene_logfc <- edges_focused %>%
  transmute(Node_ID = Source, Edge_Gene_Log2FC = if ("Gene_Log2FC" %in% names(edges_focused)) as.numeric(Gene_Log2FC) else NA_real_) %>%
  distinct(Node_ID, .keep_all = TRUE)

nodes_focused <- data.frame(Node_ID = active_ids, stringsAsFactors = FALSE) %>%
  left_join(mutate(nodes_raw, Node_ID = str_trim(as.character(Node_ID))), by = "Node_ID") %>%
  left_join(gene_logfc, by = "Node_ID") %>%
  left_join(gene_degree, by = c("Node_ID" = "Source")) %>%
  left_join(metab_degree, by = c("Node_ID" = "Target")) %>%
  mutate(
    Node_Type = ifelse(Node_ID %in% edges_focused$Source, "Hub_Gene", "Target_Metabolite"),
    Label = ifelse(!is.na(Label) & Label != "" & Label != "<NA>", Label, Node_ID),
    Plot_Label = ifelse(Node_Type == "Target_Metabolite", short_label(Node_ID), Label),
    Log2FC_Display = case_when(
      Node_Type == "Hub_Gene" & "Log2FC" %in% names(.) ~ as.numeric(Log2FC),
      Node_Type == "Hub_Gene" ~ as.numeric(Edge_Gene_Log2FC),
      TRUE ~ 0
    ),
    Log2FC_Display = ifelse(is.na(Log2FC_Display), 0, Log2FC_Display),
    Gene_Target_Count = ifelse(is.na(Gene_Target_Count), 0, Gene_Target_Count),
    Metab_Gene_Count = ifelse(is.na(Metab_Gene_Count), 0, Metab_Gene_Count)
  )

compact <- tolower(mode) %in% c("compact", "reference17", "compactreference17")
no_leader <- tolower(mode) == "noleader"

nodes_focused <- nodes_focused %>%
  mutate(
    Node_Size = if (compact) {
      case_when(
        Node_Type == "Hub_Gene" ~ 3.15 + pmin(Gene_Target_Count, 5) * 0.23,
        Node_Type == "Target_Metabolite" ~ 4.65 + sqrt(Metab_Gene_Count) * 0.28,
        TRUE ~ 3.6
      )
    } else {
      case_when(
        Node_Type == "Hub_Gene" ~ 5.8 + pmin(Gene_Target_Count, 5) * 0.30,
        Node_Type == "Target_Metabolite" ~ 8.5 + sqrt(Metab_Gene_Count) * 0.35,
        TRUE ~ 5.0
      )
    }
  )

net_graph <- tbl_graph(nodes = nodes_focused, edges = edges_focused, directed = FALSE)

edge_range <- if (compact) c(0.12, 0.55) else c(0.55, 2.15)
edge_alpha <- if (compact) 0.42 else 0.52
label_size <- if (compact) 2.35 else 4.2
seg_size <- if (compact) 0.11 else 0.28
seg_alpha <- if (no_leader) 0 else if (compact) 0.48 else 0.72
min_seg <- if (compact) 0.28 else if (tolower(mode) == "smartleader") 0.75 else 0
metab_shape <- if (compact) 22 else 23
suffix <- if (compact) "Compact" else if (no_leader) "NoLeader" else "SmartLeader"

set.seed(123)
p_net <- ggraph(net_graph, layout = "fr", niter = if (compact) 5000 else 1200, start.temp = if (compact) 1.8 else 2.5) +
  geom_edge_link(aes(edge_width = Abs_Weight, edge_color = Direction), alpha = edge_alpha, lineend = "round") +
  scale_edge_color_manual(values = c("Positive" = "#C0392B", "Negative" = "#2980B9"), name = "Correlation Type") +
  scale_edge_width_continuous(range = edge_range, name = "Pearson |r|") +
  geom_node_point(aes(filter = Node_Type == "Hub_Gene", fill = Log2FC_Display, size = Node_Size), shape = 21, color = "black", stroke = if (compact) 0.24 else 0.45) +
  scale_fill_gradient2(low = "#0881A3", mid = "#F7F7F7", high = "#911F27", midpoint = 0, name = "Gene Log2FC") +
  geom_node_point(aes(filter = Node_Type == "Target_Metabolite", size = Node_Size), shape = metab_shape, fill = "#276231", color = "black", stroke = if (compact) 0.32 else 0.55) +
  scale_size_identity() +
  geom_node_text(
    aes(label = Plot_Label),
    repel = TRUE,
    size = label_size,
    fontface = "bold.italic",
    max.overlaps = Inf,
    box.padding = if (compact) 0.34 else 0.75,
    point.padding = if (compact) 0.12 else 0.35,
    min.segment.length = min_seg,
    segment.color = "black",
    segment.size = seg_size,
    segment.alpha = seg_alpha
  ) +
  theme_graph(base_family = "Arial") +
  theme(
    legend.position = "right",
    plot.margin = margin(14, 14, 12, 14),
    text = element_text(family = "Arial"),
    legend.title = element_text(size = if (compact) 6.2 else 9, face = "bold"),
    legend.text = element_text(size = if (compact) 5.5 else 8),
    plot.title = element_text(size = if (compact) 9.5 else 15, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = if (compact) 6.5 else 10, color = "#4D4D4D", hjust = 0)
  ) +
  labs(
    title = if (compact) "Metabolites (Dark Green Squares) vs. Hub Genes (Red/Blue Circles)" else "Focused gene-metabolite network",
    subtitle = paste0(prefix, " focused network | Top20 per metabolite | Pearson |r| >= 0.90, P < 0.05")
  )

base_name <- file.path(output_dir, paste0(prefix, "_Focused_Network_OldLayout_", suffix, "_R0.90"))
w <- if (compact) 8.8 else 10
h <- if (compact) 8.2 else 9

grDevices::cairo_pdf(paste0(base_name, ".pdf"), width = w, height = h, family = "Arial")
print(p_net)
dev.off()
svglite::svglite(paste0(base_name, ".svg"), width = w, height = h)
print(p_net)
dev.off()
ragg::agg_tiff(paste0(base_name, ".tiff"), width = w, height = h, units = "in", res = 600, compression = "lzw")
print(p_net)
dev.off()
ragg::agg_png(paste0(base_name, "_preview.png"), width = w, height = h, units = "in", res = 180)
print(p_net)
dev.off()

write.csv(edges_focused, file.path(output_dir, paste0(prefix, "_Focused_Network_OldLayout_", suffix, "_Edges_R0.90.csv")), row.names = FALSE)
write.csv(nodes_focused, file.path(output_dir, paste0(prefix, "_Focused_Network_OldLayout_", suffix, "_Nodes_R0.90.csv")), row.names = FALSE)
cat("Saved: ", base_name, "\n", sep = "")
