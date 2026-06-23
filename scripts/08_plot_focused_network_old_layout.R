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
output_dir <- if (length(args) >= 1) args[1] else getwd()
mode <- if (length(args) >= 2) args[2] else "SmartLeader"

showtext_auto()
while (!is.null(dev.list())) dev.off()

edge_file <- file.path(output_dir, "GH_Target_Axis_Evidence_Edges_R0.90.csv")
node_file <- file.path(output_dir, "GH_Target_Axis_Nodes_R0.90.csv")
if (!file.exists(edge_file)) stop("Missing edge file: ", edge_file)
if (!file.exists(node_file)) stop("Missing node file: ", node_file)

edges_raw <- read.csv(edge_file, stringsAsFactors = FALSE, check.names = FALSE)
nodes_raw <- read.csv(node_file, stringsAsFactors = FALSE, check.names = FALSE)

target_metabs <- sort(unique(edges_raw$Target))
edges_focused <- edges_raw %>%
  mutate(
    Source = str_trim(as.character(Source)),
    Target = str_trim(as.character(Target)),
    Weight = as.numeric(Weight),
    Direction = ifelse(Weight > 0, "Positive", "Negative")
  ) %>%
  filter(abs(Weight) >= 0.90)

active_ids <- unique(c(edges_focused$Source, edges_focused$Target))
nodes_focused <- data.frame(Node_ID = active_ids, stringsAsFactors = FALSE) %>%
  left_join(mutate(nodes_raw, Node_ID = str_trim(as.character(Node_ID))), by = "Node_ID") %>%
  mutate(
    Node_Type = ifelse(str_detect(Node_ID, "^Os"), "Hub_Gene", "Target_Metabolite"),
    Label = ifelse(!is.na(Label) & Label != "", Label, Node_ID),
    Log2FC_Display = ifelse(Node_Type == "Hub_Gene", as.numeric(Log2FC), 0),
    Log2FC_Display = ifelse(is.na(Log2FC_Display), 0, Log2FC_Display)
  )

net_graph <- tbl_graph(nodes = nodes_focused, edges = edges_focused, directed = FALSE)

segment_alpha <- ifelse(tolower(mode) == "noleader", 0, 0.72)
min_segment <- ifelse(tolower(mode) == "smartleader", 0.75, 0)
suffix <- ifelse(tolower(mode) == "noleader", "NoLeader", "SmartLeader")

set.seed(123)
p_net <- ggraph(net_graph, layout = "fr") +
  geom_edge_link(aes(edge_width = abs(Weight), edge_color = Direction), alpha = 0.52) +
  scale_edge_color_manual(values = c("Positive" = "#C0392B", "Negative" = "#2980B9"), name = "Correlation") +
  scale_edge_width_continuous(range = c(0.55, 2.15), guide = "none") +
  geom_node_point(aes(filter = Node_Type == "Hub_Gene", fill = Log2FC_Display), shape = 21, color = "black", stroke = 0.45, size = 6.8) +
  scale_fill_gradient2(low = "#0881A3", mid = "#F7F7F7", high = "#911F27", midpoint = 0, name = "Gene log2FC") +
  geom_node_point(aes(filter = Node_Type == "Target_Metabolite"), shape = 23, fill = "#276231", color = "black", stroke = 0.55, size = 10) +
  geom_node_text(
    aes(label = Label),
    repel = TRUE,
    size = 4.2,
    fontface = "bold.italic",
    max.overlaps = Inf,
    box.padding = 0.75,
    point.padding = 0.35,
    segment.color = "black",
    segment.size = 0.28,
    segment.alpha = segment_alpha,
    min.segment.length = min_segment
  ) +
  theme_graph(base_family = "Arial") +
  theme(
    legend.position = "right",
    plot.margin = margin(18, 18, 18, 18),
    text = element_text(family = "Arial"),
    plot.title = element_text(size = 15, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 10, color = "#4D4D4D", hjust = 0),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8)
  ) +
  labs(
    title = "Focused gene-metabolite network",
    subtitle = paste0("Pearson |R| >= 0.90, P < 0.05 | ", paste(target_metabs, collapse = " / "))
  )

base_name <- file.path(output_dir, paste0("GH_Focused_Network_OldLayout_", suffix, "_R0.90"))

grDevices::cairo_pdf(paste0(base_name, ".pdf"), width = 10, height = 9, family = "Arial")
print(p_net)
dev.off()
svglite::svglite(paste0(base_name, ".svg"), width = 10, height = 9)
print(p_net)
dev.off()
ragg::agg_tiff(paste0(base_name, ".tiff"), width = 10, height = 9, units = "in", res = 600, compression = "lzw")
print(p_net)
dev.off()
ragg::agg_png(paste0(base_name, "_preview.png"), width = 10, height = 9, units = "in", res = 180)
print(p_net)
dev.off()

write.csv(edges_focused, file.path(output_dir, paste0("GH_Focused_Network_OldLayout_", suffix, "_Edges_R0.90.csv")), row.names = FALSE)
write.csv(nodes_focused, file.path(output_dir, paste0("GH_Focused_Network_OldLayout_", suffix, "_Nodes_R0.90.csv")), row.names = FALSE)
cat("Saved: ", base_name, "\n", sep = "")
