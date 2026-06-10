# =======================================================
# 🌟 KEGG 通路富集 + 气泡图 (水稻专属 dosa)
# 适配自 新代码/kegg绘图.txt
# 用法: Rscript 06_KEGG_bubble.R <数据目录> [organism_code]
# 默认 organism: dosa (Oryza sativa japonica)
# =======================================================
rm(list=ls())
graphics.off()

# 读命令行参数
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  work_dir <- args[1]
} else {
  work_dir <- getwd()
}
# 优先命令行, 其次环境变量 METAB_KEGG_ORG (run_all.R 设置), 默认 dosa
org_code <- if (length(args) > 1) args[2] else Sys.getenv("METAB_KEGG_ORG", unset = "dosa")
setwd(work_dir)
message(">>> 📂 数据目录: ", work_dir)
message(">>> 🧬 KEGG 物种: ", org_code)

# =======================================================
# 1. 环境准备
# =======================================================
if (!requireNamespace("KEGGREST", quietly = TRUE))
  BiocManager::install("KEGGREST", update = FALSE)
if (!requireNamespace("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler", update = FALSE)

library(clusterProfiler)
library(tidyverse)
library(ggplot2)
library(scales)
library(stringr)
library(KEGGREST)

# =======================================================
# 2. 构建物种专属 KEGG 背景数据库
# =======================================================
message(sprintf(">>> 🌐 正在获取 %s 的 KEGG 通路数据...", org_code))

pathways_raw <- keggList("pathway", org_code)
valid_ids <- str_extract(names(pathways_raw), "\\d{5}")

map_cpd_link <- keggLink("compound", "pathway")
term2gene <- data.frame(
  Pathway = str_extract(names(map_cpd_link), "\\d{5}"),
  Compound = gsub("cpd:", "", as.character(map_cpd_link)),
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(Pathway) & Pathway %in% valid_ids) %>%
  select(Pathway, Compound)

term2name <- data.frame(
  Pathway = str_extract(names(pathways_raw), "\\d{5}"),
  Name = as.character(pathways_raw),
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(Pathway)) %>%
  mutate(Name = gsub(sprintf(" - %s.*",
    ifelse(org_code == "dosa", "Oryza sativa", ""), ""), "", Name))

message(sprintf("   ✅ 数据库: %d 个通路, %d 个化合物-通路关联",
                nrow(term2name), nrow(term2gene)))

# =======================================================
# 3. 气泡图主函数
# =======================================================
run_kegg_bubble <- function(treat_group) {

  target_file <- sprintf("04_KEGG_Target_%s_vs_WT.csv", treat_group)

  if (!file.exists(target_file)) {
    message(sprintf("⚠️ 未找到: %s，跳过", target_file))
    return(NULL)
  }

  message(sprintf("\n>>> 🚀 KEGG 富集: %s vs WT", treat_group))

  df <- read.csv(target_file, stringsAsFactors = FALSE)

  id_col <- intersect(c("KEGG_ID", "KEGG Compound ID"), colnames(df))[1]
  if (is.na(id_col)) stop("❌ 表格中没有 KEGG_ID 列！")

  target_cpds <- unique(unlist(strsplit(as.character(df[[id_col]]), "[;,]\\s*")))
  target_cpds <- gsub("^cpd:", "", trimws(target_cpds))
  target_cpds <- target_cpds[target_cpds != "" & !is.na(target_cpds)]

  cat(sprintf("   输入 KEGG ID: %d 个\n", length(target_cpds)))

  # 富集分析 (pvalueCutoff = 1 获取全量结果，不做硬阈值)
  enrich_res <- enricher(
    gene       = target_cpds,
    pvalueCutoff = 1,
    TERM2GENE  = term2gene,
    TERM2NAME  = term2name
  )

  if (is.null(enrich_res) || nrow(enrich_res@result) == 0) {
    message(sprintf("   ⚠️ %s 无富集结果", treat_group))
    return(NULL)
  }

  # 全量数据 + Rich Factor
  full_res <- enrich_res@result %>%
    mutate(
      M = as.numeric(sub("/.*", "", BgRatio)),
      RichFactor = Count / M
    ) %>%
    arrange(p.adjust)

  # 保存全量通路表
  csv_out <- sprintf("05_KEGG_All_Pathways_%s_vs_WT.csv", treat_group)
  write.csv(full_res, csv_out, row.names = FALSE)
  message(sprintf("   📄 全量通路表: %s (%d 条)", csv_out, nrow(full_res)))

  # 取 Top 10 画气泡图
  plot_df <- full_res %>% head(10)

  if (nrow(plot_df) == 0) {
    message("   ⚠️ 无数据可绘图")
    return(full_res)
  }

  # 文字折行 + Y轴排序
  plot_df$Description <- str_wrap(plot_df$Description, width = 25)
  plot_df$Description <- factor(plot_df$Description, levels = rev(plot_df$Description))

  # 绘制气泡图
  p <- ggplot(plot_df, aes(x = RichFactor, y = Description)) +
    geom_segment(aes(x = 0, xend = RichFactor, y = Description, yend = Description),
                 color = "grey85", linewidth = 0.5) +
    geom_point(aes(size = Count, fill = p.adjust),
               shape = 21, color = "black", stroke = 0.6, alpha = 0.9) +
    scale_fill_gradient(low = "#E64B35", high = "#4DBBD5", name = "p.adjust") +
    scale_size_continuous(range = c(5, 12), name = "Count") +
    guides(
      fill = guide_colorbar(order = 1, frame.colour = "black", ticks.colour = "black"),
      size = guide_legend(order = 2)
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.15, 0.30))) +
    labs(title = paste0("KEGG: ", treat_group, " vs WT"),
         x = "Rich Factor", y = NULL) +
    theme_bw() +
    theme(
      text = element_text(color = "black"),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.5),
      panel.grid.minor = element_line(color = "grey96", linewidth = 0.3),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
      axis.text.y = element_text(size = 12, color = "black", face = "italic"),
      axis.text.x = element_text(size = 13, color = "black"),
      axis.title.x = element_text(size = 15, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    )

  # 保存 PDF
  pdf_out <- sprintf("05_KEGG_Bubble_%s_vs_WT.pdf", treat_group)
  ggsave(pdf_out, p, width = 5.25, height = 7.2, device = cairo_pdf)
  message(sprintf("   📊 气泡图: %s", pdf_out))

  return(full_res)
}

# =======================================================
# 4. 批量处理
# =======================================================
target_files <- list.files(pattern = "^04_KEGG_Target_.*_vs_WT\\.csv$")
if (length(target_files) == 0) {
  stop("❌ 未找到 04_KEGG_Target_*.csv，请先运行 KEGG 靶点提取。")
}

cat(sprintf("\n>>> 发现 %d 个靶点文件\n", length(target_files)))
for (f in target_files) cat(sprintf("    %s\n", f))

for (f in target_files) {
  m <- regmatches(f, regexec("04_KEGG_Target_(.+)_vs_WT\\.csv", f))[[1]]
  if (length(m) >= 2) {
    run_kegg_bubble(m[2])
  }
}

message("\n>>> 🏆 KEGG 通路富集分析完成！")
