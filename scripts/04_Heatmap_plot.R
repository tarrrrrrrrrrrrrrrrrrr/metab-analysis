# =======================================================
# 🌟 差异代谢物聚类热图
# KO_2 vs WT / KO_6 vs WT 显著差异代谢物
# 完全适配自 新代码/热图.txt 的配色与风格
# =======================================================
rm(list=ls())
graphics.off()

required_pkgs <- c("tidyverse", "ComplexHeatmap", "circlize")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

args <- commandArgs(trailingOnly = TRUE)
work_dir <- if (length(args) > 0) args[1] else getwd()
setwd(work_dir)
message(">>> 📂 数据目录: ", work_dir)

# =======================================================
# 1. 读取差异代谢物结果和丰度矩阵
# =======================================================
res_ko2 <- read.csv("OPLSDA_Diff_Results_KO_2_vs_WT.csv", stringsAsFactors = FALSE)
res_ko6 <- read.csv("OPLSDA_Diff_Results_KO_6_vs_WT.csv", stringsAsFactors = FALSE)

# 合并所有显著差异代谢物（去重）
sig_ko2 <- res_ko2$Name[res_ko2$Status != "NS"]
sig_ko6 <- res_ko6$Name[res_ko6$Status != "NS"]
all_sig_names <- unique(c(sig_ko2, sig_ko6))
message(sprintf("KO_2 差异: %d, KO_6 差异: %d, 合并(去重): %d",
                length(sig_ko2), length(sig_ko6), length(all_sig_names)))

if (length(all_sig_names) == 0) {
  stop("❌ 无显著差异代谢物，无法绘制热图。")
}

# 读取丰度矩阵
abund <- read.table("metab_abund_named.txt", header = TRUE,
                    sep = "\t", check.names = FALSE, row.names = 1)

# 提取差异代谢物 + 去除QC只保留生物学样本
bio_cols <- grep("^(WT|KO_2|KO_6)_", colnames(abund), value = TRUE)
abund_bio <- abund[, bio_cols]

# 取交集
target_names <- all_sig_names[all_sig_names %in% rownames(abund_bio)]
message(sprintf("丰度矩阵中匹配到: %d 个差异代谢物", length(target_names)))

if (length(target_names) < 2) {
  stop("❌ 匹配到的差异代谢物数量不足（<2），无法聚类。")
}

vst_target <- abund_bio[target_names, , drop = FALSE]

# =======================================================
# 2. Z-score 标准化 + 热图
# =======================================================
z_matrix <- t(apply(vst_target, 1, function(x) (x - mean(x)) / sd(x)))
z_matrix <- na.omit(z_matrix)

# NPG 配色：深蓝-白-砖红
col_fun <- colorRamp2(c(-2, 0, 2), c("#3C5488", "white", "#E64B35"))

# 分组注释
group_info <- factor(
  case_when(
    grepl("^WT",   colnames(z_matrix)) ~ "WT",
    grepl("^KO_2", colnames(z_matrix)) ~ "KO_2",
    grepl("^KO_6", colnames(z_matrix)) ~ "KO_6"
  ),
  levels = c("WT", "KO_2", "KO_6")
)

top_anno <- HeatmapAnnotation(
  Group = group_info,
  col = list(Group = c("WT" = "#3C5488", "KO_2" = "#E64B35", "KO_6" = "#FFCC99")),
  annotation_legend_param = list(title_gp = gpar(fontsize = 10))
)

# 确定是否显示行名
show_names <- nrow(z_matrix) <= 60

cairo_pdf("03_Heatmap_Diff_Metabolites.pdf", width = 8.0, height = 6.0)

ht <- Heatmap(z_matrix,
  name = "Z-score",
  col = col_fun,
  top_annotation = top_anno,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_dend_width = unit(2.5, "cm"),
  show_row_names = show_names,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 7),
  column_title = "Differential Metabolites Heatmap",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  column_names_gp = gpar(fontsize = 10),
  row_title = paste0("Diff Metabolites (n=", nrow(z_matrix), ")"),
  row_title_gp = gpar(fontsize = 12),
  km = ifelse(nrow(z_matrix) >= 10, min(5, floor(nrow(z_matrix)/2)), 1),
  border = TRUE,
  heatmap_legend_param = list(title_gp = gpar(fontsize = 10))
)

draw(ht)
dev.off()

message(sprintf("   ✅ 热图已保存 (共 %d 个差异代谢物)", nrow(z_matrix)))

# =======================================================
# 3. 分别画 KO_2 和 KO_6 各自差异代谢物热图
# =======================================================
for (treat in c("KO_2", "KO_6")) {
  sig_names <- if (treat == "KO_2") sig_ko2 else sig_ko6
  sig_names <- sig_names[sig_names %in% rownames(abund_bio)]

  if (length(sig_names) < 3) next

  # 仅取该对比组的样本
  gcols <- grep(paste0("^(WT_|", treat, "_)"), colnames(abund_bio), value = TRUE)
  sub_mat <- abund_bio[sig_names, gcols, drop = FALSE]
  z_sub <- t(apply(sub_mat, 1, function(x) (x - mean(x)) / sd(x)))
  z_sub <- na.omit(z_sub)

  sub_group <- factor(
    ifelse(grepl("^WT", colnames(z_sub)), "WT", treat),
    levels = c("WT", treat)
  )

  sub_anno <- HeatmapAnnotation(
    Group = sub_group,
    col = list(Group = setNames(c("#3C5488", if(treat=="KO_2") "#E64B35" else "#FFCC99"),
                                 c("WT", treat)))
  )

  show_n <- nrow(z_sub) <= 80
  pdf_out <- sprintf("03_Heatmap_%s_vs_WT.pdf", treat)

  cairo_pdf(pdf_out, width = 6.5, height = 6.0)
  ht_sub <- Heatmap(z_sub,
    name = "Z-score",
    col = col_fun,
    top_annotation = sub_anno,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    row_dend_width = unit(2, "cm"),
    show_row_names = show_n,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 7),
    column_title = paste0(treat, " vs WT (n=", nrow(z_sub), ")"),
    column_title_gp = gpar(fontsize = 13, fontface = "bold"),
    column_names_gp = gpar(fontsize = 10),
    km = ifelse(nrow(z_sub) >= 10, min(4, floor(nrow(z_sub)/2)), 1),
    border = TRUE
  )
  draw(ht_sub)
  dev.off()
  message(sprintf("   ✅ %s 单独热图已保存", treat))
}

message("\n>>> 🏆 热图任务完成！")
