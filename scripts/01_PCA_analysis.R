# =======================================================
# 🌟 PCA 分析：WT vs KO_2 / WT vs KO_6（含 QC 质控）
# 数据源：metab_abund_named.txt（tab分隔，237代谢物 × 21样本）
# 完全适配自 新代码/pca.txt
# =======================================================
rm(list=ls())
graphics.off()

library(tidyverse)
library(ggplot2)

# 🚨 工作路径：命令行参数或当前目录
args <- commandArgs(trailingOnly = TRUE)
work_dir <- if (length(args) > 0) args[1] else getwd()
setwd(work_dir)
message(">>> 📂 数据目录: ", work_dir)

# ==========================================
# 封装 PCA 绘图函数（适配 WT/KO 命名）
# ==========================================
plot_split_pca <- function(treat_group, target_color, target_shape) {

  # 读取数据（tab分隔，行名为代谢物名）
  raw_data <- read.table("metab_abund_named.txt", header = TRUE,
                         sep = "\t", check.names = FALSE, row.names = 1)

  # 智能匹配样本列：CK → WT, 处理组 → KO_2/KO_6, 含 QC
  pattern <- sprintf("^WT_|^%s_|^QC", treat_group)
  sample_cols <- grep(pattern, colnames(raw_data), value = TRUE)

  expr_matrix <- raw_data[, sample_cols, drop = FALSE]
  expr_matrix <- apply(expr_matrix, 2, as.numeric)
  rownames(expr_matrix) <- rownames(raw_data)

  # 自动识别分组标签
  groups <- case_when(
    grepl("^WT", sample_cols)     ~ "WT",
    grepl(paste0("^", treat_group), sample_cols) ~ treat_group,
    grepl("^QC", sample_cols)     ~ "QC",
    TRUE                          ~ "Other"
  )

  # 防止 log2(0)
  expr_matrix[is.na(expr_matrix)] <- 0
  min_val <- min(expr_matrix[expr_matrix > 0], na.rm = TRUE)
  expr_matrix[expr_matrix == 0] <- min_val / 2

  # 对数转换
  expr_matrix_log <- log2(expr_matrix)

  # 去除零方差行（避免 prcomp 报错）
  row_vars <- apply(expr_matrix_log, 1, var, na.rm = TRUE)
  expr_matrix_log <- expr_matrix_log[row_vars > 1e-10, , drop = FALSE]

  # PCA
  pca_res <- prcomp(t(expr_matrix_log), scale. = TRUE)
  pca_scores <- as.data.frame(pca_res$x[, 1:2])
  colnames(pca_scores) <- c("PC1", "PC2")
  pca_scores$Group <- factor(groups, levels = c("WT", treat_group, "QC"))

  # 贡献率
  vars <- pca_res$sdev^2
  percentVar <- round(100 * vars / sum(vars))

  # NPG 配色：WT 深蓝，QC 翠绿，KO 专属色
  my_colors <- c("WT" = "#3C5488", "QC" = "#00A087")
  my_colors[treat_group] <- target_color

  my_shapes <- c("WT" = 16, "QC" = 18)
  my_shapes[treat_group] <- target_shape

  # 绘图
  p <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Group, shape = Group)) +
    stat_ellipse(level = 0.95, linetype = "solid", linewidth = 0.8,
                 alpha = 0.6, show.legend = FALSE) +
    geom_point(size = 1.67, alpha = 0.9) +
    xlab(paste0("PC1 (", percentVar[1], "%)")) +
    ylab(paste0("PC2 (", percentVar[2], "%)")) +
    labs(title = paste(treat_group, "vs WT")) +
    scale_color_manual(values = my_colors) +
    scale_shape_manual(values = my_shapes) +
    scale_x_continuous(expand = expansion(mult = 0.25)) +
    scale_y_continuous(expand = expansion(mult = 0.25)) +
    theme_bw() +
    theme(
      text = element_text(color = "black"),
      aspect.ratio = 1,
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
      axis.text = element_text(size = 13, color = "black"),
      axis.title = element_text(size = 15, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 15, face = "bold",
                                margin = margin(b = 15)),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 13),
      legend.key.size = unit(1.2, "lines"),
      legend.background = element_blank(),
      plot.margin = margin(t = 15, r = 15, b = 15, l = 15)
    )

  return(p)
}

# ==========================================
# 批量执行：KO_2 和 KO_6
# ==========================================
message("\n>>> 🚀 正在 PCA 降维分析 [KO_2 vs WT] (含 QC) ...")
p_ko2 <- plot_split_pca("KO_2", target_color = "#E64B35", target_shape = 17)
if (!is.null(p_ko2)) {
  ggsave("01_PCA_KO_2_vs_WT.pdf", p_ko2, width = 6.0, height = 5.5, device = cairo_pdf)
  message("   ✅ KO_2 vs WT PCA 图已保存")
}

message("\n>>> 🚀 正在 PCA 降维分析 [KO_6 vs WT] (含 QC) ...")
p_ko6 <- plot_split_pca("KO_6", target_color = "#FFCC99", target_shape = 15)
if (!is.null(p_ko6)) {
  ggsave("01_PCA_KO_6_vs_WT.pdf", p_ko6, width = 6.0, height = 5.5, device = cairo_pdf)
  message("   ✅ KO_6 vs WT PCA 图已保存")
}

# ==========================================
# 补充：三组 PCA（WT + KO_2 + KO_6 + QC）
# ==========================================
raw_data <- read.table("metab_abund_named.txt", header = TRUE,
                       sep = "\t", check.names = FALSE, row.names = 1)

sample_cols <- grep("^WT_|^KO_2_|^KO_6_|^QC", colnames(raw_data), value = TRUE)
expr_matrix <- raw_data[, sample_cols]
expr_matrix <- apply(expr_matrix, 2, as.numeric)
rownames(expr_matrix) <- rownames(raw_data)

groups_all <- case_when(
  grepl("^WT", sample_cols)   ~ "WT",
  grepl("^KO_2", sample_cols) ~ "KO_2",
  grepl("^KO_6", sample_cols) ~ "KO_6",
  grepl("^QC", sample_cols)   ~ "QC"
)

expr_matrix[is.na(expr_matrix)] <- 0
min_val <- min(expr_matrix[expr_matrix > 0], na.rm = TRUE)
expr_matrix[expr_matrix == 0] <- min_val / 2
expr_matrix_log <- log2(expr_matrix)

# 去除零方差行
row_vars_all <- apply(expr_matrix_log, 1, var, na.rm = TRUE)
expr_matrix_log <- expr_matrix_log[row_vars_all > 1e-10, , drop = FALSE]

pca_res <- prcomp(t(expr_matrix_log), scale. = TRUE)
pca_scores <- as.data.frame(pca_res$x[, 1:2])
colnames(pca_scores) <- c("PC1", "PC2")
pca_scores$Group <- factor(groups_all, levels = c("WT", "KO_2", "KO_6", "QC"))

vars <- pca_res$sdev^2
percentVar <- round(100 * vars / sum(vars))

my_colors_all <- c("WT" = "#3C5488", "KO_2" = "#E64B35",
                    "KO_6" = "#FFCC99", "QC" = "#00A087")
my_shapes_all <- c("WT" = 16, "KO_2" = 17, "KO_6" = 15, "QC" = 18)

p_all <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Group, shape = Group)) +
  stat_ellipse(level = 0.95, linetype = "solid", linewidth = 0.8,
               alpha = 0.6, show.legend = FALSE) +
  geom_point(size = 2, alpha = 0.9) +
  xlab(paste0("PC1 (", percentVar[1], "%)")) +
  ylab(paste0("PC2 (", percentVar[2], "%)")) +
  labs(title = "All Groups PCA") +
  scale_color_manual(values = my_colors_all) +
  scale_shape_manual(values = my_shapes_all) +
  scale_x_continuous(expand = expansion(mult = 0.25)) +
  scale_y_continuous(expand = expansion(mult = 0.25)) +
  theme_bw() +
  theme(
    text = element_text(color = "black"),
    aspect.ratio = 1,
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    axis.text = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 15, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1.2, "lines")
  )

ggsave("01_PCA_All_Groups.pdf", p_all, width = 7, height = 5.5, device = cairo_pdf)
message("   ✅ 三组 PCA 图已保存")

message("\n>>> 🏆 PCA 分析完成！3张PDF已输出至 qijng 文件夹。")
