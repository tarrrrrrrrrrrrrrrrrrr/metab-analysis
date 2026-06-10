# =======================================================
# 🌟 火山图：对称纯净版（无标签）
# KO_2 vs WT / KO_6 vs WT
# 完全适配自 新代码/火山图.txt
# =======================================================
rm(list=ls())
graphics.off()

library(tidyverse)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
work_dir <- if (length(args) > 0) args[1] else getwd()
setwd(work_dir)
message(">>> 📂 数据目录: ", work_dir)

#' 火山图绘制函数
#' @param treat_group 处理组名称 ("KO_2", "KO_6")
#' @param target_color 上调点高亮色
#' @param fc_threshold 倍数阈值 (默认 1.5)
#' @param pval_cut 显著性阈值 (默认 0.05)
#' @param y_metric Y轴: "Pvalue" 或 "FDR"
run_volcano <- function(treat_group, target_color,
                        fc_threshold = 1.5, pval_cut = 0.05,
                        y_metric = "Pvalue") {

  res_file <- sprintf("OPLSDA_Diff_Results_%s_vs_WT.csv", treat_group)

  if (!file.exists(res_file)) return(message("❌ 未找到: ", res_file))
  message(sprintf("\n>>> 🚀 火山图: %s vs WT ...", treat_group))

  # 读取
  plot_df <- read.csv(res_file, stringsAsFactors = FALSE)

  # 数据清洗 + Z-order 排序
  plot_df <- plot_df %>%
    mutate(
      Status = factor(Status, levels = c("NS", "Down", "Up")),
      Y_val = -log10(.data[[y_metric]])
    ) %>%
    arrange(Status)

  # 参考线
  lfc_cut <- log10(fc_threshold)

  # 对称轴延展
  max_abs_logFC <- max(abs(plot_df$Log10FC), na.rm = TRUE)
  if (max_abs_logFC < lfc_cut) max_abs_logFC <- lfc_cut * 1.5
  x_limit <- max_abs_logFC * 1.1
  y_limit <- max(plot_df$Y_val[is.finite(plot_df$Y_val)], na.rm = TRUE) * 1.05

  # 绘图
  p <- ggplot(plot_df, aes(x = Log10FC, y = Y_val, color = Status)) +
    geom_point(alpha = 0.8, size = 1.8) +

    # NPG 配色：NS 灰，Down 深蓝，Up 处理组色
    scale_color_manual(
      values = c("NS" = "#e0e0e0", "Down" = "#3C5488", "Up" = target_color),
      drop = FALSE
    ) +

    # 阈值线
    geom_vline(xintercept = c(-lfc_cut, lfc_cut),
               linetype = "dashed", color = "grey40", linewidth = 0.5) +
    geom_hline(yintercept = -log10(pval_cut),
               linetype = "dashed", color = "grey40", linewidth = 0.5) +

    labs(x = expression(Log[10]*" (Fold Change)"),
         y = bquote(-Log[10]*" ("*.(y_metric)*")"),
         title = paste(treat_group, "vs WT")) +

    theme_bw() +
    theme(
      panel.grid = element_blank(),
      aspect.ratio = 1,
      text = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      axis.title = element_text(size = 13, face = "bold"),
      axis.text = element_text(size = 11, color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
      legend.position = "right",
      legend.title = element_blank(),
      legend.background = element_blank(),
      legend.key = element_blank()
    ) +
    scale_x_continuous(limits = c(-x_limit, x_limit)) +
    scale_y_continuous(limits = c(0, y_limit), expand = expansion(mult = c(0, 0.05)))

  # 保存
  pdf_out <- sprintf("02_Volcano_%s_vs_WT.pdf", treat_group)
  ggsave(pdf_out, p, width = 7, height = 5.5, device = cairo_pdf)
  message(sprintf("   ✅ 已保存: %s", pdf_out))
}

# ==========================================
# 批量执行
# ==========================================
run_volcano("KO_2", target_color = "#E64B35")   # 砖红
run_volcano("KO_6", target_color = "#FFCC99")   # 暖砂金

message("\n>>> 🏆 火山图完成！")
