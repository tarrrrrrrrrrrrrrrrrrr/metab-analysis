# =======================================================
# 🌟 OPLS-DA 建模 + VIP + 差异筛选
# KO_2 vs WT / KO_6 vs WT
# 数据为 Log10 尺度丰度矩阵
# 完全适配自 新代码/OPLS-DA 建模与提取 VIP.txt
# =======================================================
rm(list=ls())
graphics.off()

if (!requireNamespace("ropls", quietly = TRUE))
  BiocManager::install("ropls", update = FALSE)
library(ropls)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
work_dir <- if (length(args) > 0) args[1] else getwd()
setwd(work_dir)
message(">>> 📂 数据目录: ", getwd())

# =======================================================
# 1. 读取数据
# =======================================================
# 丰度矩阵（已命名）
abund <- read.table("metab_abund_named.txt", header = TRUE,
                    sep = "\t", check.names = FALSE, row.names = 1)

# 代谢物注释信息（含 KEGG Compound ID、Metabolite 名称）
# 与丰度矩阵行一一对应（同顺序：metab_0, metab_1, ... metab_236）
desc <- read.table("1.Preprocess/pos/metab_desc.txt", header = TRUE,
                   sep = "\t", check.names = FALSE, comment.char = "")

# =======================================================
# 2. OPLS-DA 核心函数
# =======================================================
run_oplsda <- function(treat_group) {

  message(sprintf("\n>>> 🚀 OPLS-DA: %s vs WT", treat_group))

  ck_cols <- grep("^WT_", colnames(abund), value = TRUE)
  tr_cols <- grep(paste0("^", treat_group, "_"), colnames(abund), value = TRUE)

  n_ck <- length(ck_cols)
  n_tr <- length(tr_cols)
  n_total <- n_ck + n_tr

  target_cols <- c(ck_cols, tr_cols)

  # 表达矩阵（保留行名）
  expr_matrix <- as.matrix(abund[, target_cols])
  expr_matrix[is.na(expr_matrix)] <- 0

  # 分组向量
  y_group <- factor(rep(c("WT", treat_group), times = c(n_ck, n_tr)),
                    levels = c("WT", treat_group))

  # 输出 PDF
  pdf_out <- sprintf("OPLSDA_Report_%s_vs_WT.pdf", treat_group)
  pdf(pdf_out, width = 8, height = 8)

  opls_mod <- suppressMessages(
    ropls::opls(x = t(expr_matrix), y = y_group, predI = 1, orthoI = 1,
                permI = 200, crossvalI = min(n_total, 7))
  )
  dev.off()

  # 统计计算
  # VIP: getVipVn 返回建模所用变量（已排除近零方差）的命名向量
  vips_named <- ropls::getVipVn(opls_mod)

  # 初始化为 NA，然后匹配填入
  vips_full <- rep(NA_real_, nrow(expr_matrix))
  names(vips_full) <- rownames(expr_matrix)
  common_names <- intersect(names(vips_full), names(vips_named))
  vips_full[common_names] <- vips_named[common_names]
  # 未进模型的代谢物 VIP 填充为 0（近零方差，无区分力）
  vips_full[is.na(vips_full)] <- 0

  # P值：log 空间 t 检验
  pvals <- apply(expr_matrix, 1, function(x) {
    tryCatch({ t.test(x[1:n_ck], x[(n_ck+1):n_total])$p.value },
             error = function(e) { 1 })
  })

  # Log10FC = 处理组均值 - WT均值
  log10fc <- rowMeans(expr_matrix[, (n_ck+1):n_total, drop = FALSE]) -
             rowMeans(expr_matrix[, 1:n_ck, drop = FALSE])

  log2fc <- log10fc / log10(2)

  # 阈值：1.5 倍差异 → log10(1.5) ≈ 0.176
  fc_threshold <- log10(1.5)

  # 获取代谢物名称
  metab_names <- rownames(expr_matrix)

  # 整理结果表（desc 与 abund 行一一对应，直接用位置合并注释）
  res_df <- data.frame(
    Name         = metab_names,
    metab_id     = desc$metab_id,
    Metabolite   = desc$Metabolite,
    VIP          = vips_full[1:length(metab_names)],
    Pvalue       = pvals[1:length(metab_names)],
    Log10FC      = log10fc[1:length(metab_names)],
    log2FC       = log2fc[1:length(metab_names)],
    KEGG_ID      = desc[["KEGG Compound ID"]],
    Mode         = desc$Mode,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      FDR = p.adjust(Pvalue, method = "BH"),
      Status = case_when(
        VIP > 1 & Pvalue < 0.05 & Log10FC >= fc_threshold  ~ "Up",
        VIP > 1 & Pvalue < 0.05 & Log10FC <= -fc_threshold ~ "Down",
        TRUE ~ "NS"
      )
    ) %>%
    arrange(desc(VIP)) %>%
    select(Name, metab_id, VIP, Pvalue, FDR, Log10FC, log2FC, Status,
           Metabolite, KEGG_ID, Mode)

  # 保存
  csv_out <- sprintf("OPLSDA_Diff_Results_%s_vs_WT.csv", treat_group)
  write.csv(res_df, csv_out, row.names = FALSE)

  message(sprintf("   ✅ %s 完成！显著: %d (Up: %d, Down: %d)",
                  treat_group,
                  sum(res_df$Status != "NS"),
                  sum(res_df$Status == "Up"),
                  sum(res_df$Status == "Down")))

  return(res_df)
}

# =======================================================
# 3. 执行两组对比
# =======================================================
res_ko2 <- run_oplsda("KO_2")
res_ko6 <- run_oplsda("KO_6")

# =======================================================
# 4. 汇总统计
# =======================================================
cat("\n\n========== 差异代谢物汇总 ==========\n")

cat("\nKO_2 vs WT:\n")
cat(sprintf("  总计: %d 代谢物\n", nrow(res_ko2)))
cat(sprintf("  上调: %d\n", sum(res_ko2$Status == "Up")))
cat(sprintf("  下调: %d\n", sum(res_ko2$Status == "Down")))
cat(sprintf("  不显著: %d\n", sum(res_ko2$Status == "NS")))

cat("\nKO_6 vs WT:\n")
cat(sprintf("  总计: %d 代谢物\n", nrow(res_ko6)))
cat(sprintf("  上调: %d\n", sum(res_ko6$Status == "Up")))
cat(sprintf("  下调: %d\n", sum(res_ko6$Status == "Down")))
cat(sprintf("  不显著: %d\n", sum(res_ko6$Status == "NS")))

message("\n>>> 🏆 OPLS-DA 任务结束！")
