# =======================================================
# 🌟 KEGG 靶点提取：从 OPLS-DA 差异结果中提取 KEGG ID
# 适配自 新代码/kegg通路分析.txt
# =======================================================
rm(list=ls())
graphics.off()

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
work_dir <- if (length(args) > 0) args[1] else getwd()
setwd(work_dir)
message(">>> 📂 数据目录: ", work_dir)

# =======================================================
# 自动发现对比组（扫描 OPLSDA_Diff_Results_*_vs_*.csv）
# =======================================================
csv_files <- list.files(pattern = "^OPLSDA_Diff_Results_.*_vs_.*\\.csv$")
if (length(csv_files) == 0) {
  stop("❌ 未找到 OPLSDA_Diff_Results_*.csv，请先运行 OPLS-DA 分析。")
}

# 从文件名提取对比组
parse_comp <- function(fname) {
  m <- regmatches(fname, regexec("OPLSDA_Diff_Results_(.+)_vs_(.+)\\.csv", fname))[[1]]
  if (length(m) == 3) {
    return(list(treat = m[2], ctrl = m[3]))
  }
  NULL
}

cat(sprintf(">>> 发现 %d 个差异结果文件\n", length(csv_files)))
for (f in csv_files) cat(sprintf("    %s\n", f))

# =======================================================
# 核心函数：提取并去重 KEGG ID
# =======================================================
extract_kegg_targets <- function(treat_group) {

  # 找对应 CSV（正离子模式 — 仅 pos 数据）
  res_file <- sprintf("OPLSDA_Diff_Results_%s_vs_WT.csv", treat_group)

  if (!file.exists(res_file)) {
    message(sprintf("⚠️ 未找到: %s，跳过", res_file))
    return(NULL)
  }

  res <- read.csv(res_file, stringsAsFactors = FALSE)

  # 筛选显著差异代谢物
  df_sig <- res %>%
    filter(Status %in% c("Up", "Down"))

  if (nrow(df_sig) == 0) {
    message(sprintf("⚠️ %s 无显著差异代谢物", treat_group))
    return(NULL)
  }

  # 提取 KEGG ID 并去重
  # 列名可能是 KEGG_ID 或 KEGG Compound ID
  kegg_col <- intersect(c("KEGG_ID", "KEGG Compound ID"), colnames(df_sig))[1]

  if (is.na(kegg_col)) {
    message(sprintf("⚠️ %s 结果表中没有 KEGG 列", treat_group))
    return(NULL)
  }

  df_clean <- df_sig %>%
    rename(KEGG_ID = all_of(kegg_col)) %>%
    filter(!is.na(KEGG_ID) & trimws(KEGG_ID) != "" & trimws(KEGG_ID) != "-") %>%
    mutate(KEGG_ID = sub("[;,].*", "", trimws(KEGG_ID))) %>%
    arrange(desc(VIP)) %>%
    distinct(KEGG_ID, .keep_all = TRUE) %>%
    select(Name, metab_id, KEGG_ID, VIP, Pvalue, FDR, Log10FC, log2FC, Status, Metabolite)

  out_file <- sprintf("04_KEGG_Target_%s_vs_WT.csv", treat_group)
  write.csv(df_clean, out_file, row.names = FALSE)

  message(sprintf("   ✅ %s: %d 个 KEGG ID（%d 上调, %d 下调）",
                  treat_group, nrow(df_clean),
                  sum(df_clean$Status == "Up"),
                  sum(df_clean$Status == "Down")))

  return(df_clean)
}

# =======================================================
# 批量处理所有对比组
# =======================================================
comps <- sapply(seq_along(csv_files), function(i) parse_comp(csv_files[[i]]), simplify = FALSE)

cat("\n========== KEGG 靶点提取 ==========\n")
for (comp in comps) {
  extract_kegg_targets(comp$treat)
}

message("\n>>> 🏆 KEGG 靶点提取完成！")
