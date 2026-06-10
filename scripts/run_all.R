# =======================================================
# 🌟 代谢组学 R 一键分析（主入口）
# 用法: Rscript run_all.R <数据目录> [organism_code]
# organism_code: KEGG 物种代码，默认 dosa (水稻)
# =======================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  work_dir <- args[1]
} else {
  work_dir <- getwd()
}
org_code <- if (length(args) > 1) args[2] else "dosa"
setwd(work_dir)
cat(sprintf(">>> 📂 数据目录: %s\n", work_dir))
cat(sprintf(">>> 🧬 KEGG 物种: %s\n", org_code))

script_dir <- dirname(normalizePath(sys.frame(1)$ofile))
if (is.na(script_dir) || script_dir == ".") {
  script_dir <- dirname(normalizePath(
    sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))
  ))
}

cat(sprintf(">>> 📜 脚本目录: %s\n", script_dir))

# Step 1: OPLS-DA（必须最先跑，生成 CSV）
cat("\n\n############ Step 1/6: OPLS-DA 差异分析 ############\n")
source(file.path(script_dir, "02_OPLSDA_analysis.R"))

# Step 2: PCA
cat("\n\n############ Step 2/6: PCA 分析 ############\n")
source(file.path(script_dir, "01_PCA_analysis.R"))

# Step 3: Volcano
cat("\n\n############ Step 3/6: 火山图 ############\n")
source(file.path(script_dir, "03_Volcano_plot.R"))

# Step 4: Heatmap
cat("\n\n############ Step 4/6: 差异代谢物热图 ############\n")
source(file.path(script_dir, "04_Heatmap_plot.R"))

# Step 5: KEGG 靶点提取
cat("\n\n############ Step 5/6: KEGG 靶点提取 ############\n")
source(file.path(script_dir, "05_KEGG_targets.R"))

# Step 6: KEGG 通路富集 + 气泡图
cat("\n\n############ Step 6/6: KEGG 通路富集 ############\n")
# Pass org_code to the bubble script via a temp env var
Sys.setenv(METAB_KEGG_ORG = org_code)
source(file.path(script_dir, "06_KEGG_bubble.R"))

cat("\n\n========================================\n")
cat(">>> 🏆 全部分析完成！ <<<\n")
cat("========================================\n")
