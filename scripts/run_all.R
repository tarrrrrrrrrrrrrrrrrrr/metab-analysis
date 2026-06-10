# =======================================================
# 🌟 代谢组学 R 一键分析（主入口）
# 用法: Rscript run_all.R <数据目录>
# =======================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  work_dir <- args[1]
} else {
  work_dir <- getwd()
}
setwd(work_dir)
cat(sprintf(">>> 📂 数据目录: %s\n", work_dir))

script_dir <- dirname(normalizePath(sys.frame(1)$ofile))
# 如果上述不适用，使用脚本所在目录
if (is.na(script_dir) || script_dir == ".") {
  script_dir <- dirname(normalizePath(
    sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))
  ))
}

cat(sprintf(">>> 📜 脚本目录: %s\n", script_dir))

# Step 1: OPLS-DA（必须最先跑，生成 CSV）
cat("\n\n############ Step 1/4: OPLS-DA 差异分析 ############\n")
source(file.path(script_dir, "02_OPLSDA_analysis.R"))

# Step 2: PCA
cat("\n\n############ Step 2/4: PCA 分析 ############\n")
source(file.path(script_dir, "01_PCA_analysis.R"))

# Step 3: Volcano
cat("\n\n############ Step 3/4: 火山图 ############\n")
source(file.path(script_dir, "03_Volcano_plot.R"))

# Step 4: Heatmap
cat("\n\n############ Step 4/4: 差异代谢物热图 ############\n")
source(file.path(script_dir, "04_Heatmap_plot.R"))

cat("\n\n========================================\n")
cat(">>> 🏆 全部分析完成！ <<<\n")
cat("========================================\n")
