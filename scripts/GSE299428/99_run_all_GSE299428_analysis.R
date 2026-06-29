# GSE299428分析流水线一键运行脚本
#
# 默认顺序运行参考模板同款全套分析：
# 00 样本聚类热图 -> 01 limma/edgeR差异分析 -> 02 显著差异基因交集 ->
# 03 火山图 -> 04 多组火山图 -> 05 Top DEG热图 ->
# 06 GSEA运算 -> 07 GSEA绘图 -> 08 TF富集 -> 09 TF结果整合。
#
# 注意：01号及之后脚本需要先在
# data/ngs/GSE299428/data_prepare/GSE299428_clinical_edit.csv
# 中增加至少一个analysis_开头的实验设计列。


# 0. 可修改配置 ---------------------------------------------------------------

SCRIPTS_TO_RUN <- c(
  "scripts/GSE299428/00_sample_clustering_heatmap.R",
  "scripts/GSE299428/01_limma_differential_expression.R",
  "scripts/GSE299428/02_intersect_significant_genes.R",
  "scripts/GSE299428/03_volcano_plot.R",
  "scripts/GSE299428/04_multiple_volcano_plot.R",
  "scripts/GSE299428/05_top_deg_gene_heatmap.R",
  "scripts/GSE299428/06_gsea_analysis.R",
  "scripts/GSE299428/07_gsea_plot.R",
  "scripts/GSE299428/08_tf_enrichment_analysis.R",
  "scripts/GSE299428/09_integrate_tf_enrichment_results.R"
)

options(width = 200)


# 1. 运行 ----------------------------------------------------------------------

cat("\nRunning GSE299428 analysis pipeline...\n")
cat("Scripts: ", length(SCRIPTS_TO_RUN), "\n\n", sep = "")

rscript_bin <- file.path(R.home("bin"), "Rscript")

for (script_file in SCRIPTS_TO_RUN) {
  if (!file.exists(script_file)) {
    stop("Missing script: ", script_file)
  }

  cat("\n============================================================\n")
  cat("Running: ", script_file, "\n", sep = "")
  cat("============================================================\n")

  exit_status <- system2(rscript_bin, args = script_file)

  if (!identical(exit_status, 0L)) {
    stop("Script failed: ", script_file, " (exit status ", exit_status, ")")
  }
}

cat("\nGSE299428 analysis pipeline finished.\n")
