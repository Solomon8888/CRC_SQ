# GSE211859 differential expression analysis with limma
#
# microarray: expression matrix -> optional normalizeBetweenArrays -> limma
# ngs: count matrix -> edgeR filtering/TMM -> voom -> limma


# 0. 可修改配置 ---------------------------------------------------------------

DATASET_ID <- "GSE211859"
DATA_TYPE <- "microarray"  # 可选："microarray" 或 "ngs"

SE_RDS_FILE <- "data/microarray/GSE211859/data_prepare/GSE211859_se_raw.rds"
CLINICAL_FILE <- "data/microarray/GSE211859/data_prepare/GSE211859_clinical_edit.csv"
FUNCTION_FILE <- "scripts/functions/limma_de_functions.R"

# 只配置实验组；对照组会从clinical文件的group列中自动识别。
EXPERIMENT_GROUP <- "OXLP"

# 显著差异筛选阈值
# P_VALUE_COLUMN可选："P.Value" 或 "adj.P.Val"
P_VALUE_COLUMN <- "adj.P.Val"
P_VALUE_CUTOFF <- 0.05
LOGFC_CUTOFF <- 1

OUTPUT_DIR <- file.path("results", DATA_TYPE, DATASET_ID, "tables")


# 1. 加载包和函数 --------------------------------------------------------------

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(limma)
})

if (DATA_TYPE == "ngs") {
  suppressPackageStartupMessages(library(edgeR))
}

source(FUNCTION_FILE)


# 2. 读取数据 -----------------------------------------------------------------

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

se <- readRDS(SE_RDS_FILE)
clinical_data <- read.csv(
  CLINICAL_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

class(se)
dim(se)
names(assays(se))
head(clinical_data)
table(clinical_data$group)


# 3. 检查样本信息 --------------------------------------------------------------

stopifnot(inherits(se, "SummarizedExperiment"))
stopifnot(DATA_TYPE %in% c("microarray", "ngs"))
stopifnot(all(c("Sample_ID", "group") %in% colnames(clinical_data)))
stopifnot(!any(duplicated(clinical_data$Sample_ID)))

exprSet <- get_assay_matrix(se, DATA_TYPE)

dim(exprSet)
head(exprSet[, 1:min(4, ncol(exprSet))])
summary(as.vector(exprSet))

stopifnot(is.numeric(exprSet))

missing_samples <- setdiff(colnames(exprSet), clinical_data$Sample_ID)
stopifnot(length(missing_samples) == 0)

sample_info <- clinical_data[
  match(colnames(exprSet), clinical_data$Sample_ID),
  ,
  drop = FALSE
]

head(sample_info)
table(sample_info$group)
stopifnot(all(sample_info$Sample_ID == colnames(exprSet)))


# 4. 设置分组和比较 ------------------------------------------------------------

group_names <- unique(sample_info$group)
control_group <- setdiff(group_names, EXPERIMENT_GROUP)

stopifnot(EXPERIMENT_GROUP %in% group_names)
stopifnot(length(control_group) == 1)

group_list <- factor(
  sample_info$group,
  levels = c(control_group, EXPERIMENT_GROUP)
)

table(group_list)

design <- model.matrix(~ 0 + group_list)
colnames(design) <- make.names(levels(group_list))
rownames(design) <- colnames(exprSet)

design

contrast_name <- paste0(EXPERIMENT_GROUP, "_vs_", control_group)
contrast_formula <- paste0(
  make.names(EXPERIMENT_GROUP),
  " - ",
  make.names(control_group)
)

contrast.matrix <- makeContrasts(
  contrasts = contrast_formula,
  levels = design
)
colnames(contrast.matrix) <- contrast_name

contrast.matrix


# 5. 准备基因注释 --------------------------------------------------------------

feature_id <- rownames(exprSet)

if (is.null(feature_id)) {
  feature_id <- paste0("Feature_", seq_len(nrow(exprSet)))
  rownames(exprSet) <- feature_id
}

gene_annotation <- data.frame(
  Feature_ID = feature_id,
  as.data.frame(rowData(se), stringsAsFactors = FALSE),
  check.names = FALSE
)
rownames(gene_annotation) <- rownames(exprSet)

head(gene_annotation)


# 6. 根据数据类型准备limma输入 -------------------------------------------------

if (DATA_TYPE == "microarray") {
  analysis_data <- prepare_microarray_data(exprSet)
  limma_input <- analysis_data$data
  genes_for_output <- gene_annotation
}

if (DATA_TYPE == "ngs") {
  analysis_data <- prepare_ngs_data(
    counts = exprSet,
    gene_annotation = gene_annotation,
    group_list = group_list,
    design = design
  )
  limma_input <- analysis_data$data
  genes_for_output <- limma_input$genes
}

if (DATA_TYPE == "microarray") {
  get_distribution_diagnostics(limma_input)
  analysis_data$log2_transformed
  analysis_data$normalized_between_arrays
}


# 7. limma差异分析 -------------------------------------------------------------

fit <- lmFit(limma_input, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

diff_results <- topTable(
  fit2,
  coef = contrast_name,
  number = Inf,
  adjust.method = "BH",
  sort.by = "P",
  genelist = genes_for_output
)

dim(diff_results)
head(diff_results)
colnames(diff_results)


# 8. 筛选显著差异基因 ----------------------------------------------------------

stopifnot(P_VALUE_COLUMN %in% colnames(diff_results))

up_index <- diff_results$logFC > LOGFC_CUTOFF &
  diff_results[[P_VALUE_COLUMN]] < P_VALUE_CUTOFF

down_index <- diff_results$logFC < -LOGFC_CUTOFF &
  diff_results[[P_VALUE_COLUMN]] < P_VALUE_CUTOFF

significant_results <- diff_results[up_index | down_index, , drop = FALSE]

de_summary <- data.frame(
  Up = sum(up_index),
  Down = sum(down_index),
  Not_Significant = nrow(diff_results) - sum(up_index) - sum(down_index)
)

de_summary
dim(significant_results)
head(significant_results)


# 9. 保存结果 -----------------------------------------------------------------

all_results_file <- file.path(
  OUTPUT_DIR,
  paste0(DATASET_ID, "_limma_", contrast_name, "_all_genes.csv")
)

significant_results_file <- file.path(
  OUTPUT_DIR,
  paste0(DATASET_ID, "_limma_", contrast_name, "_significant_genes.csv")
)

summary_file <- file.path(
  OUTPUT_DIR,
  paste0(DATASET_ID, "_limma_", contrast_name, "_summary.csv")
)

summary_table <- data.frame(
  Dataset = DATASET_ID,
  Data_Type = DATA_TYPE,
  Contrast = contrast_name,
  Control_Group = control_group,
  Experiment_Group = EXPERIMENT_GROUP,
  Total_Genes_Input = nrow(exprSet),
  Total_Genes_Analyzed = nrow(diff_results),
  Up = de_summary$Up,
  Down = de_summary$Down,
  Not_Significant = de_summary$Not_Significant,
  P_Value_Column = P_VALUE_COLUMN,
  P_Value_Cutoff = P_VALUE_CUTOFF,
  LogFC_Cutoff = LOGFC_CUTOFF,
  Microarray_Log2_Transformed = analysis_data$log2_transformed,
  Microarray_Normalized_Between_Arrays = analysis_data$normalized_between_arrays,
  Microarray_Final_Median_Spread = analysis_data$median_spread,
  Microarray_Final_IQR_Spread = analysis_data$iqr_spread,
  NGS_Filtered_Genes = analysis_data$filtered_genes,
  stringsAsFactors = FALSE
)

summary_table

write.csv(diff_results, all_results_file, row.names = FALSE)
write.csv(significant_results, significant_results_file, row.names = FALSE)
write.csv(summary_table, summary_file, row.names = FALSE)
