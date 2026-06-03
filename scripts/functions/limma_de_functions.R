# Common helper functions for limma differential expression analysis.

get_assay_matrix <- function(se, data_type) {
  assay_names <- names(SummarizedExperiment::assays(se))

  if (data_type == "ngs" && "counts" %in% assay_names) {
    return(as.matrix(SummarizedExperiment::assay(se, "counts")))
  }

  if ("expr" %in% assay_names) {
    return(as.matrix(SummarizedExperiment::assay(se, "expr")))
  }

  as.matrix(SummarizedExperiment::assay(se, 1))
}

needs_log2_transform <- function(exprSet) {
  q99 <- as.numeric(quantile(exprSet, 0.99, na.rm = TRUE))
  max_value <- max(exprSet, na.rm = TRUE)

  q99 > 50 || max_value > 100
}

get_distribution_diagnostics <- function(exprSet) {
  sample_median <- apply(exprSet, 2, median, na.rm = TRUE)
  sample_iqr <- apply(exprSet, 2, IQR, na.rm = TRUE)

  data.frame(
    Median_Spread = diff(range(sample_median)),
    IQR_Spread = diff(range(sample_iqr))
  )
}

prepare_microarray_data <- function(exprSet) {
  log2_transformed <- FALSE
  normalized_between_arrays <- FALSE

  if (needs_log2_transform(exprSet)) {
    if (min(exprSet, na.rm = TRUE) < 0) {
      stop("Expression values look unlogged, but contain negative values.")
    }

    exprSet <- log2(exprSet + 1)
    log2_transformed <- TRUE
  }

  diagnostics <- get_distribution_diagnostics(exprSet)

  if (diagnostics$Median_Spread > 0.5 || diagnostics$IQR_Spread > 0.5) {
    exprSet <- limma::normalizeBetweenArrays(exprSet, method = "quantile")
    normalized_between_arrays <- TRUE
    diagnostics <- get_distribution_diagnostics(exprSet)
  }

  list(
    data = exprSet,
    log2_transformed = log2_transformed,
    normalized_between_arrays = normalized_between_arrays,
    median_spread = diagnostics$Median_Spread,
    iqr_spread = diagnostics$IQR_Spread,
    filtered_genes = NA_integer_
  )
}

prepare_ngs_data <- function(counts, gene_annotation, group_list, design) {
  y <- edgeR::DGEList(
    counts = counts,
    group = group_list,
    genes = gene_annotation
  )

  keep <- edgeR::filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y)

  v <- limma::voom(y, design = design, plot = FALSE)

  list(
    data = v,
    log2_transformed = NA,
    normalized_between_arrays = NA,
    median_spread = NA,
    iqr_spread = NA,
    filtered_genes = sum(!keep)
  )
}
