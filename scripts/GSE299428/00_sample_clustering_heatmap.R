# GSE299428样本聚类热图
#
# 从当前数据集的SummarizedExperiment对象中提取全部样本，
# 使用TPM表达量计算全样本相关性，并绘制层次聚类热图。
# 样本显示名称来自临床信息表中的Title列。


# 0. 可修改配置 ---------------------------------------------------------------

# 当前脚本只服务于GSE299428这个NGS数据集；目录结构也依赖这两个字段。
DATASET_ID <- "GSE299428"
DATA_TYPE <- "ngs"

# 输入文件：SE对象提供TPM矩阵和基因注释，临床表提供样本Title和样本分组。
SE_RDS_FILE <- "data/ngs/GSE299428/data_prepare/GSE299428_se_raw.rds"
CLINICAL_FILE <- "data/ngs/GSE299428/data_prepare/GSE299428_clinical_edit.csv"

# FUNCTION_FILE用于基因类型筛选；
# PLOTTING_FUNCTION_FILE保存跨绘图脚本共用的风格配置和基础函数。
FUNCTION_FILE <- "scripts/functions/limma_de_functions.R"
PLOTTING_FUNCTION_FILE <- "scripts/functions/plotting_common_functions.R"
PARALLEL_FUNCTION_FILE <- "scripts/functions/parallel_runtime_functions.R"

# 输出根目录。最终图片保存到plots/sample_clustering_heatmap/<gene_biotype>/ALL_SAMPLES/。
RESULT_ROOT <- file.path("results", DATA_TYPE, DATASET_ID)
PLOT_ROOT <- file.path(RESULT_ROOT, "plots", "sample_clustering_heatmap")

# 重跑时清理当前热图输出目录里的旧图片，避免新旧命名混在一起。
CLEAN_PLOT_OUTPUT_DIR <- TRUE

# 基因类型筛选。可选："coding", "protein", "protein_coding", "non_coding", "all"。
# 这里沿用01号差异分析脚本的基因类型筛选逻辑。
GENE_BIOTYPE_FILTER <- "coding"

# SE对象中TPM矩阵的assay名称。聚类热图使用TPM，不使用原始count。
TPM_ASSAY_NAME <- "tpm"

# 样本筛选列和值。GSE299428的样本聚类默认使用全部样本，
# 以便观察PAR/LD/MD/HD和HCT116/SW480在全局表达空间中的自然聚类趋势。
# 若后续确实需要画某个子集，可把SAMPLE_GROUPS改成对应列里的具体取值。
SAMPLE_GROUP_COLUMN <- "dose_group"
ALL_SAMPLE_GROUP_VALUE <- "__ALL_SAMPLES__"
SAMPLE_GROUPS <- c(
  ALL_SAMPLES = ALL_SAMPLE_GROUP_VALUE
)

# 图中展示的样本名来源。若Title为空，会自动回退到Sample_ID。
SAMPLE_LABEL_COLUMN <- "Title"

# 全样本热图顶部和左侧的样本注释条。样本名本身已经展示剂量和重复编号；
# 这里额外用颜色标记剂量组和细胞系，便于快速判断聚类趋势。
SAMPLE_ANNOTATION_COLUMNS <- c(
  Dose_group = "dose_group",
  Cell_line = "cell_line"
)
DOSE_GROUP_COLORS <- c(
  PAR = "#4D4D4D",
  LD = "#2C7BB6",
  MD = "#F28E2B",
  HD = "#D7191C"
)

# 样本相关性和层次聚类方法。聚类只基于TPM表达模式，不按细胞系分组强行排序。
# 全样本质控聚类默认使用高变异基因，避免低信息量基因稀释样本间距离。
CORRELATION_METHOD <- "pearson"
CLUSTERING_METHOD <- "average"
TOP_VARIABLE_GENE_N <- 5000

# 样本名较长时按45字符自动换行；列名垂直显示以避免24个长样本名互相重叠。
ROW_LABEL_WIDTH <- 45
COLUMN_LABEL_WIDTH <- 45
COLUMN_NAMES_ROT <- 90

# 热图主体中每个小格子的边长。行列样本数一致时，热图主体保持正方形。
CELL_SIZE_MM <- 10.5

# 聚类树和样本注释条的空间。聚类树只反映TPM表达相关性，不按分组强行排序。
ROW_DEND_WIDTH_MM <- 28
COLUMN_DEND_HEIGHT_MM <- 24
SAMPLE_ANNOTATION_SIZE_MM <- 4
DENDROGRAM_LINE_WIDTH_SCALE <- 0.95

# 加粗倍率。字体使用bold，线条适度加粗；全样本图中字号不能过大，否则长样本名会重叠。
BOLDNESS_MULTIPLIER <- 2.2
BASE_CELL_BORDER_WIDTH <- 0.55
BASE_SAMPLE_FONT_SIZE <- 3.9
BASE_ANNOTATION_FONT_SIZE <- 4.0
BASE_LEGEND_FONT_SIZE <- 4.0

CELL_BORDER_WIDTH <- BASE_CELL_BORDER_WIDTH * BOLDNESS_MULTIPLIER
DENDROGRAM_LINE_WIDTH <- CELL_BORDER_WIDTH * DENDROGRAM_LINE_WIDTH_SCALE

# 样本名与热图之间的空格，避免粗体样本名贴到热图本体。
LABEL_HEATMAP_GAP_SPACES <- 2

# 红白蓝渐变配色。需要换色时只改这里即可。
HEATMAP_COLOR_LOW <- "#0d0dbb7f"   # 深蓝
HEATMAP_COLOR_MID <- "#FFFFFF"     # 白色
HEATMAP_COLOR_HIGH <- "#cd0e0e"    # 鲜红
CORRELATION_COLOR_MIN <- 0.80
CORRELATION_COLOR_MAX <- 1.00

# 图片大小会根据样本数量和样本名长度自动调整；上下限用于避免图片过小或过大。
MIN_PDF_WIDTH <- 8.0
MIN_PDF_HEIGHT <- 8.0
MAX_PDF_WIDTH <- 34.0
MAX_PDF_HEIGHT <- 30.0

SAMPLE_FONT_SIZE <- BASE_SAMPLE_FONT_SIZE * BOLDNESS_MULTIPLIER
COLUMN_FONT_SCALE <- 1.00
COLUMN_SAMPLE_FONT_SIZE <- SAMPLE_FONT_SIZE * COLUMN_FONT_SCALE
ANNOTATION_FONT_SIZE <- BASE_ANNOTATION_FONT_SIZE * BOLDNESS_MULTIPLIER
LEGEND_FONT_SIZE <- BASE_LEGEND_FONT_SIZE * BOLDNESS_MULTIPLIER

# 图例放在整体图片左侧，并与热图主体保持适度间距。
# 这些参数只控制整体排版，不改变热图本体的聚类结果。
OUTER_MARGIN_INCH <- 0.35
LEGEND_LEFT_WIDTH_INCH <- 1.65
LEGEND_HEATMAP_GAP_INCH <- 0.45
LEGEND_INNER_MARGIN_INCH <- 0.05
LEGEND_TOP_EXTRA_MM <- 26
LEGEND_ITEM_GAP_MM <- 5
LEGEND_GROUP_GAP_MM <- 18
CORRELATION_LEGEND_HEIGHT_MM <- 50
LEGEND_TOP_MARGIN_INCH <- (
  COLUMN_DEND_HEIGHT_MM + SAMPLE_ANNOTATION_SIZE_MM + LEGEND_TOP_EXTRA_MM
) / 25.4

LABEL_SPACE_PADDING_INCH <- 0.28
ROW_LABEL_MIN_SPACE_INCH <- 2.3
COLUMN_LABEL_MIN_SPACE_INCH <- 2.7
LABEL_LINE_HEIGHT <- 1.08

options(width = 200)


# 1. 加载包和函数 --------------------------------------------------------------

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(RColorBrewer)
})

source(FUNCTION_FILE)
source(PLOTTING_FUNCTION_FILE)
source(PARALLEL_FUNCTION_FILE)

SCRIPT_START_TIME <- start_runtime_timer()


# 2. 读取数据 ------------------------------------------------------------------

se <- readRDS(SE_RDS_FILE)
clinical_data <- read.csv(
  CLINICAL_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(inherits(se, "SummarizedExperiment"))
stopifnot("Sample_ID" %in% colnames(clinical_data))
stopifnot(SAMPLE_LABEL_COLUMN %in% colnames(clinical_data))
stopifnot(SAMPLE_GROUP_COLUMN %in% colnames(clinical_data))
stopifnot(all(unname(SAMPLE_ANNOTATION_COLUMNS) %in% colnames(clinical_data)))
stopifnot(!any(duplicated(clinical_data$Sample_ID)))

missing_samples <- setdiff(colnames(se), clinical_data$Sample_ID)
stopifnot(length(missing_samples) == 0)

sample_info_all <- clinical_data[
  match(colnames(se), clinical_data$Sample_ID),
  ,
  drop = FALSE
]
rownames(sample_info_all) <- sample_info_all$Sample_ID
stopifnot(all(sample_info_all$Sample_ID == colnames(se)))


# 3. 准备TPM表达矩阵 -----------------------------------------------------------

assay_names <- names(SummarizedExperiment::assays(se))
stopifnot(TPM_ASSAY_NAME %in% assay_names)

tpm_all <- as.matrix(SummarizedExperiment::assay(se, TPM_ASSAY_NAME))
stopifnot(is.numeric(tpm_all))

feature_id <- rownames(tpm_all)
if (is.null(feature_id)) {
  feature_id <- paste0("Feature_", seq_len(nrow(tpm_all)))
  rownames(tpm_all) <- feature_id
}

gene_annotation <- data.frame(
  Feature_ID = feature_id,
  as.data.frame(rowData(se), stringsAsFactors = FALSE),
  check.names = FALSE
)
rownames(gene_annotation) <- rownames(tpm_all)

gene_biotype_filter <- trimws(tolower(GENE_BIOTYPE_FILTER))
gene_biotype_filter <- gsub("-", "_", gene_biotype_filter)
if (gene_biotype_filter %in% c("protein", "protein_coding")) {
  gene_biotype_filter <- "coding"
}

gene_filter <- filter_genes_by_biotype(
  exprSet = tpm_all,
  gene_annotation = gene_annotation,
  biotype_filter = gene_biotype_filter
)

# 输出目录用基因类型区分；图片文件名只保留结果类型。
plot_filter_name <- sanitize_file_name(gene_filter$filter)
PLOT_DIR <- file.path(PLOT_ROOT, plot_filter_name)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

if (CLEAN_PLOT_OUTPUT_DIR) {
  unlink(list.files(
    PLOT_DIR,
    pattern = "[.](pdf|png)$",
    recursive = TRUE,
    full.names = TRUE
  ))

  legacy_plot_files <- list.files(
    PLOT_ROOT,
    pattern = "[.](pdf|png)$",
    full.names = TRUE
  )
  if (length(legacy_plot_files) > 0) {
    unlink(legacy_plot_files)
  }
}

tpm_filtered <- gene_filter$exprSet


# 4. 样本集合聚类热图绘制函数 --------------------------------------------------

wrap_sample_label <- function(x, width = 45) {
  # 使用词边界优先的45字符换行；只有单个词超长时才强制断词。
  x <- as.character(x)

  vapply(x, function(label) {
    wrapped <- strwrap(
      label,
      width = width,
      simplify = FALSE
    )[[1]]

    if (length(wrapped) == 0) {
      return(label)
    }

    paste(wrapped, collapse = "\n")
  }, character(1))
}

measure_multiline_label_size <- function(labels, fontsize, line_height = 1.08) {
  # 根据换行后的实际文字宽度/高度估算标签占位，避免长样本名和热图主体互相挤压。
  label_parts <- strsplit(labels, "\n", fixed = TRUE)

  label_widths <- vapply(label_parts, function(parts) {
    max(grid::convertWidth(grid::stringWidth(parts), "in", valueOnly = TRUE))
  }, numeric(1))

  label_heights <- vapply(label_parts, function(parts) {
    length(parts) * fontsize / 72 * line_height
  }, numeric(1))

  list(
    width = max(label_widths),
    height = max(label_heights),
    max_lines = max(lengths(label_parts))
  )
}

select_top_variable_genes <- function(expr_matrix, top_n = 5000) {
  # 样本聚类常用高变异基因以突出样本间主要结构；不使用分组信息，不改变无监督性质。
  log_expr <- log2(expr_matrix + 1)
  gene_sd <- apply(log_expr, 1, sd, na.rm = TRUE)
  valid_index <- is.finite(gene_sd) & gene_sd > 0
  valid_gene_ids <- rownames(expr_matrix)[valid_index]

  if (length(valid_gene_ids) <= 1) {
    stop("Too few non-zero-variance genes for sample clustering.")
  }

  if (is.finite(top_n) && length(valid_gene_ids) > top_n) {
    ordered_gene_ids <- names(sort(gene_sd[valid_index], decreasing = TRUE))
    selected_gene_ids <- ordered_gene_ids[seq_len(top_n)]
  } else {
    selected_gene_ids <- valid_gene_ids
  }

  list(
    expr_matrix = expr_matrix[selected_gene_ids, , drop = FALSE],
    nonzero_gene_count = length(valid_gene_ids),
    selected_gene_count = length(selected_gene_ids),
    top_variable_gene_n = ifelse(is.finite(top_n), top_n, NA_integer_)
  )
}

draw_sample_clustering_heatmap <- function(sample_group_name, sample_group_value) {
  # 对一个样本集合完成筛选、相关性计算、层次聚类、排版和保存。
  # 聚类只基于TPM表达谱得到的样本相关性，不按细胞系分组强行排序。
  if (identical(sample_group_value, ALL_SAMPLE_GROUP_VALUE)) {
    sample_index <- rep(TRUE, nrow(sample_info_all))
  } else {
    sample_status <- trimws(as.character(sample_info_all[[SAMPLE_GROUP_COLUMN]]))
    sample_status[is.na(sample_status)] <- ""
    sample_index <- sample_status == sample_group_value
  }

  sample_info <- sample_info_all[sample_index, , drop = FALSE]
  stopifnot(nrow(sample_info) >= 2)

  sample_labels <- get_display_labels(
    sample_info = sample_info,
    label_column = SAMPLE_LABEL_COLUMN
  )

  tpm <- tpm_filtered[, sample_info$Sample_ID, drop = FALSE]
  variable_gene_result <- select_top_variable_genes(
    expr_matrix = tpm,
    top_n = TOP_VARIABLE_GENE_N
  )
  tpm_for_clustering <- variable_gene_result$expr_matrix

  correlation_result <- prepare_sample_correlation(
    expr_matrix = tpm_for_clustering,
    correlation_method = CORRELATION_METHOD,
    clustering_method = CLUSTERING_METHOD
  )
  expr_for_correlation <- correlation_result$expr_for_correlation
  cor_matrix <- correlation_result$cor_matrix
  sample_hclust <- correlation_result$sample_hclust

  plot_row_labels <- wrap_sample_label(sample_labels, ROW_LABEL_WIDTH)
  plot_column_labels <- wrap_sample_label(sample_labels, COLUMN_LABEL_WIDTH)

  label_padding <- paste(rep(" ", LABEL_HEATMAP_GAP_SPACES), collapse = "")
  plot_row_labels <- paste0(label_padding, plot_row_labels)
  plot_column_labels <- paste0(label_padding, plot_column_labels)

  rownames(cor_matrix) <- sample_info$Sample_ID
  colnames(cor_matrix) <- sample_info$Sample_ID

  annotation_data <- data.frame(
    .sample_id = sample_info$Sample_ID,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(annotation_data) <- sample_info$Sample_ID
  annotation_color_list <- list()

  for (annotation_name in names(SAMPLE_ANNOTATION_COLUMNS)) {
    annotation_column <- SAMPLE_ANNOTATION_COLUMNS[[annotation_name]]
    annotation_values <- trimws(as.character(sample_info[[annotation_column]]))
    annotation_values[is.na(annotation_values) | annotation_values == ""] <- "Unknown"
    annotation_data[[annotation_name]] <- annotation_values

    if (identical(annotation_column, "dose_group")) {
      annotation_levels <- unique(annotation_values)
      known_levels <- intersect(names(DOSE_GROUP_COLORS), annotation_levels)
      extra_levels <- setdiff(annotation_levels, names(DOSE_GROUP_COLORS))
      annotation_palette <- DOSE_GROUP_COLORS[known_levels]

      if (length(extra_levels) > 0) {
        annotation_palette <- c(
          annotation_palette,
          get_named_brewer_palette(extra_levels, palette = "Set3")
        )
      }

      annotation_color_list[[annotation_name]] <- annotation_palette[annotation_levels]
    } else {
      annotation_color_list[[annotation_name]] <- get_named_brewer_palette(
        annotation_values
      )
    }
  }

  annotation_data$.sample_id <- NULL
  annotation_track_size_inch <- ncol(annotation_data) * SAMPLE_ANNOTATION_SIZE_MM / 25.4

  top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    df = annotation_data,
    col = annotation_color_list,
    show_legend = FALSE,
    show_annotation_name = FALSE,
    annotation_name_gp = grid::gpar(
      fontsize = ANNOTATION_FONT_SIZE,
      fontface = TEXT_FONT_FACE,
      fontfamily = TEXT_FONT_FAMILY,
      col = TEXT_COLOR
    ),
    simple_anno_size = grid::unit(SAMPLE_ANNOTATION_SIZE_MM, "mm"),
    gp = grid::gpar(col = "black", lwd = CELL_BORDER_WIDTH),
    border = TRUE
  )

  left_annotation <- ComplexHeatmap::rowAnnotation(
    df = annotation_data,
    col = annotation_color_list,
    show_annotation_name = FALSE,
    show_legend = FALSE,
    simple_anno_size = grid::unit(SAMPLE_ANNOTATION_SIZE_MM, "mm"),
    gp = grid::gpar(col = "black", lwd = CELL_BORDER_WIDTH),
    border = TRUE
  )

  n_samples <- ncol(cor_matrix)
  heatmap_body_inch <- n_samples * CELL_SIZE_MM / 25.4
  row_label_size <- measure_multiline_label_size(
    labels = plot_row_labels,
    fontsize = SAMPLE_FONT_SIZE,
    line_height = LABEL_LINE_HEIGHT
  )
  column_label_size <- measure_multiline_label_size(
    labels = plot_column_labels,
    fontsize = COLUMN_SAMPLE_FONT_SIZE,
    line_height = LABEL_LINE_HEIGHT
  )

  column_rotation_rad <- COLUMN_NAMES_ROT * pi / 180
  rotated_column_label_height <- column_label_size$width * abs(sin(column_rotation_rad)) +
    column_label_size$height * abs(cos(column_rotation_rad))

  row_label_space <- max(
    ROW_LABEL_MIN_SPACE_INCH,
    row_label_size$width + LABEL_SPACE_PADDING_INCH
  )
  col_label_space <- max(
    COLUMN_LABEL_MIN_SPACE_INCH,
    rotated_column_label_height + LABEL_SPACE_PADDING_INCH
  )

  heatmap_panel_width <- heatmap_body_inch + row_label_space +
    ROW_DEND_WIDTH_MM / 25.4 + annotation_track_size_inch + 0.95
  heatmap_panel_height <- heatmap_body_inch + col_label_space +
    COLUMN_DEND_HEIGHT_MM / 25.4 + annotation_track_size_inch + 0.95

  pdf_width <- OUTER_MARGIN_INCH * 2 +
    LEGEND_LEFT_WIDTH_INCH + LEGEND_HEATMAP_GAP_INCH + heatmap_panel_width
  pdf_height <- OUTER_MARGIN_INCH * 2 + heatmap_panel_height

  pdf_width <- min(max(pdf_width, MIN_PDF_WIDTH), MAX_PDF_WIDTH)
  pdf_height <- min(max(pdf_height, MIN_PDF_HEIGHT), MAX_PDF_HEIGHT)

  output_dir <- file.path(PLOT_DIR, sanitize_file_name(sample_group_name))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_file <- file.path(output_dir, "heatmap.pdf")

  color_min <- CORRELATION_COLOR_MIN
  color_max <- CORRELATION_COLOR_MAX
  color_mid <- (color_min + color_max) / 2

  heatmap_colors <- circlize::colorRamp2(
    c(color_min, color_mid, color_max),
    c(HEATMAP_COLOR_LOW, HEATMAP_COLOR_MID, HEATMAP_COLOR_HIGH)
  )

  heatmap_body_size <- grid::unit(n_samples * CELL_SIZE_MM, "mm")

  ht <- ComplexHeatmap::Heatmap(
    cor_matrix,
    name = "Correlation",
    col = heatmap_colors,
    cluster_rows = sample_hclust,
    cluster_columns = sample_hclust,
    top_annotation = top_annotation,
    left_annotation = left_annotation,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_labels = plot_row_labels,
    column_labels = plot_column_labels,
    row_names_side = "right",
    column_names_side = "bottom",
    row_names_max_width = grid::unit(row_label_space, "in"),
    column_names_max_height = grid::unit(col_label_space, "in"),
    row_names_gp = grid::gpar(
      fontsize = SAMPLE_FONT_SIZE,
      fontface = TEXT_FONT_FACE,
      fontfamily = TEXT_FONT_FAMILY,
      col = TEXT_COLOR
    ),
    column_names_gp = grid::gpar(
      fontsize = COLUMN_SAMPLE_FONT_SIZE,
      fontface = TEXT_FONT_FACE,
      fontfamily = TEXT_FONT_FAMILY,
      col = TEXT_COLOR
    ),
    row_names_centered = TRUE,
    column_names_rot = COLUMN_NAMES_ROT,
    column_names_centered = TRUE,
    rect_gp = grid::gpar(
      col = "black",
      lwd = CELL_BORDER_WIDTH
    ),
    border = TRUE,
    border_gp = grid::gpar(col = "black", lwd = CELL_BORDER_WIDTH),
    row_dend_width = grid::unit(ROW_DEND_WIDTH_MM, "mm"),
    column_dend_height = grid::unit(COLUMN_DEND_HEIGHT_MM, "mm"),
    row_dend_gp = grid::gpar(col = "black", lwd = DENDROGRAM_LINE_WIDTH),
    column_dend_gp = grid::gpar(col = "black", lwd = DENDROGRAM_LINE_WIDTH),
    width = heatmap_body_size,
    height = heatmap_body_size,
    show_heatmap_legend = FALSE
  )

  annotation_legends <- lapply(names(annotation_color_list), function(annotation_name) {
    annotation_palette <- annotation_color_list[[annotation_name]]

    ComplexHeatmap::Legend(
      title = annotation_name,
      at = names(annotation_palette),
      type = "grid",
      legend_gp = grid::gpar(
        fill = annotation_palette,
        col = "black",
        lwd = CELL_BORDER_WIDTH
      ),
      labels_gp = grid::gpar(
        fontsize = LEGEND_FONT_SIZE,
        fontface = TEXT_FONT_FACE,
        fontfamily = TEXT_FONT_FAMILY,
        col = TEXT_COLOR
      ),
      title_gp = grid::gpar(
        fontsize = LEGEND_FONT_SIZE,
        fontface = TEXT_FONT_FACE,
        fontfamily = TEXT_FONT_FAMILY,
        col = TEXT_COLOR
      ),
      grid_height = grid::unit(5, "mm"),
      grid_width = grid::unit(5, "mm"),
      gap = grid::unit(LEGEND_ITEM_GAP_MM, "mm"),
      row_gap = grid::unit(LEGEND_ITEM_GAP_MM, "mm"),
      title_gap = grid::unit(LEGEND_ITEM_GAP_MM, "mm")
    )
  })

  correlation_legend <- ComplexHeatmap::Legend(
    title = "Correlation",
    col_fun = heatmap_colors,
    at = round(seq(color_min, color_max, length.out = 5), 2),
    labels_gp = grid::gpar(
      fontsize = LEGEND_FONT_SIZE,
      fontface = TEXT_FONT_FACE,
      fontfamily = TEXT_FONT_FAMILY,
      col = TEXT_COLOR
    ),
    title_gp = grid::gpar(
      fontsize = LEGEND_FONT_SIZE,
      fontface = TEXT_FONT_FACE,
      fontfamily = TEXT_FONT_FAMILY,
      col = TEXT_COLOR
    ),
    legend_height = grid::unit(CORRELATION_LEGEND_HEIGHT_MM, "mm"),
    grid_width = grid::unit(6, "mm"),
    title_gap = grid::unit(LEGEND_ITEM_GAP_MM, "mm")
  )

  legend_pack <- do.call(
    ComplexHeatmap::packLegend,
    c(
      annotation_legends,
      list(
        correlation_legend,
        direction = "vertical",
        gap = grid::unit(LEGEND_GROUP_GAP_MM, "mm")
      )
    )
  )

  draw_heatmap_output <- function() {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(
      layout = grid::grid.layout(
        nrow = 3,
        ncol = 5,
        widths = grid::unit.c(
          grid::unit(OUTER_MARGIN_INCH, "in"),
          grid::unit(LEGEND_LEFT_WIDTH_INCH, "in"),
          grid::unit(LEGEND_HEATMAP_GAP_INCH, "in"),
          grid::unit(heatmap_panel_width, "in"),
          grid::unit(OUTER_MARGIN_INCH, "in")
        ),
        heights = grid::unit.c(
          grid::unit(OUTER_MARGIN_INCH, "in"),
          grid::unit(heatmap_panel_height, "in"),
          grid::unit(OUTER_MARGIN_INCH, "in")
        )
      )
    ))

    grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
    ComplexHeatmap::draw(
      legend_pack,
      x = grid::unit(LEGEND_INNER_MARGIN_INCH, "in"),
      y = grid::unit(1, "npc") - grid::unit(LEGEND_TOP_MARGIN_INCH, "in"),
      just = c("left", "top")
    )
    grid::popViewport()

    grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 4))
    ComplexHeatmap::draw(
      ht,
      newpage = FALSE,
      show_heatmap_legend = FALSE,
      show_annotation_legend = FALSE
    )
    grid::popViewport(2)
  }

  output_files <- save_grid_pdf_png(
    pdf_file = pdf_file,
    width = pdf_width,
    height = pdf_height,
    draw_fun = draw_heatmap_output
  )

  data.frame(
    Sample_Group = sample_group_name,
    Group_Value = ifelse(
      identical(sample_group_value, ALL_SAMPLE_GROUP_VALUE),
      "ALL_SAMPLES",
      sample_group_value
    ),
    Samples = nrow(sample_info),
    Genes_After_Biotype_Filter = gene_filter$selected_gene_count,
    Genes_After_Removing_Zero_Variance = variable_gene_result$nonzero_gene_count,
    Genes_Used_For_Clustering = nrow(expr_for_correlation),
    Top_Variable_Gene_N = variable_gene_result$top_variable_gene_n,
    Correlation_Method = CORRELATION_METHOD,
    Clustering_Method = CLUSTERING_METHOD,
    Column_Name_Rotation = COLUMN_NAMES_ROT,
    PDF_Width = round(pdf_width, 2),
    PDF_Height = round(pdf_height, 2),
    PDF_File = output_files$pdf_file,
    PNG_File = output_files$png_file,
    stringsAsFactors = FALSE
  )
}


# 5. 绘制全样本聚类热图 -------------------------------------------------------

cat("\nRunning sample TPM clustering heatmap generation...\n")
parallel_strategy <- setup_parallel_strategy(
  total_tasks = length(SAMPLE_GROUPS),
  inner_label = "Heatmap inner workers",
  nested_label = "Nested workers"
)

sample_group_names <- names(SAMPLE_GROUPS)
summary_list <- run_indexed_tasks_with_progress(
  total_tasks = length(sample_group_names),
  workers = parallel_strategy$task_workers,
  task_function = function(i) {
    sample_group_name <- sample_group_names[i]

    draw_sample_clustering_heatmap(
      sample_group_name = sample_group_name,
      sample_group_value = SAMPLE_GROUPS[[sample_group_name]]
    )
  }
)
stop_on_parallel_errors(summary_list, task_ids = sample_group_names, label = "sample heatmap tasks")

summary_table <- do.call(rbind, summary_list)
rownames(summary_table) <- NULL

summary_file <- file.path(PLOT_DIR, "sample_clustering_heatmap_summary.csv")
write.csv(summary_table, summary_file, row.names = FALSE, quote = TRUE, na = "")


# 6. 输出信息 ------------------------------------------------------------------

cat("\nSample TPM clustering heatmaps finished.\n")
cat("Gene biotype filter: ", gene_filter$filter, "\n", sep = "")
cat("TPM assay: ", TPM_ASSAY_NAME, "\n", sep = "")
cat("\nHeatmap summary:\n")
print(
  summary_table[
    ,
    c(
      "Sample_Group", "Samples",
      "Genes_After_Biotype_Filter", "Genes_After_Removing_Zero_Variance",
      "Genes_Used_For_Clustering", "Correlation_Method", "Clustering_Method",
      "PDF_Width", "PDF_Height"
    )
  ],
  row.names = FALSE
)

cat("\nSummary file:\n")
cat(summary_file, "\n", sep = "")

cat("\nPDF files:\n")
cat(paste(summary_table$PDF_File, collapse = "\n"), "\n", sep = "")
cat("\nPNG files:\n")
cat(paste(summary_table$PNG_File, collapse = "\n"), "\n", sep = "")
print_runtime_summary(SCRIPT_START_TIME, label = "Total runtime")
