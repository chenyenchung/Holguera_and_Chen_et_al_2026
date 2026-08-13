#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
if (is.na(script_arg)) {
  stop("Cannot determine the path of this script")
}
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

renv::load(repo_root)
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(R.utils))

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

cam_file <- if (is.null(argvs$cam)) {
  file.path(repo_root, "data", "P15_CAM.csv")
} else {
  argvs$cam
}
ann_file <- if (is.null(argvs$ann)) {
  file.path(repo_root, "data", "visual_neurons_anno.csv")
} else {
  argvs$ann
}
output_file <- if (is.null(argvs$output)) {
  file.path(repo_root, "int", "P15_CSM_neuronal_types.xlsx")
} else {
  argvs$output
}

required_annotation_columns <- c(
  "cell_type", "ozel2021_cluster", "Confident_annotation", "Notch",
  "ntype", "temporal_label", "broad_temp", "newly_ann", "subsystem",
  "func"
)

validate_file <- function(path, label) {
  if (length(path) != 1 || is.na(path) || !nzchar(path) || !file.exists(path)) {
    stop(sprintf("%s file not found: %s", label, path))
  }
}

parse_binary <- function(values, column_name) {
  if (is.logical(values)) {
    if (anyNA(values)) {
      stop(sprintf("CAM column '%s' contains missing values", column_name))
    }
    return(values)
  }

  normalized <- toupper(trimws(as.character(values)))
  valid <- normalized %in% c("TRUE", "FALSE", "T", "F", "1", "0")
  if (any(!valid)) {
    invalid <- unique(normalized[!valid])
    stop(sprintf(
      "CAM column '%s' contains non-binary values: %s",
      column_name,
      paste(head(invalid, 5), collapse = ", ")
    ))
  }
  normalized %in% c("TRUE", "T", "1")
}

normalize_cluster <- function(values) {
  values <- trimws(as.character(values))
  values[values %in% c("", "NA", "NaN")] <- NA_character_
  values <- sub("\\.0+$", "", values)
  values
}

write_sheet <- function(wb, sheet_name, values, widths = "auto") {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, values, keepNA = TRUE, na.string = "")
  freezePane(wb, sheet_name, firstRow = TRUE, firstCol = TRUE)
  if (ncol(values) > 0) {
    setColWidths(wb, sheet_name, cols = seq_len(ncol(values)), widths = widths)
    if (nrow(values) > 0) {
      addFilter(wb, sheet_name, rows = 1, cols = seq_len(ncol(values)))
    }
  }
}

validate_file(cam_file, "CAM expression")
validate_file(ann_file, "Annotation")

# The source CSV intentionally has a blank first header for its row-name
# column, so fread's automatic header detection would otherwise treat the
# header as a data row.
cam <- fread(cam_file, header = TRUE, check.names = FALSE)
if (ncol(cam) < 2 || nrow(cam) == 0) {
  stop("CAM expression matrix must contain at least one CSM and one cluster column")
}

setnames(cam, 1, "CSM")
cam[, CSM := trimws(as.character(CSM))]
if (anyNA(cam$CSM) || any(!nzchar(cam$CSM))) {
  stop("CAM expression matrix contains missing or empty CSM names")
}
if (anyDuplicated(cam$CSM)) {
  duplicates <- unique(cam$CSM[duplicated(cam$CSM)])
  stop(sprintf("CAM expression matrix contains duplicate CSM names: %s",
               paste(head(duplicates, 5), collapse = ", ")))
}

cluster_columns <- names(cam)[-1]
normalized_cluster_columns <- normalize_cluster(cluster_columns)
if (anyNA(normalized_cluster_columns) || any(!grepl("^[0-9]+$", normalized_cluster_columns))) {
  stop("All CAM expression columns after the first column must be numeric Ozel cluster IDs")
}
if (anyDuplicated(normalized_cluster_columns)) {
  stop("CAM expression matrix contains duplicate Ozel cluster columns")
}
setnames(cam, cluster_columns, normalized_cluster_columns)
cluster_columns <- normalized_cluster_columns

for (column_name in cluster_columns) {
  set(cam, j = column_name, value = parse_binary(cam[[column_name]], column_name))
}

ann <- fread(ann_file, header = TRUE, check.names = FALSE)
missing_annotation_columns <- setdiff(required_annotation_columns, names(ann))
if (length(missing_annotation_columns) > 0) {
  stop(sprintf(
    "Annotation file is missing required columns: %s",
    paste(missing_annotation_columns, collapse = ", ")
  ))
}
ann[, cell_type := trimws(as.character(cell_type))]
if (anyNA(ann$cell_type) || any(!nzchar(ann$cell_type))) {
  stop("Annotation file contains missing or empty cell_type values")
}
if (anyDuplicated(ann$cell_type)) {
  duplicates <- unique(ann$cell_type[duplicated(ann$cell_type)])
  stop(sprintf("Annotation file contains duplicate cell types: %s",
               paste(head(duplicates, 5), collapse = ", ")))
}
ann[, ozel2021_cluster := normalize_cluster(ozel2021_cluster)]
ann[, .annotation_order := .I]

mapping_audit <- ann[
  Confident_annotation == "Y" & !is.na(ozel2021_cluster),
  c(
    list(
      cell_type = cell_type,
      ozel2021_cluster = ozel2021_cluster,
      annotation_order = .annotation_order
    ),
    mget(setdiff(required_annotation_columns, c(
      "cell_type", "ozel2021_cluster"
    )))
  )
]
mapping_audit[, in_CAM_matrix := ozel2021_cluster %in% cluster_columns]
mapping_audit[, types_per_cluster := .N, by = ozel2021_cluster]
mapping_audit[, shared_cluster := types_per_cluster > 1]
setorder(mapping_audit, annotation_order)

type_mapping <- mapping_audit[in_CAM_matrix == TRUE]
if (nrow(type_mapping) == 0) {
  stop("No confidently annotated neuronal types map to CAM expression columns")
}
if (anyDuplicated(type_mapping$cell_type)) {
  stop("A confidently annotated neuronal type maps to more than one Ozel cluster")
}

cam_matrix <- as.matrix(cam[, ..cluster_columns])
storage.mode(cam_matrix) <- "logical"
rownames(cam_matrix) <- cam$CSM
positive_index <- which(cam_matrix, arr.ind = TRUE)

positive_clusters <- data.table(
  CSM = rownames(cam_matrix)[positive_index[, "row"]],
  ozel2021_cluster = colnames(cam_matrix)[positive_index[, "col"]],
  CSM_order = positive_index[, "row"],
  cluster_order = positive_index[, "col"]
)

long_columns <- c(
  "cell_type", "ozel2021_cluster", "Notch", "ntype", "temporal_label",
  "broad_temp", "newly_ann", "subsystem", "func", "Confident_annotation",
  "annotation_order"
)
expression_long <- merge(
  positive_clusters,
  type_mapping[, ..long_columns],
  by = "ozel2021_cluster",
  allow.cartesian = TRUE,
  sort = FALSE
)
setorder(expression_long, CSM_order, annotation_order)
expression_long[, c("CSM_order", "cluster_order", "annotation_order") := NULL]
setcolorder(expression_long, c(
  "CSM", "cell_type", "ozel2021_cluster", "Notch", "ntype",
  "temporal_label", "broad_temp", "newly_ann", "subsystem", "func",
  "Confident_annotation"
))

expression_summary <- expression_long[, .(
  n_expressing_neuronal_types = uniqueN(cell_type),
  expressing_neuronal_types = paste(unique(cell_type), collapse = "; ")
), by = CSM]
csm_summary <- data.table(CSM = cam$CSM, CSM_order = seq_len(nrow(cam)))
csm_summary <- merge(csm_summary, expression_summary, by = "CSM", all.x = TRUE, sort = FALSE)
setorder(csm_summary, CSM_order)
csm_summary[is.na(n_expressing_neuronal_types), n_expressing_neuronal_types := 0L]
csm_summary[is.na(expressing_neuronal_types), expressing_neuronal_types := ""]
csm_summary[, CSM_order := NULL]

expression_matrix <- data.table(CSM = cam$CSM)
for (mapping_row in seq_len(nrow(type_mapping))) {
  cell_type <- type_mapping$cell_type[mapping_row]
  cluster_id <- type_mapping$ozel2021_cluster[mapping_row]
  expression_matrix[, (cell_type) := ifelse(cam[[cluster_id]], "On", "Off")]
}

type_metadata_columns <- c(
  "cell_type", "ozel2021_cluster", "Notch", "ntype", "temporal_label",
  "broad_temp", "newly_ann", "subsystem", "func", "Confident_annotation",
  "types_per_cluster", "shared_cluster", "annotation_order"
)
type_metadata <- type_mapping[, ..type_metadata_columns]
setorder(type_metadata, annotation_order)
type_metadata[, annotation_order := NULL]

mapping_audit_output <- copy(mapping_audit)
mapping_audit_output[, annotation_order := NULL]

git_commit <- tryCatch(
  system2(
    "git",
    c("-C", repo_root, "rev-parse", "HEAD"),
    stdout = TRUE,
    stderr = FALSE
  )[1],
  error = function(e) "unavailable"
)
if (length(git_commit) == 0 || is.na(git_commit) || !nzchar(git_commit)) {
  git_commit <- "unavailable"
}

metadata <- data.table(
  field = c(
    "workbook_schema_version", "expression_stage", "CAM_input",
    "annotation_input", "output", "confidence_filter",
    "shared_cluster_rule", "matrix_values", "n_CSMs", "n_CAM_clusters",
    "n_confident_neuronal_types", "n_positive_CSM_type_pairs",
    "n_shared_confident_clusters", "generated_at_UTC", "script", "git_commit"
  ),
  value = as.character(c(
    "1.0.0", "P15", normalizePath(cam_file), normalizePath(ann_file),
    normalizePath(output_file, mustWork = FALSE),
    "Confident_annotation == Y and non-missing ozel2021_cluster",
    "Each positive cluster call expands to every confidently annotated neuronal type mapped to that cluster",
    "On/Off", nrow(cam), length(cluster_columns), nrow(type_mapping),
    nrow(expression_long),
    uniqueN(type_mapping[shared_cluster == TRUE, ozel2021_cluster]),
    format(Sys.time(), tz = "UTC", usetz = TRUE),
    file.path("src", basename(script_path)), git_commit
  ))
)

output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(output_dir)) {
  stop(sprintf("Could not create output directory: %s", output_dir))
}

wb <- createWorkbook(creator = "flyem CSM exporter")
write_sheet(wb, "CSM_summary", csm_summary, widths = c(20, 18, 80))
write_sheet(wb, "Expression_long", expression_long)
write_sheet(wb, "Expression_matrix", expression_matrix, widths = 15)
write_sheet(wb, "Type_metadata", type_metadata)
write_sheet(wb, "Mapping_audit", mapping_audit_output)
write_sheet(wb, "Metadata", metadata, widths = c(28, 100))
saveWorkbook(wb, output_file, overwrite = TRUE)

cat(sprintf("Wrote %s\n", normalizePath(output_file)))
cat(sprintf("CSMs: %d\n", nrow(cam)))
cat(sprintf("Confident neuronal types: %d\n", nrow(type_mapping)))
cat(sprintf("Positive CSM-type pairs: %d\n", nrow(expression_long)))
cat(sprintf(
  "Shared confident clusters: %d\n",
  uniqueN(type_mapping[shared_cluster == TRUE, ozel2021_cluster])
))
