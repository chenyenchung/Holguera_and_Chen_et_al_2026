#!/usr/bin/env Rscript
# Run from the repository root with the installed R 4.5.1 environment.
# This exports existing results; it does not rerun tests or adjust P values.
suppressPackageStartupMessages(library(openxlsx))

outdir <- "int/20260922_deliver"
overwrite <- "--overwrite" %in% commandArgs(trailingOnly = TRUE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
depth_file <- "int/stats/deep_superficial/combined_broad_depth_results.xlsx"
enrichment_file <- "int/stats/functional_enrichment/functional_enrichment_pvalues_corrected.xlsx"
spatial_file <- "int/ds_plot/spatial_origin/spatial_origin_neuronal_types.xlsx"
template_file <- "int/Supp. Table 5-TS and CAM.xlsx"
manifest <- list()
column_manifest <- list()
outputs <- list()
source_files <- character()
neuropil_labels <- c(ME_R = "Medulla", LO_R = "Lobula", LOP_R = "Lobula plate")

read_sheet <- function(path, sheet) {
  source_files <<- union(source_files, path)
  read.xlsx(path, sheet = sheet, check.names = FALSE, sep.names = " ", skipEmptyCols = FALSE)
}
read_csv <- function(path) {
  source_files <<- union(source_files, path)
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
}
select_columns <- function(d, map) {
  stopifnot(all(unname(map) %in% names(d)))
  result <- d[, unname(map), drop = FALSE]
  names(result) <- names(map)
  result
}
add_result <- function(sheets, name, d, map, source, source_sheet = "", selection = "All source rows") {
  stopifnot(nchar(name) <= 31, !name %in% names(sheets))
  if ("neuropil" %in% names(d)) {
    d <- d[!is.na(d$neuropil) & d$neuropil %in% names(neuropil_labels), , drop = FALSE]
    d$neuropil <- unname(neuropil_labels[d$neuropil])
    selection <- paste(selection, "; right hemisphere only (ME_R, LO_R, LOP_R); labels = Medulla, Lobula, Lobula plate")
  }
  sheets[[name]] <- select_columns(d, map)
  attr(sheets[[name]], "provenance") <- list(source = source, source_sheet = source_sheet, selection = selection, map = map)
  sheets
}

# Names and order match the supplied Table 5 example; identity/context varies.
depth_map <- c(
  "Observed Bias Ratio" = "observed_bias_ratio",
  "p-value" = "p_value_exceeds_threshold",
  "q-value (FDR)" = "p_value_fdr",
  "Synapse Type" = "syn_type",
  "Observed Distance to Reference (Superficial)" = "observed_sup_d",
  "Observed Distance to Reference (Deep)" = "observed_deep_d",
  "Observed Distance Difference" = "observed_distance_diff",
  "Observed Reference Difference" = "observed_delta_thres_base",
  "Bootstrapped Distance Difference (Median)" = "bootstrap_distance_diff_median",
  "Bootstrapped Distance Difference (2.5%)" = "bootstrap_distance_diff_lower",
  "Bootstrapped Distance Difference (97.5%)" = "bootstrap_distance_diff_upper",
  "Bootstrapped Reference Difference (Median)" = "bootstrap_delta_thres_base_median",
  "Bootstrapped Reference Difference (2.5%)" = "bootstrap_delta_thres_base_lower",
  "Bootstrapped Reference Difference (97.5%)" = "bootstrap_delta_thres_base_upper",
  "Bootstrapped Bias Ratio (Median)" = "bootstrap_bias_ratio_median",
  "Bootstrapped Bias Ratio (2.5%)" = "bootstrap_bias_ratio_lower",
  "Bootstrapped Bias Ratio (97.5%)" = "bootstrap_bias_ratio_upper"
)
add_depth <- function(sheets, name, source_sheet, identity) {
  d <- read_sheet(depth_file, source_sheet)
  stopifnot(all(d$conf_level == 95), all(d$n_bootstrap == 1000), all(d$coefficient == 0.5))
  d$notch_label <- ifelse(d$split == "all", "All", sub("_.*$", "", d$split))
  d$class_label <- ifelse(d$split == "all", "All", sub("^[^_]*_", "", d$split))
  stopifnot(all(d$notch_label %in% c("All", "Notch On", "Notch Off")))
  map <- c(setNames("types_of_interest", identity), "Notch Status" = "notch_label",
           "Neuronal Class" = "class_label", "Neuropil" = "neuropil", depth_map)
  add_result(sheets, name, d, map, depth_file, source_sheet)
}

common_notes <- c(
  "Data are re-exported from existing analysis results. No statistical tests were rerun and no P or Q values were recomputed.",
  "Depth columns follow the supplied Table 5 example. p-value is p_value_exceeds_threshold; q-value (FDR) is the corresponding existing p_value_fdr. The alternative p_value_bias_ratio is not exported.",
  "Depth bias is (distance to deep reference minus distance to superficial reference) / reference difference. Positive values indicate superficial bias; negative values indicate deep bias.",
  "Observed Reference Difference is observed_delta_thres_base, before multiplication by the 0.5 depth-test coefficient. Bootstrap bounds are the 2.5th and 97.5th percentiles (95% interval; 1,000 resamples in these depth tables).",
  "Distance columns retain the original rotated-coordinate units. No unit conversion or numerical rounding is applied to stored results.",
  "Zero bootstrap P values are preserved as reported by the analysis; they mean no non-exceeding resamples were counted, not an exact population probability of zero.",
  "Neuropil-specific results use only the right hemisphere: ME_R is labeled Medulla, LO_R Lobula, and LOP_R Lobula plate. pre/post denote presynaptic/postsynaptic compartments. intrinsic is the source annotation for interneurons.",
  "Known/new/combined and putative analysis groups remain separate. Left-hemisphere rows are excluded. Tables of cell-type annotations or enrichment without a neuropil field retain their original biological scope. A context value of All means the source test pooled that dimension.",
  "Blank cells mean unavailable or inapplicable values, not zero. Number formatting affects display only. Documentation and exact source-column mappings accompany the workbook."
)

save_table <- function(number, filename, sheets, notes) {
  directory <- file.path(outdir, paste0("Supp_Table_", number))
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(directory, filename)
  if (file.exists(path) && !overwrite) stop("Refusing to overwrite ", path, "; use --overwrite for regenerated delivery files")
  wb <- createWorkbook(creator = "FlyEM supplementary table export")
  header <- createStyle(fontName = "Calibri", fontSize = 11, textDecoration = "bold", fgFill = "#E2EFDA", wrapText = TRUE, valign = "center")
  numeric_style <- createStyle(numFmt = "0.######")
  probability_style <- createStyle(numFmt = "0.000E+00")
  percentage_style <- createStyle(numFmt = "0.0")
  for (sheet in names(sheets)) {
    d <- sheets[[sheet]]
    provenance <- attr(d, "provenance")
    attr(d, "provenance") <- NULL
    addWorksheet(wb, sheet)
    writeData(wb, sheet, d, headerStyle = header, withFilter = nrow(d) > 0, keepNA = FALSE)
    # Excel has no numeric infinity. Keep finite entries numeric and write
    # only infinite estimates as explicit text rather than Excel errors.
    for (col in which(vapply(d, is.numeric, logical(1)))) {
      for (row in which(is.infinite(d[[col]]))) {
        writeData(wb, sheet, as.character(d[[col]][row]), startRow = row + 1,
                  startCol = col, colNames = FALSE)
      }
    }
    freezePane(wb, sheet, firstRow = TRUE)
    setRowHeights(wb, sheet, 1, 66)
    widths <- ifelse(nchar(names(d)) > 32, 23, ifelse(nchar(names(d)) > 18, 21, 18))
    widths[grepl("Explanation|Test Used|Alternative|Source", names(d))] <- 38
    setColWidths(wb, sheet, seq_len(ncol(d)), widths)
    if (nrow(d)) {
      numeric_cols <- which(vapply(d, is.numeric, logical(1)))
      if (length(numeric_cols)) addStyle(wb, sheet, numeric_style, rows = 2:(nrow(d)+1), cols = numeric_cols, gridExpand = TRUE, stack = TRUE)
      pcols <- which(names(d) %in% c("p-value", "q-value (FDR)"))
      if (length(pcols)) addStyle(wb, sheet, probability_style, rows = 2:(nrow(d)+1), cols = pcols, gridExpand = TRUE, stack = TRUE)
      pct <- which(names(d) == "Percentage")
      if (length(pct)) addStyle(wb, sheet, percentage_style, rows = 2:(nrow(d)+1), cols = pct, gridExpand = TRUE, stack = TRUE)
    }
    manifest[[length(manifest)+1]] <<- data.frame(Table = number, Workbook = filename, Sheet = sheet, Rows = nrow(d), Source = provenance$source, Source_sheet = provenance$source_sheet, Selection = provenance$selection)
    column_manifest[[length(column_manifest)+1]] <<- data.frame(Table = number, Sheet = sheet, Export_column = names(provenance$map), Source_column = unname(provenance$map))
    sheets[[sheet]] <- d
  }
  saveWorkbook(wb, path, overwrite = overwrite)
  stopifnot(identical(getSheetNames(path), names(sheets)))
  for (sheet in names(sheets)) {
    actual <- read.xlsx(path, sheet = sheet, check.names = FALSE, sep.names = " ", skipEmptyCols = FALSE)
    expected <- sheets[[sheet]]
    stopifnot(identical(names(actual), names(expected)), nrow(actual) == nrow(expected))
    for (col in names(expected)) {
      a <- actual[[col]]; e <- expected[[col]]
      stopifnot(identical(is.na(a), is.na(e)))
      ok <- !is.na(e)
      if (is.numeric(e)) {
        stopifnot(is.numeric(a) || !any(ok) || any(is.infinite(e)))
        stopifnot(isTRUE(all.equal(as.numeric(a[ok]), as.numeric(e[ok]), tolerance = 1e-13)))
      } else stopifnot(identical(as.character(a[ok]), as.character(e[ok])))
    }
  }
  writeLines(c(paste0("# Supplementary Table ", number), "", paste0("Workbook: ", filename), "", paste0("- ", c(common_notes, notes))), file.path(directory, "README.md"))
  outputs[[as.character(number)]] <<- list(path = path, sheets = sheets)
  cat("Verified", path, ":", length(sheets), "sheets,", sum(vapply(sheets, nrow, integer(1))), "data rows\n")
}

sheets <- list()
for (pair in list(c("Temporal Known", "Temporal Known"), c("Temporal New", "temporal_new"), c("Temporal Combined", "temporal_all"), c("Early-Late Known", "broad_known"), c("Early-Late New", "broad_new"))) {
  sheets <- add_depth(sheets, pair[1], pair[2], "Temporal Origin")
}
depth_fdr_note <- "Aggregate depth Q values preserve Benjamini-Hochberg correction within each source preset, neuropil, synapse type, and biological split. They are not re-adjusted across exported sheets."
save_table(1, "Supp_Table_1_Temporal_origin.xlsx", sheets, depth_fdr_note)

sheets <- list()
for (i in 1:9) sheets <- add_depth(sheets, paste("Putative OPC", i), paste0("type_putative_", i), "Cell Type")
save_table(3, "Supp_Table_3_Putative_OPC.xlsx", sheets, c(depth_fdr_note, "The nine putative-OPC preset groups are retained separately to preserve the source analysis and correction families."))

sheets <- list()
for (pair in list(c("Subsystem Known", "Subsystem Known"), c("Subsystem New", "subsystem_new"), c("Subsystem Putative", "subsystem_putative"))) sheets <- add_depth(sheets, pair[1], pair[2], "Functional Subsystem")
enrichment_map <- c("Functional Subsystem" = "Functional_Subsystem", "Notch Status" = "Notch_Status", "Total Cell Types" = "N_Total", "Subsystem Cell Types" = "N_Subsystem", "Test Used" = "Test_Used", "Statistic" = "Statistic", "Direction" = "Direction", "Null Expectation" = "Null_Expectation", "Alternative" = "Alternative", "p-value" = "P_value_raw", "q-value (FDR)" = "P_value_FDR", "Untestable Explanation" = "Skip_Reason")
enrichment_names <- c(Temporal = "Temporal Fisher", Cochran_Armitage = "Cochran-Armitage", Wald_Wolfowitz = "Wald-Wolfowitz", Broad_Temporal = "Early-Late by Temporal ID", Notch = "Notch Fisher", Broad_Temp = "Early-Late by Annotation")
for (original in names(enrichment_names)) {
  d <- read_sheet(enrichment_file, original)
  sheets <- add_result(sheets, enrichment_names[[original]], d, enrichment_map, enrichment_file, original)
}
per_window_file <- "int/Sup_tbl_3_per_window_func_enrichment.csv"
sheets <- add_result(sheets, "Per-window Fisher", read_csv(per_window_file), c("Temporal Origin" = "temporal", "Functional Subsystem" = "subsystem", "Odds Ratio" = "or", "p-value" = "p.val", "q-value (FDR)" = "adj.p"), per_window_file)
save_table(4, "Supp_Table_4_Functional_subsystems.xlsx", sheets, c(depth_fdr_note,
  "Temporal Fisher, Cochran-Armitage, and Wald-Wolfowitz each retain a separate FDR family pooling subsystem hypotheses across Notch On and Off. Other Fisher tests retain the pooled Existing_Fisher family; per-window Fisher tests retain their own BH correction.",
  "Cochran-Armitage uses the exact conditional test ordered by temporal_id. Wald-Wolfowitz uses the exact tie-averaged test for contiguous concentration. The reported statistic, null expectation, alternative, and direction are retained where applicable.",
  "Early-Late by Temporal ID and Early-Late by Annotation preserve both original analyses (Broad_Temporal and Broad_Temp); the former uses temporal_id < el_cut, the latter broad_temp labels. Their current numerical results coincide, but neither source test is removed.",
  "Per-window Fisher uses confidently annotated cell types, comparing membership in each temporal window against subsystem membership. Its odds ratio is preserved; an infinite estimate is represented by the text Inf, not a missing value."))

# Editorial annotations come only from the author-supplied workbook.
template <- lapply(c("TS-Medulla", "TS-Lobula", "CAM-Medulla", "CAM-Lobula"), function(s) read_sheet(template_file, s))
names(template) <- c("TS-Medulla", "TS-Lobula", "CAM-Medulla", "CAM-Lobula")
concentric <- rbind(template[["TS-Medulla"]][, c("Selector", "Concentric gene")], template[["TS-Lobula"]][, c("Selector", "Concentric gene")])
concentric_values <- split(concentric[["Concentric gene"]], concentric$Selector)
stopifnot(all(vapply(concentric_values, function(x) length(unique(x[!is.na(x)])) == 1L, logical(1))))
concentric_lookup <- vapply(concentric_values, function(x) unique(x[!is.na(x)])[1], character(1))
row_key <- function(gene, notch, syn) paste(gene, notch, syn, sep = "|")
sheets <- list(); summary_rows <- list(); lop_audit <- character()
for (molecule in c("TS", "CAM")) {
  path <- if (molecule == "TS") "int/selector_test/combined_selector_depth.csv" else "int/cam_test/combined_cam_depth.csv"
  all <- read_csv(path)
  tested <- all[is.na(all$skip_reason) | trimws(all$skip_reason) == "", , drop = FALSE]
  stopifnot(all(is.finite(tested$p_value_exceeds_threshold)), all(is.finite(tested$p_value_fdr)), all(tested$conf_level == 95), all(tested$n_bootstrap == 1000))
  unique_key <- paste(tested$neuropil, row_key(tested$types_of_interest, tested$notch_category, tested$syn_type))
  stopifnot(!anyDuplicated(unique_key))
  lop_genes <- unique(tested$types_of_interest[tested$neuropil == "LOP_R"])
  other_genes <- unique(tested$types_of_interest[tested$neuropil %in% c("ME_R", "LO_R")])
  stopifnot(length(setdiff(lop_genes, other_genes)) == 0)
  lop_audit <- c(lop_audit, sprintf("%s: %d tested lobula-plate genes; none exclusive to lobula plate.", molecule, length(lop_genes)))
  for (np in c("ME_R", "LO_R", "LOP_R")) {
    label <- c(ME_R = "Medulla", LO_R = "Lobula", LOP_R = "Lobula plate")[[np]]
    d <- tested[tested$neuropil == np, , drop = FALSE]
    d$Visualized <- "No"
    if (np != "LOP_R") {
      t <- template[[paste(molecule, label, sep = "-")]]
      id <- if (molecule == "TS") "Selector" else "CAM"
      keys <- row_key(t[[id]], t[["Notch Status"]], t[["Synapse Type"]])
      stopifnot(!anyDuplicated(keys))
      matched <- match(row_key(d$types_of_interest, d$notch_category, d$syn_type), keys)
      stopifnot(!anyNA(matched), length(matched) == nrow(t))
      d$Visualized <- t$Visualized[matched]
    }
    map <- c(setNames("types_of_interest", if (molecule == "TS") "Selector" else "CAM"), "Notch Status" = "notch_category", depth_map, "Visualized" = "Visualized")
    if (molecule == "TS") {
      d$concentric_gene <- unname(concentric_lookup[d$types_of_interest])
      stopifnot(!anyNA(d$concentric_gene))
      map <- c(map, "Concentric gene" = "concentric_gene")
    }
    sheets <- add_result(sheets, paste(molecule, label, sep = "-"), d, map, path, selection = paste("Successful tests; neuropil =", np, "; editorial annotations from supplied Table 5"))
    for (notch in c("Notch On", "Notch Off")) for (syn in c("pre", "post")) {
      sub <- d[d$notch_category == notch & d$syn_type == syn, , drop = FALSE]
      n <- nrow(sub); significant <- sum(sub$p_value_fdr < 0.05)
      summary_rows[[length(summary_rows)+1]] <- data.frame(molecule_class = if (molecule == "TS") "Selector TF" else "CAM", neuropil = np, notch_status = notch, synapse_type = syn, significant_genes = significant, tested_genes = n, percentage = 100 * significant / n)
    }
  }
}
summary <- do.call(rbind, summary_rows)
summary_map <- c("Molecule Class" = "molecule_class", "Neuropil" = "neuropil", "Notch Status" = "notch_status", "Synapse Type" = "synapse_type", "Significant Genes" = "significant_genes", "Tested Genes" = "tested_genes", "Percentage" = "percentage")
sheets <- add_result(sheets, "Depth association percentages", summary, summary_map, "int/selector_test/combined_selector_depth.csv; int/cam_test/combined_cam_depth.csv", selection = "Derived counts from successful tests; significance = q < 0.05; percentage = 100 * significant / tested")
sheets <- sheets[c("Depth association percentages", setdiff(names(sheets), "Depth association percentages"))]
save_table(5, "Supp_Table_5_TS_and_CAM.xlsx", sheets, c(
  "All six molecular depth sheets concern the right hemisphere. Only successful tests are exported. The complete source CSVs retain skipped analyses for audit.",
  "Existing selector/CAM Q values use BH correction within neuropil x synapse type x Notch status after combining all batches. Percentages are recomputed from successful tests using Q < 0.05 and displayed to one decimal place.",
  "Visualized values in medulla and lobula are copied by molecule, neuropil, Notch status, and synapse type from the supplied Table 5. Every lobula-plate Visualized value is No, as instructed by the author.",
  "Concentric gene annotations apply to selectors and are transferred by gene from the supplied medulla/lobula sheets. Every exported selector has an unambiguous annotation; no values are inferred from the statistics.", lop_audit))

sheets <- list()
lookup_map <- c("Cell Type" = "cell_type", "Temporal Origin" = "temporal_origin", "Spatial Origin" = "spatial_origin", "Notch Status" = "notch_status", "Neuronal Class" = "neuron_class", "Newly Annotated" = "newly_annotated", "Ozel 2021 Cluster" = "ozel2021_cluster", "Ventral Vsx" = "vVsx", "Dorsal Vsx" = "dVsx", "Ventral Optix" = "vOptix", "Dorsal Optix" = "dOptix", "Ventral Dpp" = "vDpp", "Dorsal Dpp" = "dDpp")
lookup <- read_sheet(spatial_file, "neuronal_types")
stopifnot(nrow(lookup) == 55, all(lookup$confident_annotation == "Y"))
sheets <- add_result(sheets, "Spatial origin by cell type", lookup, lookup_map, spatial_file, "neuronal_types")
d <- read_sheet(spatial_file, "figure_data")
d <- d[!is.na(d$included_in_figure) & d$included_in_figure == TRUE, , drop = FALSE]
figure_map <- c("Cell Type" = "cell_type", "Temporal Origin" = "temporal_origin", "Spatial Origin" = "spatial_origin", "Notch Status" = "notch_status", "Neuronal Class" = "neuron_class", "Neuropil" = "neuropil", "Synapse Type" = "synapse_type", "Bootstrapped Bias Ratio (Median)" = "bootstrap_bias_ratio_median", "Bootstrapped Bias Ratio (2.5%)" = "bootstrap_bias_ratio_lower", "Bootstrapped Bias Ratio (97.5%)" = "bootstrap_bias_ratio_upper", "Number of Neurons" = "n_neurons_interest", "Number of Synapses" = "n_syn_interest")
sheets <- add_result(sheets, "Within-temporal figure data", d, figure_map, spatial_file, "figure_data", "included_in_figure = TRUE; descriptive data, without significance columns")
for (pair in list(c("Supporting Spatial All", "spatial_all"), c("Supporting Spatial Notch", "spatial_notch"), c("Supporting Spatial Hth", "spatial_Hth"))) sheets <- add_depth(sheets, pair[1], pair[2], "Spatial Origin")
heatmap_file <- "int/origin_heatmap/origin_heatmap_counts.csv"
sheets <- add_result(sheets, "Spatial origin heatmap counts", read_csv(heatmap_file), c("Neuropil" = "neuropil", "Synapse Type" = "syn_type", "Temporal Origin" = "temporal_origin", "Spatial Origin" = "spatial_origin", "Synapse Count" = "count", "Log10 Synapse Count" = "count_log10"), heatmap_file)
save_table(6, "Supp_Table_6_Spatial_patterning.xlsx", sheets, c(depth_fdr_note,
  "Spatial origin by cell type includes 55 confidently annotated types and their biological annotations; technical sorting and constant confidence fields are omitted. Ventral/dorsal spatial-domain calls and source cluster identifiers are retained as biological data.",
  "Within-temporal figure data contains only rows actually included in revised Supplementary Figure 18. This is descriptive: P values, Q values, and significance flags are deliberately omitted. It is stratified by temporal origin, Notch status, neuronal class, and neuropil.",
  "Supporting Spatial All, Supporting Spatial Notch, and Supporting Spatial Hth retain aggregate statistical analyses. The former aggregate spatial figure was replaced by the within-temporal figure; these supporting analyses are not the inferential basis of that replacement.",
  "Hth results are pooled across Notch status and neuronal class, restricted to the Hth temporal window and right medulla.",
  "Heatmap counts are synapse counts grouped by temporal and spatial origin, not neuronal-type counts. Zero counts have blank log10 values because log10(0) is undefined. The accompanying origin_heatmap.pdf is copied unchanged from the existing analysis."))
pdf_source <- "int/origin_heatmap/origin_heatmap.pdf"
source_files <- union(source_files, pdf_source)
stopifnot(file.copy(pdf_source, file.path(outdir, "Supp_Table_6", "origin_heatmap.pdf"), overwrite = overwrite))

write.csv(do.call(rbind, manifest), file.path(outdir, "SOURCE_MANIFEST.csv"), row.names = FALSE, na = "")
write.csv(do.call(rbind, column_manifest), file.path(outdir, "COLUMN_MAPPING.csv"), row.names = FALSE, na = "")
writeLines(sort(source_files), file.path(outdir, "SOURCE_FILES.txt"))
writeLines(c(
  "# Supplementary-table delivery — 2026-09-22", "",
  "Five manuscript-ready Excel workbooks are organized in Supp_Table_1, Supp_Table_3, Supp_Table_4, Supp_Table_5, and Supp_Table_6. Each workbook contains only its relevant result tabs. Table 6 also contains the existing spatial-origin heatmap PDF.", "",
  "The supplied int/Supp. Table 5-TS and CAM.xlsx determines the depth-statistic labels, order, and editorial annotations. Existing analysis outputs determine numerical results. Originals are unchanged.", "",
  "Each subdirectory README defines its statistical context. SOURCE_MANIFEST.csv maps every delivered tab to the input file, source sheet, row count, and selection. COLUMN_MAPPING.csv records each exported column's source field. notch_label and class_label are derived from the source split field; Visualized and concentric_gene come from the author-supplied workbook (lobula-plate Visualized = No). Percentage-summary fields are derived counts as described in the Table 5 README.", "",
  "Depth P values consistently use p_value_exceeds_threshold to match the existing FDR Q values. Statistics and stored numerical precision are preserved. No new tests or FDR corrections were run. Table 5 percentages use exact counts, with display rounded to one decimal place.", "",
  "All neuropil-specific tables are restricted to ME_R, LO_R, and LOP_R and display these as Medulla, Lobula, and Lobula plate. This matches the right-hemisphere selector/CAM and spatial data in the August 15/21 deliveries; bilateral spatial plots were separate supporting material. Enrichment and cell-type annotation tables without a neuropil field retain their original scope. Known/new/combined analyses remain distinct.", "",
  "Table 6 distinguishes the revised descriptive figure data from supporting aggregate spatial tests. The two early/late Fisher variants in Table 4 remain separate and are labeled by their source grouping definition.", "",
  "Validation: all output workbooks were reopened and checked column-by-column against the selected source data at numerical tolerance 1e-13, including missing-value masks. Editorial joins must be unique and complete. Source hashes and delivery hashes accompany this package.", "",
  paste0("- ", lop_audit), "",
  "Reproduction: module load r/4.5.1; Rscript --vanilla src/export_supplementary_tables_20260922.R. Add --overwrite only to regenerate this delivery. Then run python src/verify_supplementary_tables_20260922.py to independently validate the exports and generate checksums. No input files are modified."
), file.path(outdir, "README.md"))
