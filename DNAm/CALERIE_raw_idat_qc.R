#!/usr/bin/env Rscript

# CALERIE_raw_idat_qc.R
#
# CALERIE-specific raw IDAT QC script.
#
# Edit only the CONFIGURATION section below, then run:
#   Rscript CALERIE_raw_idat_qc.R
#
# This script is intentionally chunked. It does not load all 600+ samples into
# one RGset at once. Each chunk is loaded as an RGset, QC metrics are extracted,
# then the chunk RGset is removed from memory.

# ==============================================================================
# CONFIGURATION - EDIT THESE PATHS
# ==============================================================================

# Full CALERIE sample sheet CSV.
sample_sheet_file <- Sys.getenv("CALERIE_SAMPLE_SHEET", "Preprocessing_cell_counts/Data/idats_samplesheet/Calerie_SampleSheet_DL_FINAL.csv")

# Folder containing IDAT files named like Slide_Array_Red.idat and Slide_Array_Grn.idat.
idat_dir <- Sys.getenv("CALERIE_IDAT_DIR", "Preprocessing_cell_counts/Data/idats_samplesheet/IDAT_FILES")

# Output folder. The script will create this if needed.
output_dir <- Sys.getenv("CALERIE_OUTPUT_DIR", "Code/FAST_DNAmQC/CALERIE_raw_idat_qc_output")

# Number of samples loaded at one time. If memory errors occur, lower to 24 or 48.
chunk_size <- as.integer(Sys.getenv("CALERIE_CHUNK_SIZE", "608"))

# QC thresholds requested for collaborator-facing QC.
bead_mean_threshold <- 2
detection_p_threshold <- 0.05
control_pca_variance_threshold <- 0.90

# Save each chunk RGset? Usually FALSE. Set TRUE only if you explicitly need RDS objects.
save_chunk_rgsets <- FALSE

# ==============================================================================
# PACKAGE CHECKS
# ==============================================================================

require_package <- function(package, install_hint = NULL) {
  if (!requireNamespace(package, quietly = TRUE)) {
    msg <- paste0("Required package is not installed: ", package)
    if (!is.null(install_hint)) {
      msg <- paste0(msg, "\n", install_hint)
    }
    stop(msg, call. = FALSE)
  }
}

require_package("minfi", "Install with: BiocManager::install(\"minfi\")")
require_package("Biobase", "Install with: BiocManager::install(\"Biobase\")")

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

stop_if_unedited <- function(path, label) {
  if (is.null(path) || !nzchar(path) || grepl("^EDIT_ME", path)) {
    stop("Edit ", label, " in the CONFIGURATION section before running this script.", call. = FALSE)
  }
}

read_calerie_sample_sheet <- function(path) {
  stop_if_unedited(path, "sample_sheet_file")
  if (!file.exists(path)) {
    stop("sample_sheet_file does not exist: ", path, call. = FALSE)
  }

  lines <- readLines(path, warn = FALSE)
  data_line <- grep("^\\[Data\\]", trimws(lines))

  if (length(data_line)) {
    csv_text <- paste(lines[(data_line[[1]] + 1):length(lines)], collapse = "\n")
    samples <- utils::read.csv(
      text = csv_text,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      colClasses = "character",
      na.strings = c("", "NA")
    )
  } else {
    samples <- utils::read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      colClasses = "character",
      na.strings = c("", "NA")
    )
  }

  samples <- samples[rowSums(!is.na(samples)) > 0, , drop = FALSE]
  rownames(samples) <- NULL
  samples
}

add_calerie_barcode <- function(samples) {
  if ("Barcode" %in% names(samples)) {
    samples$Barcode <- as.character(samples$Barcode)
    return(samples)
  }

  if ("barcode" %in% names(samples)) {
    samples$Barcode <- as.character(samples$barcode)
    return(samples)
  }

  if (all(c("Slide", "Array") %in% names(samples))) {
    samples$Barcode <- paste(samples$Slide, samples$Array, sep = "_")
    return(samples)
  }

  # Fallback for Illumina-style sample sheets if needed.
  if (all(c("Sentrix_ID", "Sentrix_Position") %in% names(samples))) {
    samples$Barcode <- paste(samples$Sentrix_ID, samples$Sentrix_Position, sep = "_")
    return(samples)
  }

  stop(
    "Could not create Barcode. The sample sheet must contain either Barcode, ",
    "Slide + Array, or Sentrix_ID + Sentrix_Position.",
    call. = FALSE
  )
}

validate_idat_files <- function(samples, idat_dir) {
  stop_if_unedited(idat_dir, "idat_dir")
  if (!dir.exists(idat_dir)) {
    stop("idat_dir does not exist: ", idat_dir, call. = FALSE)
  }

  samples$Basename <- file.path(idat_dir, samples$Barcode)
  samples$Red_IDAT <- paste0(samples$Basename, "_Red.idat")
  samples$Green_IDAT <- paste0(samples$Basename, "_Grn.idat")
  samples$red_idat_exists <- file.exists(samples$Red_IDAT)
  samples$green_idat_exists <- file.exists(samples$Green_IDAT)
  samples$idat_pair_exists <- samples$red_idat_exists & samples$green_idat_exists
  samples
}

ensure_manifest_package <- function(rgset) {
  array_name <- minfi::annotation(rgset)[["array"]]
  if (is.null(array_name) || is.na(array_name) || !nzchar(array_name)) {
    return(invisible(TRUE))
  }

  manifest_package <- paste0(array_name, "manifest")
  if (!requireNamespace(manifest_package, quietly = TRUE)) {
    stop(
      "Required array manifest package is not installed: ", manifest_package, "\n",
      "Install with: BiocManager::install(\"", manifest_package, "\")",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

extract_nbeads_matrix <- function(rgset) {
  if ("getNBeads" %in% getNamespaceExports("minfi")) {
    nbeads <- try(minfi::getNBeads(rgset), silent = TRUE)
    if (!inherits(nbeads, "try-error")) {
      return(nbeads)
    }
  }

  assay_elements <- Biobase::assayDataElementNames(rgset)
  if ("NBeads" %in% assay_elements) {
    return(Biobase::assayDataElement(rgset, "NBeads"))
  }

  if (requireNamespace("wateRmelon", quietly = TRUE)) {
    return(wateRmelon::beadcount(rgset))
  }

  stop(
    "Could not extract bead counts. Load IDATs with extended = TRUE, or install wateRmelon.",
    call. = FALSE
  )
}

load_rgset_chunk <- function(samples_chunk) {
  rgset <- minfi::read.metharray(
    basenames = samples_chunk$Basename,
    extended = TRUE,
    verbose = TRUE
  )
  colnames(rgset) <- samples_chunk$Barcode
  rgset
}

calculate_control_beta_chunk <- function(rgset) {
  ensure_manifest_package(rgset)
  control_probes <- minfi::getControlAddress(rgset)
  control_probes <- intersect(control_probes, rownames(minfi::getRed(rgset)))

  red_intensities <- minfi::getRed(rgset)[control_probes, , drop = FALSE]
  green_intensities <- minfi::getGreen(rgset)[control_probes, , drop = FALSE]
  control_beta <- red_intensities / (red_intensities + green_intensities + 100)
  control_beta
}

calculate_bead_qc_chunk <- function(rgset) {
  nbeads <- extract_nbeads_matrix(rgset)
  mean_bead_count <- colMeans(nbeads, na.rm = TRUE)

  data.frame(
    Barcode = names(mean_bead_count),
    mean_bead_count = as.numeric(mean_bead_count),
    bead_mean_threshold = bead_mean_threshold,
    passed_mean_bead_count = as.numeric(mean_bead_count) > bead_mean_threshold,
    stringsAsFactors = FALSE
  )
}

calculate_detection_qc_chunk <- function(rgset) {
  ensure_manifest_package(rgset)
  det_p <- minfi::detectionP(rgset)
  mean_detection_p_value <- colMeans(det_p, na.rm = TRUE)

  data.frame(
    Barcode = names(mean_detection_p_value),
    mean_detection_p_value = as.numeric(mean_detection_p_value),
    detection_p_threshold = detection_p_threshold,
    passed_mean_detection_p_value = as.numeric(mean_detection_p_value) < detection_p_threshold,
    stringsAsFactors = FALSE
  )
}

bind_control_beta <- function(existing, new) {
  if (is.null(existing)) {
    return(new)
  }

  common_probes <- intersect(rownames(existing), rownames(new))
  if (length(common_probes) < 2) {
    stop("Control probe rows do not overlap across chunks.", call. = FALSE)
  }

  cbind(existing[common_probes, , drop = FALSE], new[common_probes, , drop = FALSE])
}

run_control_pca <- function(control_beta) {
  keep_probe <- apply(control_beta, 1, function(x) all(is.finite(x)) && stats::sd(x) > 0)
  control_beta <- control_beta[keep_probe, , drop = FALSE]

  if (nrow(control_beta) < 2) {
    stop("Fewer than two usable control probes remained after filtering.", call. = FALSE)
  }

  pca <- stats::prcomp(t(control_beta), center = TRUE, scale. = TRUE)
  variance_explained <- pca$sdev^2 / sum(pca$sdev^2)
  cumulative_variance_explained <- cumsum(variance_explained)
  n_pcs_to_threshold <- which(cumulative_variance_explained >= control_pca_variance_threshold)[[1]]

  scores <- as.data.frame(pca$x[, seq_len(n_pcs_to_threshold), drop = FALSE])
  scores <- data.frame(Barcode = rownames(scores), scores, row.names = NULL, check.names = FALSE)

  variance <- data.frame(
    pc = paste0("PC", seq_along(variance_explained)),
    variance_explained = variance_explained,
    cumulative_variance_explained = cumulative_variance_explained,
    stringsAsFactors = FALSE
  )

  summary <- data.frame(
    metric = c("n_samples", "n_control_probes_used", "variance_threshold", "n_pcs_to_threshold"),
    value = c(ncol(control_beta), nrow(control_beta), control_pca_variance_threshold, n_pcs_to_threshold),
    stringsAsFactors = FALSE
  )

  list(scores = scores, variance = variance, summary = summary, pca = pca)
}

# ==============================================================================
# RUN QC
# ==============================================================================

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "01_manifest_rgset"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "02_control_probe_pca"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "03_bead_qc"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "04_detection_pval_qc"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "restricted_objects"), recursive = TRUE, showWarnings = FALSE)

cat("=== CALERIE Raw IDAT QC ===\n")
cat("Sample sheet:", sample_sheet_file, "\n")
cat("IDAT dir:", idat_dir, "\n")
cat("Output dir:", output_dir, "\n")
cat("Chunk size:", chunk_size, "\n\n")

samples <- read_calerie_sample_sheet(sample_sheet_file)
samples <- add_calerie_barcode(samples)
samples <- validate_idat_files(samples, idat_dir)

utils::write.csv(
  samples,
  file.path(output_dir, "01_manifest_rgset", "raw_sample_manifest_validated.csv"),
  row.names = FALSE
)

if (!all(samples$idat_pair_exists)) {
  missing_idats <- samples[!samples$idat_pair_exists, , drop = FALSE]
  utils::write.csv(
    missing_idats,
    file.path(output_dir, "01_manifest_rgset", "missing_idats.csv"),
    row.names = FALSE
  )
  stop("Missing IDAT files detected. See 01_manifest_rgset/missing_idats.csv", call. = FALSE)
}

if (anyDuplicated(samples$Barcode)) {
  duplicated_barcodes <- unique(samples$Barcode[duplicated(samples$Barcode)])
  stop("Duplicate Barcode values found: ", paste(duplicated_barcodes, collapse = ", "), call. = FALSE)
}

n_samples <- nrow(samples)
chunk_ids <- split(seq_len(n_samples), ceiling(seq_len(n_samples) / chunk_size))

control_beta_all <- NULL
bead_qc_all <- list()
detection_qc_all <- list()

for (i in seq_along(chunk_ids)) {
  cat("\n=== Chunk ", i, " of ", length(chunk_ids), " ===\n", sep = "")
  idx <- chunk_ids[[i]]
  samples_chunk <- samples[idx, , drop = FALSE]
  cat("Samples:", nrow(samples_chunk), "\n")

  rgset <- load_rgset_chunk(samples_chunk)

  if (save_chunk_rgsets) {
    saveRDS(rgset, file.path(output_dir, "restricted_objects", paste0("RGset_chunk_", i, ".rds")))
  }

  cat("Calculating control-probe values\n")
  control_beta_chunk <- calculate_control_beta_chunk(rgset)
  control_beta_all <- bind_control_beta(control_beta_all, control_beta_chunk)

  cat("Calculating bead QC\n")
  bead_qc_all[[i]] <- calculate_bead_qc_chunk(rgset)

  cat("Calculating detection p-value QC\n")
  detection_qc_all[[i]] <- calculate_detection_qc_chunk(rgset)

  rm(rgset, control_beta_chunk)
  gc(verbose = FALSE)
}

cat("\n=== Combining QC Results ===\n")

bead_qc <- do.call(rbind, bead_qc_all)
detection_qc <- do.call(rbind, detection_qc_all)

bead_qc <- merge(samples, bead_qc, by = "Barcode", all.x = TRUE, sort = FALSE)
detection_qc <- merge(samples, detection_qc, by = "Barcode", all.x = TRUE, sort = FALSE)

bead_summary <- data.frame(
  metric = c(
    "n_samples",
    "bead_mean_threshold",
    "n_samples_mean_bead_count_gt_threshold",
    "percent_samples_mean_bead_count_gt_threshold"
  ),
  value = c(
    nrow(bead_qc),
    bead_mean_threshold,
    sum(bead_qc$mean_bead_count > bead_mean_threshold, na.rm = TRUE),
    mean(bead_qc$mean_bead_count > bead_mean_threshold, na.rm = TRUE) * 100
  ),
  stringsAsFactors = FALSE
)

detection_summary <- data.frame(
  metric = c(
    "n_samples",
    "detection_p_threshold",
    "n_samples_mean_detection_p_lt_threshold",
    "percent_samples_mean_detection_p_lt_threshold"
  ),
  value = c(
    nrow(detection_qc),
    detection_p_threshold,
    sum(detection_qc$mean_detection_p_value < detection_p_threshold, na.rm = TRUE),
    mean(detection_qc$mean_detection_p_value < detection_p_threshold, na.rm = TRUE) * 100
  ),
  stringsAsFactors = FALSE
)

cat("Running control-probe PCA on combined control-probe matrix\n")
control_pca <- run_control_pca(control_beta_all)

utils::write.csv(
  control_pca$scores,
  file.path(output_dir, "02_control_probe_pca", "control_probe_pca_scores.csv"),
  row.names = FALSE
)
utils::write.csv(
  control_pca$variance,
  file.path(output_dir, "02_control_probe_pca", "control_probe_pca_variance_explained.csv"),
  row.names = FALSE
)
utils::write.csv(
  control_pca$summary,
  file.path(output_dir, "02_control_probe_pca", "control_probe_pca_summary.csv"),
  row.names = FALSE
)
saveRDS(control_pca$pca, file.path(output_dir, "restricted_objects", "control_probe_pca.rds"))

utils::write.csv(
  bead_qc,
  file.path(output_dir, "03_bead_qc", "bead_qc_per_sample.csv"),
  row.names = FALSE
)
utils::write.csv(
  bead_summary,
  file.path(output_dir, "03_bead_qc", "bead_qc_summary.csv"),
  row.names = FALSE
)

utils::write.csv(
  detection_qc,
  file.path(output_dir, "04_detection_pval_qc", "detection_pval_per_sample.csv"),
  row.names = FALSE
)
utils::write.csv(
  detection_summary,
  file.path(output_dir, "04_detection_pval_qc", "detection_pval_summary.csv"),
  row.names = FALSE
)

raw_qc_summary <- data.frame(
  metric = c(
    "n_samples_total",
    "chunk_size",
    "n_control_probes_used",
    "control_pca_variance_threshold",
    "n_control_pcs_to_threshold",
    "bead_mean_threshold",
    "n_samples_mean_bead_count_gt_threshold",
    "percent_samples_mean_bead_count_gt_threshold",
    "detection_p_threshold",
    "n_samples_mean_detection_p_lt_threshold",
    "percent_samples_mean_detection_p_lt_threshold"
  ),
  value = c(
    n_samples,
    chunk_size,
    control_pca$summary$value[control_pca$summary$metric == "n_control_probes_used"],
    control_pca_variance_threshold,
    control_pca$summary$value[control_pca$summary$metric == "n_pcs_to_threshold"],
    bead_mean_threshold,
    bead_summary$value[bead_summary$metric == "n_samples_mean_bead_count_gt_threshold"],
    bead_summary$value[bead_summary$metric == "percent_samples_mean_bead_count_gt_threshold"],
    detection_p_threshold,
    detection_summary$value[detection_summary$metric == "n_samples_mean_detection_p_lt_threshold"],
    detection_summary$value[detection_summary$metric == "percent_samples_mean_detection_p_lt_threshold"]
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  raw_qc_summary,
  file.path(output_dir, "raw_qc_summary.csv"),
  row.names = FALSE
)

cat("\n=== CALERIE Raw IDAT QC Complete ===\n")
cat("Samples processed:", n_samples, "\n")
cat("Control probes used:", control_pca$summary$value[control_pca$summary$metric == "n_control_probes_used"], "\n")
cat("PCs to 90% variance:", control_pca$summary$value[control_pca$summary$metric == "n_pcs_to_threshold"], "\n")
cat("Samples with mean bead count >", bead_mean_threshold, ":", bead_summary$value[bead_summary$metric == "n_samples_mean_bead_count_gt_threshold"], "\n")
cat("Samples with mean detection p-value <", detection_p_threshold, ":", detection_summary$value[detection_summary$metric == "n_samples_mean_detection_p_lt_threshold"], "\n")
cat("Output directory:", output_dir, "\n")

