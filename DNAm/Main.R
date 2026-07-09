#!/usr/bin/env Rscript

# Main.R - DNAm QC Pipeline Entry Point
# Run this script to execute the full DNAm QC pipeline.

# ==============================================================================
# SETUP
# ==============================================================================

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
file_idx <- grep(file_arg, script_args, fixed = TRUE)
pipeline_dir <- if (length(file_idx)) {
  dirname(normalizePath(sub(file_arg, "", script_args[file_idx[[1]]])))
} else {
  getwd()
}

source(file.path(pipeline_dir, "dnam_io_helpers.R"))
source(file.path(pipeline_dir, "dnam_load_rgset.R"))
source(file.path(pipeline_dir, "dnam_control_probe_pca.R"))
source(file.path(pipeline_dir, "dnam_bead_qc.R"))
source(file.path(pipeline_dir, "dnam_detection_pval_qc.R"))
source(file.path(pipeline_dir, "dnam_replicate_correlations.R"))
source(file.path(pipeline_dir, "dnam_uniformity_check.R"))

# ==============================================================================
# CONFIGURATION - EDIT THIS SECTION
# ==============================================================================

# Raw IDAT inputs
sample_sheet_file <- "path/to/sample_sheet.csv"
idat_dir <- "path/to/idats"

# Processed beta matrix input. Required only for replicate correlations and uniformity.
beta_file <- "path/to/processed_beta_matrix.rds"

# Optional explicit replicate-pair file with columns sample_1 and sample_2.
# If NULL, replicate groups are inferred from replicate_group_cols in the sample sheet.
pair_file <- NULL

# Output directory
output_dir <- "output"

# Column names
sample_id_col <- "barcode"
basename_col <- "barcode"
sentrix_id_col <- "Sentrix_ID"
sentrix_position_col <- "Sentrix_Position"
replicate_group_cols <- c("Participant_ID", "Time_Point")

# QC thresholds
bead_mean_threshold <- 2
detection_p_threshold <- 0.05
control_pca_variance_threshold <- 0.90
uniformity_p_threshold <- 0.05

# Which steps to run
RUN_RAW_IDAT_QC <- TRUE
RUN_PROCESSED_BETA_QC <- TRUE
RUN_REPLICATE_CORRELATIONS <- TRUE
RUN_UNIFORMITY_CHECK <- TRUE

# Privacy and restricted output settings
MASK_IDS <- FALSE
SAVE_RGSET <- FALSE
SAVE_CONTROL_PCA_OBJECT <- TRUE

# ==============================================================================
# CREATE OUTPUT DIRECTORIES
# ==============================================================================

dir.create(file.path(output_dir, "01_manifest_rgset"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "02_control_probe_pca"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "03_bead_qc"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "04_detection_pval_qc"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "05_replicate_correlations"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "06_uniformity_check"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "restricted_objects"), recursive = TRUE, showWarnings = FALSE)

mask_suffix <- ifelse(MASK_IDS, "_masked", "")

# ==============================================================================
# STEP 1: LOAD SAMPLE SHEET, VALIDATE IDATS, AND CREATE RGSET
# ==============================================================================

if (RUN_RAW_IDAT_QC) {
  cat("=== Step 1: Sample Sheet, IDAT Validation, and RGset Loading ===\n")

  sample_sheet <- dnam_read_sample_sheet(sample_sheet_file)
  sample_sheet <- dnam_add_idat_basenames(
    sample_sheet = sample_sheet,
    idat_dir = idat_dir,
    basename_col = basename_col,
    sentrix_id_col = sentrix_id_col,
    sentrix_position_col = sentrix_position_col
  )
  sample_sheet <- dnam_validate_idat_pairs(sample_sheet)

  manifest_to_save <- if (MASK_IDS) dnam_mask_columns(sample_sheet, c(sample_id_col)) else sample_sheet
  utils::write.csv(
    manifest_to_save,
    file.path(output_dir, "01_manifest_rgset", paste0("raw_sample_manifest_validated", mask_suffix, ".csv")),
    row.names = FALSE
  )

  if (!all(sample_sheet$idat_pair_exists)) {
    missing_idats <- sample_sheet[!sample_sheet$idat_pair_exists, , drop = FALSE]
    missing_to_save <- if (MASK_IDS) dnam_mask_columns(missing_idats, c(sample_id_col)) else missing_idats
    utils::write.csv(
      missing_to_save,
      file.path(output_dir, "01_manifest_rgset", paste0("missing_idats", mask_suffix, ".csv")),
      row.names = FALSE
    )
    stop("One or more samples are missing red/green IDAT files. See missing_idats.csv.", call. = FALSE)
  }

  rgset <- dnam_load_rgset(sample_sheet, sample_id_col = sample_id_col, extended = TRUE)

  rgset_summary <- data.frame(
    metric = c("n_samples", "n_idat_pairs", "sample_id_col", "basename_col"),
    value = c(ncol(rgset), sum(sample_sheet$idat_pair_exists), sample_id_col, basename_col),
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    rgset_summary,
    file.path(output_dir, "01_manifest_rgset", "rgset_summary.csv"),
    row.names = FALSE
  )

  if (SAVE_RGSET) {
    saveRDS(rgset, file.path(output_dir, "restricted_objects", "RGset.rds"))
  }

  # ==============================================================================
  # STEP 2: CONTROL-PROBE PCA
  # ==============================================================================

  cat("\n=== Step 2: Control-Probe PCA ===\n")

  control_pca <- dnam_control_probe_pca(
    rgset,
    variance_threshold = control_pca_variance_threshold,
    pseudocount = 100
  )

  pca_scores <- if (MASK_IDS) dnam_mask_columns(control_pca$scores, "sample_id") else control_pca$scores
  utils::write.csv(
    pca_scores,
    file.path(output_dir, "02_control_probe_pca", paste0("control_probe_pca_scores", mask_suffix, ".csv")),
    row.names = FALSE
  )
  utils::write.csv(
    control_pca$variance_explained,
    file.path(output_dir, "02_control_probe_pca", "control_probe_pca_variance_explained.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    control_pca$summary,
    file.path(output_dir, "02_control_probe_pca", "control_probe_pca_summary.csv"),
    row.names = FALSE
  )

  if (SAVE_CONTROL_PCA_OBJECT) {
    saveRDS(control_pca$pca_object, file.path(output_dir, "restricted_objects", "control_probe_pca.rds"))
  }

  cat("PCs to threshold:", control_pca$summary$value[control_pca$summary$metric == "n_pcs_to_threshold"], "\n")

  # ==============================================================================
  # STEP 3: BEAD QC
  # ==============================================================================

  cat("\n=== Step 3: Bead QC ===\n")

  bead_qc <- dnam_bead_qc(rgset, mean_bead_threshold = bead_mean_threshold)
  bead_per_sample <- if (MASK_IDS) dnam_mask_columns(bead_qc$per_sample, "sample_id") else bead_qc$per_sample
  utils::write.csv(
    bead_per_sample,
    file.path(output_dir, "03_bead_qc", paste0("bead_qc_per_sample", mask_suffix, ".csv")),
    row.names = FALSE
  )
  utils::write.csv(
    bead_qc$summary,
    file.path(output_dir, "03_bead_qc", "bead_qc_summary.csv"),
    row.names = FALSE
  )
  cat("Samples passing bead threshold:", bead_qc$summary$value[bead_qc$summary$metric == "n_samples_mean_bead_count_gt_threshold"], "\n")

  # ==============================================================================
  # STEP 4: DETECTION P-VALUE QC
  # ==============================================================================

  cat("\n=== Step 4: Detection P-Value QC ===\n")

  detection_qc <- dnam_detection_pval_qc(rgset, detection_p_threshold = detection_p_threshold)
  detection_per_sample <- if (MASK_IDS) dnam_mask_columns(detection_qc$per_sample, "sample_id") else detection_qc$per_sample
  utils::write.csv(
    detection_per_sample,
    file.path(output_dir, "04_detection_pval_qc", paste0("detection_pval_per_sample", mask_suffix, ".csv")),
    row.names = FALSE
  )
  utils::write.csv(
    detection_qc$summary,
    file.path(output_dir, "04_detection_pval_qc", "detection_pval_summary.csv"),
    row.names = FALSE
  )
  cat("Samples passing detection threshold:", detection_qc$summary$value[detection_qc$summary$metric == "n_samples_mean_detection_p_lt_threshold"], "\n")
}

# ==============================================================================
# STEP 5: PROCESSED BETA MATRIX QC
# ==============================================================================

if (RUN_PROCESSED_BETA_QC) {
  cat("\n=== Step 5: Processed Beta Matrix Loading ===\n")

  beta_mat <- dnam_read_beta_matrix(beta_file)
  cat("Beta matrix dimensions:", nrow(beta_mat), "features x", ncol(beta_mat), "samples\n")

  pheno <- NULL
  if (file.exists(sample_sheet_file)) {
    pheno <- dnam_read_sample_sheet(sample_sheet_file)
  }

  if (RUN_REPLICATE_CORRELATIONS) {
    cat("\n=== Step 6: Replicate Correlations ===\n")

    pair_data <- NULL
    if (!is.null(pair_file)) {
      dnam_check_file(pair_file, "pair_file")
      pair_data <- utils::read.csv(pair_file, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
    }

    replicate_qc <- dnam_replicate_correlations(
      beta_mat = beta_mat,
      pheno = pheno,
      sample_id_col = sample_id_col,
      replicate_group_cols = replicate_group_cols,
      pair_data = pair_data
    )

    replicate_results <- if (MASK_IDS) dnam_mask_columns(replicate_qc$results, c("sample_1", "sample_2")) else replicate_qc$results
    utils::write.csv(
      replicate_results,
      file.path(output_dir, "05_replicate_correlations", paste0("replicate_correlations", mask_suffix, ".csv")),
      row.names = FALSE
    )
    utils::write.csv(
      replicate_qc$summary,
      file.path(output_dir, "05_replicate_correlations", "replicate_correlations_summary.csv"),
      row.names = FALSE
    )
    cat("Replicate pairs:", replicate_qc$summary$value[replicate_qc$summary$metric == "n_pairs"], "\n")
  }

  if (RUN_UNIFORMITY_CHECK) {
    cat("\n=== Step 7: Uniformity Check ===\n")

    uniformity_qc <- dnam_uniformity_check(beta_mat, p_threshold = uniformity_p_threshold)
    uniformity_results <- if (MASK_IDS) dnam_mask_columns(uniformity_qc$results, "sample_id") else uniformity_qc$results
    utils::write.csv(
      uniformity_results,
      file.path(output_dir, "06_uniformity_check", paste0("uniformity_check", mask_suffix, ".csv")),
      row.names = FALSE
    )
    utils::write.csv(
      uniformity_qc$summary,
      file.path(output_dir, "06_uniformity_check", "uniformity_summary.csv"),
      row.names = FALSE
    )
    cat("Non-unimodal samples:", uniformity_qc$summary$value[uniformity_qc$summary$metric == "n_non_unimodal"], "\n")
  }
}

# ==============================================================================
# COMPLETE
# ==============================================================================

cat("\n=== Pipeline Complete ===\n")
cat("Output directory:", output_dir, "\n")
cat("Results saved to:\n")
cat("  - 01_manifest_rgset/\n")
cat("  - 02_control_probe_pca/\n")
cat("  - 03_bead_qc/\n")
cat("  - 04_detection_pval_qc/\n")
cat("  - 05_replicate_correlations/\n")
cat("  - 06_uniformity_check/\n")
cat("  - restricted_objects/\n")
