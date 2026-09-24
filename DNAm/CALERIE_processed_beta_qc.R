#!/usr/bin/env Rscript

# CALERIE_processed_beta_qc.R
#
# CALERIE-specific processed beta matrix QC script.
#
# Edit only the CONFIGURATION section below, then run:
#   Rscript CALERIE_processed_beta_qc.R
#
# This script is intentionally separate from CALERIE_raw_idat_qc.R because raw
# IDAT QC can take a long time. Run this script after you have a processed beta
# matrix with sample columns named by Barcode.

# ==============================================================================
# CONFIGURATION - EDIT THESE PATHS
# ==============================================================================

# Full CALERIE sample sheet CSV. This matches CALERIE_raw_idat_qc.R.
sample_sheet_file <- Sys.getenv("CALERIE_SAMPLE_SHEET", "Preprocessing_cell_counts/Data/idats_samplesheet/Calerie_SampleSheet_DL_FINAL.csv")

# Processed beta matrix. Strongly recommend .rds for large CALERIE data.
# Expected shape: CpGs/features as rows, samples as columns.
beta_file <- Sys.getenv("CALERIE_BETA_FILE", "EDIT_ME/path/to/CALERIE_processed_betas.rds")

# Output folder. This follows the same output root used by CALERIE_raw_idat_qc.R.
output_dir <- Sys.getenv("CALERIE_PROCESSED_QC_OUTPUT_DIR", "Code/FAST_DNAmQC/CALERIE_processed_beta_qc_output")

# CALERIE replicate definition.
# Barcode is the unique sample key and should match beta matrix column names.
# Replicates are samples with the same person column AND same timepoint column.
person_col <- Sys.getenv("CALERIE_PERSON_COL", "Participant_ID")
timepoint_col <- Sys.getenv("CALERIE_TIMEPOINT_COL", "Time_Point")

# Optional explicit replicate pair file. Leave as NULL for CALERIE default
# grouping by person_col + timepoint_col. If used, it must contain sample_1 and
# sample_2 columns with Barcode values.
pair_file <- Sys.getenv("CALERIE_PAIR_FILE", "")
if (!nzchar(pair_file)) pair_file <- NULL

# QC settings.
correlation_method <- "pearson"
uniformity_p_threshold <- 0.05

# Which steps to run.
RUN_REPLICATE_CORRELATIONS <- TRUE
RUN_UNIFORMITY_CHECK <- TRUE

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

require_package("diptest", "Install with: install.packages(\"diptest\")")

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

  if (all(c("Sentrix_ID", "Sentrix_Position") %in% names(samples))) {
    samples$Barcode <- paste(samples$Sentrix_ID, samples$Sentrix_Position, sep = "_")
    return(samples)
  }

  stop(
    "Could not create Barcode. The sample sheet must contain either Barcode, ",
    "barcode, Slide + Array, or Sentrix_ID + Sentrix_Position.",
    call. = FALSE
  )
}

read_beta_matrix <- function(path) {
  stop_if_unedited(path, "beta_file")
  if (!file.exists(path)) {
    stop("beta_file does not exist: ", path, call. = FALSE)
  }

  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    beta_mat <- readRDS(path)
  } else {
    warning("Reading a large beta matrix from CSV can be slow and memory-intensive. Prefer .rds when possible.")
    beta_mat <- utils::read.csv(path, row.names = 1, check.names = FALSE)
  }

  beta_mat <- as.matrix(beta_mat)
  storage.mode(beta_mat) <- "double"

  if (is.null(colnames(beta_mat)) || any(!nzchar(colnames(beta_mat)))) {
    stop("Beta matrix must have sample IDs as column names. These should match Barcode.", call. = FALSE)
  }

  beta_mat
}

validate_beta_sample_match <- function(samples, beta_mat) {
  samples_in_beta <- samples$Barcode %in% colnames(beta_mat)
  beta_in_samples <- colnames(beta_mat) %in% samples$Barcode

  match_summary <- data.frame(
    metric = c(
      "n_sample_sheet_rows",
      "n_beta_matrix_columns",
      "n_sample_sheet_barcodes_in_beta",
      "n_beta_columns_in_sample_sheet"
    ),
    value = c(
      nrow(samples),
      ncol(beta_mat),
      sum(samples_in_beta),
      sum(beta_in_samples)
    ),
    stringsAsFactors = FALSE
  )

  list(
    match_summary = match_summary,
    missing_from_beta = samples[!samples_in_beta, , drop = FALSE],
    missing_from_sample_sheet = data.frame(
      Barcode = colnames(beta_mat)[!beta_in_samples],
      stringsAsFactors = FALSE
    )
  )
}

calculate_replicate_correlations <- function(beta_mat,
                                             samples,
                                             person_col,
                                             timepoint_col,
                                             pair_file = NULL,
                                             method = "pearson") {
  if (!is.null(pair_file)) {
    if (!file.exists(pair_file)) {
      stop("pair_file does not exist: ", pair_file, call. = FALSE)
    }
    pairs <- utils::read.csv(pair_file, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
    required_cols <- c("sample_1", "sample_2")
    missing_cols <- setdiff(required_cols, names(pairs))
    if (length(missing_cols)) {
      stop("pair_file is missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
    }
    if (!"replicate_group" %in% names(pairs)) {
      pairs$replicate_group <- NA_character_
    }
  } else {
    missing_cols <- setdiff(c("Barcode", person_col, timepoint_col), names(samples))
    if (length(missing_cols)) {
      stop("Sample sheet is missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
    }

    samples_for_reps <- samples[samples$Barcode %in% colnames(beta_mat), , drop = FALSE]
    group_key <- paste(samples_for_reps[[person_col]], samples_for_reps[[timepoint_col]], sep = "||")
    split_barcodes <- split(samples_for_reps$Barcode, group_key)
    split_barcodes <- split_barcodes[lengths(split_barcodes) > 1]

    if (!length(split_barcodes)) {
      pairs <- data.frame(
        replicate_group = character(),
        sample_1 = character(),
        sample_2 = character(),
        stringsAsFactors = FALSE
      )
    } else {
      pairs <- do.call(rbind, lapply(names(split_barcodes), function(group_name) {
        pair_matrix <- utils::combn(split_barcodes[[group_name]], 2)
        data.frame(
          replicate_group = group_name,
          sample_1 = pair_matrix[1, ],
          sample_2 = pair_matrix[2, ],
          stringsAsFactors = FALSE
        )
      }))
    }
  }

  pairs <- pairs[pairs$sample_1 %in% colnames(beta_mat) & pairs$sample_2 %in% colnames(beta_mat), , drop = FALSE]

  if (!nrow(pairs)) {
    results <- data.frame(
      replicate_group = character(),
      sample_1 = character(),
      sample_2 = character(),
      correlation = numeric(),
      stringsAsFactors = FALSE
    )
  } else {
    pairs$correlation <- vapply(seq_len(nrow(pairs)), function(i) {
      stats::cor(
        beta_mat[, pairs$sample_1[[i]]],
        beta_mat[, pairs$sample_2[[i]]],
        use = "pairwise.complete.obs",
        method = method
      )
    }, numeric(1))
    results <- pairs[, c("replicate_group", "sample_1", "sample_2", "correlation"), drop = FALSE]
  }

  summary <- data.frame(
    metric = c("n_pairs", "mean_correlation", "median_correlation", "min_correlation"),
    value = c(
      nrow(results),
      if (nrow(results)) mean(results$correlation, na.rm = TRUE) else NA_real_,
      if (nrow(results)) stats::median(results$correlation, na.rm = TRUE) else NA_real_,
      if (nrow(results)) min(results$correlation, na.rm = TRUE) else NA_real_
    ),
    stringsAsFactors = FALSE
  )

  list(results = results, summary = summary)
}

calculate_uniformity_check <- function(beta_mat, p_threshold = 0.05) {
  dip_p_values <- apply(beta_mat, 2, function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 2) {
      return(NA_real_)
    }

    p_value <- NA_real_
    suppressMessages(utils::capture.output({
      p_value <- diptest::dip.test(x)$p.value
    }))
    p_value
  })

  results <- data.frame(
    Barcode = names(dip_p_values),
    dip_p_value = as.numeric(dip_p_values),
    p_threshold = p_threshold,
    non_unimodal = as.numeric(dip_p_values) < p_threshold,
    stringsAsFactors = FALSE
  )

  summary <- data.frame(
    metric = c("n_samples", "p_threshold", "n_non_unimodal", "percent_non_unimodal"),
    value = c(
      nrow(results),
      p_threshold,
      sum(results$non_unimodal, na.rm = TRUE),
      mean(results$non_unimodal, na.rm = TRUE) * 100
    ),
    stringsAsFactors = FALSE
  )

  list(results = results, summary = summary)
}

# ==============================================================================
# RUN QC
# ==============================================================================

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "01_beta_sample_matching"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "02_replicate_correlations"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "03_uniformity_check"), recursive = TRUE, showWarnings = FALSE)

cat("=== CALERIE Processed Beta QC ===\n")
cat("Sample sheet:", sample_sheet_file, "\n")
cat("Beta file:", beta_file, "\n")
cat("Output dir:", output_dir, "\n")
cat("Person column:", person_col, "\n")
cat("Timepoint column:", timepoint_col, "\n\n")

samples <- read_calerie_sample_sheet(sample_sheet_file)
samples <- add_calerie_barcode(samples)

if (anyDuplicated(samples$Barcode)) {
  duplicated_barcodes <- unique(samples$Barcode[duplicated(samples$Barcode)])
  stop("Duplicate Barcode values found: ", paste(duplicated_barcodes, collapse = ", "), call. = FALSE)
}

beta_mat <- read_beta_matrix(beta_file)

cat("Beta matrix dimensions:", nrow(beta_mat), "features x", ncol(beta_mat), "samples\n")

match_info <- validate_beta_sample_match(samples, beta_mat)
utils::write.csv(
  match_info$match_summary,
  file.path(output_dir, "01_beta_sample_matching", "beta_sample_match_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  match_info$missing_from_beta,
  file.path(output_dir, "01_beta_sample_matching", "sample_sheet_barcodes_missing_from_beta.csv"),
  row.names = FALSE
)
utils::write.csv(
  match_info$missing_from_sample_sheet,
  file.path(output_dir, "01_beta_sample_matching", "beta_columns_missing_from_sample_sheet.csv"),
  row.names = FALSE
)

cat("Sample sheet barcodes in beta:", match_info$match_summary$value[match_info$match_summary$metric == "n_sample_sheet_barcodes_in_beta"], "\n")
cat("Beta columns in sample sheet:", match_info$match_summary$value[match_info$match_summary$metric == "n_beta_columns_in_sample_sheet"], "\n")

processed_qc_summary <- match_info$match_summary

if (RUN_REPLICATE_CORRELATIONS) {
  cat("\n=== Replicate Correlations ===\n")
  replicate_qc <- calculate_replicate_correlations(
    beta_mat = beta_mat,
    samples = samples,
    person_col = person_col,
    timepoint_col = timepoint_col,
    pair_file = pair_file,
    method = correlation_method
  )

  utils::write.csv(
    replicate_qc$results,
    file.path(output_dir, "02_replicate_correlations", "replicate_correlations.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    replicate_qc$summary,
    file.path(output_dir, "02_replicate_correlations", "replicate_correlations_summary.csv"),
    row.names = FALSE
  )

  processed_qc_summary <- rbind(
    processed_qc_summary,
    data.frame(metric = paste0("replicate_", replicate_qc$summary$metric), value = replicate_qc$summary$value)
  )

  cat("Replicate pairs:", replicate_qc$summary$value[replicate_qc$summary$metric == "n_pairs"], "\n")
}

if (RUN_UNIFORMITY_CHECK) {
  cat("\n=== Uniformity Check ===\n")
  uniformity_qc <- calculate_uniformity_check(beta_mat, p_threshold = uniformity_p_threshold)

  utils::write.csv(
    uniformity_qc$results,
    file.path(output_dir, "03_uniformity_check", "uniformity_check.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    uniformity_qc$summary,
    file.path(output_dir, "03_uniformity_check", "uniformity_summary.csv"),
    row.names = FALSE
  )

  processed_qc_summary <- rbind(
    processed_qc_summary,
    data.frame(metric = paste0("uniformity_", uniformity_qc$summary$metric), value = uniformity_qc$summary$value)
  )

  cat("Non-unimodal samples:", uniformity_qc$summary$value[uniformity_qc$summary$metric == "n_non_unimodal"], "\n")
}

utils::write.csv(
  processed_qc_summary,
  file.path(output_dir, "processed_beta_qc_summary.csv"),
  row.names = FALSE
)

cat("\n=== CALERIE Processed Beta QC Complete ===\n")
cat("Output directory:", output_dir, "\n")
