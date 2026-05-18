olink_sample_qc <- function(df,
                             sample_col = "SAMPLE_ID",
                             value_col = "NPX",
                             assay_type_col = "AssayType",
                             assay_keep = "assay",
                             cutoff = 3,
                             max_prune = 0,
                             output_dir = "output",
                             mask_sample_ids = FALSE) {
  
  if (!is.data.frame(df)) {
    stop("df must be a data.frame")
  }
  
  if (!sample_col %in% names(df)) {
    stop(paste0("Column '", sample_col, "' not found in df"))
  }
  
  if (!value_col %in% names(df)) {
    stop(paste0("Column '", value_col, "' not found in df"))
  }
  
  if (!assay_type_col %in% names(df)) {
    stop(paste0("Column '", assay_type_col, "' not found in df"))
  }
  
  if (!is.numeric(cutoff) || length(cutoff) != 1 || is.na(cutoff)) {
    stop("cutoff must be a single numeric value")
  }
  
  if (!is.numeric(max_prune) || max_prune < 0 || max_prune != floor(max_prune)) {
    stop("max_prune must be a non-negative integer")
  }
  
  compute_stats <- function(dat) {
    split_vals <- split(dat[[value_col]], dat[[sample_col]])
    sample_ids <- names(split_vals)
    sample_median <- vapply(split_vals, median, numeric(1), na.rm = TRUE)
    sample_iqr <- vapply(split_vals, IQR, numeric(1), na.rm = TRUE)
    
    median_mean <- mean(sample_median, na.rm = TRUE)
    median_sd <- sd(sample_median, na.rm = TRUE)
    iqr_mean <- mean(sample_iqr, na.rm = TRUE)
    iqr_sd <- sd(sample_iqr, na.rm = TRUE)
    
    median_low <- median_mean - cutoff * median_sd
    median_high <- median_mean + cutoff * median_sd
    iqr_low <- iqr_mean - cutoff * iqr_sd
    iqr_high <- iqr_mean + cutoff * iqr_sd
    
    data.frame(
      sample_id = sample_ids,
      sample_median = sample_median,
      sample_iqr = sample_iqr,
      global_median = median_mean,
      global_iqr = iqr_mean,
      median_z = (sample_median - median_mean) / median_sd,
      iqr_z = (sample_iqr - iqr_mean) / iqr_sd,
      median_fail = sample_median < median_low | sample_median > median_high,
      iqr_fail = sample_iqr < iqr_low | sample_iqr > iqr_high,
      stringsAsFactors = FALSE
    )
  }
  
  dat <- df[
    !is.na(df[[assay_type_col]]) &
      df[[assay_type_col]] == assay_keep &
      !is.na(df[[sample_col]]) &
      !is.na(df[[value_col]]),
    c(sample_col, value_col)
  ]
  
  if (nrow(dat) == 0) {
    stop("No non-missing assay data available after filtering")
  }
  
  removed_samples <- character(0)
  iter <- 0
  
  sample_stats <- compute_stats(dat)
  sample_stats$any_fail <- sample_stats$median_fail | sample_stats$iqr_fail
  
  total_samples <- nrow(sample_stats)
  n_fail_any <- sum(sample_stats$any_fail, na.rm = TRUE)
  pct_fail_any <- (n_fail_any / total_samples) * 100
  
  if (pct_fail_any > 15 && max_prune == 0) {
    warning(">15% of samples fail QC (", round(pct_fail_any, 1), "%). ",
            "To prune failing samples, re-run with max_prune > 0")
    message("\nTo prune samples, re-run with: max_prune = N ",
            "(where N is the number of samples to remove)")
  }
  
  while (iter < max_prune && n_fail_any > 0 && total_samples > 1) {
    iter <- iter + 1
    
    sample_stats$combined_deviation <- abs(sample_stats$sample_median - mean(sample_stats$sample_median)) +
      abs(sample_stats$sample_iqr - mean(sample_stats$sample_iqr))
    
    max_dev <- max(sample_stats$combined_deviation, na.rm = TRUE)
    worst_idx <- which(sample_stats$combined_deviation == max_dev)
    
    if (length(worst_idx) > 1) {
      worst_idx <- sample(worst_idx, 1)
    }
    
    removed <- as.character(sample_stats$sample_id[worst_idx])
    removed_samples <- c(removed_samples, removed)
    
    dat <- dat[dat[[sample_col]] != removed, ]
    sample_stats <- compute_stats(dat)
    sample_stats$any_fail <- sample_stats$median_fail | sample_stats$iqr_fail
    
    total_samples <- nrow(sample_stats)
    n_fail_any <- sum(sample_stats$any_fail, na.rm = TRUE)
    pct_fail_any <- (n_fail_any / total_samples) * 100
  }
  
  n_pass <- sum(!sample_stats$any_fail, na.rm = TRUE)
  pct_pass <- (n_pass / total_samples) * 100
  
  summary_df <- data.frame(
    total_samples = total_samples,
    n_fail_median = sum(sample_stats$median_fail, na.rm = TRUE),
    pct_fail_median = sum(sample_stats$median_fail, na.rm = TRUE) / total_samples * 100,
    n_fail_iqr = sum(sample_stats$iqr_fail, na.rm = TRUE),
    pct_fail_iqr = sum(sample_stats$iqr_fail, na.rm = TRUE) / total_samples * 100,
    n_fail_any = n_fail_any,
    pct_fail_any = pct_fail_any,
    n_pass = n_pass,
    pct_pass = pct_pass,
    n_removed = length(removed_samples)
  )
  
  summary_file <- file.path(output_dir, "olink_sample_qc_summary.csv")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(summary_df, file = summary_file, row.names = FALSE)
  
  passed_samples <- sample_stats$sample_id[!sample_stats$any_fail]
  failed_samples <- sample_stats$sample_id[sample_stats$any_fail]
  
  out_per_sample <- data.frame(
    sample_id = sample_stats$sample_id,
    sample_median = sample_stats$sample_median,
    sample_iqr = sample_stats$sample_iqr,
    global_median = sample_stats$global_median,
    global_iqr = sample_stats$global_iqr,
    median_z = sample_stats$median_z,
    iqr_z = sample_stats$iqr_z,
    median_fail = sample_stats$median_fail,
    iqr_fail = sample_stats$iqr_fail,
    any_fail = sample_stats$any_fail,
    stringsAsFactors = FALSE
  )
  
  if (mask_sample_ids && nrow(out_per_sample) > 0) {
    unique_ids <- unique(out_per_sample$sample_id)
    id_map <- setNames(paste0("Sample_", seq_along(unique_ids)), unique_ids)
    out_per_sample$sample_id <- id_map[out_per_sample$sample_id]
  }
  
  if (length(removed_samples) > 0 && max_prune > 0) {
    message("\n=== Olink Sample QC: Samples Removed ===")
    message("Number of samples removed: ", length(removed_samples))
    if (mask_sample_ids) {
      message("Sample IDs have been masked in saved outputs.")
    } else {
      message("Removed samples: ", paste(removed_samples, collapse = ", "))
    }
    message("\nSamples remaining: ", total_samples)
    message("To use only passing samples in subsequent steps, run:")
    message("  passed_sample_ids <- qc_result$passed_samples")
    message("==========================================\n")
  }
  
  list(
    summary = summary_df,
    per_sample = out_per_sample,
    passed_samples = passed_samples,
    failed_samples = failed_samples,
    removed_samples = removed_samples,
    summary_file = summary_file
  )
}
