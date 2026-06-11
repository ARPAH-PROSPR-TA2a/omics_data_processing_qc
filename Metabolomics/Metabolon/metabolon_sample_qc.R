metabolon_sample_qc <- function(phenotype_df,
                                metabolite_df,
                                sample_id_col = "SampleID",
                                matrix_id_col = NULL,
                                cutoff = 3,
                                prune = 0,
                                target_sample_prop = 0.80,
                                mask_sample_ids = FALSE,
                                verbose = TRUE) {
  if (!is.data.frame(phenotype_df)) stop("phenotype_df must be a data.frame")
  if (!(is.data.frame(metabolite_df) || is.matrix(metabolite_df))) stop("metabolite_df must be a data.frame or matrix")
  if (!sample_id_col %in% names(phenotype_df)) stop("Missing sample ID column in phenotype_df: ", sample_id_col)
  if (!is.numeric(cutoff) || length(cutoff) != 1 || is.na(cutoff)) stop("cutoff must be a single numeric value")
  if (!is.numeric(prune) || prune < 0 || prune != floor(prune)) stop("prune must be a non-negative integer")
  if (!is.numeric(target_sample_prop) || length(target_sample_prop) != 1 || is.na(target_sample_prop)) stop("target_sample_prop must be a single numeric value")

  get_matrix_ids <- function(x) {
    if (!is.null(matrix_id_col) && is.data.frame(x) && matrix_id_col %in% names(x)) {
      ids <- as.character(x[[matrix_id_col]])
      x <- x[, setdiff(names(x), matrix_id_col), drop = FALSE]
      return(list(ids = ids, mat = as.matrix(x)))
    }
    if (is.null(rownames(x))) stop("metabolite_df must have rownames or a matrix_id_col")
    list(ids = as.character(rownames(x)), mat = as.matrix(x))
  }

  mat_info <- get_matrix_ids(metabolite_df)
  matrix_ids <- mat_info$ids
  X0 <- mat_info$mat

  pheno_ids <- as.character(phenotype_df[[sample_id_col]])
  if (anyDuplicated(pheno_ids)) warning("Duplicate phenotype sample IDs detected in ", sample_id_col)
  if (anyDuplicated(matrix_ids)) warning("Duplicate metabolite sample IDs detected in matrix IDs")

  in_matrix_not_pheno <- setdiff(matrix_ids, pheno_ids)
  in_pheno_not_matrix <- setdiff(pheno_ids, matrix_ids)
  if (length(in_matrix_not_pheno) > 0) {
    warning("A: ", length(in_matrix_not_pheno), " metabolite matrix sample(s) not found in phenotype data. They will be ignored.")
  }
  if (length(in_pheno_not_matrix) > 0) {
    warning("B: ", length(in_pheno_not_matrix), " phenotype sample(s) not found in metabolite matrix. They will be ignored.")
  }

  matched_ids <- intersect(pheno_ids, matrix_ids)
  if (length(matched_ids) == 0) stop("No overlapping sample IDs between phenotype data and metabolite matrix")

  pheno_idx <- match(matched_ids, pheno_ids)
  matrix_idx <- match(matched_ids, matrix_ids)
  pheno_matched <- phenotype_df[pheno_idx, , drop = FALSE]
  rownames(pheno_matched) <- matched_ids
  X <- X0[matrix_idx, , drop = FALSE]
  rownames(X) <- matched_ids

  compute_stats <- function(mat) {
    sample_median <- apply(mat, 1, median, na.rm = TRUE)
    sample_iqr <- apply(mat, 1, stats::IQR, na.rm = TRUE)
    median_mean <- mean(sample_median, na.rm = TRUE)
    median_sd <- stats::sd(sample_median, na.rm = TRUE)
    iqr_mean <- mean(sample_iqr, na.rm = TRUE)
    iqr_sd <- stats::sd(sample_iqr, na.rm = TRUE)
    median_low <- median_mean - cutoff * median_sd
    median_high <- median_mean + cutoff * median_sd
    iqr_low <- iqr_mean - cutoff * iqr_sd
    iqr_high <- iqr_mean + cutoff * iqr_sd
    data.frame(
      sample_id = rownames(mat),
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

  removed_samples <- character(0)
  iter <- 0
  sample_stats <- compute_stats(X)
  sample_stats$any_fail <- sample_stats$median_fail | sample_stats$iqr_fail

  total_samples <- nrow(sample_stats)
  n_pass <- sum(!sample_stats$any_fail, na.rm = TRUE)
  pct_pass <- 100 * n_pass / total_samples
  pct_fail_any <- 100 - pct_pass

  if (pct_pass < 80 && prune == 0) {
    warning("<80% of matched samples pass QC (", round(pct_pass, 1), "%). Re-run with prune > 0 to remove worst samples and define the analysis dataset.")
    if (verbose) {
      message("\n==================================================")
      message("WARNING: LESS THAN 80% OF SAMPLES PASS QC")
      message("Only ", round(pct_pass, 1), "% of matched samples passed this QC step.")
      message("Prune a small number of low-quality samples, re-run QC, and save the filtered analysis dataset before continuing.")
      message("The filtered dataset is the analysis dataset for downstream steps.")
      message("==================================================\n")
    }
  }

  while (iter < prune && pct_pass < 80 && total_samples > 1) {
    iter <- iter + 1
    sample_stats$combined_deviation <- abs(sample_stats$sample_median - mean(sample_stats$sample_median, na.rm = TRUE)) +
      abs(sample_stats$sample_iqr - mean(sample_stats$sample_iqr, na.rm = TRUE))
    max_dev <- max(sample_stats$combined_deviation, na.rm = TRUE)
    worst_idx <- which(sample_stats$combined_deviation == max_dev)
    if (length(worst_idx) > 1) worst_idx <- sample(worst_idx, 1)
    removed <- sample_stats$sample_id[worst_idx]
    removed_samples <- c(removed_samples, removed)
    X <- X[rownames(X) != removed, , drop = FALSE]
    pheno_matched <- pheno_matched[rownames(pheno_matched) != removed, , drop = FALSE]
    sample_stats <- compute_stats(X)
    sample_stats$any_fail <- sample_stats$median_fail | sample_stats$iqr_fail
    total_samples <- nrow(sample_stats)
    n_pass <- sum(!sample_stats$any_fail, na.rm = TRUE)
    pct_pass <- 100 * n_pass / total_samples
    pct_fail_any <- 100 - pct_pass
  }

  summary_df <- data.frame(
    total_samples = total_samples,
    n_pass = n_pass,
    pct_pass = pct_pass,
    n_fail_any = sum(sample_stats$any_fail, na.rm = TRUE),
    pct_fail_any = pct_fail_any,
    n_fail_median = sum(sample_stats$median_fail, na.rm = TRUE),
    n_fail_iqr = sum(sample_stats$iqr_fail, na.rm = TRUE),
    n_removed = length(removed_samples),
    stringsAsFactors = FALSE
  )

  per_sample <- data.frame(
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

  if (mask_sample_ids && nrow(per_sample) > 0) {
    id_map <- setNames(paste0("Sample_", seq_along(per_sample$sample_id)), per_sample$sample_id)
    per_sample$sample_id <- id_map[per_sample$sample_id]
  }

  if (length(removed_samples) > 0 && verbose) {
    message("\n=== Metabolon Sample QC: Samples Removed ===")
    message("Number of samples removed: ", length(removed_samples))
    message("Removed samples: ", paste(removed_samples, collapse = ", "))
    message("Samples remaining: ", total_samples)
    message("To use only passing samples in subsequent steps, run:")
    message("  analysis_ids <- qc_result$analysis_sample_ids")
    message("===========================================\n")
  }

  list(
    summary = summary_df,
    per_sample = per_sample,
    matched_sample_ids = matched_ids,
    passed_samples = sample_stats$sample_id[!sample_stats$any_fail],
    failed_samples = sample_stats$sample_id[sample_stats$any_fail],
    removed_samples = removed_samples,
    analysis_sample_ids = sample_stats$sample_id[!sample_stats$any_fail],
    analysis_matrix = X,
    analysis_pheno = pheno_matched
  )
}
