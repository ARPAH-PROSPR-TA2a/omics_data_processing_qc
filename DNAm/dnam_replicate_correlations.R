dnam_replicate_correlations <- function(beta_mat,
                                        pheno = NULL,
                                        sample_id_col = NULL,
                                        replicate_group_cols = NULL,
                                        pair_data = NULL,
                                        sample_1_col = "sample_1",
                                        sample_2_col = "sample_2",
                                        method = "pearson") {
  beta_mat <- as.matrix(beta_mat)

  if (!is.null(pair_data)) {
    required_cols <- c(sample_1_col, sample_2_col)
    missing_cols <- setdiff(required_cols, names(pair_data))
    if (length(missing_cols)) {
      stop("Pair file is missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
    }

    pairs <- data.frame(
      replicate_group = if ("replicate_group" %in% names(pair_data)) as.character(pair_data$replicate_group) else NA_character_,
      sample_1 = as.character(pair_data[[sample_1_col]]),
      sample_2 = as.character(pair_data[[sample_2_col]]),
      stringsAsFactors = FALSE
    )
  } else {
    if (is.null(pheno) || is.null(sample_id_col) || is.null(replicate_group_cols)) {
      stop("Provide either pair_data or pheno + sample_id_col + replicate_group_cols.", call. = FALSE)
    }

    missing_cols <- setdiff(c(sample_id_col, replicate_group_cols), names(pheno))
    if (length(missing_cols)) {
      stop("Phenotype data is missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
    }

    pheno <- pheno[pheno[[sample_id_col]] %in% colnames(beta_mat), , drop = FALSE]
    group_key <- do.call(paste, c(pheno[replicate_group_cols], sep = "||"))
    split_samples <- split(as.character(pheno[[sample_id_col]]), group_key)
    split_samples <- split_samples[lengths(split_samples) > 1]

    if (!length(split_samples)) {
      pairs <- data.frame(replicate_group = character(), sample_1 = character(), sample_2 = character())
    } else {
      pairs <- do.call(rbind, lapply(names(split_samples), function(group_name) {
        pair_matrix <- utils::combn(split_samples[[group_name]], 2)
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
    results <- pairs
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
