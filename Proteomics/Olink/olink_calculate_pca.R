olink_calculate_pca <- function(df,
                                  sample_col = "SAMPLE_ID",
                                  olink_id_col = "OlinkID",
                                  value_col = "NPX",
                                  assay_type_col = "AssayType",
                                  assay_keep = "assay",
                                  sample_type_col = "SampleType",
                                  sample_type_keep = "SAMPLE",
                                  impute_missing = TRUE,
                                  missingness_threshold = 0.10) {
  
  if (!is.data.frame(df)) {
    stop("df must be a data.frame")
  }
  
  required_cols <- c(sample_col, olink_id_col, value_col)
  for (col in required_cols) {
    if (!col %in% names(df)) {
      stop(paste0("Column '", col, "' not found in df"))
    }
  }
  
  has_sample_type <- sample_type_col %in% names(df)
  
  dat <- df[
    !is.na(df[[assay_type_col]]) &
      df[[assay_type_col]] == assay_keep &
      !is.na(df[[sample_col]]) &
      !is.na(df[[olink_id_col]]) &
      !is.na(df[[value_col]]) &
      (!has_sample_type | df[[sample_type_col]] == sample_type_keep),
    c(sample_col, olink_id_col, value_col)
  ]
  
  if (nrow(dat) == 0) {
    stop("No non-missing assay data available after filtering")
  }
  
  wide <- tidyr::pivot_wider(
    dat,
    id_cols = NULL,
    names_from = olink_id_col,
    values_from = value_col
  )
  
  sample_ids <- wide[[sample_col]]
  wide_matrix <- as.matrix(wide[, -1, drop = FALSE])
  rownames(wide_matrix) <- sample_ids
  
  n_samples <- nrow(wide_matrix)
  n_proteins_total <- ncol(wide_matrix)
  
  if (n_proteins_total < 3) {
    stop("Need at least 3 proteins for PCA")
  }
  if (n_samples < 3) {
    stop("Need at least 3 samples for PCA")
  }
  
  missingness_threshold_adj <- missingness_threshold
  if (n_samples <= 88) {
    missingness_threshold_adj <- 0.05
    message("Using 5% missingness threshold for small dataset (n <= 88)")
  }
  
  percent_missing <- colMeans(is.na(wide_matrix))
  proteins_to_keep <- percent_missing <= missingness_threshold_adj
  
  if (any(!proteins_to_keep)) {
    n_removed <- sum(!proteins_to_keep)
    message(n_removed, " protein(s) removed due to high missingness (>", 
            missingness_threshold_adj * 100, "%)")
    wide_matrix <- wide_matrix[, proteins_to_keep, drop = FALSE]
  }
  
  n_proteins_remaining <- ncol(wide_matrix)
  if (n_proteins_remaining < 3) {
    stop("Too few proteins remain after filtering for PCA (need >= 3)")
  }
  
  if (impute_missing && any(is.na(wide_matrix))) {
    for (j in seq_len(ncol(wide_matrix))) {
      if (any(is.na(wide_matrix[, j]))) {
        wide_matrix[is.na(wide_matrix[, j]), j] <- median(wide_matrix[, j], na.rm = TRUE)
      }
    }
  }
  
  if (any(colSums(is.na(wide_matrix)) > 0)) {
    stop("Missing values remain after imputation. Check data quality.")
  }
  
  pca <- stats::prcomp(wide_matrix, scale. = TRUE, center = TRUE)
  
  variance_explained <- pca$sdev^2 / sum(pca$sdev^2)
  
  n_pcs_to_return <- min(5, ncol(pca$x))
  
  loadings_cols <- c("protein", paste0("PC", 1:n_pcs_to_return))
  loadings <- data.frame(
    protein = rownames(pca$rotation),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(n_pcs_to_return)) {
    loadings[[paste0("PC", i)]] <- pca$rotation[, i]
  }
  
  scores <- data.frame(
    SampleID = sample_ids,
    stringsAsFactors = FALSE
  )
  for (i in seq_len(n_pcs_to_return)) {
    scores[[paste0("PC", i)]] <- pca$x[, i]
  }
  
  list(
    scores = scores,
    variance_explained = variance_explained,
    loadings = loadings,
    n_proteins = n_proteins_remaining,
    n_proteins_total = n_proteins_total,
    n_samples = n_samples,
    pca_object = pca,
    proteins_removed = sum(!proteins_to_keep),
    missingness_threshold_used = missingness_threshold_adj
  )
}
