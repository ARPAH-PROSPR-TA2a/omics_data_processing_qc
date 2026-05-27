metabolon_calculate_pca <- function(mat,
                                    sample_id_col = NULL,
                                    log2_transform = FALSE,
                                    missingness_threshold = 0.10,
                                    mask_sample_ids = FALSE) {
  if (!(is.data.frame(mat) || is.matrix(mat))) stop("mat must be a data.frame or matrix")
  X <- as.matrix(mat)

  sample_ids <- if (!is.null(sample_id_col) && sample_id_col %in% colnames(mat)) {
    ids <- as.character(mat[[sample_id_col]])
    X <- X[, setdiff(colnames(X), sample_id_col), drop = FALSE]
    ids
  } else if (!is.null(rownames(X))) {
    rownames(X)
  } else {
    stop("Matrix must have rownames or a sample_id_col")
  }

  rownames(X) <- sample_ids
  storage.mode(X) <- "double"
  if (log2_transform) X <- log2(X)

  n_samples <- nrow(X)
  n_metabolites_total <- ncol(X)
  if (n_samples < 3) stop("Need at least 3 samples for PCA")
  if (n_metabolites_total < 3) stop("Need at least 3 metabolites for PCA")

  threshold_adj <- missingness_threshold
  if (n_samples <= 88) {
    threshold_adj <- 0.05
    message("Using 5% missingness threshold for small dataset (n <= 88)")
  }

  percent_missing <- colMeans(is.na(X))
  keep <- percent_missing <= threshold_adj
  if (any(!keep)) {
    message(sum(!keep), " metabolite(s) removed due to high missingness (>", threshold_adj * 100, "%)")
    X <- X[, keep, drop = FALSE]
  }

  if (ncol(X) < 3) stop("Too few metabolites remain after filtering for PCA")

  if (any(is.na(X))) {
    for (j in seq_len(ncol(X))) {
      if (any(is.na(X[, j]))) {
        X[is.na(X[, j]), j] <- median(X[, j], na.rm = TRUE)
      }
    }
  }
  if (any(is.na(X))) stop("Missing values remain after imputation")

  zero_var <- apply(X, 2, stats::var, na.rm = TRUE)
  keep_var <- is.finite(zero_var) & zero_var > 0
  if (any(!keep_var)) {
    X <- X[, keep_var, drop = FALSE]
  }
  if (ncol(X) < 3) stop("Too few variable metabolites remain after filtering for PCA")

  pca <- stats::prcomp(X, center = TRUE, scale. = TRUE)
  variance_explained <- pca$sdev^2 / sum(pca$sdev^2)
  n_pcs_to_return <- min(5, ncol(pca$x))

  loadings <- data.frame(metabolite = rownames(pca$rotation), stringsAsFactors = FALSE)
  for (i in seq_len(n_pcs_to_return)) loadings[[paste0("PC", i)]] <- pca$rotation[, i]

  scores <- data.frame(SampleID = rownames(X), stringsAsFactors = FALSE)
  for (i in seq_len(n_pcs_to_return)) scores[[paste0("PC", i)]] <- pca$x[, i]

  masked_scores <- scores
  if (mask_sample_ids && nrow(masked_scores) > 0) {
    id_map <- setNames(paste0("Sample_", seq_along(masked_scores$SampleID)), masked_scores$SampleID)
    masked_scores$SampleID <- id_map[masked_scores$SampleID]
  }

  list(
    scores = scores,
    masked_scores = masked_scores,
    variance_explained = variance_explained,
    loadings = loadings,
    pca_object = pca,
    n_samples = nrow(X),
    n_metabolites = ncol(X),
    n_metabolites_total = n_metabolites_total,
    metabolites_removed = sum(!keep),
    missingness_threshold_used = threshold_adj
  )
}
