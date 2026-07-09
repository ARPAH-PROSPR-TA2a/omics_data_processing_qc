dnam_ensure_manifest_package <- function(rgset) {
  dnam_require_package("minfi")

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

dnam_control_probe_pca <- function(rgset,
                                   variance_threshold = 0.90,
                                   pseudocount = 100) {
  dnam_require_package("minfi")
  dnam_ensure_manifest_package(rgset)

  control_probes <- minfi::getControlAddress(rgset)
  control_probes <- intersect(control_probes, rownames(minfi::getRed(rgset)))
  if (!length(control_probes)) {
    stop("No control probes were found in the RGset.", call. = FALSE)
  }

  red_intensities <- minfi::getRed(rgset)[control_probes, , drop = FALSE]
  green_intensities <- minfi::getGreen(rgset)[control_probes, , drop = FALSE]
  control_beta <- red_intensities / (red_intensities + green_intensities + pseudocount)

  keep_probe <- apply(control_beta, 1, function(x) all(is.finite(x)) && stats::sd(x) > 0)
  control_beta <- control_beta[keep_probe, , drop = FALSE]
  if (nrow(control_beta) < 2) {
    stop("Fewer than two usable control probes remained after filtering.", call. = FALSE)
  }

  pca_object <- stats::prcomp(t(control_beta), center = TRUE, scale. = TRUE)
  variance_explained <- pca_object$sdev^2 / sum(pca_object$sdev^2)
  cumulative_variance_explained <- cumsum(variance_explained)
  n_pcs_to_threshold <- which(cumulative_variance_explained >= variance_threshold)[[1]]

  scores <- as.data.frame(pca_object$x[, seq_len(n_pcs_to_threshold), drop = FALSE])
  scores <- data.frame(sample_id = rownames(scores), scores, row.names = NULL, check.names = FALSE)

  variance <- data.frame(
    pc = paste0("PC", seq_along(variance_explained)),
    variance_explained = variance_explained,
    cumulative_variance_explained = cumulative_variance_explained,
    stringsAsFactors = FALSE
  )

  summary <- data.frame(
    metric = c("n_samples", "n_control_probes_used", "variance_threshold", "n_pcs_to_threshold"),
    value = c(ncol(rgset), nrow(control_beta), variance_threshold, n_pcs_to_threshold),
    stringsAsFactors = FALSE
  )

  list(
    scores = scores,
    variance_explained = variance,
    summary = summary,
    pca_object = pca_object
  )
}
