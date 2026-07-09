dnam_detection_pval_qc <- function(rgset, detection_p_threshold = 0.05) {
  dnam_require_package("minfi")
  dnam_ensure_manifest_package(rgset)

  det_p <- minfi::detectionP(rgset)
  mean_detection_p_value <- colMeans(det_p, na.rm = TRUE)

  per_sample <- data.frame(
    sample_id = names(mean_detection_p_value),
    mean_detection_p_value = as.numeric(mean_detection_p_value),
    detection_p_threshold = detection_p_threshold,
    passed_mean_detection_p_value = as.numeric(mean_detection_p_value) < detection_p_threshold,
    stringsAsFactors = FALSE
  )

  summary <- data.frame(
    metric = c(
      "n_samples",
      "detection_p_threshold",
      "n_samples_mean_detection_p_lt_threshold",
      "percent_samples_mean_detection_p_lt_threshold"
    ),
    value = c(
      nrow(per_sample),
      detection_p_threshold,
      sum(per_sample$passed_mean_detection_p_value, na.rm = TRUE),
      mean(per_sample$passed_mean_detection_p_value, na.rm = TRUE) * 100
    ),
    stringsAsFactors = FALSE
  )

  list(per_sample = per_sample, summary = summary)
}
