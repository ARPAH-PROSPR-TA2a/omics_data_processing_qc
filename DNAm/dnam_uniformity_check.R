dnam_uniformity_check <- function(beta_mat, p_threshold = 0.05) {
  dnam_require_package(
    "diptest",
    "Install with: install.packages(\"diptest\")"
  )

  beta_mat <- as.matrix(beta_mat)
  storage.mode(beta_mat) <- "double"

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
    sample_id = names(dip_p_values),
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
