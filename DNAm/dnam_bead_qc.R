dnam_extract_nbeads_matrix <- function(rgset) {
  dnam_require_package("Biobase")

  if ("getNBeads" %in% getNamespaceExports("minfi")) {
    nbeads <- try(minfi::getNBeads(rgset), silent = TRUE)
    if (!inherits(nbeads, "try-error")) {
      return(nbeads)
    }
  }

  assay_elements <- Biobase::assayDataElementNames(rgset)
  if ("NBeads" %in% assay_elements) {
    return(Biobase::assayDataElement(rgset, "NBeads"))
  }

  if (requireNamespace("wateRmelon", quietly = TRUE)) {
    return(wateRmelon::beadcount(rgset))
  }

  stop(
    "Could not extract bead counts. Load IDATs with extended = TRUE, or install wateRmelon.",
    call. = FALSE
  )
}

dnam_bead_qc <- function(rgset, mean_bead_threshold = 2) {
  nbeads <- dnam_extract_nbeads_matrix(rgset)
  mean_bead_count <- colMeans(nbeads, na.rm = TRUE)

  per_sample <- data.frame(
    sample_id = names(mean_bead_count),
    mean_bead_count = as.numeric(mean_bead_count),
    bead_mean_threshold = mean_bead_threshold,
    passed_mean_bead_count = as.numeric(mean_bead_count) > mean_bead_threshold,
    stringsAsFactors = FALSE
  )

  summary <- data.frame(
    metric = c(
      "n_samples",
      "bead_mean_threshold",
      "n_samples_mean_bead_count_gt_threshold",
      "percent_samples_mean_bead_count_gt_threshold"
    ),
    value = c(
      nrow(per_sample),
      mean_bead_threshold,
      sum(per_sample$passed_mean_bead_count, na.rm = TRUE),
      mean(per_sample$passed_mean_bead_count, na.rm = TRUE) * 100
    ),
    stringsAsFactors = FALSE
  )

  list(per_sample = per_sample, summary = summary)
}
