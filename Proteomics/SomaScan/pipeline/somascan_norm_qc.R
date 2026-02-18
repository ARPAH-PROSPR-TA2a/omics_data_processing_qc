somascan_norm_qc <- function(dat,
                             lower = 0.4,
                             upper = 2.5,
                             sample_type_col = "SampleType",
                             sample_id_col   = "SampleId",
                             sample_label    = "Sample",
                             mask_sample_ids = FALSE) {
  
  required_cols <- c("NormScale_20", "NormScale_0_5", "NormScale_0_005")
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  samples <- dat[dat[[sample_type_col]] == sample_label, , drop = FALSE]
  if (nrow(samples) == 0) stop("No Sample rows found.")
  
  # Flag per dilution
  fail_20     <- samples$NormScale_20     < lower | samples$NormScale_20     > upper
  fail_0_5    <- samples$NormScale_0_5    < lower | samples$NormScale_0_5    > upper
  fail_0_005  <- samples$NormScale_0_005  < lower | samples$NormScale_0_005  > upper
  
  # Any dilution failure
  fail_any <- fail_20 | fail_0_5 | fail_0_005
  
  per_sample <- data.frame(
    SampleId = samples[[sample_id_col]],
    NormScale_20_fail    = fail_20,
    NormScale_0_5_fail   = fail_0_5,
    NormScale_0_005_fail = fail_0_005,
    fail_any = fail_any,
    stringsAsFactors = FALSE
  )
  
  if (mask_sample_ids) {
    unique_ids <- unique(per_sample$SampleId)
    id_map <- setNames(paste0("Sample_", seq_along(unique_ids)), unique_ids)
    per_sample$SampleId <- id_map[per_sample$SampleId]
  }
  
  list(
    n_samples = nrow(samples),
    fail_counts = list(
      NormScale_20     = sum(fail_20, na.rm = TRUE),
      NormScale_0_5    = sum(fail_0_5, na.rm = TRUE),
      NormScale_0_005  = sum(fail_0_005, na.rm = TRUE),
      any_dilution     = sum(fail_any, na.rm = TRUE)
    ),
    pass_counts = list(
      NormScale_20     = sum(!fail_20, na.rm = TRUE),
      NormScale_0_5    = sum(!fail_0_5, na.rm = TRUE),
      NormScale_0_005  = sum(!fail_0_005, na.rm = TRUE),
      all_dilutions    = sum(!fail_any, na.rm = TRUE)
    ),
    per_sample = per_sample
  )
}