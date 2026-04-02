olink_sample_qc <- function(df,
                           sample_id_col = "SAMPLE_ID",
                           sample_type_col = "SampleType",
                           sample_type_keep = "SAMPLE",
                           assay_type_col = "AssayType",
                           count_col = "Count",
                           assay_label = "assay",
                           ext_label = "ext_ctrl",
                           inc_label = "inc_ctrl",
                           amp_label = "amp_ctrl",
                           mask_sample_id = FALSE) {
  
  needed_cols <- c(sample_id_col, sample_type_col, assay_type_col, count_col)
  missing_cols <- needed_cols[!needed_cols %in% names(df)]
  
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  dat <- df[df[[sample_type_col]] == sample_type_keep,
            c(sample_id_col, assay_type_col, count_col)]
  
  if (nrow(dat) == 0) {
    stop("No rows remain after filtering to the requested sample type")
  }
  
  qc_counts <- dat %>%
    dplyr::group_by(.data[[sample_id_col]], .data[[assay_type_col]]) %>%
    dplyr::summarise(total_count = sum(.data[[count_col]], na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from = tidyselect::all_of(assay_type_col),
      values_from = total_count,
      values_fill = 0
    )
  
  needed_types <- c(assay_label, ext_label, inc_label, amp_label)
  missing_types <- needed_types[!needed_types %in% names(qc_counts)]
  
  for (nm in missing_types) {
    qc_counts[[nm]] <- 0
  }
  
  qc_counts <- qc_counts %>%
    dplyr::mutate(
      assay_pass = .data[[assay_label]] > 10000,
      ext_pass   = .data[[ext_label]] > 500,
      inc_pass   = .data[[inc_label]] > 500,
      amp_pass   = .data[[amp_label]] > 500,
      all_pass   = assay_pass & ext_pass & inc_pass & amp_pass
    )
  
  qc_summary <- qc_counts %>%
    dplyr::summarise(
      n_samples  = dplyr::n(),
      assay_fail = sum(!assay_pass, na.rm = TRUE),
      ext_fail   = sum(!ext_pass, na.rm = TRUE),
      inc_fail   = sum(!inc_pass, na.rm = TRUE),
      amp_fail   = sum(!amp_pass, na.rm = TRUE),
      any_fail   = sum(!all_pass, na.rm = TRUE)
    )
  
  if (mask_sample_id) {
    original_ids <- qc_counts[[sample_id_col]]
    masked_ids <- paste0("MASKED_SAMPLE_", sprintf("%04d", seq_along(original_ids)))
    
    id_map <- data.frame(
      original_sample_id = original_ids,
      masked_sample_id = masked_ids,
      stringsAsFactors = FALSE
    )
    
    qc_counts[[sample_id_col]] <- masked_ids
    passed_samples <- qc_counts[[sample_id_col]][qc_counts$all_pass]
  } else {
    id_map <- NULL
    passed_samples <- qc_counts[[sample_id_col]][qc_counts$all_pass]
  }
  
  return(list(
    qc_counts = qc_counts,
    qc_summary = qc_summary,
    passed_samples = passed_samples,
    id_map = id_map
  ))
}
