plate_control_qc <- function(df,
                                   sample_type_col = "SampleType",
                                   sample_type_keep = "PLATE_CONTROL",
                                   block_col = "Block",
                                   plate_col = "PlateID",
                                   replicate_col = "WellID",
                                   assay_type_col = "AssayType",
                                   count_col = "Count",
                                   npx_col = "NPX",
                                   assay_label = "assay",
                                   ext_label = "ext_ctrl",
                                   inc_label = "inc_ctrl",
                                   amp_label = "amp_ctrl") {
  
  if (!is.data.frame(df)) {
    stop("df must be a data.frame")
  }
  
  needed_cols <- c(sample_type_col, block_col, plate_col,
                   replicate_col, assay_type_col, count_col, npx_col)
  missing_cols <- needed_cols[!needed_cols %in% names(df)]
  
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # Keep only requested sample type
  dat <- df[
    !is.na(df[[sample_type_col]]) &
      df[[sample_type_col]] == sample_type_keep,
    c(block_col, plate_col, replicate_col, assay_type_col, count_col)
  ]
  
  if (nrow(dat) == 0) {
    stop("No rows remain after filtering to the requested sample type")
  }
  
  # Summarize count within block + plate + replicate + assay type
  qc_counts <- dat %>%
    dplyr::group_by(
      .data[[block_col]],
      .data[[plate_col]],
      .data[[replicate_col]],
      .data[[assay_type_col]]
    ) %>%
    dplyr::summarise(total_count = sum(.data[[count_col]], na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from = tidyselect::all_of(assay_type_col),
      values_from = total_count,
      values_fill = 0
    )
  
  # Make sure expected assay/control columns exist
  needed_types <- c(assay_label, ext_label, inc_label, amp_label)
  missing_types <- needed_types[!needed_types %in% names(qc_counts)]
  
  for (nm in missing_types) {
    qc_counts[[nm]] <- 0
  }
  
  # Fixed QC thresholds
  qc_counts <- qc_counts %>%
    dplyr::mutate(
      assay_pass = .data[[assay_label]] > 10000,
      ext_pass   = .data[[ext_label]] > 1000,
      inc_pass   = .data[[inc_label]] > 1000,
      amp_pass   = .data[[amp_label]] > 500,
      all_pass   = assay_pass & ext_pass & inc_pass & amp_pass
    )
  
  # Plate-level summary: at least 3 replicates must pass
  qc_summary <- qc_counts %>%
    dplyr::group_by(.data[[block_col]], .data[[plate_col]]) %>%
    dplyr::summarise(
      n_plate_controls = dplyr::n(),
      n_pass = sum(all_pass, na.rm = TRUE),
      n_fail = sum(!all_pass, na.rm = TRUE),
      plate_pass = n_pass >= 3,
      .groups = "drop"
    )
  
  # Failed plates as user-friendly labels
  failed_plate_table <- qc_summary %>%
    dplyr::filter(!plate_pass) %>%
    dplyr::mutate(
      failed_plate = paste0("Block_", .data[[block_col]], "_Plate_", .data[[plate_col]])
    )
  
  failed_plates <- as.character(failed_plate_table$failed_plate)
  
  failed_plates_message <- if (length(failed_plates) == 0) {
    "No plates failed"
  } else {
    failed_plates
  }
  
  # Check whether failed plates have any non-missing NPX values in the original data
  if (nrow(failed_plate_table) == 0) {
    failed_plate_npx_check <- data.frame(
      failed_plate = character(0),
      has_reported_npx = logical(0),
      n_reported_npx_rows = integer(0),
      stringsAsFactors = FALSE
    )
    
    failed_plate_npx_message <- "No failed plates, so no NPX check needed"
    
  } else {
    failed_npx <- df %>%
      dplyr::inner_join(
        failed_plate_table %>% dplyr::select(dplyr::all_of(c(block_col, plate_col)), failed_plate),
        by = setNames(c(block_col, plate_col), c(block_col, plate_col))
      ) %>%
      dplyr::group_by(failed_plate) %>%
      dplyr::summarise(
        n_reported_npx_rows = sum(!is.na(.data[[npx_col]])),
        has_reported_npx = n_reported_npx_rows > 0,
        .groups = "drop"
      )
    
    failed_plate_npx_check <- failed_plate_table %>%
      dplyr::select(failed_plate) %>%
      dplyr::left_join(failed_npx, by = "failed_plate")
    
    failed_plate_npx_check$n_reported_npx_rows[is.na(failed_plate_npx_check$n_reported_npx_rows)] <- 0
    failed_plate_npx_check$has_reported_npx[is.na(failed_plate_npx_check$has_reported_npx)] <- FALSE
    
    if (any(failed_plate_npx_check$has_reported_npx)) {
      failed_plate_npx_message <- "Some failed plates had reported NPX values"
    } else {
      failed_plate_npx_message <- "No failed plates had reported NPX values"
    }
  }
  
  return(list(
    qc_counts = qc_counts,
    qc_summary = qc_summary,
    failed_plates = failed_plates_message,
    failed_plate_npx_check = failed_plate_npx_check,
    failed_plate_npx_message = failed_plate_npx_message
  ))
}