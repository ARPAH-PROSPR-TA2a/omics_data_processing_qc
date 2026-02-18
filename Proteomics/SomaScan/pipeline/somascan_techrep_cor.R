# # QC-only technical replicate correlations
# qc_cor <- somascan_techrep_cor(sl)
# 
# # QC + a few biological duplicated SampleIds
# qc_plus_bio <- somascan_techrep_cor(sl, replicate_ids = c("S123", "S456"))

somascan_techrep_cor <- function(dat,
                                  replicate_ids = NULL,
                                  use_qc = TRUE,
                                  sample_type_col = "SampleType",
                                  qc_label = "QC",
                                  sample_id_col = "SampleId",
                                  plate_id_col = "PlateId",
                                  seq_prefix = "seq.",
                                  method = "pearson",
                                  log2_transform = TRUE,
                                  include_only_baseline = TRUE,
                                  time_col = NULL,
                                  baseline_values = NULL,
                                  mask_sample_ids = FALSE) {
  
  stopifnot(is.data.frame(dat))
  if (!all(c(sample_id_col, plate_id_col) %in% names(dat))) {
    stop("Missing required columns: ", sample_id_col, " and/or ", plate_id_col)
  }

  if (include_only_baseline) {
    if (is.null(time_col) || is.null(baseline_values)) {
      stop("include_only_baseline = TRUE requires both 'time_col' and 'baseline_values' to be specified. ",
           "Set include_only_baseline = FALSE to skip time filtering.")
    }
    if (!(time_col %in% names(dat))) {
      stop("time_col '", time_col, "' not found in data")
    }
    dat <- dat[dat[[time_col]] %in% baseline_values, , drop = FALSE]
    if (nrow(dat) == 0) {
      stop("No samples remaining after filtering to baseline (time_col='", time_col, 
           "', baseline_values=", paste(baseline_values, collapse=", "), ")")
    }
  } else if (!is.null(time_col) && !is.null(baseline_values)) {
    message("Note: time_col and baseline_values provided but include_only_baseline = FALSE. ",
            "Ignoring time filtering.")
  }

  analyte_cols <- grep(paste0("^", seq_prefix), names(dat), value = TRUE)
  if (length(analyte_cols) == 0) stop("No analyte columns found (expected columns starting with '", seq_prefix, "').")
  
  # ---- choose rows to include ----
  keep_idx <- rep(FALSE, nrow(dat))
  
  if (use_qc) {
    if (!(sample_type_col %in% names(dat))) stop("Missing column: ", sample_type_col)
    keep_idx <- keep_idx | (dat[[sample_type_col]] == qc_label)
  }
  
  if (!is.null(replicate_ids)) {
    keep_idx <- keep_idx | (dat[[sample_id_col]] %in% replicate_ids)
  }
  
  d <- dat[keep_idx, , drop = FALSE]
  if (nrow(d) == 0) stop("No rows selected. Check use_qc / qc_label / replicate_ids.")
  
  # Only groups with >=2 occurrences can form replicate pairs
  tab <- table(d[[sample_id_col]])
  rep_ids <- names(tab)[tab >= 2]
  if (length(rep_ids) == 0) {
    return(list(
      n_rows_used = nrow(d),
      n_ids_with_reps = 0,
      results = data.frame()
    ))
  }
  
  # Prepare RFU matrix
  X <- as.matrix(d[, analyte_cols, drop = FALSE])
  storage.mode(X) <- "double"
  if (log2_transform) X <- log2(X)
  
  # Helper: correlate two rows with missing handling
  cor_row <- function(i, j) {
    stats::cor(X[i, ], X[j, ], method = method, use = "pairwise.complete.obs")
  }
  
  # Build pairwise results for each SampleId
  out_list <- vector("list", length(rep_ids))
  names(out_list) <- rep_ids
  
  for (sid in rep_ids) {
    idx <- which(d[[sample_id_col]] == sid)
    cmb <- utils::combn(idx, 2)
    
    # Preallocate
    m <- ncol(cmb)
    res <- data.frame(
      SampleId = rep(sid, m),
      row_i = cmb[1, ],
      row_j = cmb[2, ],
      PlateId_i = as.character(d[[plate_id_col]][cmb[1, ]]),
      PlateId_j = as.character(d[[plate_id_col]][cmb[2, ]]),
      same_plate = d[[plate_id_col]][cmb[1, ]] == d[[plate_id_col]][cmb[2, ]],
      r = NA_real_,
      stringsAsFactors = FALSE
    )
    
    # Compute correlations
    for (k in seq_len(m)) {
      res$r[k] <- cor_row(res$row_i[k], res$row_j[k])
    }
    
    out_list[[sid]] <- res
  }
  
  results <- do.call(rbind, out_list)
  
  # Drop internal row indices (optional; keep if useful for debugging)
  results$row_i <- NULL
  results$row_j <- NULL
  
  if (mask_sample_ids && nrow(results) > 0) {
    unique_ids <- unique(c(results$SampleId))
    id_map <- setNames(paste0("Sample_", seq_along(unique_ids)), unique_ids)
    results$SampleId <- id_map[results$SampleId]
  }
  
  list(
    n_rows_used = nrow(d),
    n_ids_with_reps = length(rep_ids),
    results = results
  )
}


