somascan_techrep_cor <- function(dat,
                                 subject_id_col,
                                 time_col,
                                 replicate_ids = NULL,
                                 use_qc = TRUE,
                                 sample_type_col = "SampleType",
                                 qc_label = "QC",
                                 plate_id_col = "PlateId",
                                 seq_prefix = "seq.",
                                 method = "pearson",
                                 log2_transform = TRUE,
                                 mask_sample_ids = TRUE) {

  stopifnot(is.data.frame(dat))
  required_cols <- c(subject_id_col, time_col, plate_id_col)
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  analyte_cols <- grep(paste0("^", seq_prefix), names(dat), value = TRUE)
  if (length(analyte_cols) == 0) {
    stop("No analyte columns found (expected columns starting with '", seq_prefix, "').")
  }

  # Keep rows that can be assigned to a subject-timepoint replicate group.
  d <- dat[!is.na(dat[[subject_id_col]]) & !is.na(dat[[time_col]]), , drop = FALSE]

  if (!is.null(replicate_ids)) {
    d <- d[d[[subject_id_col]] %in% replicate_ids, , drop = FALSE]
  }

  if (!use_qc && sample_type_col %in% names(d)) {
    d <- d[d[[sample_type_col]] != qc_label, , drop = FALSE]
  }

  if (nrow(d) == 0) {
    stop("No rows selected. Check subject_id_col, time_col, replicate_ids, and use_qc settings.")
  }

  group_key <- paste(as.character(d[[subject_id_col]]), as.character(d[[time_col]]), sep = "__TIME__")
  tab <- table(group_key)
  rep_keys <- names(tab)[tab >= 2]
  if (length(rep_keys) == 0) {
    return(list(
      n_rows_used = nrow(d),
      n_ids_with_reps = 0,
      results = data.frame()
    ))
  }

  X <- as.matrix(d[, analyte_cols, drop = FALSE])
  storage.mode(X) <- "double"
  if (log2_transform) X <- log2(X)

  cor_row <- function(i, j) {
    stats::cor(X[i, ], X[j, ], method = method, use = "pairwise.complete.obs")
  }

  out_list <- vector("list", length(rep_keys))
  names(out_list) <- rep_keys

  for (g in rep_keys) {
    idx <- which(group_key == g)
    cmb <- utils::combn(idx, 2)
    m <- ncol(cmb)
    subject_id <- as.character(d[[subject_id_col]][idx[1]])
    timepoint <- as.character(d[[time_col]][idx[1]])

    res <- data.frame(
      SubjectId = rep(subject_id, m),
      Timepoint = rep(timepoint, m),
      row_i = cmb[1, ],
      row_j = cmb[2, ],
      PlateId_i = as.character(d[[plate_id_col]][cmb[1, ]]),
      PlateId_j = as.character(d[[plate_id_col]][cmb[2, ]]),
      same_plate = d[[plate_id_col]][cmb[1, ]] == d[[plate_id_col]][cmb[2, ]],
      r = NA_real_,
      stringsAsFactors = FALSE
    )

    for (k in seq_len(m)) {
      res$r[k] <- cor_row(res$row_i[k], res$row_j[k])
    }

    out_list[[g]] <- res
  }

  results <- do.call(rbind, out_list)
  results$row_i <- NULL
  results$row_j <- NULL

  if (mask_sample_ids && nrow(results) > 0) {
    unique_ids <- unique(results$SubjectId)
    id_map <- setNames(paste0("Sample_", seq_along(unique_ids)), unique_ids)
    results$SubjectId <- id_map[results$SubjectId]
  }

  list(
    n_rows_used = nrow(d),
    n_ids_with_reps = length(rep_keys),
    results = results
  )
}
