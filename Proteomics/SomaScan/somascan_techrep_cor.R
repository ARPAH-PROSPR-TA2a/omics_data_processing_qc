somascan_techrep_cor <- function(dat,
                                 subject_id_col,
                                 time_col,
                                 replicate_ids = NULL,
                                 use_qc = TRUE,
                                 sample_type_col = "SampleType",
                                 qc_label = "QC",
                                 plate_id_col = "PlateId",
                                 sample_id_col = "SampleId",
                                 seq_prefix = "seq.",
                                 method = "pearson",
                                 log2_transform = TRUE,
                                 mask_sample_ids = TRUE,
                                 id_map = NULL) {

  stopifnot(is.data.frame(dat))
  required_cols <- c(subject_id_col, time_col, plate_id_col, sample_id_col)
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

  # Normalize subject/time values so whitespace or factor-label differences
  # do not silently split replicate groups.
  subj <- trimws(as.character(d[[subject_id_col]]))
  time <- trimws(as.character(d[[time_col]]))
  group_key <- paste(subj, time, sep = "__TIME__")
  tab <- table(group_key)

  # Only groups with >= 2 members are technical replicates.
  rep_keys <- names(tab)[tab >= 2]

  # Per-group replicate table (matched pairs only; singletons are excluded).
  rep_member_ids <- vapply(rep_keys, function(g) {
    paste(d[[sample_id_col]][group_key == g], collapse = "|")
  }, character(1))

  groups <- data.frame(
    SubjectId = sub("__TIME__.*$", "", rep_keys),
    Timepoint = sub("^.*__TIME__", "", rep_keys),
    n_members = as.integer(tab[rep_keys]),
    SampleIds = rep_member_ids,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # Masking helpers. Member (sample-well) IDs map through the LOD id_map when
  # provided so labels agree with per_sample_lod; otherwise a local map is used.
  if (mask_sample_ids && is.null(id_map)) {
    samp_ids <- unique(as.character(d[[sample_id_col]]))
    id_map <- setNames(paste0("Sample_", seq_along(samp_ids)), samp_ids)
  }
  map_or_raw <- function(ids) {
    if (!mask_sample_ids) return(ids)
    mapped <- id_map[ids]
    mapped[is.na(mapped)] <- ids[is.na(mapped)]
    mapped
  }

  # Subject IDs are masked with a distinct prefix (Subj_#) so they never
  # collide with the per-sample Sample_# labels used elsewhere.
  subj_map <- NULL
  if (mask_sample_ids) {
    unique_subj <- unique(subj)
    subj_map <- setNames(paste0("Subj_", seq_along(unique_subj)), unique_subj)
  }
  mask_subj <- function(ids) {
    if (!mask_sample_ids) return(ids)
    mapped <- subj_map[ids]
    mapped[is.na(mapped)] <- ids[is.na(mapped)]
    mapped
  }

  results <- data.frame(
    SubjectId = character(0),
    Timepoint = character(0),
    Sample_i = character(0),
    Sample_j = character(0),
    PlateId_i = character(0),
    PlateId_j = character(0),
    same_plate = logical(0),
    r = numeric(0),
    stringsAsFactors = FALSE
  )

  if (length(rep_keys) > 0) {
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

      res <- data.frame(
        SubjectId = rep(subj[idx[1]], m),
        Timepoint = rep(time[idx[1]], m),
        Sample_i = d[[sample_id_col]][cmb[1, ]],
        Sample_j = d[[sample_id_col]][cmb[2, ]],
        PlateId_i = as.character(d[[plate_id_col]][cmb[1, ]]),
        PlateId_j = as.character(d[[plate_id_col]][cmb[2, ]]),
        same_plate = d[[plate_id_col]][cmb[1, ]] == d[[plate_id_col]][cmb[2, ]],
        r = NA_real_,
        stringsAsFactors = FALSE
      )

      for (k in seq_len(m)) {
        res$r[k] <- cor_row(cmb[1, k], cmb[2, k])
      }

      out_list[[g]] <- res
    }

    results <- do.call(rbind, out_list)
  }

  if (mask_sample_ids) {
    if (nrow(results) > 0) {
      results$Sample_i <- map_or_raw(results$Sample_i)
      results$Sample_j <- map_or_raw(results$Sample_j)
      results$SubjectId <- mask_subj(results$SubjectId)
    }
    groups$SubjectId <- mask_subj(groups$SubjectId)
    groups$SampleIds <- vapply(
      strsplit(groups$SampleIds, "|", fixed = TRUE),
      function(x) paste(map_or_raw(x), collapse = "|"),
      character(1)
    )
  }

  list(
    n_samples_analyzed = nrow(d),
    n_selected_rows = length(group_key),
    n_replicate_groups = length(rep_keys),
    n_ids_with_reps = length(rep_keys),
    n_replicate_rows = sum(as.integer(tab[rep_keys])),
    n_replicate_pairs = nrow(results),
    n_rows_used = nrow(d),
    groups = groups,
    results = results
  )
}