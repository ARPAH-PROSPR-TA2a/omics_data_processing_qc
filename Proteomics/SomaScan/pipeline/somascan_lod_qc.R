somascan_lod_qc <- function(dat,
                            snr_thresh = 2,
                            target_protein_prop = 0.80,
                            target_sample_prop = 0.80,
                            max_prune = 0,
                            lod_k = 3,
                            sample_type_col = "SampleType",
                            sample_id_col = "SampleId",
                            sample_label = "Sample",
                            buffer_label = "Buffer",
                            seq_prefix = "seq.",
                            mask_sample_ids = FALSE,
                            verbose = FALSE) {
  
  stopifnot(is.data.frame(dat))
  if (!all(c(sample_type_col, sample_id_col) %in% names(dat))) {
    stop("Missing required columns: ", sample_type_col, " and/or ", sample_id_col)
  }
  
  # Identify analyte columns
  analyte_cols <- grep(paste0("^", gsub("\\.", "\\\\.", seq_prefix)), names(dat), value = TRUE)
  if (length(analyte_cols) == 0) stop("No analyte columns found (expected columns starting with '", seq_prefix, "').")
  
  # Split
  buffers <- dat[dat[[sample_type_col]] == buffer_label, , drop = FALSE]
  samples <- dat[dat[[sample_type_col]] == sample_label, , drop = FALSE]
  
  if (nrow(buffers) == 0) stop("No Buffer rows found (", sample_type_col, " == '", buffer_label, "').")
  if (nrow(samples) == 0) stop("No Sample rows found (", sample_type_col, " == '", sample_label, "').")
  
  # Compute LOD per analyte from buffers
  buffer_mat <- as.matrix(buffers[, analyte_cols, drop = FALSE])
  storage.mode(buffer_mat) <- "double"
  
    lod <- apply(buffer_mat, 2, function(x) {
      med <- stats::median(x, na.rm = TRUE)
      md  <- stats::mad(x, constant = 1, na.rm = TRUE) # MAD (not scaled)
      med + lod_k * md
    })
  
  # Guard: avoid division by 0 or negative LOD
  lod[!is.finite(lod) | lod <= 0] <- NA_real_
  
  # Helper to score samples
  score_samples <- function(samples_df) {
    smat <- as.matrix(samples_df[, analyte_cols, drop = FALSE])
    storage.mode(smat) <- "double"
    
    # SNR per sample x analyte
    snr <- sweep(smat, 2, lod, "/")
    
    # For each sample: proportion of analytes with SNR > threshold
    prot_pass_prop <- rowMeans(snr > snr_thresh, na.rm = TRUE)
    
    # Sample passes if ≥ target_protein_prop
    sample_pass <- prot_pass_prop >= target_protein_prop
    
    list(
      prot_pass_prop = prot_pass_prop,
      sample_pass = sample_pass,
      pct_samples_passing = mean(sample_pass, na.rm = TRUE)
    )
  }
  
  removed <- character(0)
  iter <- 0
  
  res <- score_samples(samples)
  
  # Prune worst samples until target achieved or max_prune reached
  while (iter < max_prune && is.finite(res$pct_samples_passing) &&
         res$pct_samples_passing < target_sample_prop && nrow(samples) > 1) {
    iter <- iter + 1
    
    # Worst = smallest proportion of proteins passing
    worst_idx <- which.min(res$prot_pass_prop)
    worst_id  <- as.character(samples[[sample_id_col]][worst_idx])
    
    removed <- c(removed, worst_id)
    samples <- samples[-worst_idx, , drop = FALSE]
    
    if (verbose) {
      message("Pruned: ", worst_id, " (prot_pass_prop=", round(res$prot_pass_prop[worst_idx], 3), ")")
    }
    
    res <- score_samples(samples)
  }
  
  # Print message if samples were pruned
  if (length(removed) > 0) {
    message("\n=== LOD QC: Samples Removed ===")
    message("Number of samples removed: ", length(removed))
    message("To use only passing samples in subsequent steps, run:")
    message("  passed_sample_ids <- lod_result$per_sample$SampleId[lod_result$per_sample$sample_pass]")
    message("================================\n")
  }
  
  # Final per-sample table
  out_samples <- data.frame(
    SampleId = as.character(samples[[sample_id_col]]),
    proteins_pass_prop = res$prot_pass_prop,
    sample_pass = res$sample_pass,
    stringsAsFactors = FALSE
  )
  
  # Apply masking if requested
  if (mask_sample_ids) {
    unique_ids <- unique(out_samples$SampleId)
    id_map <- setNames(paste0("Sample_", seq_along(unique_ids)), unique_ids)
    out_samples$SampleId <- id_map[out_samples$SampleId]
  }
  
  list(
    pct_samples_with_ge80pct_proteins = res$pct_samples_passing,
    target_sample_prop = target_sample_prop,
    target_protein_prop = target_protein_prop,
    snr_thresh = snr_thresh,
    lod_k = lod_k,
    n_samples_initial = sum(dat[[sample_type_col]] == sample_label, na.rm = TRUE),
    n_samples_final = nrow(samples),
    n_pruned = length(removed),
    pruned_sample_ids = removed,
    per_sample = out_samples,
    lod_per_analyte = lod
  )
}
