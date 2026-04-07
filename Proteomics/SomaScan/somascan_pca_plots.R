somascan_pca_plots <- function(dat,
                               sample_type_col = "SampleType",
                               sample_label = "Sample",
                               sample_id_col = "SampleId",
                               seq_prefix = "seq.",
                               id_cols = c("SampleId", "PlateId", "SlideId"),
                               color_vars = c("PlateId", "SlideId", "PlatePosition", "Subarray",
                                              "ScannerID", "Sex", "Age", "SampleGroup"),
                               pcs = c(1, 2),
                               log2_transform = TRUE,
                               center = TRUE,
                               scale. = FALSE,
                               drop_na_color = TRUE,
                               mask_sample_ids = FALSE) {
  stopifnot(is.data.frame(dat))
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install ggplot2 first: install.packages('ggplot2')")
  }
  if (!requireNamespace("rlang", quietly = TRUE)) {
    stop("Install rlang first: install.packages('rlang')")
  }
  
  # subset to biological samples
  if (sample_type_col %in% names(dat)) {
    dat <- dat[dat[[sample_type_col]] == sample_label, , drop = FALSE]
  }
  if (nrow(dat) < 3) stop("Not enough rows after filtering to samples for PCA.")
  
  analyte_cols <- grep(paste0("^", seq_prefix), names(dat), value = TRUE)
  if (length(analyte_cols) == 0) stop("No analyte columns found (expected columns starting with '", seq_prefix, "').")
  
  X <- as.matrix(dat[, analyte_cols, drop = FALSE])
  storage.mode(X) <- "double"
  if (log2_transform) X <- log2(X)
  
  v <- apply(X, 2, stats::var, na.rm = TRUE)
  keep <- is.finite(v) & v > 0
  X <- X[, keep, drop = FALSE]
  if (ncol(X) < 3) stop("Too few variable analytes for PCA after filtering.")
  
  pca <- stats::prcomp(X, center = center, scale. = scale.)
  
  ve <- (pca$sdev^2) / sum(pca$sdev^2)
  ve_lbl <- function(k) paste0("PC", k, " (", sprintf("%.1f", 100 * ve[k]), "%)")
  
  pcx <- pcs[1]; pcy <- pcs[2]
  xcol <- paste0("PC", pcx)
  ycol <- paste0("PC", pcy)
  
  scores <- as.data.frame(pca$x)
  scores <- scores[, c(xcol, ycol), drop = FALSE]
  
  keep_meta <- intersect(unique(c(id_cols, color_vars)), names(dat))
  meta <- dat[, keep_meta, drop = FALSE]
  
  if (mask_sample_ids && sample_id_col %in% names(meta)) {
    unique_ids <- unique(meta[[sample_id_col]])
    id_map <- setNames(paste0("Sample_", seq_along(unique_ids)), unique_ids)
    meta[[sample_id_col]] <- id_map[as.character(meta[[sample_id_col]])]
  }
  
  plot_df <- cbind(meta, scores)
  
  used_color_vars <- intersect(color_vars, names(plot_df))
  
  plots <- setNames(vector("list", length(used_color_vars)), used_color_vars)
  
  for (cv in used_color_vars) {
    dfp <- plot_df
    if (drop_na_color) dfp <- dfp[!is.na(dfp[[cv]]), , drop = FALSE]
    
    n_levels <- length(unique(dfp[[cv]]))
    omit_legend <- n_levels > 50
    
    p <- ggplot2::ggplot(dfp, ggplot2::aes(x = .data[[xcol]], y = .data[[ycol]], color = .data[[cv]])) +
      ggplot2::geom_point(alpha = 0.75, size = 2) +
      ggplot2::labs(
        title = paste0("PCA colored by ", cv),
        x = ve_lbl(pcx),
        y = ve_lbl(pcy)
      ) +
      ggplot2::theme_classic()
    
    if (omit_legend) {
      p <- p + ggplot2::labs(subtitle = paste0(cv)) +
        ggplot2::theme(legend.position = "none")
    } else {
      p <- p + ggplot2::labs(color = cv)
    }
    
    plots[[cv]] <- p
  }
  
  list(
    pca = pca,
    variance_explained = ve,
    scores = plot_df,
    plots = plots
  )
}
