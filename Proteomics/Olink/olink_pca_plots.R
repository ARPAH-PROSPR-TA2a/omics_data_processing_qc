olink_pca_plots <- function(pca_result,
                            metadata,
                            sample_col = "SAMPLE_ID",
                            color_vars = c("PlateID", "SampleQC"),
                            pcs = c(1, 2)) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install ggplot2 first: install.packages('ggplot2')")
  }
  
  if (!is.list(pca_result) || !all(c("scores", "variance_explained") %in% names(pca_result))) {
    stop("pca_result must be output from olink_calculate_pca()")
  }
  
  if (!is.data.frame(metadata)) {
    stop("metadata must be a data.frame")
  }
  
  if (!sample_col %in% names(metadata)) {
    stop(paste0("Column '", sample_col, "' not found in metadata"))
  }
  
  scores <- pca_result$scores
  ve <- pca_result$variance_explained
  
  pcx <- pcs[1]
  pcy <- pcs[2]
  xcol <- paste0("PC", pcx)
  ycol <- paste0("PC", pcy)
  
  if (!xcol %in% names(scores)) {
    stop(paste0("PC", pcx, " not found in PCA scores. Only ", length(ve), " PCs available."))
  }
  if (!ycol %in% names(scores)) {
    stop(paste0("PC", pcy, " not found in PCA scores. Only ", length(ve), " PCs available."))
  }
  
  plot_df <- merge(scores, metadata, by.x = "SampleID", by.y = sample_col, all.x = TRUE)
  
  available_colors <- intersect(color_vars, names(plot_df))
  
  if (length(available_colors) == 0) {
    stop("None of the requested color_vars found in metadata. Available: ", 
         paste(names(plot_df), collapse = ", "))
  }
  
  plots <- setNames(vector("list", length(available_colors)), available_colors)
  
  for (cv in available_colors) {
    dfp <- plot_df
    dfp <- dfp[!is.na(dfp[[cv]]), , drop = FALSE]
    
    n_levels <- length(unique(dfp[[cv]]))
    is_numeric <- is.numeric(dfp[[cv]])
    omit_legend <- n_levels > 50 && !is_numeric
    
    p <- ggplot2::ggplot(dfp, ggplot2::aes(x = .data[[xcol]], y = .data[[ycol]], 
                                            color = .data[[cv]])) +
      ggplot2::geom_point(alpha = 0.75, size = 2.5) +
      ggplot2::xlab(paste0("PC", pcx, " (", sprintf("%.1f", 100 * ve[pcx]), "%)")) +
      ggplot2::ylab(paste0("PC", pcy, " (", sprintf("%.1f", 100 * ve[pcy]), "%)")) +
      ggplot2::labs(title = paste0("PCA colored by ", cv)) +
      ggplot2::theme_classic()
    
    if (omit_legend) {
      p <- p + ggplot2::labs(subtitle = paste0(cv)) +
        ggplot2::theme(legend.position = "none")
    } else if (!is_numeric) {
      p <- p + ggplot2::labs(color = cv)
    }
    
    plots[[cv]] <- p
  }
  
  list(
    plots = plots,
    scores = plot_df,
    variance_explained = ve
  )
}
