metabolon_pca_plots <- function(pca_result,
                                phenotype_df,
                                sample_id_col = "SubjectID",
                                color_vars = c("PlateID", "Batch", "Sex", "Age"),
                                pcs = c(1, 2)) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2 first: install.packages('ggplot2')")
  if (!is.list(pca_result) || !all(c("scores", "variance_explained") %in% names(pca_result))) {
    stop("pca_result must come from metabolon_calculate_pca()")
  }
  if (!is.data.frame(phenotype_df)) stop("phenotype_df must be a data.frame")
  if (!sample_id_col %in% names(phenotype_df)) stop("Column '", sample_id_col, "' not found in phenotype_df")

  scores <- pca_result$scores
  ve <- pca_result$variance_explained
  pcx <- pcs[1]; pcy <- pcs[2]
  xcol <- paste0("PC", pcx)
  ycol <- paste0("PC", pcy)
  if (!xcol %in% names(scores)) stop("", xcol, " not found in PCA scores")
  if (!ycol %in% names(scores)) stop("", ycol, " not found in PCA scores")

  plot_df <- merge(scores, phenotype_df, by.x = "SampleID", by.y = sample_id_col, all.x = TRUE)
  available_colors <- intersect(color_vars, names(plot_df))
  if (length(available_colors) == 0) stop("None of the requested color_vars found in phenotype data")

  plots <- setNames(vector("list", length(available_colors)), available_colors)
  for (cv in available_colors) {
    dfp <- plot_df[!is.na(plot_df[[cv]]), , drop = FALSE]
    n_levels <- length(unique(dfp[[cv]]))
    is_numeric <- is.numeric(dfp[[cv]])
    omit_legend <- !is_numeric && n_levels > 50

    p <- ggplot2::ggplot(dfp, ggplot2::aes(x = .data[[xcol]], y = .data[[ycol]], color = .data[[cv]])) +
      ggplot2::geom_point(alpha = 0.75, size = 2.5) +
      ggplot2::xlab(paste0("PC", pcx, " (", sprintf("%.1f", 100 * ve[pcx]), "%)")) +
      ggplot2::ylab(paste0("PC", pcy, " (", sprintf("%.1f", 100 * ve[pcy]), "%)")) +
      ggplot2::labs(title = paste0("PCA colored by ", cv)) +
      ggplot2::theme_classic()

    if (omit_legend) p <- p + ggplot2::theme(legend.position = "none")
    else if (!is_numeric) p <- p + ggplot2::labs(color = cv)

    plots[[cv]] <- p
  }

  list(plots = plots, scores = plot_df, variance_explained = ve)
}
