#!/usr/bin/env Rscript

# Main.R - Olink QC Pipeline Entry Point
# Run this script to execute the full QC pipeline

# ==============================================================================
# SETUP
# ==============================================================================

# Load functions
source("olink_sample_qc.R")
source("olink_calculate_pca.R")
source("olink_pca_plots.R")

# Load data (replace with your Olink CSV file)
# df <- read.csv("your_data.csv", stringsAsFactors = FALSE)

# For testing, use example data:
# df <- read.csv("CARDIA_OlinkHT_NPX_c2_20240617.csv", stringsAsFactors = FALSE)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Output directory
output_dir <- "output"

# Privacy settings (set to TRUE to mask sample IDs)
MASK_IDS <- TRUE

# Column names (adjust if your data uses different names)
SAMPLE_COL <- "SAMPLE_ID"
PLATE_COL <- "PlateID"
METADATA_COLS <- c("SAMPLE_ID", "Investigator_ID", "SampleType", "WellID", 
                    "PlateID", "OlinkID", "Assay", "AssayType", "SampleQC")

# PCA coloring variables
PCA_COLOR_VARS <- c("PlateID", "SampleQC")

# ==============================================================================
# CREATE OUTPUT DIRECTORIES
# ==============================================================================

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# STEP 1: Sample QC
# ==============================================================================

cat("=== Step 1: Sample QC ===\n")

qc_result <- olink_sample_qc(df,
                              sample_col = SAMPLE_COL,
                              value_col = "NPX",
                              assay_type_col = "AssayType",
                              assay_keep = "assay",
                              cutoff = 3,
                              max_prune = 0,
                              mask_sample_ids = MASK_IDS)

cat("Samples passing:", qc_result$summary$n_pass, "/", qc_result$summary$total_samples, "\n")
cat("Summary saved to:", qc_result$summary_file, "\n")

qc_pass_pct <- 100 * qc_result$summary$n_pass / qc_result$summary$total_samples
if (qc_pass_pct < 80) {
  cat("\n==================================================\n")
  cat("WARNING: LESS THAN 80% OF SAMPLES PASS QC\n")
  cat("Only ", round(qc_pass_pct, 1), "% of samples passed this QC step.\n", sep = "")
  cat("Prune a small number of low-quality samples, re-run QC, and save the filtered analysis dataset before continuing.\n")
  cat("The filtered dataset is the analysis dataset for downstream steps.\n")
  cat("==================================================\n")
}

# Save per-sample results
qc_suffix <- ifelse(MASK_IDS, "_masked", "")
write.csv(qc_result$per_sample,
          file.path(output_dir, paste0("sample_qc_per_sample", qc_suffix, ".csv")),
          row.names = FALSE)

# Extract passing sample IDs
passed_sample_ids <- qc_result$passed_samples

# ==============================================================================
# STEP 2: Filter to passing samples
# ==============================================================================

cat("\n=== Filtering to passing samples ===\n")
df_filtered <- df[df[[SAMPLE_COL]] %in% passed_sample_ids, ]
cat("Samples after filtering:", length(unique(df_filtered[[SAMPLE_COL]])), "\n")

# ==============================================================================
# STEP 3: PCA Calculation
# ==============================================================================

cat("\n=== Step 2: PCA Calculation ===\n")

pca_result <- olink_calculate_pca(df_filtered,
                                    sample_col = SAMPLE_COL,
                                    olink_id_col = "OlinkID",
                                    value_col = "NPX",
                                    assay_type_col = "AssayType",
                                    assay_keep = "assay",
                                    impute_missing = TRUE,
                                    mask_sample_ids = MASK_IDS)

cat("Proteins used:", pca_result$n_proteins, "/", pca_result$n_proteins_total, "\n")
cat("Samples used:", pca_result$n_samples, "\n")
cat("Variance explained by PC1:", round(pca_result$variance_explained[1] * 100, 1), "%\n")
cat("Variance explained by PC2:", round(pca_result$variance_explained[2] * 100, 1), "%\n")

# Save PCA outputs
pca_result_to_save <- pca_result
if (MASK_IDS && !is.null(pca_result$masked_scores)) {
  pca_result_to_save$scores <- pca_result$masked_scores
}
saveRDS(pca_result_to_save, file.path(output_dir, "pca_result.rds"))

write.csv(data.frame(
  PC = paste0("PC", seq_along(pca_result$variance_explained)),
  variance_explained = pca_result$variance_explained
), file.path(output_dir, "variance_explained.csv"), row.names = FALSE)

write.csv(pca_result$loadings,
          file.path(output_dir, paste0("pca_loadings", qc_suffix, ".csv")),
          row.names = FALSE)

# ==============================================================================
# STEP 4: PCA Plots
# ==============================================================================

cat("\n=== Step 3: PCA Plots ===\n")

# Extract metadata for plotting
metadata <- df_filtered[, intersect(METADATA_COLS, names(df_filtered))]
metadata <- metadata[!duplicated(metadata[[SAMPLE_COL]]), ]

# Generate plots
pca_plots <- olink_pca_plots(pca_result,
                              metadata,
                              sample_col = SAMPLE_COL,
                              color_vars = PCA_COLOR_VARS, 
							  show_legend = TRUE)
							  # If your batch labels are causing the legend to cover your plot, set show_legend = FALSE.

cat("Generated", length(pca_plots$plots), "plots\n")

# Save individual plots
dir.create(file.path(output_dir, "plots"), showWarnings = FALSE)
for (name in names(pca_plots$plots)) {
  ggplot2::ggsave(
    file.path(output_dir, "plots", paste0("pca_", name, ".png")),
    pca_plots$plots[[name]],
    width = 7,
    height = 5
  )
}

# Save combined scores
scores_to_save <- if (MASK_IDS && !is.null(pca_result$masked_scores)) pca_result$masked_scores else pca_result$scores
write.csv(scores_to_save,
          file.path(output_dir, paste0("pca_scores", qc_suffix, ".csv")),
          row.names = FALSE)

# ==============================================================================
# COMPLETE
# ==============================================================================

cat("\n=== Pipeline Complete ===\n")
cat("Output directory:", output_dir, "\n")
cat("Results saved to:\n")
cat("  - sample_qc_per_sample.csv\n")
cat("  - pca_result.rds\n")
cat("  - variance_explained.csv\n")
cat("  - pca_loadings.csv\n")
cat("  - pca_scores.csv\n")
cat("  - plots/\n")
