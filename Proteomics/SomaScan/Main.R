#!/usr/bin/env Rscript

# Main.R - SomaScan QC Pipeline Entry Point
# Run this script to execute the full QC pipeline

# ==============================================================================
# SETUP
# ==============================================================================

# Load functions
source("somascan_lod_qc.R")
source("somascan_norm_qc.R")
source("somascan_techrep_cor.R")
source("somascan_pca_plots.R")

# Load data (replace with your ADAT file)
# dat <- read.table("your_data.adat", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# OR for tab-delimited:
# dat <- read.delim("your_data.adat", stringsAsFactors = FALSE)

# For testing, use example data:
# install.packages(SomaDataIO)
# library(SomaDataIO)
# dat <- read_adat("example_data_v4.1_plasma.adat")

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Output directory
output_dir <- "output"

# Privacy settings (set to TRUE to mask sample IDs)
MASK_IDS <- TRUE

# Run technical replicate correlation (set to FALSE if dataset has no technical replicates)
RUN_TECHREP <- TRUE

# Timepoint column for baseline filtering (set NULL if not applicable)
TIME_COL <- "Followup"        # Column name containing timepoint info
BASELINE_VALUES <- c(0, "Baseline", "BL")  # Values representing baseline

# PCA coloring variables
PCA_COLOR_VARS <- c("PlateId", "SlideId", "Sex", "Age", "SampleGroup")

# ==============================================================================
# CREATE OUTPUT DIRECTORIES
# ==============================================================================

dir.create(file.path(output_dir, "01_lod"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "02_norm"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "03_techrep"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "04_pca", "plots"), recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# STEP 1: LOD QC
# ==============================================================================

cat("=== Step 1: Limit of Detection QC ===\n")

# Store original SampleIds BEFORE passing to lod_qc (mask may alter them)
original_sample_ids <- dat$SampleId

lod_result <- somascan_lod_qc(dat,
                              snr_thresh = 2,
                              target_protein_prop = 0.80,
                              target_sample_prop = 0.80,
                              max_prune = 0,
                              mask_sample_ids = MASK_IDS)

cat("Samples passing LOD:", round(lod_result$pct_samples_with_ge80pct_proteins * 100, 1), "%\n")

# Save outputs
lod_suffix <- ifelse(MASK_IDS, "_masked", "")
write.csv(data.frame(
  metric = c("pct_samples_passing", "n_samples_initial", "n_samples_final", "n_pruned"),
  value = c(lod_result$pct_samples_with_ge80pct_proteins,
            lod_result$n_samples_initial,
            lod_result$n_samples_final,
            lod_result$n_pruned)
), file.path(output_dir, "01_lod", paste0("lod_summary", lod_suffix, ".csv")), row.names = FALSE)

saveRDS(lod_result$lod_per_analyte, 
        file.path(output_dir, "01_lod", paste0("lod_per_analyte", lod_suffix, ".rds")))

write.csv(lod_result$per_sample,
          file.path(output_dir, "01_lod", paste0("per_sample_lod", lod_suffix, ".csv")), 
          row.names = FALSE)

# ==============================================================================
# Extract passing sample indices robustly (handles masked vs original IDs)
# ==============================================================================

# lod_result$per_sample rows correspond 1-to-1 with the input dat rows,
# so we can use positional indexing regardless of whether IDs were masked.
passing_indices <- which(lod_result$per_sample$sample_pass)

# ==============================================================================
# STEP 2: Filter to passing samples
# ==============================================================================

cat("\n=== Filtering to passing samples ===\n")

dat_filtered <- dat[passing_indices, ]

cat("Samples after filtering:", nrow(dat_filtered), "\n")

# ==============================================================================
# STEP 3: Normalization QC
# ==============================================================================

cat("\n=== Step 2: Normalization QC ===\n")

norm_result <- somascan_norm_qc(dat_filtered,
                                lower = 0.4,
                                upper = 2.5,
                                mask_sample_ids = MASK_IDS)

cat("Samples passing all normalizations:", norm_result$pass_counts$all_dilutions, "/", 
    norm_result$n_samples, "\n")

# Save outputs
write.csv(data.frame(
  dilution = c("NormScale_20", "NormScale_0_5", "NormScale_0_005", "all_dilutions"),
  pass_count = c(norm_result$pass_counts$NormScale_20,
                 norm_result$pass_counts$NormScale_0_5,
                 norm_result$pass_counts$NormScale_0_005,
                 norm_result$pass_counts$all_dilutions),
  fail_count = c(norm_result$fail_counts$NormScale_20,
                 norm_result$fail_counts$NormScale_0_5,
                 norm_result$fail_counts$NormScale_0_005,
                 norm_result$fail_counts$any_dilution)
), file.path(output_dir, "02_norm", paste0("norm_summary", lod_suffix, ".csv")), row.names = FALSE)

write.csv(norm_result$per_sample,
          file.path(output_dir, "02_norm", paste0("per_sample_norm", lod_suffix, ".csv")), 
          row.names = FALSE)

# ==============================================================================
# STEP 4: Technical Replicate Correlations
# ==============================================================================

if (RUN_TECHREP) {
  
  cat("\n=== Step 3: Technical Replicate Correlations ===\n")
  
  techrep_result <- somascan_techrep_cor(dat_filtered,
                                           use_qc = TRUE,
                                           include_only_baseline = TRUE,
                                           time_col = TIME_COL,
                                           baseline_values = BASELINE_VALUES,
                                           mask_sample_ids = MASK_IDS)
  
  cat("Sample IDs with replicates:", techrep_result$n_ids_with_reps, "\n")
  
  if (nrow(techrep_result$results) > 0) {
    cat("Mean correlation:", round(mean(techrep_result$results$r, na.rm = TRUE), 3), "\n")
    
    write.csv(data.frame(
      n_rows = techrep_result$n_rows_used,
      n_ids_with_reps = techrep_result$n_ids_with_reps,
      mean_r = mean(techrep_result$results$r, na.rm = TRUE),
      median_r = median(techrep_result$results$r, na.rm = TRUE)
    ), file.path(output_dir, "03_techrep", "techrep_summary.csv"), row.names = FALSE)
    
    write.csv(techrep_result$results,
              file.path(output_dir, "03_techrep", paste0("pairwise_correlations", lod_suffix, ".csv")),
              row.names = FALSE)
  }
  
} else {
  cat("\n=== Step 3: Technical Replicate Correlations ===\n")
  cat("Skipped (RUN_TECHREP = FALSE)\n")
}

# ==============================================================================
# STEP 5: PCA Plots
# ==============================================================================

cat("\n=== Step 4: PCA Plots ===\n")

pca_result <- somascan_pca_plots(dat_filtered,
                                  color_vars = PCA_COLOR_VARS,
                                  mask_sample_ids = MASK_IDS)

cat("Variance explained by PC1:", round(pca_result$variance_explained[1] * 100, 1), "%\n")
cat("Variance explained by PC2:", round(pca_result$variance_explained[2] * 100, 1), "%\n")

# Save outputs
write.csv(pca_result$scores,
          file.path(output_dir, "04_pca", paste0("pca_scores", lod_suffix, ".csv")),
          row.names = FALSE)

write.csv(data.frame(
  PC = paste0("PC", 1:length(pca_result$variance_explained)),
  variance_explained = pca_result$variance_explained
), file.path(output_dir, "04_pca", "variance_explained.csv"), row.names = FALSE)

# Save individual plots
for (name in names(pca_result$plots)) {
  ggsave(file.path(output_dir, "04_pca", "plots", paste0("pca_", name, ".png")),
         pca_result$plots[[name]], width = 7, height = 5)
}

# ==============================================================================
# COMPLETE
# ==============================================================================

cat("\n=== Pipeline Complete ===\n")
cat("Output directory:", output_dir, "\n")
cat("Results saved to:\n")
cat("  - 01_lod/\n")
cat("  - 02_norm/\n")
cat("  - 03_techrep/\n")
cat("  - 04_pca/\n")
