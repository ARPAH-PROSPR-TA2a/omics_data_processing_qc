#!/usr/bin/env Rscript
library(tidyverse)
# Main.R - Metabolon QC Pipeline Entry Point

source("Code/Metabolon/pipeline/metabolon_sample_qc.R")
source("Code/Metabolon/pipeline/metabolon_calculate_pca.R")
source("Code/Metabolon/pipeline/metabolon_pca_plots.R")

# Replace these with your real objects
obt <-read_csv("Data/Metabolon RScript for Duke analyses/OBT.csv")
phenotype_df <- obt %>% select(PARENT_SAMPLE_NAME:PARAM_VOLUME_EXTRACTED)
metabolite_mat <- obt %>% select(-c(PARENT_SAMPLE_NAME:PARAM_VOLUME_EXTRACTED))

rownames(metabolite_mat) <-obt$PARENT_SAMPLE_ID

dim(phenotype_df)
dim(metabolite_mat)

output_dir <- "Output/QC_test"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)

processed_metabolomics_data_dir_name <- "processed_metabolomics_data"
processed_metabolomics_data_dir <- file.path(dirname(output_dir), processed_metabolomics_data_dir_name)
dir.create(processed_metabolomics_data_dir, recursive = TRUE, showWarnings = FALSE)

MASK_IDS <- TRUE

# The full PCA RDS can contain participant-specific structure. Save it only when needed,
# and keep it outside the standard shareable output directory.
SAVE_FULL_PCA_OBJECT <- TRUE
full_pca_output_dir_name <- "restricted_pca_output"
full_pca_output_dir <- file.path(dirname(output_dir), full_pca_output_dir_name)

# Column in phenotype_df that uniquely identifies each sample row
SAMPLE_ID_COL <- "PARENT_SAMPLE_ID"

qc_result <- metabolon_sample_qc(
  phenotype_df = phenotype_df,
  metabolite_df = metabolite_mat,
  sample_id_col = SAMPLE_ID_COL,
  cutoff = 3,
  prune = 0,
  target_sample_prop = 0.80,
  mask_sample_ids = MASK_IDS
)

write.csv(qc_result$summary, file.path(output_dir, "metabolon_sample_qc_summary.csv"), row.names = FALSE)
write.csv(qc_result$per_sample, file.path(output_dir, paste0("metabolon_per_sample", ifelse(MASK_IDS, "_masked", ""), ".csv")), row.names = FALSE)

if (qc_result$summary$pct_pass < 80) {
  cat("\n==================================================\n")
  cat("WARNING: LESS THAN 80% OF SAMPLES PASS QC\n")
  cat("Only ", round(qc_result$summary$pct_pass, 1), "% of samples passed this QC step.\n", sep = "")
  cat("Re-run with prune > 0 to remove the worst samples and save the filtered analysis dataset before downstream steps.\n")
  cat("==================================================\n")
}

analysis_mat <- qc_result$analysis_matrix
analysis_pheno <- qc_result$analysis_pheno

write.csv(analysis_pheno, file.path(processed_metabolomics_data_dir, "metabolon_analysis_pheno.csv"), row.names = FALSE)
saveRDS(analysis_mat, file.path(processed_metabolomics_data_dir, "metabolon_analysis_matrix.rds"))

# STEP 2: Calculate PCA and save outputs
pca_result <- metabolon_calculate_pca(
  analysis_mat,
  log2_transform = FALSE,
  mask_sample_ids = MASK_IDS
)


pca_summary <- data.frame(
  n_metabolites = pca_result$n_metabolites,
  n_metabolites_total = pca_result$n_metabolites_total,
  n_samples = pca_result$n_samples,
  metabolites_removed = pca_result$metabolites_removed,
  missingness_threshold_used = pca_result$missingness_threshold_used,
  stringsAsFactors = FALSE
)
write.csv(pca_summary, file.path(output_dir, "metabolon_pca_summary.csv"), row.names = FALSE)


pca_result_to_save <- pca_result
if (MASK_IDS && !is.null(pca_result$masked_scores)) {
  pca_result_to_save$scores <- pca_result$masked_scores
}

if (SAVE_FULL_PCA_OBJECT) {
  dir.create(full_pca_output_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(pca_result_to_save, file.path(full_pca_output_dir, "metabolon_pca_result.rds"))
}
write.csv(
  data.frame(
    PC = paste0("PC", seq_along(pca_result$variance_explained)),
    variance_explained = pca_result$variance_explained
  ),
  file.path(output_dir, "metabolon_variance_explained.csv"),
  row.names = FALSE
)
write.csv(
  pca_result$loadings,
  file.path(output_dir, paste0("metabolon_pca_loadings", ifelse(MASK_IDS, "_masked", ""), ".csv")),
  row.names = FALSE
)

scores_to_save <- if (MASK_IDS && !is.null(pca_result$masked_scores)) pca_result$masked_scores else pca_result$scores
write.csv(
  scores_to_save,
  file.path(output_dir, paste0("metabolon_pca_scores", ifelse(MASK_IDS, "_masked", ""), ".csv")),
  row.names = FALSE
)

# STEP 3: Generate PCA plots and save
pca_plots <- metabolon_pca_plots(
  pca_result = pca_result,
  phenotype_df = analysis_pheno,
  sample_id_col = SAMPLE_ID_COL,
  color_vars = c("PARAM_COMMENT", "PARAM_TP", "PARAM_GENDER", "PARAM_AGE", "PARAM_BMI")
)

for (name in names(pca_plots$plots)) {
  ggplot2::ggsave(
    file.path(output_dir, "plots", paste0("metabolon_pca_", name, ".png")),
    pca_plots$plots[[name]],
    width = 7,
    height = 5
  )
}

cat("Metabolon pipeline complete\n")
