#!/usr/bin/env Rscript

# Main.R - Olink QC Pipeline Entry Point
# Run this script to execute the QC pipeline

# ==============================================================================
# SETUP
# ==============================================================================

# Load all functions
source("olink_sample_qc.R")
source("plate_control_qc.R")
source("save_sample_qc_outputs.R")
source("save_plate_control_qc_outputs.R")

# Load data (replace with your Olink CSV file)
# df <- read.csv("your_data.csv", stringsAsFactors = FALSE)

# For testing, use example data:
# df <- read.csv("CARDIA_OlinkHT_NPX_c2_20240617.csv", stringsAsFactors = FALSE)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Study name for output files
STUDY_NAME <- "MyStudy"

# Privacy settings (set to TRUE to mask sample IDs)
MASK_SAMPLE_IDS <- FALSE

# ==============================================================================
# STEP 1: Sample QC
# ==============================================================================

cat("=== Step 1: Sample QC ===\n")

sample_qc_result <- olink_sample_qc(df,
                                     sample_id_col = "SAMPLE_ID",
                                     sample_type_col = "SampleType",
                                     sample_type_keep = "SAMPLE",
                                     assay_type_col = "AssayType",
                                     count_col = "Count",
                                     mask_sample_id = MASK_SAMPLE_IDS)

cat("Total samples:", sample_qc_result$qc_summary$n_samples, "\n")
cat("Samples passing all QC:", sample_qc_result$qc_summary$n_samples - sample_qc_result$qc_summary$any_fail, "\n")

# Save outputs
sample_output <- save_sample_qc_outputs(sample_qc_result, STUDY_NAME)
cat("Output saved to:", sample_output$output_dir, "\n")

# ==============================================================================
# STEP 2: Plate Control QC
# ==============================================================================

cat("\n=== Step 2: Plate Control QC ===\n")

plate_qc_result <- plate_control_qc(df,
                                    sample_type_col = "SampleType",
                                    sample_type_keep = "PLATE_CONTROL",
                                    block_col = "Block",
                                    plate_col = "PlateID",
                                    replicate_col = "Investigator_ID",
                                    assay_type_col = "AssayType",
                                    count_col = "Count",
                                    npx_col = "NPX")

cat("Plates assessed:", nrow(plate_qc_result$qc_summary), "\n")
cat("Plates passing:", sum(plate_qc_result$qc_summary$plate_pass), "\n")
cat("Plates failing:", sum(!plate_qc_result$qc_summary$plate_pass), "\n")

# Save outputs
plate_output <- save_plate_control_qc_outputs(plate_qc_result, STUDY_NAME)
cat("Output saved to:", plate_output$output_dir, "\n")

# ==============================================================================
# COMPLETE
# ==============================================================================

cat("\n=== Pipeline Complete ===\n")
cat("Output directory:", paste0("QC_", STUDY_NAME), "\n")
cat("Results saved to:\n")
cat("  - QC_", STUDY_NAME, "/", STUDY_NAME, "_qc_counts.csv\n", sep = "")
cat("  - QC_", STUDY_NAME, "/", STUDY_NAME, "_qc_summary.csv\n", sep = "")
cat("  - QC_", STUDY_NAME, "/", STUDY_NAME, "_plate_control_qc_counts.csv\n", sep = "")
cat("  - QC_", STUDY_NAME, "/", STUDY_NAME, "_plate_control_qc_summary.csv\n", sep = "")
cat("  - QC_", STUDY_NAME, "/", STUDY_NAME, "_failed_plate_npx_check.csv\n", sep = "")
cat("  - QC_", STUDY_NAME, "/", STUDY_NAME, "_failed_plates.txt\n", sep = "")
