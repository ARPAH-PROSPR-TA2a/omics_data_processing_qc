save_plate_control_qc_outputs <- function(qc_res, study_name, base_dir = ".") {
  
  study_clean <- gsub("[^A-Za-z0-9_]", "_", study_name)
  
  out_dir <- file.path(base_dir, paste0("QC_", study_clean))
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  counts_path <- file.path(out_dir, paste0(study_clean, "_plate_control_qc_counts.csv"))
  summary_path <- file.path(out_dir, paste0(study_clean, "_plate_control_qc_summary.csv"))
  failed_npx_path <- file.path(out_dir, paste0(study_clean, "_failed_plate_npx_check.csv"))
  failed_plates_path <- file.path(out_dir, paste0(study_clean, "_failed_plates.txt"))
  
  write.csv(qc_res$qc_counts, counts_path, row.names = FALSE)
  write.csv(qc_res$qc_summary, summary_path, row.names = FALSE)
  write.csv(qc_res$failed_plate_npx_check, failed_npx_path, row.names = FALSE)
  
  if (length(qc_res$failed_plates) == 1 && grepl("^No plates failed$", qc_res$failed_plates)) {
    writeLines("No plates failed", failed_plates_path)
  } else {
    writeLines(qc_res$failed_plates, failed_plates_path)
  }
  
  return(list(
    output_dir = out_dir,
    qc_counts_file = counts_path,
    qc_summary_file = summary_path,
    failed_plate_npx_check_file = failed_npx_path,
    failed_plates_file = failed_plates_path
  ))
}
