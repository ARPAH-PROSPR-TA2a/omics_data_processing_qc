save_sample_qc_outputs <- function(qc_res, study_name, base_dir = ".") {
  
  # clean study name (avoid spaces / weird chars)
  study_clean <- gsub("[^A-Za-z0-9_]", "_", study_name)
  
  # create folder
  out_dir <- file.path(base_dir, paste0("QC_", study_clean))
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # file paths
  counts_path  <- file.path(out_dir, paste0(study_clean, "_qc_counts.csv"))
  summary_path <- file.path(out_dir, paste0(study_clean, "_qc_summary.csv"))
  
  # save
  write.csv(qc_res$qc_counts, counts_path, row.names = FALSE)
  write.csv(qc_res$qc_summary, summary_path, row.names = FALSE)
  
  # return paths (useful)
  return(list(
    output_dir = out_dir,
    qc_counts_file = counts_path,
    qc_summary_file = summary_path
  ))
}