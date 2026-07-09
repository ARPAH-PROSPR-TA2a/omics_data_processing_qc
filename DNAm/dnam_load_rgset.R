dnam_load_rgset <- function(sample_sheet,
                            sample_id_col = "barcode",
                            extended = TRUE,
                            verbose = TRUE) {
  dnam_require_package(
    "minfi",
    "Install with: BiocManager::install(\"minfi\")"
  )

  if (!"Basename" %in% names(sample_sheet)) {
    stop("sample_sheet must contain a Basename column.", call. = FALSE)
  }
  if (!sample_id_col %in% names(sample_sheet)) {
    stop("Sample sheet is missing sample_id_col: ", sample_id_col, call. = FALSE)
  }

  rgset <- minfi::read.metharray(
    basenames = sample_sheet$Basename,
    extended = extended,
    verbose = verbose
  )

  colnames(rgset) <- make.unique(as.character(sample_sheet[[sample_id_col]]))
  rgset
}
