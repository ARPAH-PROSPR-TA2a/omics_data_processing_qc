dnam_require_package <- function(package, install_hint = NULL) {
  if (!requireNamespace(package, quietly = TRUE)) {
    msg <- paste0("Required package is not installed: ", package)
    if (!is.null(install_hint)) {
      msg <- paste0(msg, "\n", install_hint)
    }
    stop(msg, call. = FALSE)
  }
}

dnam_check_path_configured <- function(path, label) {
  if (is.null(path) || !nzchar(path) || grepl("^path/to/", path)) {
    stop("Configure ", label, " in Main.R before running this step.", call. = FALSE)
  }
}

dnam_check_file <- function(path, label) {
  dnam_check_path_configured(path, label)
  if (!file.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }
}

dnam_check_dir <- function(path, label) {
  dnam_check_path_configured(path, label)
  if (!dir.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }
}

dnam_read_sample_sheet <- function(sample_sheet_file) {
  dnam_check_file(sample_sheet_file, "sample_sheet_file")

  lines <- readLines(sample_sheet_file, warn = FALSE)
  data_line <- grep("^\\[Data\\]", trimws(lines))

  if (length(data_line)) {
    csv_text <- paste(lines[(data_line[[1]] + 1):length(lines)], collapse = "\n")
    sample_sheet <- utils::read.csv(
      text = csv_text,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      colClasses = "character",
      na.strings = c("", "NA")
    )
  } else {
    sample_sheet <- utils::read.csv(
      sample_sheet_file,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      colClasses = "character",
      na.strings = c("", "NA")
    )
  }

  sample_sheet <- sample_sheet[rowSums(!is.na(sample_sheet)) > 0, , drop = FALSE]
  rownames(sample_sheet) <- NULL
  sample_sheet
}

dnam_add_idat_basenames <- function(sample_sheet,
                                    idat_dir,
                                    basename_col = NULL,
                                    sentrix_id_col = "Sentrix_ID",
                                    sentrix_position_col = "Sentrix_Position") {
  dnam_check_dir(idat_dir, "idat_dir")

  if (!is.null(basename_col) && nzchar(basename_col)) {
    if (!basename_col %in% names(sample_sheet)) {
      stop("Sample sheet is missing basename_col: ", basename_col, call. = FALSE)
    }
    sample_sheet$Basename <- file.path(idat_dir, sample_sheet[[basename_col]])
    return(sample_sheet)
  }

  required_cols <- c(sentrix_id_col, sentrix_position_col)
  missing_cols <- setdiff(required_cols, names(sample_sheet))
  if (length(missing_cols)) {
    stop("Sample sheet is missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  sample_sheet$Basename <- file.path(
    idat_dir,
    paste0(sample_sheet[[sentrix_id_col]], "_", sample_sheet[[sentrix_position_col]])
  )
  sample_sheet
}

dnam_validate_idat_pairs <- function(sample_sheet) {
  if (!"Basename" %in% names(sample_sheet)) {
    stop("Sample sheet must contain a Basename column. Run dnam_add_idat_basenames() first.", call. = FALSE)
  }

  sample_sheet$Red_IDAT <- paste0(sample_sheet$Basename, "_Red.idat")
  sample_sheet$Green_IDAT <- paste0(sample_sheet$Basename, "_Grn.idat")
  sample_sheet$red_idat_exists <- file.exists(sample_sheet$Red_IDAT)
  sample_sheet$green_idat_exists <- file.exists(sample_sheet$Green_IDAT)
  sample_sheet$idat_pair_exists <- sample_sheet$red_idat_exists & sample_sheet$green_idat_exists
  sample_sheet
}

dnam_read_beta_matrix <- function(beta_file) {
  dnam_check_file(beta_file, "beta_file")

  if (grepl("\\.rds$", beta_file, ignore.case = TRUE)) {
    beta_mat <- readRDS(beta_file)
  } else {
    beta_mat <- utils::read.csv(beta_file, row.names = 1, check.names = FALSE)
  }

  beta_mat <- as.matrix(beta_mat)
  storage.mode(beta_mat) <- "double"
  beta_mat
}

dnam_mask_columns <- function(dat, columns, prefix = "Sample_") {
  columns <- intersect(columns, names(dat))
  if (!length(columns)) {
    return(dat)
  }

  ids <- unique(unlist(dat[columns], use.names = FALSE))
  ids <- ids[!is.na(ids)]
  id_map <- setNames(paste0(prefix, seq_along(ids)), ids)

  for (col in columns) {
    dat[[col]] <- unname(id_map[as.character(dat[[col]])])
  }

  dat
}
