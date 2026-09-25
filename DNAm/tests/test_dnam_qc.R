#!/usr/bin/env Rscript

# Privacy regression tests for the DNAm QC pipeline.
#
# Every collaborator-facing output must contain barcodes and QC metrics only.
# Sample-level metadata (Sample_Name, Participant_ID, Time_Point, plate, well,
# absolute IDAT paths, replicate group) must be absent unless the matching
# opt-in flag is set, and all opt-in flags must default to FALSE.
#
# Usage:
#   Rscript tests/test_dnam_qc.R --data_dir=/path/to/DNAm_QC/Data
#   Rscript tests/test_dnam_qc.R --data_dir=... --scratch_dir=/path/to/DNAm_QC
#
# DNAM_TEST_DATA_DIR is used when --data_dir is not passed.
# DNAM_TEST_BETA_ROWS caps the beta rows used for speed (0 uses the full matrix).
#
# Required in --data_dir: a sample sheet CSV, IDATs/, and a processed beta
# matrix .rds.

`%||%` <- function(x, y) if (is.null(x)) y else x

args <- list()
for (arg in commandArgs(trailingOnly = TRUE)) {
  if (!grepl("^--", arg)) next
  parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
  args[[parts[[1]]]] <- if (length(parts) > 1) paste(parts[-1], collapse = "=") else TRUE
}

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
tests_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
} else {
  getwd()
}
pipeline_dir <- normalizePath(file.path(tests_dir, ".."))

data_dir <- if (!is.null(args$data_dir)) args$data_dir else Sys.getenv("DNAM_TEST_DATA_DIR", "")
if (!nzchar(data_dir) || !dir.exists(data_dir)) {
  stop("Pass --data_dir=<DNAm_QC/Data> or set DNAM_TEST_DATA_DIR.", call. = FALSE)
}
scratch_dir <- if (!is.null(args$scratch_dir)) args$scratch_dir else data_dir
scratch_code_dir <- file.path(scratch_dir, "Code")

sample_sheet_candidates <- list.files(
  data_dir,
  pattern = "\\.csv$",
  full.names = TRUE,
  ignore.case = TRUE
)
sample_sheet_candidates <- sample_sheet_candidates[!grepl("idat", sample_sheet_candidates, ignore.case = TRUE)]
sample_sheet_file <- sample_sheet_candidates[[1]]
beta_file <- list.files(data_dir, pattern = "\\.rds$", full.names = TRUE, ignore.case = TRUE)[[1]]
idat_dir <- file.path(data_dir, "IDATs")
if (is.na(sample_sheet_file)) stop("No sample sheet CSV found in ", data_dir, call. = FALSE)
if (is.na(beta_file)) stop("No processed beta matrix .rds found in ", data_dir, call. = FALSE)
if (!dir.exists(idat_dir)) stop("IDAT directory not found: ", idat_dir, call. = FALSE)

# Columns that must never reach a shared output without an explicit opt-in.
forbidden_columns <- c(
  "Sample_Name", "Sample_Group", "Sample_Well", "Sample_Plate", "Pool_ID",
  "Participant_ID", "Time_Point", "REP", "Sentrix_ID", "Sentrix_Position",
  "Array_Batch", "Basename", "Red_IDAT", "Green_IDAT", "replicate_group"
)

results <- list()
failures <- character()

record <- function(name, passed, detail = "") {
  results[[name]] <<- list(passed = passed, detail = detail)
  if (!passed) failures <<- c(failures, name)
  cat(if (passed) "PASS  " else "FAIL  ", name, if (nzchar(detail)) paste0(" - ", detail) else "", "\n", sep = "")
}

run_script <- function(script, args = character(), env = character()) {
  out <- suppressWarnings(system2(
    "Rscript",
    c(shQuote(script), args),
    stdout = TRUE,
    stderr = TRUE,
    env = env
  ))
  list(status = attr(out, "status") %||% 0L, output = out)
}

read_header <- function(path) {
  names(utils::read.csv(path, nrows = 0, check.names = FALSE))
}

assert_no_forbidden_columns <- function(dir, label) {
  csvs <- list.files(dir, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  if (!length(csvs)) {
    record(paste0(label, ": outputs exist"), FALSE, "no CSV files written")
    return(invisible(NULL))
  }
  record(paste0(label, ": outputs exist"), TRUE, paste(length(csvs), "CSV files"))

  leaked <- character()
  for (csv in csvs) {
    hits <- intersect(read_header(csv), forbidden_columns)
    if (length(hits)) {
      leaked <- c(leaked, paste0(basename(csv), " -> ", paste(hits, collapse = ",")))
    }
  }
  record(
    paste0(label, ": no sample-level metadata columns"),
    !length(leaked),
    if (length(leaked)) paste(leaked, collapse = " | ") else ""
  )
}

assert_absent <- function(paths, label) {
  present <- paths[file.exists(paths)]
  record(label, !length(present), if (length(present)) paste(basename(present), collapse = ", ") else "")
}

assert_present <- function(path, label) {
  record(label, file.exists(path), if (!file.exists(path)) "file not written" else "")
}

assert_header <- function(path, expected, label) {
  if (!file.exists(path)) {
    record(label, FALSE, "file not written")
    return(invisible(NULL))
  }
  actual <- read_header(path)
  record(label, identical(actual, expected), paste(actual, collapse = ","))
}

assert_rows <- function(path, label, min_rows = 1L) {
  if (!file.exists(path)) {
    record(label, FALSE, "file not written")
    return(invisible(NULL))
  }
  n <- nrow(utils::read.csv(path, check.names = FALSE))
  record(label, n >= min_rows, paste(n, "rows"))
}

# A small beta matrix keeps the uniformity dip test fast while still covering
# every output path. Set DNAM_TEST_BETA_ROWS=0 to use the full matrix.
beta_rows <- as.integer(Sys.getenv("DNAM_TEST_BETA_ROWS", "2000"))
beta_matrix <- readRDS(beta_file)
beta_matrix <- as.matrix(beta_matrix)
if (beta_rows > 0 && beta_rows < nrow(beta_matrix)) {
  beta_matrix <- beta_matrix[seq_len(beta_rows), , drop = FALSE]
}
test_beta_file <- file.path(tempdir(), "test_beta_matrix.rds")
saveRDS(beta_matrix, test_beta_file)

# ==============================================================================
# 1. CALERIE raw IDAT QC, default flags
# ==============================================================================

raw_out_default <- file.path(tempdir(), "raw_default")
unlink(raw_out_default, recursive = TRUE)
raw_run <- run_script(
  file.path(pipeline_dir, "CALERIE", "CALERIE_raw_idat_qc.R"),
  env = c(
    paste0("CALERIE_SAMPLE_SHEET=", sample_sheet_file),
    paste0("CALERIE_IDAT_DIR=", idat_dir),
    paste0("CALERIE_OUTPUT_DIR=", raw_out_default),
    "CALERIE_CHUNK_SIZE=608"
  )
)
record("raw QC runs with default flags", raw_run$status == 0,
       if (raw_run$status != 0) paste(tail(raw_run$output, 3), collapse = " | ") else "")

assert_no_forbidden_columns(raw_out_default, "raw default")
assert_absent(
  c(
    file.path(raw_out_default, "01_manifest_rgset", "raw_sample_manifest_validated.csv"),
    file.path(raw_out_default, "01_manifest_rgset", "missing_idats.csv"),
    file.path(raw_out_default, "restricted_objects", "control_probe_pca.rds")
  ),
  "raw default: manifest, missing_idats and PCA object not written"
)
assert_absent(
  list.files(file.path(raw_out_default, "restricted_objects"), full.names = TRUE),
  "raw default: restricted_objects is empty"
)
assert_header(
  file.path(raw_out_default, "03_bead_qc", "bead_qc_per_sample.csv"),
  c("Barcode", "mean_bead_count", "bead_mean_threshold", "passed_mean_bead_count"),
  "raw default: bead_qc_per_sample columns"
)
assert_header(
  file.path(raw_out_default, "04_detection_pval_qc", "detection_pval_per_sample.csv"),
  c("Barcode", "mean_detection_p_value", "detection_p_threshold", "passed_mean_detection_p_value"),
  "raw default: detection_pval_per_sample columns"
)
scores_path <- file.path(raw_out_default, "02_control_probe_pca", "control_probe_pca_scores.csv")
scores_header <- if (file.exists(scores_path)) read_header(scores_path) else character()
record(
  "raw default: control_probe_pca_scores is Barcode plus PC columns",
  length(scores_header) > 1 && scores_header[[1]] == "Barcode" && all(grepl("^PC[0-9]+$", scores_header[-1])),
  paste(scores_header, collapse = ",")
)
assert_header(
  file.path(raw_out_default, "raw_qc_summary.csv"),
  c("metric", "value"),
  "raw default: raw_qc_summary is aggregate only"
)

# ==============================================================================
# 2. CALERIE raw IDAT QC, opt-in flags enabled
# ==============================================================================

raw_out_optin <- file.path(tempdir(), "raw_optin")
unlink(raw_out_optin, recursive = TRUE)
raw_optin_run <- run_script(
  file.path(pipeline_dir, "CALERIE", "CALERIE_raw_idat_qc.R"),
  env = c(
    paste0("CALERIE_SAMPLE_SHEET=", sample_sheet_file),
    paste0("CALERIE_IDAT_DIR=", idat_dir),
    paste0("CALERIE_OUTPUT_DIR=", raw_out_optin),
    "CALERIE_CHUNK_SIZE=608",
    "CALERIE_WRITE_RAW_SAMPLE_MANIFEST=TRUE",
    "CALERIE_SAVE_CONTROL_PCA_OBJECT=TRUE"
  )
)
record("raw QC runs with opt-in flags", raw_optin_run$status == 0,
       if (raw_optin_run$status != 0) paste(tail(raw_optin_run$output, 3), collapse = " | ") else "")

assert_present(
  file.path(raw_out_optin, "01_manifest_rgset", "raw_sample_manifest_validated.csv"),
  "raw opt-in: manifest written"
)
assert_present(
  file.path(raw_out_optin, "restricted_objects", "control_probe_pca.rds"),
  "raw opt-in: PCA object written"
)
assert_header(
  file.path(raw_out_optin, "03_bead_qc", "bead_qc_per_sample.csv"),
  c("Barcode", "mean_bead_count", "bead_mean_threshold", "passed_mean_bead_count"),
  "raw opt-in: per-sample bead output stays barcode-only"
)
assert_header(
  file.path(raw_out_optin, "04_detection_pval_qc", "detection_pval_per_sample.csv"),
  c("Barcode", "mean_detection_p_value", "detection_p_threshold", "passed_mean_detection_p_value"),
  "raw opt-in: per-sample detection output stays barcode-only"
)

# ==============================================================================
# 3. CALERIE processed beta QC, default flags and opt-in flags
# ==============================================================================

# The example sample sheet has no duplicate Participant_ID + Time_Point group, so
# build a copy where row 2 joins row 1's group, and add one barcode that is absent
# from the beta matrix so the missing-barcode output has content.
sheet <- utils::read.csv(sample_sheet_file, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("Participant_ID", "Time_Point") %in% names(sheet)))
sheet[2, "Participant_ID"] <- sheet[1, "Participant_ID"]
sheet[2, "Time_Point"] <- sheet[1, "Time_Point"]

absent_row <- sheet[1, , drop = FALSE]
barcode_col <- intersect(c("barcode", "Barcode"), names(sheet))[[1]]
absent_row[[barcode_col]] <- "TEST_BARCODE_NOT_IN_BETA"
if ("Sample_Name" %in% names(absent_row)) absent_row$Sample_Name <- "TEST_SAMPLE_NOT_IN_BETA"
sheet <- rbind(sheet, absent_row)

test_sheet_file <- file.path(tempdir(), "test_samplesheet.csv")
utils::write.csv(sheet, test_sheet_file, row.names = FALSE)

run_processed <- function(out_dir, env = character()) {
  unlink(out_dir, recursive = TRUE)
  run_script(
    file.path(pipeline_dir, "CALERIE", "CALERIE_processed_beta_qc.R"),
    env = c(
      paste0("CALERIE_SAMPLE_SHEET=", test_sheet_file),
      paste0("CALERIE_BETA_FILE=", test_beta_file),
      paste0("CALERIE_PROCESSED_QC_OUTPUT_DIR=", out_dir),
      env
    )
  )
}

proc_default <- file.path(tempdir(), "proc_default")
proc_default_run <- run_processed(proc_default)
record("processed QC runs with default flags", proc_default_run$status == 0,
       if (proc_default_run$status != 0) paste(tail(proc_default_run$output, 3), collapse = " | ") else "")

assert_no_forbidden_columns(proc_default, "processed default")
assert_header(
  file.path(proc_default, "01_beta_sample_matching", "sample_sheet_barcodes_missing_from_beta.csv"),
  "Barcode",
  "processed default: sample_sheet_barcodes_missing_from_beta is barcode-only"
)
assert_rows(
  file.path(proc_default, "01_beta_sample_matching", "sample_sheet_barcodes_missing_from_beta.csv"),
  "processed default: missing barcode is reported"
)
replicate_path <- file.path(proc_default, "02_replicate_correlations", "replicate_correlations.csv")
assert_header(
  replicate_path,
  c("sample_1", "sample_2", "correlation"),
  "processed default: replicate_correlations drops replicate_group"
)
assert_rows(replicate_path, "processed default: replicate pairs were produced")
assert_header(
  file.path(proc_default, "03_uniformity_check", "uniformity_check.csv"),
  c("Barcode", "dip_p_value", "p_threshold", "non_unimodal"),
  "processed default: uniformity output unchanged"
)
assert_header(
  file.path(proc_default, "processed_beta_qc_summary.csv"),
  c("metric", "value"),
  "processed default: processed_beta_qc_summary is aggregate only"
)

proc_optin <- file.path(tempdir(), "proc_optin")
proc_optin_run <- run_processed(
  proc_optin,
  env = c("CALERIE_INCLUDE_SAMPLE_NAME=TRUE", "CALERIE_INCLUDE_REPLICATE_GROUP=TRUE")
)
record("processed QC runs with opt-in flags", proc_optin_run$status == 0,
       if (proc_optin_run$status != 0) paste(tail(proc_optin_run$output, 3), collapse = " | ") else "")
assert_header(
  file.path(proc_optin, "01_beta_sample_matching", "sample_sheet_barcodes_missing_from_beta.csv"),
  c("Barcode", "Sample_Name"),
  "processed opt-in: Sample_Name restored"
)
assert_header(
  file.path(proc_optin, "02_replicate_correlations", "replicate_correlations.csv"),
  c("replicate_group", "sample_1", "sample_2", "correlation"),
  "processed opt-in: replicate_group restored"
)

# ==============================================================================
# 4. Shared helper column contract
# ==============================================================================

source(file.path(pipeline_dir, "dnam_replicate_correlations.R"))
beta_stub <- matrix(
  seq(0.05, 0.8, length.out = 16),
  nrow = 4,
  dimnames = list(NULL, c("A", "B", "C", "D"))
)
pheno_stub <- data.frame(
  barcode = c("A", "B", "C", "D"),
  Participant_ID = c("1", "1", "2", "2"),
  Time_Point = c("base", "base", "24", "24"),
  stringsAsFactors = FALSE
)
replicate_stub <- function(...) {
  dnam_replicate_correlations(
    beta_stub,
    pheno = pheno_stub,
    sample_id_col = "barcode",
    replicate_group_cols = c("Participant_ID", "Time_Point"),
    ...
  )
}
record(
  "helper: default omits replicate_group",
  identical(colnames(replicate_stub()$results), c("sample_1", "sample_2", "correlation")),
  paste(colnames(replicate_stub()$results), collapse = ",")
)
record(
  "helper: opt-in includes replicate_group",
  identical(
    colnames(replicate_stub(include_replicate_group = TRUE)$results),
    c("replicate_group", "sample_1", "sample_2", "correlation")
  )
)
record(
  "helper: produces one pair per group",
  nrow(replicate_stub()$results) == 2
)
empty_helper <- dnam_replicate_correlations(
  beta_stub,
  pheno = pheno_stub[0, , drop = FALSE],
  sample_id_col = "barcode",
  replicate_group_cols = c("Participant_ID", "Time_Point")
)
record(
  "helper: empty result matches default column contract",
  identical(colnames(empty_helper$results), c("sample_1", "sample_2", "correlation")) &&
    nrow(empty_helper$results) == 0
)
pair_helper <- dnam_replicate_correlations(
  beta_stub,
  pair_data = data.frame(
    sample_1 = c("A", "C"),
    sample_2 = c("B", "D"),
    stringsAsFactors = FALSE
  )
)
record(
  "helper: pair file path omits replicate_group",
  identical(colnames(pair_helper$results), c("sample_1", "sample_2", "correlation"))
)

# ==============================================================================
# 5. Main.R configuration defaults
# ==============================================================================

main_lines <- readLines(file.path(pipeline_dir, "Main.R"), warn = FALSE)
default_of <- function(name) {
  hit <- grep(paste0("^", name, " <- "), main_lines, value = TRUE)
  if (!length(hit)) return(NA_character_)
  trimws(sub(paste0("^", name, " <- "), "", hit[[1]]))
}
flag_defaults <- c(
  MASK_IDS = "FALSE",
  SAVE_RGSET = "FALSE",
  SAVE_CONTROL_PCA_OBJECT = "FALSE",
  WRITE_RAW_SAMPLE_MANIFEST = "FALSE",
  INCLUDE_REPLICATE_GROUP = "FALSE"
)
record(
  "Main.R: privacy flags default to FALSE",
  all(vapply(names(flag_defaults), function(nm) identical(default_of(nm), unname(flag_defaults[[nm]])), logical(1))),
  paste(
    sprintf("%s=%s", names(flag_defaults), vapply(names(flag_defaults), default_of, character(1))),
    collapse = " "
  )
)
record(
  "Main.R: replicate correlations opt out of the group column",
  any(grepl("include_replicate_group = INCLUDE_REPLICATE_GROUP", main_lines, fixed = TRUE))
)

# ==============================================================================
# 6. Scratch CLI scripts share the same defaults
# ==============================================================================

if (dir.exists(scratch_code_dir)) {
  scratch_raw_path <- file.path(scratch_code_dir, "run_raw_idat_qc.R")
  scratch_processed_path <- file.path(scratch_code_dir, "run_processed_beta_qc.R")

  if (file.exists(scratch_raw_path)) {
    src <- readLines(scratch_raw_path, warn = FALSE)
    flag_line <- function(name) {
      hit <- grep(name, src, fixed = TRUE, value = TRUE)
      hit <- hit[grepl("as_logical_arg", hit, fixed = TRUE)]
      hit <- hit[grepl("default = FALSE", hit, fixed = TRUE)]
      length(hit) > 0
    }
    record("scratch run_raw_idat_qc.R: manifest flag defaults to FALSE", flag_line("write_raw_sample_manifest"))
    record("scratch run_raw_idat_qc.R: PCA object flag defaults to FALSE", flag_line("save_control_pca_object"))
    record(
      "scratch run_raw_idat_qc.R: no combined sample-sheet metrics file",
      !any(grepl("raw_qc_sample_metrics", src, fixed = TRUE))
    )
    record(
      "scratch run_raw_idat_qc.R: manifest and PCA object are gated",
      any(grepl("if (write_raw_sample_manifest)", src, fixed = TRUE)) &&
        any(grepl("if (save_control_pca_object)", src, fixed = TRUE))
    )
  }

  if (file.exists(scratch_processed_path)) {
    src <- readLines(scratch_processed_path, warn = FALSE)
    hit <- grep("include_replicate_group", src, fixed = TRUE, value = TRUE)
    hit <- hit[grepl("as_logical_arg", hit, fixed = TRUE) & grepl("default = FALSE", hit, fixed = TRUE)]
    record("scratch run_processed_beta_qc.R: replicate group flag defaults to FALSE", length(hit) > 0)
  }

  if (file.exists(scratch_raw_path)) {
    scratch_out <- file.path(tempdir(), "scratch_raw")
    unlink(scratch_out, recursive = TRUE)
    scratch_run <- run_script(
      scratch_raw_path,
      args = c(
        paste0("--sample_sheet=", shQuote(sample_sheet_file)),
        paste0("--idat_dir=", shQuote(idat_dir)),
        paste0("--output_dir=", shQuote(scratch_out))
      )
    )
    record(
      "scratch run_raw_idat_qc.R: runs with default flags",
      scratch_run$status == 0,
      if (scratch_run$status != 0) paste(tail(scratch_run$output, 3), collapse = " | ") else ""
    )
    if (scratch_run$status == 0) {
      assert_absent(
        c(
          file.path(scratch_out, "raw_sample_manifest_validated.csv"),
          file.path(scratch_out, "missing_idats.csv"),
          file.path(scratch_out, "control_probe_pca.rds")
        ),
        "scratch raw default: manifest, missing_idats and PCA object not written"
      )
      assert_no_forbidden_columns(scratch_out, "scratch raw default")
      assert_header(
        file.path(scratch_out, "bead_qc_per_sample.csv"),
        c("barcode", "mean_bead_count", "bead_mean_threshold", "passed_mean_bead_count"),
        "scratch raw default: bead_qc_per_sample columns"
      )
      assert_header(
        file.path(scratch_out, "detection_pval_per_sample.csv"),
        c("barcode", "mean_detection_p_value", "detection_p_threshold", "passed_mean_detection_p_value"),
        "scratch raw default: detection_pval_per_sample columns"
      )
    }
  }

  if (file.exists(scratch_processed_path)) {
    scratch_proc_out <- file.path(tempdir(), "scratch_proc")
    unlink(scratch_proc_out, recursive = TRUE)
    scratch_proc_run <- run_script(
      scratch_processed_path,
      args = c(
        paste0("--beta_file=", shQuote(test_beta_file)),
        paste0("--pheno_file=", shQuote(test_sheet_file)),
        paste0("--output_dir=", shQuote(scratch_proc_out))
      )
    )
    record(
      "scratch run_processed_beta_qc.R: runs with default flags",
      scratch_proc_run$status == 0,
      if (scratch_proc_run$status != 0) paste(tail(scratch_proc_run$output, 3), collapse = " | ") else ""
    )
    if (scratch_proc_run$status == 0) {
      assert_header(
        file.path(scratch_proc_out, "replicate_correlations.csv"),
        c("sample_1", "sample_2", "correlation"),
        "scratch processed default: replicate_correlations drops replicate_group"
      )
      assert_rows(
        file.path(scratch_proc_out, "replicate_correlations.csv"),
        "scratch processed default: replicate pairs were produced"
      )
      assert_header(
        file.path(scratch_proc_out, "uniformity_check.csv"),
        c("sample_id", "dip_p_value", "p_threshold", "non_unimodal"),
        "scratch processed default: uniformity output unchanged"
      )
    }
  }
} else {
  cat("SKIP  scratch CLI checks - no Code directory at ", scratch_code_dir, "\n", sep = "")
}

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n", strrep("-", 60), "\n", sep = "")
cat(
  "Tests run: ", length(results),
  "  Passed: ", sum(vapply(results, function(r) r$passed, logical(1))),
  "  Failed: ", length(failures), "\n",
  sep = ""
)

if (length(failures)) {
  cat("\nFailed tests:\n")
  for (name in failures) {
    cat("  - ", name, ": ", results[[name]]$detail, "\n", sep = "")
  }
  quit(status = 1L)
}

cat("All tests passed.\n")
