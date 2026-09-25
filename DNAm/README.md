# DNAm QC Pipeline

Quality control pipeline for DNA methylation array data. The pipeline is designed for collaborator-facing QC using raw IDAT files and, when available, a processed beta matrix.

The workflow consists of six independent functions that can be run in sequence:

1. **RGset Loading** - Reads sample sheet metadata, validates IDAT files, and creates an RGset
2. **Control-Probe PCA** - Calculates principal components from control probes
3. **Bead QC** - Summarizes average bead count per sample
4. **Detection P-Value QC** - Summarizes average detection p-values per sample
5. **Replicate Correlations** - Calculates correlations between replicate samples in processed beta data
6. **Uniformity Check** - Tests beta-value distributions for non-unimodality

---

## Installation

1. Install R, version 4.0 or higher recommended.

2. Install required packages:

```r
install.packages("diptest")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("minfi", "Biobase", "IlluminaHumanMethylationEPICmanifest"))
```

3. Install the correct array manifest package for your array type if different from EPIC. Examples:

```r
BiocManager::install("IlluminaHumanMethylation450kmanifest")
BiocManager::install("IlluminaHumanMethylationEPICv2manifest")
```

4. Load the pipeline functions:

```r
source("dnam_io_helpers.R")
source("dnam_load_rgset.R")
source("dnam_control_probe_pca.R")
source("dnam_bead_qc.R")
source("dnam_detection_pval_qc.R")
source("dnam_replicate_correlations.R")
source("dnam_uniformity_check.R")
```

---

## Data Format

### Raw IDAT Inputs

The raw QC steps require:

| Input | Description |
|-------|-------------|
| Sample sheet | CSV file with one row per sample |
| IDAT directory | Directory containing paired red and green IDAT files |

The sample sheet can be either an Illumina sample sheet with a `[Data]` section or a plain CSV.

Preferred sample sheet columns:

| Column | Description | Example |
|--------|-------------|---------|
| `barcode` | IDAT basename and sample ID | 203546000006_R01C01 |
| `Sample_Name` | Optional biological/sample label | S001_base |
| `Sentrix_ID` | Array barcode, used if no basename column exists | 203546000006 |
| `Sentrix_Position` | Array position, used if no basename column exists | R01C01 |

IDAT files should be named like:

```text
203546000006_R01C01_Red.idat
203546000006_R01C01_Grn.idat
```

If your sample sheet does not have a `barcode` column, set `basename_col <- NULL` in `Main.R` and the pipeline will construct basenames from `Sentrix_ID` and `Sentrix_Position`.

### Processed Beta Matrix

The processed beta QC steps require a matrix-like object with:

| Dimension | Contents |
|-----------|----------|
| Rows | CpGs/features |
| Columns | Samples |

The beta matrix can be an `.rds` object or a CSV with feature IDs in the first column. Column names should match the sample ID column in the metadata, usually `barcode`.

---

## Functions

### 1. `dnam_load_rgset`

**Purpose:** Load raw IDATs into an RGset.

```r
rgset <- dnam_load_rgset(sample_sheet,
                         sample_id_col = "barcode",
                         extended = TRUE)
```

**Key behavior:** `extended = TRUE` is used so bead counts are available for downstream QC.

---

### 2. `dnam_control_probe_pca`

**Purpose:** Calculate principal components from control probes.

```r
result <- dnam_control_probe_pca(rgset,
                                 variance_threshold = 0.90,
                                 pseudocount = 100)
```

**Method:** Control-probe signal is calculated as:

```r
red / (red + green + 100)
```

**Returns:**

- `scores`: Sample PC scores through the number of PCs needed to reach the variance threshold
- `variance_explained`: Per-PC variance explained and cumulative variance explained
- `summary`: Number of samples, probes, and PCs used
- `pca_object`: Full `prcomp` object

---

### 3. `dnam_bead_qc`

**Purpose:** Calculate average bead count per sample.

```r
result <- dnam_bead_qc(rgset, mean_bead_threshold = 2)
```

**Returns:**

- `per_sample`: Mean bead count and pass/fail flag per sample
- `summary`: Number and percentage of samples with mean bead count greater than the threshold

---

### 4. `dnam_detection_pval_qc`

**Purpose:** Calculate average detection p-value per sample.

```r
result <- dnam_detection_pval_qc(rgset, detection_p_threshold = 0.05)
```

**Returns:**

- `per_sample`: Mean detection p-value and pass/fail flag per sample
- `summary`: Number and percentage of samples with mean detection p-value below the threshold

---

### 5. `dnam_replicate_correlations`

**Purpose:** Calculate correlations between replicate samples in a processed beta matrix.

Replicates can be defined in two ways.

Using metadata grouping columns:

```r
result <- dnam_replicate_correlations(beta_mat,
                                      pheno = sample_sheet,
                                      sample_id_col = "barcode",
                                      replicate_group_cols = c("Participant_ID", "Time_Point"))
```

Using an explicit pair file:

```r
result <- dnam_replicate_correlations(beta_mat,
                                      pair_data = pair_data,
                                      sample_1_col = "sample_1",
                                      sample_2_col = "sample_2")
```

**Returns:**

- `results`: `sample_1`, `sample_2` and `correlation`. The `replicate_group`
  column, which combines the grouping columns and can contain
  `Participant_ID`, is only returned when `include_replicate_group = TRUE`.
- `summary`: Number of pairs plus mean, median, and minimum correlation

If no replicate pairs are found, the function returns an empty results table with
the same columns and a valid summary.

---

### 6. `dnam_uniformity_check`

**Purpose:** Run Hartigan's dip test on each sample's beta-value distribution.

```r
result <- dnam_uniformity_check(beta_mat, p_threshold = 0.05)
```

**Returns:**

- `results`: Dip-test p-value and non-unimodal flag per sample
- `summary`: Number and percentage of samples flagged as non-unimodal

---

## Running the Full Pipeline

Edit the configuration section at the top of `Main.R`, then run:

```bash
Rscript Main.R
```

Key user-controlled settings:

```r
sample_sheet_file <- "path/to/sample_sheet.csv"
idat_dir <- "path/to/idats"
beta_file <- "path/to/processed_beta_matrix.rds"
output_dir <- "output"

sample_id_col <- "barcode"
basename_col <- "barcode"
replicate_group_cols <- c("Participant_ID", "Time_Point")

RUN_RAW_IDAT_QC <- TRUE
RUN_PROCESSED_BETA_QC <- TRUE
RUN_REPLICATE_CORRELATIONS <- TRUE
RUN_UNIFORMITY_CHECK <- TRUE

MASK_IDS <- FALSE
SAVE_RGSET <- FALSE
SAVE_CONTROL_PCA_OBJECT <- FALSE
WRITE_RAW_SAMPLE_MANIFEST <- FALSE
INCLUDE_REPLICATE_GROUP <- FALSE
```

Every privacy setting defaults to `FALSE`, so a default run shares barcodes and
QC metrics only. Set a flag to `TRUE` only when you have a reason to produce a
restricted artifact.

| Setting | Controls | Default output when TRUE |
|---------|----------|--------------------------|
| `WRITE_RAW_SAMPLE_MANIFEST` | `raw_sample_manifest_validated.csv` and `missing_idats.csv`, which contain the full sample sheet | `01_manifest_rgset/` |
| `INCLUDE_REPLICATE_GROUP` | adds the `Participant_ID` + `Time_Point` group key to `replicate_correlations.csv` | `05_replicate_correlations/` |
| `SAVE_CONTROL_PCA_OBJECT` | `restricted_objects/control_probe_pca.rds` | `restricted_objects/` |
| `SAVE_RGSET` | `restricted_objects/RGset.rds` | `restricted_objects/` |
| `MASK_IDS` | replaces sample IDs in outputs with `Sample_1`, `Sample_2`, ... | `_masked` filename suffix |

## CALERIE Raw IDAT QC Shortcut

For the CALERIE raw IDAT QC workflow, use the focused script:

```bash
Rscript CALERIE/CALERIE_raw_idat_qc.R
```

Edit only these lines at the top of `CALERIE/CALERIE_raw_idat_qc.R`:

```r
sample_sheet_file <- "path/to/CALERIE_samplesheet.csv"
idat_dir <- "path/to/IDATs"
output_dir <- "CALERIE_raw_idat_qc_output"
chunk_size <- 48
```

The script creates `Barcode` using:

```r
Barcode = paste(Slide, Array, sep = "_")
```

If the sample sheet already contains `Barcode` or `barcode`, that existing column is used instead. The expected IDAT filenames are:

```text
Barcode_Red.idat
Barcode_Grn.idat
```

The CALERIE shortcut loads one chunk of samples at a time, extracts control-probe values, bead QC, and detection p-value QC, removes the chunk RGset from memory, and runs control-probe PCA after all chunks are complete. This is the recommended script for cloud runs with hundreds of samples.

Both CALERIE scripts take their configuration from environment variables, so no
editing is needed if you prefer that. Every optional, sample-level output is off
by default:

| Environment variable | Default | Controls |
|-----------------------|---------|----------|
| `CALERIE_SAMPLE_SHEET` | - | Sample sheet path |
| `CALERIE_IDAT_DIR` | - | IDAT folder path |
| `CALERIE_OUTPUT_DIR` | - | Raw QC output folder |
| `CALERIE_PROCESSED_QC_OUTPUT_DIR` | - | Processed beta QC output folder |
| `CALERIE_BETA_FILE` | `EDIT_ME/...` | Processed beta matrix path |
| `CALERIE_CHUNK_SIZE` | `608` | Samples loaded per RGset |
| `CALERIE_WRITE_RAW_SAMPLE_MANIFEST` | `FALSE` | `raw_sample_manifest_validated.csv` and `missing_idats.csv` |
| `CALERIE_SAVE_CONTROL_PCA_OBJECT` | `FALSE` | `restricted_objects/control_probe_pca.rds` |
| `CALERIE_SAVE_CHUNK_RGSETS` | `FALSE` | Per-chunk `RGset_chunk_*.rds` |
| `CALERIE_INCLUDE_SAMPLE_NAME` | `FALSE` | `Sample_Name` in `sample_sheet_barcodes_missing_from_beta.csv` |
| `CALERIE_INCLUDE_REPLICATE_GROUP` | `FALSE` | `replicate_group` in `replicate_correlations.csv` |
| `CALERIE_PERSON_COL` | `Participant_ID` | Replicate grouping column |
| `CALERIE_TIMEPOINT_COL` | `Time_Point` | Replicate grouping column |
| `CALERIE_PAIR_FILE` | empty | Explicit replicate pair file |

```bash
CALERIE_SAMPLE_SHEET=/path/to/samplesheet.csv \
CALERIE_IDAT_DIR=/path/to/IDATs \
CALERIE_OUTPUT_DIR=./raw_out \
Rscript CALERIE/CALERIE_raw_idat_qc.R
```

Run the processed beta QC separately, since raw IDAT QC can take a long time:

```bash
CALERIE_SAMPLE_SHEET=/path/to/samplesheet.csv \
CALERIE_BETA_FILE=/path/to/processed_betas.rds \
CALERIE_PROCESSED_QC_OUTPUT_DIR=./processed_out \
Rscript CALERIE/CALERIE_processed_beta_qc.R
```

---

## Output Files

Outputs are organized as:

```text
output/
├── 01_manifest_rgset/
│   ├── rgset_summary.csv
│   ├── raw_sample_manifest_validated.csv   (only if WRITE_RAW_SAMPLE_MANIFEST)
│   └── missing_idats.csv                   (only if WRITE_RAW_SAMPLE_MANIFEST)
├── 02_control_probe_pca/
│   ├── control_probe_pca_scores.csv
│   ├── control_probe_pca_variance_explained.csv
│   └── control_probe_pca_summary.csv
├── 03_bead_qc/
│   ├── bead_qc_per_sample.csv
│   └── bead_qc_summary.csv
├── 04_detection_pval_qc/
│   ├── detection_pval_per_sample.csv
│   └── detection_pval_summary.csv
├── 05_replicate_correlations/
│   ├── replicate_correlations.csv
│   └── replicate_correlations_summary.csv
├── 06_uniformity_check/
│   ├── uniformity_check.csv
│   └── uniformity_summary.csv
└── restricted_objects/
    ├── RGset.rds                (only if SAVE_RGSET)
    └── control_probe_pca.rds    (only if SAVE_CONTROL_PCA_OBJECT)
```

`replicate_correlations.csv` contains `sample_1`, `sample_2` and `correlation`.
The `replicate_group` column is only present when `INCLUDE_REPLICATE_GROUP` is
`TRUE`.

The CALERIE processed-beta script adds one folder:

```text
CALERIE_processed_beta_qc_output/
├── 01_beta_sample_matching/
│   ├── beta_sample_match_summary.csv
│   ├── sample_sheet_barcodes_missing_from_beta.csv
│   └── beta_columns_missing_from_sample_sheet.csv
├── 02_replicate_correlations/
├── 03_uniformity_check/
└── processed_beta_qc_summary.csv
```

`sample_sheet_barcodes_missing_from_beta.csv` contains only `Barcode`. It also
carries `Sample_Name` when `CALERIE_INCLUDE_SAMPLE_NAME=TRUE`.

When `MASK_IDS = TRUE`, sample-level output filenames include `_masked`.

---

## Tests

`tests/test_dnam_qc.R` checks that no output carries sample-level metadata unless
the matching opt-in flag is set, and that every opt-in flag defaults to `FALSE`.
No data is stored in this repository; point the test at your own data.

```bash
Rscript tests/test_dnam_qc.R --data_dir=/path/to/Data --scratch_dir=/path/to/project
```

| Argument | Meaning |
|----------|---------|
| `--data_dir` | Folder holding the sample sheet CSV, `IDATs/`, and a processed beta matrix `.rds` |
| `--scratch_dir` | Project folder containing the CLI scripts in `Code/`, used for the scratch-script checks |
| `DNAM_TEST_DATA_DIR` | Environment variable fallback for `--data_dir` |
| `DNAM_TEST_BETA_ROWS` | Beta rows to use for the processed-beta runs. Defaults to `2000`; set `0` for the full matrix |

---

## Privacy

Default outputs contain barcodes and QC metrics only. Sample-level identifiers
such as `Sample_Name`, `Participant_ID`, `Time_Point`, plate, well, and file
paths never reach a shared output. Each of these is available only through an
explicit opt-in:

- `raw_sample_manifest_validated.csv` and `missing_idats.csv` via
  `WRITE_RAW_SAMPLE_MANIFEST` or `CALERIE_WRITE_RAW_SAMPLE_MANIFEST`
- the `replicate_group` column via `INCLUDE_REPLICATE_GROUP` or
  `CALERIE_INCLUDE_REPLICATE_GROUP`
- the `Sample_Name` column in `sample_sheet_barcodes_missing_from_beta.csv` via
  `CALERIE_INCLUDE_SAMPLE_NAME`

Full R objects such as `RGset.rds` and `control_probe_pca.rds` may contain
participant-level structure, are written to `restricted_objects/`, and are opt-in
only.

If a sample sheet contains direct identifiers or sensitive metadata, remove those
columns before sharing outputs or before running this collaborator-facing QC
pipeline.

---

## Contact

For questions or issues, contact: cpr2139@cumc.columbia.edu
