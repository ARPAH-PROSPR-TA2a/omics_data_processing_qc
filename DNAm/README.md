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

If no replicate pairs are found, the function returns an empty results table and a valid summary.

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
SAVE_CONTROL_PCA_OBJECT <- TRUE
```

---

## Output Files

Outputs are organized as:

```text
output/
├── 01_manifest_rgset/
│   ├── raw_sample_manifest_validated.csv
│   ├── missing_idats.csv
│   └── rgset_summary.csv
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
    ├── RGset.rds
    └── control_probe_pca.rds
```

When `MASK_IDS = TRUE`, sample-level output filenames include `_masked`.

---

## Privacy

Set `MASK_IDS = TRUE` in `Main.R` to mask sample IDs in standard outputs. Full R objects such as `RGset.rds` and `control_probe_pca.rds` may contain participant-level structure and should be treated as restricted outputs.

If a sample sheet contains direct identifiers or sensitive metadata, remove those columns before sharing outputs or before running this collaborator-facing QC pipeline.

---

## Contact

For questions or issues, contact: cpr2139@cumc.columbia.edu
