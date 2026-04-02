# SomaScan QC Pipeline

A quality control pipeline for SomaScan proteomics data consisting of four independent functions that can be run in sequence:

1. **LOD QC** – Identifies samples with adequate signal
2. **Normalization QC** – Checks plate normalization factors
3. **Technical Replicate Correlations** – Assesses technical variation
4. **PCA Plots** – Visualizes batch effects

---

## Installation

1. Install R (version 4.0 or higher recommended)
2. Install required packages:
   ```r
   install.packages(c("ggplot2", "rlang"))
   ```
3. Load the pipeline functions:
   ```r
   source("somascan_lod_qc.R")
   source("somascan_norm_qc.R")
   source("somascan_techrep_cor.R")
   source("somascan_pca_plots.R")
   ```

> **Dependencies:** R (>= 4.0), `ggplot2` (for PCA plots), `rlang` (for non-standard evaluation in plots). All core calculations use base R.

---

## Data Format

The pipeline expects an ADAT-style data frame with:

| Column Type | Description | Examples |
|-------------|-------------|----------|
| Sample metadata | Identifying information | SampleId, PlateId, SlideId |
| Normalization scales | | NormScale_20, NormScale_0_005, NormScale_0_5 |
| SampleType | Type of well | Sample, Buffer, Calibrator, QC |
| Timepoint (optional) | Visit/follow-up time | Followup, Visit, Timepoint |
| Biological covariates | For PCA coloring | Sex, Age, Treatment, SampleGroup |
| Analyte columns | Protein measurements | Columns starting with `seq.` |

### Example Structure

```
SampleId | SampleType | PlateId | SlideId | NormScale_20 | NormScale_0_005 | NormScale_0_5 | Followup | Sex | Age | seq.12345 | seq.12346 | ...
---------|------------|---------|---------|--------------|-----------------|---------------|----------|-----|-----|-----------|-----------| ...
S001     | Sample     | P1      | A1      | 0.7039301    | 0.7556493       | 0.6901323     | 0        | M   | 45  | 1234.5    | 5678.9    | ...
S001     | Sample     | P1      | A1      | 0.7855821    | 0.7489244       | 0.6667662     | 6        | M   | 45  | 1245.6    | 5700.1    | ...
Buffer1  | Buffer     | P1      | A1      | 1.0400000    | 1.0331248       | 1.0000000     | NA       | NA  | NA  | 100.2     | 50.3      | ...
QC_1     | QC         | P1      | A1      | NA           | NA              | NA            | NA       | NA  | NA  | 5000.0    | 3000.0    | ...
```

---

## Functions

### 1. `somascan_lod_qc`

**Purpose:** Assess limit of detection (LOD) using Buffer wells.

**Method:**
1. Extract Buffer wells from the data.
2. Calculate LOD per analyte: `median(Buffer) + k * MAD(Buffer)` (k=3 by default). MAD is used instead of SD for robustness to outliers.
3. For each sample, calculate the signal-to-noise ratio (SNR) per analyte.
4. Count the proportion of analytes passing the SNR threshold per sample.
5. Flag samples failing the target proportion.
6. Optionally prune the worst-performing samples to meet `target_sample_prop`.

> **Key thresholds:** SNR threshold of 2 and 80% protein pass rate are manufacturer-recommended defaults for SomaScan.

```r
result <- somascan_lod_qc(dat,
                          snr_thresh = 2,
                          target_protein_prop = 0.80,
                          target_sample_prop = 0.80,
                          max_prune = 0,
                          sample_type_col = "SampleType",
                          sample_id_col = "SampleId",
                          sample_label = "Sample",
                          buffer_label = "Buffer",
                          mask_sample_ids = FALSE)
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `snr_thresh` | Signal-to-noise ratio threshold | `2` |
| `target_protein_prop` | Minimum proportion of proteins passing SNR per sample | `0.80` |
| `target_sample_prop` | Minimum proportion of samples passing | `0.80` |
| `max_prune` | Number of worst-performing samples to remove | `0` |
| `mask_sample_ids` | If TRUE, replaces SampleIds with generic labels | `FALSE` |

**Returns:**
- `pct_samples_with_ge80pct_proteins`: Percentage of samples passing
- `n_samples_initial`, `n_samples_final`: Sample counts before/after
- `n_pruned`: Number of samples removed
- `per_sample`: Per-sample pass/fail data (`SampleId`, `proteins_pass_prop`, `sample_pass`)
- `lod_per_analyte`: Named vector of LOD values

**Note:** If samples are pruned, extract passing IDs with:
```r
passed_sample_ids <- result$per_sample$SampleId[result$per_sample$sample_pass]
```

---

### 2. `somascan_norm_qc`

**Purpose:** Validate normalization scale factors for the three dilutions.

**Method:**
1. Extract normalization scale columns: `NormScale_20`, `NormScale_0_5`, `NormScale_0_005`.
2. Flag samples where any scale factor falls outside `[lower, upper]`.
3. Count pass/fail per dilution.

> **Key range:** 0.4–2.5 is the standard SomaScan range, representing ~2.5x up/down regulation. Samples outside this range may have hybridization or calibration issues.

```r
result <- somascan_norm_qc(dat,
                           lower = 0.4,
                           upper = 2.5,
                           sample_type_col = "SampleType",
                           sample_id_col = "SampleId",
                           sample_label = "Sample",
                           mask_sample_ids = FALSE)
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `lower`, `upper` | Acceptable range for normalization scale factors | `0.4`, `2.5` |
| `mask_sample_ids` | If TRUE, replaces SampleIds with generic labels | `FALSE` |

**Returns:**
- `fail_counts`: Number of samples failing per dilution
- `pass_counts`: Number of samples passing per dilution
- `per_sample`: Per-sample pass/fail for each dilution

---

### 3. `somascan_techrep_cor`

**Purpose:** Measure technical reproducibility using replicate wells.

**Method:**
1. Optionally filter to baseline timepoint only.
2. Identify rows sharing the same `SampleId`.
3. Compute pairwise Pearson correlations on log2-transformed data.
4. Track whether replicates are on the same or different plates.

> **Key decisions:** Log2 transformation and Pearson correlation are standard for SomaScan. Cross-plate vs. within-plate variation is tracked explicitly.

```r
result <- somascan_techrep_cor(dat,
                               replicate_ids = NULL,
                               use_qc = TRUE,
                               sample_type_col = "SampleType",
                               qc_label = "QC",
                               sample_id_col = "SampleId",
                               plate_id_col = "PlateId",
                               include_only_baseline = TRUE,
                               time_col = "Followup",
                               baseline_values = c(0, "Baseline", "BL"),
                               mask_sample_ids = FALSE)
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `use_qc` | Include QC wells in correlation analysis | `TRUE` |
| `replicate_ids` | Optional vector of SampleIds to include (e.g., biological duplicates) | `NULL` |
| `include_only_baseline` | Filter to baseline timepoint only | `TRUE` |
| `time_col` | Column containing timepoint information | `"Followup"` |
| `baseline_values` | Values in `time_col` indicating baseline | `c(0, "Baseline", "BL")` |
| `mask_sample_ids` | If TRUE, replaces SampleIds with generic labels | `FALSE` |

> **Note:** If you have single-timepoint data, set `include_only_baseline = FALSE`.

**Returns:**
- `n_rows_used`: Number of rows analyzed
- `n_ids_with_reps`: Number of unique SampleIds with replicates
- `results`: Pairwise correlations with plate information

---

### 4. `somascan_pca_plots`

**Purpose:** Visualize sources of variation for batch effect assessment.

**Method:**
1. Filter to `Sample` rows only.
2. Log2 transform data and remove zero-variance analytes.
3. Run PCA.
4. Generate scatter plots colored by each specified variable.

**Suggested `color_vars`:**
- **Technical/Batch:** `PlateId`, `SlideId`, `PlatePosition`, `Subarray`, `ScannerID`
- **Biological:** `Sex`, `Age`, `Treatment`, `SampleGroup`

```r
result <- somascan_pca_plots(dat,
                             sample_type_col = "SampleType",
                             sample_label = "Sample",
                             color_vars = c("PlateId", "SlideId", "Sex", "Age"),
                             pcs = c(1, 2),
                             log2_transform = TRUE,
                             mask_sample_ids = FALSE)
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `color_vars` | Variables to color PCA plots by | *(required)* |
| `pcs` | Which principal components to plot | `c(1, 2)` |
| `log2_transform` | Log2 transform data before PCA | `TRUE` |
| `mask_sample_ids` | If TRUE, replaces SampleIds with generic labels | `FALSE` |

**Returns:**
- `pca`: The `prcomp` object
- `variance_explained`: Variance explained by each PC (shown in axis labels)
- `scores`: Sample scores with metadata
- `plots`: List of ggplot objects

---

## Running the Full Pipeline

See `Main.R` for a complete working example. The recommended sequence is:

```r
# Step 1: LOD QC
lod <- somascan_lod_qc(dat)
passed_ids <- lod$per_sample$SampleId[lod$per_sample$sample_pass]

# Step 2: Filter to passing samples
dat_passed <- dat[dat$SampleId %in% passed_ids, ]

# Step 3: Normalization QC
norm <- somascan_norm_qc(dat_passed)

# Step 4: Technical replicate correlations
techrep <- somascan_techrep_cor(dat_passed,
                                time_col = "Followup",
                                baseline_values = c(0, "Baseline"))

# Step 5: PCA
pca <- somascan_pca_plots(dat_passed,
                          color_vars = c("PlateId", "Sex", "Age"))
```

---

## Output Files

When using the save functions, outputs are organized as:

```
output/
├── 01_lod/
│   ├── lod_summary.csv
│   ├── lod_per_analyte.rds
│   └── per_sample_lod.csv
├── 02_norm/
│   ├── norm_summary.csv
│   └── per_sample_norm.csv
├── 03_techrep/
│   ├── techrep_summary.csv
│   └── pairwise_correlations.csv
└── 04_pca/
    ├── pca_scores.csv
    ├── variance_explained.csv
    └── plots/
```

> When `mask_sample_ids = TRUE`, filenames will include `_masked`.

---

## Privacy

Set `mask_sample_ids = TRUE` in any function to prevent individual-level sample identifiers from appearing in outputs. All `SampleId` values will be replaced with generic labels (e.g., `"Sample_1"`, `"Sample_2"`, etc.).

---

## Contact

For questions or issues, contact: cpr2139@cumc.columbia.edu
