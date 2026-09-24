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
- `row_pass`: Logical vector aligned 1:1 with the input `dat` rows, marking Sample-type rows that pass LOD. Use `which(lod_result$row_pass)` to filter the ADAT safely (never positionally index into `per_sample`, which contains only Sample-type rows).
- `lod_per_analyte`: Named vector of LOD values
- `id_map`: Named vector mapping raw `SampleId` -> `Sample_#` (only when `mask_sample_ids = TRUE`). Pass this to `somascan_techrep_cor` so replicate-well labels agree with the LOD output.

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
1. Require user-specified subject, timepoint, plate, and sample-ID columns.
2. Identify replicate rows sharing the same subject and same timepoint (subject and timepoint values are trimmed before grouping so whitespace/label differences do not split pairs).
3. Compute pairwise Pearson correlations on log2-transformed data.
4. Track whether replicates are on the same or different plates.
5. Report each replicate group (≥ 2 rows) in a separate `groups` table. Singletons (n = 1) are excluded entirely.

> **Key decisions:** Log2 transformation and Pearson correlation are standard for SomaScan. Cross-plate vs. within-plate variation is tracked explicitly.

> **Privacy and timepoint note:** Replicates are only compared within the same subject and same timepoint; different follow-up visits for the same subject are never compared. The pairwise output lists both replicate wells (`Sample_i`, `Sample_j`) plus the `SubjectId`-`Timepoint` linkage. When `mask_sample_ids = TRUE`, member well IDs are masked with the **same** `Sample_#` scheme as the LOD step (pass `id_map = lod_result$id_map`) so labels are consistent across modules, while `SubjectId` is masked separately as `Subj_#` to avoid collisions.

```r
result <- somascan_techrep_cor(dat,
                               subject_id_col = "SubjectID",
                               time_col = "Followup",
                               replicate_ids = NULL,
                               use_qc = TRUE,
                               sample_type_col = "SampleType",
                               qc_label = "QC",
                               plate_id_col = "PlateId",
                               sample_id_col = "SampleId",
                               mask_sample_ids = TRUE,
                               id_map = lod_result$id_map)
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `subject_id_col` | Column containing the subject/participant identifier | required |
| `time_col` | Column containing visit/follow-up/timepoint information | required |
| `use_qc` | Include QC wells in correlation analysis | `TRUE` |
| `replicate_ids` | Optional vector of subject IDs to include | `NULL` |
| `sample_id_col` | Column containing the well/barcode IDs of the two replicates | `"SampleId"` |
| `id_map` | Mapping from raw sample ID to masked `Sample_#` (from the LOD step) | `NULL` |
| `mask_sample_ids` | If TRUE, masks subject IDs (`Subj_#`) and well IDs (`Sample_#`) | `TRUE` |

**Returns:**
- `n_samples_analyzed`: Number of rows fed into the step (all rows with non-NA subject/timepoint) — *not* the number of replicate rows
- `n_replicate_groups`: Number of subject-timepoint groups with ≥ 2 rows
- `n_replicate_rows`: Total number of replicate rows (sum of group sizes, e.g., 2 groups × 2 = 4)
- `n_replicate_pairs`: Number of pairwise correlations computed (`nrow(results)`)
- `groups`: Matched-pair groups only (`SubjectId`, `Timepoint`, `n_members`, member `SampleIds`). Singletons (n = 1) are excluded — a group appears here only if it has ≥ 2 rows.
- `results`: Pairwise correlations with masked `Sample_i`/`Sample_j`, `SubjectId`, `Timepoint`, plate info, `same_plate`, and `r`

> **Note on old `n_rows`:** the previous summary wrote `n_rows = n_samples_analyzed`, which is the total rows analyzed, not replicate rows. Use `n_replicate_rows` for that.

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

# Step 2: Filter to passing samples (row_pass is aligned 1:1 with dat rows)
dat_passed <- dat[which(lod$row_pass), ]

# Step 3: Normalization QC
norm <- somascan_norm_qc(dat_passed)

# Step 4: Technical replicate correlations
techrep <- somascan_techrep_cor(dat_passed,
                                subject_id_col = "SubjectID",
                                time_col = "Followup")

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
│   ├── rep_groups.csv
│   └── pairwise_correlations.csv
└── 04_pca/
    ├── pca_scores.csv
    ├── variance_explained.csv
    └── plots/
```

> When `mask_sample_ids = TRUE`, filenames will include `_masked`.

---

## Privacy

Set `mask_sample_ids = TRUE` in any function to prevent individual-level identifiers from appearing in outputs. Sample/well `SampleId` values are replaced with generic labels (e.g., `"Sample_1"`, `"Sample_2"`, etc.). In the technical-replicate step, replicate well IDs reuse the LOD `Sample_#` mapping, and subject IDs are masked separately as `"Subj_1"`, `"Subj_2"`, etc., so the two schemes never collide.

---

## Contact

For questions or issues, contact: cpr2139@cumc.columbia.edu
