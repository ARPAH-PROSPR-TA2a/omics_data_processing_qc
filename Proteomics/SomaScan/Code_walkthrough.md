# Code Walkthrough

This document provides technical details on how each function works.

## Overview

The SomaScan QC pipeline consists of four independent functions that can be run in sequence:

1. **LOD QC** - Identifies samples with adequate signal
2. **Normalization QC** - Checks plate normalization factors
3. **Technical Replicate Correlations** - Assesses technical variation
4. **PCA Plots** - Visualizes batch effects

## Function Details

### somascan_lod_qc

**Purpose:** Assess limit of detection using Buffer wells.

**Method:**
1. Extract Buffer wells from the data
2. Calculate LOD per analyte: `median(Buffer) + k * MAD(Buffer)` where k=3 by default
3. For each sample, calculate signal-to-noise ratio (SNR) for each analyte
4. Count proportion of analytes passing SNR threshold per sample
5. Flag samples failing the target proportion
6. Optionally prune worst-performing samples to meet target_sample_prop

**Key decisions:**
- Uses MAD (not SD) for robustness to outliers in Buffer wells
- SNR threshold of 2 is standard for SomaScan
- 80% protein pass rate is manufacturer-recommended

**Output structure:**
```
$per_sample: data.frame with SampleId, proteins_pass_prop, sample_pass (Sample-type rows only)
$row_pass: logical vector aligned 1:1 with the input dat rows (Sample-type rows passing LOD)
$lod_per_analyte: named vector of LOD values
$n_pruned: count of removed samples
$id_map: named vector raw SampleId -> Sample_# (when mask_sample_ids = TRUE)
```
> **Filtering note:** `per_sample` contains only Sample-type rows, so its indices do NOT map to rows of the full ADAT (which also has Buffer/QC/Calibrator rows). Always filter with `dat[which(lod$row_pass), , drop = FALSE]`. Using `which(lod$per_sample$sample_pass)` as positional indices into `dat` silently drops dat rows past the last sample row — the historical cause of "missing" replicate pairs (e.g. a replicate at dat row 869 vanishing because only rows 1..858 were selected).

---

### somascan_norm_qc

**Purpose:** Validate normalization scale factors.

**Method:**
1. Extract normalization scale columns: NormScale_20, NormScale_0_5, NormScale_0_005
2. Flag samples where any scale factor falls outside [lower, upper]
3. Count pass/fail per dilution

**Key decisions:**
- Range 0.4-2.5 is standard for SomaScan (represents 2.5x up/down regulation)
- Samples outside this range may have issues with hybridization or calibration

---

### somascan_techrep_cor

**Purpose:** Measure technical reproducibility using replicate wells.

**Method:**
1. Require user-specified subject, timepoint, plate, and sample-ID columns
2. Trim subject/timepoint values and identify rows with the same subject and same timepoint
3. Compute pairwise Pearson correlations on log2-transformed data
4. Track whether replicates are on same or different plates

**Subject/timepoint grouping:**
- `subject_id_col` must be specified by the user because the correct column varies by dataset
- `time_col` must be specified by the user because visit/follow-up naming varies by dataset
- Replicate comparisons are only made within the same subject and same timepoint
- Different timepoints for the same subject are never compared
- `replicate_ids`, if supplied, are subject IDs
- Subject and timepoint values are `trimws()`ed before grouping so whitespace/label mismatches do not silently split replicate pairs

**Outputs:**
- `n_samples_analyzed` = all rows with non-NA subject/timepoint fed into the step (this was formerly mislabeled `n_rows`); *not* the replicate row count
- `n_replicate_groups` (alias `n_ids_with_reps`) = subject-timepoint groups with ≥ 2 rows
- `n_replicate_rows` = sum of sizes of those groups (e.g. 2 groups × 2 = 4)
- `n_replicate_pairs` = `nrow(results)`
- `groups` = matched-pair groups only (`SubjectId`, `Timepoint`, `n_members`, member `SampleIds`); singletons (n = 1) are excluded entirely
- `results` = one row per pair with `Sample_i`, `Sample_j`, `SubjectId`, `Timepoint`, `PlateId_i`, `PlateId_j`, `same_plate`, `r`

**Masking:**
- When `mask_sample_ids = TRUE`, replicate well IDs (Sample_i/j and group members) are mapped through the LOD step's `id_map` (pass `id_map = lod_result$id_map`) so `Sample_#` labels are identical across modules
- Subject IDs are masked with a separate `Subj_#` prefix so a subject is never labeled `Sample_#`
- If `id_map` is not supplied, a local per-sample map is generated inside the function

**Key decisions:**
- Uses log2 transformation (standard for SomaScan)
- Pearson correlation is standard metric
- Tracks plate information to identify cross-plate vs within-plate variation

---

### somascan_pca_plots

**Purpose:** Visualize sources of variation.

**Method:**
1. Filter to Sample rows only
2. Log2 transform and filter to variable analytes
3. Run PCA
4. Generate plots colored by each specified variable

**Variables for coloring:**
- **Technical/Batch:** PlateId, SlideId, PlatePosition, Subarray, ScannerID
- **Biological:** Sex, Age, Treatment, SampleGroup (user-provided)

**Key decisions:**
- Removes analytes with zero variance
- Uses classic ggplot theme for clean visualizations
- Variance explained shown in axis labels

---

## Pipeline Integration

To chain functions together:

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

## Dependencies

- R (>= 4.0)
- ggplot2 (for PCA plots)
- rlang (for non-standard evaluation in plots)

No other dependencies - all core calculations use base R functions.
