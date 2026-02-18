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
$per_sample: data.frame with SampleId, proteins_pass_prop, sample_pass
$lod_per_analyte: named vector of LOD values
$n_pruned: count of removed samples
```

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
1. Optionally filter to baseline timepoint only
2. Identify rows with same SampleId
3. Compute pairwise Pearson correlations on log2-transformed data
4. Track whether replicates are on same or different plates

**Baseline filtering:**
- Default: `include_only_baseline = TRUE` requires both `time_col` and `baseline_values`
- If you have single timepoint data, set `include_only_baseline = FALSE`

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

## Dependencies

- R (>= 4.0)
- ggplot2 (for PCA plots)
- rlang (for non-standard evaluation in plots)

No other dependencies - all core calculations use base R functions.
