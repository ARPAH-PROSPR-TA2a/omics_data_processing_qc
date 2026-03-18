# Code Walkthrough

This document provides technical details on how each function works.

## Overview

The Olink QC pipeline consists of three functions that can be run in sequence:

1. **Sample QC** - Identifies samples with abnormal median or IQR
2. **PCA Calculation** - Computes principal components from protein matrix
3. **PCA Plots** - Generates visualizations colored by batch/quality variables

## Function Details

### olink_sample_qc

**Purpose:** Identify samples with outlier median or interquartile range (IQR) of NPX values.

**Method:**
1. Filter to `assay` rows only (exclude controls)
2. For each sample, calculate:
   - Median of all NPX values
   - IQR of all NPX values
3. Calculate global thresholds:
   - `median_low/high` = mean(medians) ± cutoff × sd(medians)
   - `iqr_low/high` = mean(IQRs) ± cutoff × sd(IQRs)
4. Flag samples where median or IQR falls outside thresholds
5. If >15% fail, issue warning and suggest pruning
6. If `max_prune > 0`, iteratively remove worst samples

**Pruning algorithm:**
- For each iteration, calculate combined deviation: |median - mean| + |iqr - mean|
- Remove sample with maximum deviation
- Recalculate thresholds with remaining samples
- Repeat until `max_prune` samples removed or no failures remain

**Key decisions:**
- Uses SD-based thresholds (default: 3 SD) following standard Olink conventions
- 15% threshold is protective against over-pruning
- Combined deviation metric balances both median and IQR outliers
- Sample lists (passed/failed/removed) only printed when pruning is actually used
- Global median/IQR added to per_sample output for context

---

### olink_calculate_pca

**Purpose:** Compute principal components from the protein expression matrix.

**Method:**
1. Filter to `assay` rows AND `SampleType == "SAMPLE"` (excludes all controls)
2. Pivot from long to wide format (sample × protein using OlinkID)
3. Remove proteins with >10% missing values (or >5% for n ≤ 88)
4. Impute remaining NAs with column median
5. Run PCA with scaling and centering
6. Extract scores (PC1-PC5), loadings, and variance explained

**Missing data handling:**
- Proteins with high missingness (>10%) are removed entirely
- Remaining missing values are imputed with protein-specific median
- This approach preserves samples while ensuring complete matrix for PCA

**Key decisions:**
- Wide matrix orientation: proteins as columns, samples as rows (standard for PCA)
- Scale = TRUE, Center = TRUE (standard for proteomics data)
- 10% missingness threshold follows Olink best practices

---

### olink_pca_plots

**Purpose:** Generate publication-ready PCA visualizations.

**Method:**
1. Merge PCA scores with metadata
2. For each color variable:
   - Filter to non-missing values
   - Create scatter plot with points colored by variable
   - Format axes with variance explained percentages
3. Return list of ggplot objects

**Color handling:**
- Categorical variables: discrete color scale
- Continuous variables: gradient scale
- Handles both factor and character types

**Key decisions:**
- Separate function from calculation (cleaner modularity)
- One plot per variable (avoids cluttered multi-panel plots)
- Variance explained in axis labels (standard for PCA plots)

---

## Pipeline Integration

To chain functions together:

```r
# Step 1: Sample QC
qc <- olink_sample_qc(df, max_prune = 5)

# Step 2: Filter to passing samples
df_pass <- df[df$SAMPLE_ID %in% qc$passed_samples, ]

# Step 3: PCA calculation
pca <- olink_calculate_pca(df_pass)

# Step 4: Generate plots
metadata <- df_pass[!duplicated(df_pass$SAMPLE_ID), 
                    c("SAMPLE_ID", "PlateID", "SampleQC")]
plots <- olink_pca_plots(pca, metadata, color_vars = c("PlateID", "SampleQC"))

# View plots
plots$plots$PlateID
plots$plots$SampleQC
```

## Dependencies

- R (>= 4.0)
- ggplot2 (for PCA plots)
- tidyr (for pivot_wider in PCA calculation)

No Olink-specific packages required - all core calculations use base R functions.
