# Code Walkthrough

This document provides technical details on how each DNAm QC function works.

## Overview

The DNAm QC pipeline consists of six steps that can be run in sequence:

1. **RGset Loading** - Reads sample sheet metadata, validates IDAT pairs, and creates an RGset
2. **Control-Probe PCA** - Computes principal components from control probes
3. **Bead QC** - Calculates average bead count per sample
4. **Detection P-Value QC** - Calculates average detection p-value per sample
5. **Replicate Correlations** - Calculates pairwise correlations among replicate samples
6. **Uniformity Check** - Runs Hartigan's dip test on processed beta values

## Function Details

### dnam_read_sample_sheet and dnam_add_idat_basenames

**Purpose:** Read sample metadata and construct IDAT basenames.

**Method:**

1. Read either an Illumina sample sheet with a `[Data]` section or a plain CSV.
2. Preserve all columns as character strings to avoid scientific notation in barcode-like IDs.
3. Construct `Basename` either from a basename column such as `barcode` or from `Sentrix_ID + Sentrix_Position`.
4. Validate existence of paired `_Red.idat` and `_Grn.idat` files.

**Key decisions:**

- A single `barcode` column is preferred because it can link raw IDATs, sample metadata, and processed beta matrix columns.
- Plain CSV and Illumina sample sheet formats are both supported for portability.
- Validation happens before RGset creation so missing files fail early with a useful table.

---

### dnam_load_rgset

**Purpose:** Create an RGset from raw IDAT files.

**Method:**

1. Use `minfi::read.metharray()` with `basenames = sample_sheet$Basename`.
2. Use `extended = TRUE` so bead-count data are retained when available.
3. Set RGset column names to the configured sample ID column.

**Key decisions:**

- RGset sample names should match downstream beta matrix column names when possible.
- RGset objects are large and potentially sensitive, so saving them is optional and directed to `restricted_objects/`.

---

### dnam_control_probe_pca

**Purpose:** Quantify variation in array control-probe signal.

**Method:**

1. Use `minfi::getControlAddress()` to identify control probes.
2. Extract red and green intensities for those probes.
3. Calculate control-probe values using:

```r
red / (red + green + 100)
```

4. Remove probes with non-finite values or zero variance.
5. Run `prcomp()` on transposed control-probe values so rows are samples.
6. Return scores through the number of PCs required to reach the cumulative variance threshold.

**Key decisions:**

- The pseudocount of 100 follows the existing project code.
- PCA is centered and scaled.
- Variance explained is saved separately from PC scores for downstream reporting.

---

### dnam_bead_qc

**Purpose:** Summarize bead count quality per sample.

**Method:**

1. Extract bead counts from the extended RGset.
2. Calculate mean bead count per sample across available probes.
3. Flag samples where mean bead count is greater than the configured threshold.
4. Summarize the number and percentage passing.

**Key decisions:**

- The default threshold is `mean_bead_count > 2`.
- The function focuses on sample-level mean bead count because this is the requested collaborator-facing metric.

---

### dnam_detection_pval_qc

**Purpose:** Summarize detection p-value quality per sample.

**Method:**

1. Calculate detection p-values with `minfi::detectionP()`.
2. Calculate average detection p-value per sample.
3. Flag samples where the average detection p-value is below the configured threshold.
4. Summarize the number and percentage passing.

**Key decisions:**

- The default threshold is `mean_detection_p_value < 0.05`.
- This is a sample-level summary. Probe-level filtering is intentionally not performed in this QC-only pipeline.

---

### dnam_replicate_correlations

**Purpose:** Assess reproducibility between replicate samples using processed beta values.

**Method:**

1. Accept a beta matrix with CpGs/features as rows and samples as columns.
2. Define replicate pairs either from an explicit pair file or from metadata grouping columns.
3. Keep only sample pairs present in the beta matrix.
4. Calculate pairwise Pearson correlations by default using pairwise complete observations.
5. Return both pair-level results and a summary table.

**Key decisions:**

- Replicate definition is user-controlled because different studies define replicates differently.
- Explicit pair files take priority over metadata group inference.
- If no replicate pairs exist, an empty result table is returned rather than failing.

---

### dnam_uniformity_check

**Purpose:** Identify samples with beta-value distributions that fail a unimodality check.

**Method:**

1. For each sample, keep finite beta values.
2. Run Hartigan's dip test using `diptest::dip.test()`.
3. Flag samples where `dip_p_value < p_threshold`.
4. Return sample-level results and a summary table.

**Key decisions:**

- The default threshold is `p < 0.05`.
- The function suppresses verbose dip-test messages for large arrays so pipeline logs remain readable.

---

## Pipeline Integration

The intended sequence in `Main.R` is:

```r
# Step 1: Read sample sheet and validate IDATs
sample_sheet <- dnam_read_sample_sheet(sample_sheet_file)
sample_sheet <- dnam_add_idat_basenames(sample_sheet, idat_dir, basename_col = "barcode")
sample_sheet <- dnam_validate_idat_pairs(sample_sheet)

# Step 2: Create RGset
rgset <- dnam_load_rgset(sample_sheet, sample_id_col = "barcode")

# Step 3: Raw IDAT QC
control_pca <- dnam_control_probe_pca(rgset)
bead_qc <- dnam_bead_qc(rgset)
detection_qc <- dnam_detection_pval_qc(rgset)

# Step 4: Processed beta QC
beta_mat <- dnam_read_beta_matrix(beta_file)
replicate_qc <- dnam_replicate_correlations(beta_mat, pheno = sample_sheet,
                                            sample_id_col = "barcode",
                                            replicate_group_cols = c("Participant_ID", "Time_Point"))
uniformity_qc <- dnam_uniformity_check(beta_mat)
```

## Dependencies

- R >= 4.0
- minfi
- Biobase
- Array manifest package such as IlluminaHumanMethylationEPICmanifest
- diptest
- Optional: wateRmelon as a fallback bead-count extractor

Core tabular outputs are written with base R functions. No project-specific packages or paths are required.
