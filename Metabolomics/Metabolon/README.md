# Metabolon QC Pipeline

Quality control and PCA pipeline for Metabolon-style metabolomics data.

## Before You Run (Critical)

### 1) Split your source data into two objects
Even if your Metabolon export is one large file, split it into:

- `phenotype_df`: metadata table (one row per sample)
- `metabolite_df`: wide matrix/data frame (samples in rows, metabolites in columns)

### 2) Confirm sample ID matching
The pipeline matches samples by ID.

- `phenotype_df` must contain a sample identifier column (any name), passed as `sample_id_col`
- `metabolite_df` must have sample IDs either:
  - as rownames (preferred), or
  - in a dedicated ID column passed as `matrix_id_col`

If IDs do not match:
- **Warning A**: matrix has IDs not found in phenotype
- **Warning B**: phenotype has IDs not found in matrix

Only matched samples are used for analysis.

### 3) Confirm metabolite scale before PCA
You must decide whether PCA should apply log2 transform.

- If your metabolite values are **already log2-transformed**: use `log2_transform = FALSE`
- If your metabolite values are **raw/non-log intensities**: use `log2_transform = TRUE`

Using the wrong setting can distort PCA structure.

## Inputs

### `phenotype_df` (required)
- `data.frame`
- one row per sample
- includes sample-level metadata (timepoint, batch, sex, age, etc.)
- must include sample ID column used in `sample_id_col`

### `metabolite_df` (required)
- wide `matrix` or `data.frame`
- **rows = samples**, **columns = metabolites**
- sample IDs must be available as:
  - rownames (preferred), or
  - a column specified in `matrix_id_col`

## Transform Policy (QC vs PCA)

| Step | Function | Log2 default | Why |
|------|----------|--------------|-----|
| Sample QC | `metabolon_sample_qc` | No log2 by default | QC is based on per-sample median/IQR, which are robust summary metrics |
| PCA | `metabolon_calculate_pca` | User decides (`log2_transform`) | PCA geometry is sensitive to scale; user must match transform to data scale |

## Functions

### `metabolon_sample_qc`

Compares phenotype IDs and metabolite matrix IDs first.

- Warning A: matrix has samples not found in phenotype
- Warning B: phenotype has samples not found in matrix

Then computes per-sample median and IQR QC, with pruning if you rerun with `prune > 0`.

If fewer than 80% pass, the pipeline prints a prominent warning to prune before downstream analysis.

The example `Main.R` also saves the filtered analysis phenotype table and analysis matrix before PCA.

### `metabolon_calculate_pca`

Calculates PCA on the filtered metabolite matrix.

- removes low-variance and high-missingness metabolites
- imputes remaining NAs with median
- returns PC1-PC5
- supports `mask_sample_ids`

### `metabolon_pca_plots`

Creates PCA plots colored by phenotype covariates.

- suppresses legend when a categorical variable has more than 50 levels

## Example

```r
# Example: if your source is one large table, split it first
# phenotype_df <- source_df[, c("SampleID", "SubjectID", "Timepoint", "PlateID", "Batch", "Sex", "Age")]
# metabolite_df <- source_df[, c("SampleID", metabolite_columns)]
# rownames(metabolite_df) <- metabolite_df$SampleID
# metabolite_df$SampleID <- NULL

qc <- metabolon_sample_qc(phenotype_df, metabolite_df, sample_id_col = "SampleID")

pca <- metabolon_calculate_pca(qc$analysis_matrix)

# If matrix is not already log2 scale, use:
# pca <- metabolon_calculate_pca(qc$analysis_matrix, log2_transform = TRUE)

plots <- metabolon_pca_plots(pca, qc$analysis_pheno, sample_id_col = "SampleID")
```

## Output

Files are written to `output/` and include QC summaries, PCA scores, loadings, variance explained, and plots.
