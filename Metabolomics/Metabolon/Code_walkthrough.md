# Code Walkthrough

## Input Contract (Important)

The Metabolon pipeline expects two separate objects:

1. `phenotype_df`:
   - one row per sample
   - contains `sample_id_col` and metadata variables

2. `metabolite_df`:
   - wide matrix/data frame
   - samples in rows, metabolites in columns
   - sample IDs in rownames (or in `matrix_id_col`)

Matching is done using phenotype sample IDs vs matrix sample IDs.
Unmatched IDs trigger Warning A/Warning B; only matched samples continue.

## `metabolon_sample_qc`

1. Compare phenotype sample IDs (`sample_id_col`) to metabolite matrix IDs.
2. Emit warnings for unmatched IDs in either direction.
3. Restrict to matched samples only.
4. Compute sample medians and IQRs across metabolites.
5. Flag samples outside `cutoff` SD from the global mean of medians/IQRs.
6. If `prune > 0`, iteratively remove the worst samples until the pass rate is acceptable.
7. Return the filtered analysis phenotype table and analysis matrix for downstream PCA.

### QC scale note
QC is run on the provided scale by default (no log2 transform required).
This is intentional because median/IQR-based QC is robust to non-normality.

## `metabolon_calculate_pca`

1. Accept filtered wide matrix input.
2. Remove high-missingness metabolites.
3. Impute remaining NAs with median.
4. Remove zero-variance metabolites.
5. Run PCA and return PC1-PC5.
6. Optionally return masked sample IDs for saving.

### Transform decision (must be explicit)
- If matrix is already log2: `log2_transform = FALSE`
- If matrix is raw/non-log: `log2_transform = TRUE`

PCA is sensitive to scale, so this choice should be made deliberately for each dataset.

## `metabolon_pca_plots`

1. Merge PCA scores with phenotype metadata.
2. Plot by each requested covariate.
3. Remove the legend for categorical variables with >50 levels.

## Pipeline order

1. Run sample QC.
2. Save the filtered analysis dataset.
3. Run PCA.
4. Save PCA results and plots.
