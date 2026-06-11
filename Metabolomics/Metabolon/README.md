# Metabolon QC Pipeline

Quality control and PCA pipeline for Metabolon-style metabolomics data.

## Installation

1. Install R (version 4.0 or higher recommended)
2. Install required packages:
   ```r
   install.packages(c("ggplot2"))
   ```
3. Load the pipeline functions:
   ```r
   source("metabolon_sample_qc.R")
   source("metabolon_calculate_pca.R")
   source("metabolon_pca_plots.R")
   ```

## Data Format

The pipeline expects two objects: a phenotype table and a metabolite matrix.

### `phenotype_df`

Sample-level metadata with one row per sample.

| Column | Description | Examples |
|--------|-------------|----------|
| `SampleID` | Unique sample identifier used to match the metabolite matrix | S001, S002 |
| `SubjectID` | Participant identifier, if available | P001, P002 |
| `PlateID` | Plate identifier, if available | Plate1, Plate2 |
| `Batch` | Batch identifier, if available | Batch1, Batch2 |
| `Sex` | Sample or participant sex, if available | Female, Male |
| `Age` | Age or other numeric covariate, if available | 45, 62 |

The sample ID column can have any name, but the same name must be supplied to `sample_id_col`.

### `metabolite_df` or `metabolite_mat`

Wide metabolite abundance matrix with samples in rows and metabolites in columns.

| Element | Description | Examples |
|---------|-------------|----------|
| Rows | Samples | S001, S002 |
| Columns | Metabolites or metabolite features | glucose, lactate, X12345 |
| Values | Numeric metabolite abundances | 12.4, 8.9, NA |
| Sample IDs | Either row names or a dedicated ID column | rownames, `SampleID` |

**Note:** Sample IDs must be present either as row names or in a column passed to `matrix_id_col`.

### Example Structure

`phenotype_df`

```
SampleID | SubjectID | PlateID | Batch  | Sex    | Age
---------|-----------|---------|--------|--------|----
S001     | P001      | Plate1  | Batch1 | Female | 45
S002     | P002      | Plate1  | Batch1 | Male   | 62
S003     | P003      | Plate2  | Batch2 | Female | 50
```

`metabolite_df`

```
SampleID | glucose | lactate | citrate | X12345
---------|---------|---------|---------|-------
S001     | 12.4    | 8.9     | 3.1     | NA
S002     | 11.8    | 9.4     | 2.9     | 5.6
S003     | 18.2    | 14.0    | 4.5     | 7.1
```

## Before You Run

### 1. Split the source data into two objects

Even if the Metabolon export arrives as one large file, split it into:

- `phenotype_df`: metadata table with one row per sample
- `metabolite_df`: wide matrix or data frame with samples in rows and metabolites in columns

### 2. Confirm sample ID matching

The pipeline matches samples by ID before sample QC.

- `phenotype_df` must contain a sample identifier column passed as `sample_id_col`
- `metabolite_df` must have sample IDs as row names or in a dedicated column passed as `matrix_id_col`

If IDs do not match, the QC function reports two warning types:

- Warning A: metabolite matrix has sample IDs not found in phenotype data
- Warning B: phenotype data has sample IDs not found in the metabolite matrix

Only matched samples are used for QC, PCA, and plotting.

### 3. Confirm metabolite scale before PCA

PCA is sensitive to scale. Decide whether the PCA step should apply a log2 transform.

- If metabolite values are already log2-transformed, use `log2_transform = FALSE`
- If metabolite values are raw or non-log intensities, use `log2_transform = TRUE`

Using the wrong setting can distort PCA structure.

## Main Functions

### 1. metabolon_sample_qc

Sample-level quality control based on the median and IQR of metabolite values for each sample.

```r
result <- metabolon_sample_qc(phenotype_df = phenotype_df,
                              metabolite_df = metabolite_df,
                              sample_id_col = "SampleID",
                              matrix_id_col = NULL,
                              cutoff = 3,
                              prune = 0,
                              target_sample_prop = 0.80,
                              mask_sample_ids = FALSE,
                              verbose = TRUE)
```

**Parameters:**

- `phenotype_df`: Sample-level phenotype or metadata table
- `metabolite_df`: Wide metabolite matrix or data frame
- `sample_id_col`: Sample ID column in `phenotype_df` used for matching
- `matrix_id_col`: Optional sample ID column in `metabolite_df`; if `NULL`, row names are used
- `cutoff`: Number of standard deviations for median and IQR outlier detection (default: 3)
- `prune`: Number of worst samples to remove if fewer than 80% pass QC (default: 0)
- `target_sample_prop`: Target passing sample proportion for QC reporting (default: 0.80)
- `mask_sample_ids`: If `TRUE`, masks sample IDs in the returned per-sample QC table
- `verbose`: If `TRUE`, prints QC warnings and pruning messages

**Returns:** List containing:

- `summary`: Summary statistics including total samples, pass/fail counts, and removed sample count
- `per_sample`: Per-sample QC results with `sample_median`, `sample_iqr`, z-scores, and pass/fail flags
- `matched_sample_ids`: Sample IDs present in both phenotype and metabolite data before pruning
- `passed_samples`: Samples that pass median and IQR QC after optional pruning
- `failed_samples`: Samples that fail median or IQR QC after optional pruning
- `removed_samples`: Samples removed during pruning
- `analysis_sample_ids`: Final passing sample IDs for downstream analysis
- `analysis_matrix`: Matched and optionally pruned metabolite matrix for downstream PCA
- `analysis_pheno`: Matched and optionally pruned phenotype table for downstream plotting

**Behavior:**

- Matches phenotype rows and metabolite matrix rows by sample ID before QC
- Warns when sample IDs exist in one input object but not the other
- Computes each sample's metabolite median and IQR
- Compares those sample-level values to the global median and IQR distributions
- Flags samples using `median_fail`, `iqr_fail`, and `any_fail`
- If fewer than 80% of matched samples pass and `prune = 0`, prints a prominent warning to rerun with pruning before downstream analysis
- When `prune > 0`, iteratively removes the samples with the largest combined median/IQR deviation until either the pass rate reaches 80% or the pruning limit is reached
- Does not write files directly; `Main.R` controls output locations and file names

---

### 2. metabolon_calculate_pca

Calculate principal components from the filtered metabolite matrix.

```r
result <- metabolon_calculate_pca(mat = qc_result$analysis_matrix,
                                  sample_id_col = NULL,
                                  log2_transform = FALSE,
                                  missingness_threshold = 0.10,
                                  mask_sample_ids = FALSE)
```

**Parameters:**

- `mat`: Wide metabolite matrix or data frame with samples in rows and metabolites in columns
- `sample_id_col`: Optional sample ID column in `mat`; if `NULL`, row names are used
- `log2_transform`: If `TRUE`, applies `log2()` before PCA (default: `FALSE`)
- `missingness_threshold`: Removes metabolites with missingness above this threshold (default: 0.10)
- `mask_sample_ids`: If `TRUE`, returns a masked version of the PCA scores for saving

**Returns:** List containing:

- `scores`: Sample scores on PC1-PC5, or fewer PCs if fewer are available
- `masked_scores`: Same PCA scores with generic sample IDs when `mask_sample_ids = TRUE`
- `variance_explained`: Proportion of variance explained by each PC
- `loadings`: Metabolite loadings on PC1-PC5, or fewer PCs if fewer are available
- `pca_object`: Full `prcomp` object
- `n_samples`: Number of samples used in PCA
- `n_metabolites`: Number of metabolites used after filtering
- `n_metabolites_total`: Number of metabolites before filtering
- `metabolites_removed`: Number of metabolites removed due to missingness
- `missingness_threshold_used`: Missingness threshold applied by the function

**Behavior:**

- Requires at least 3 samples and 3 metabolites before PCA
- Uses the requested `missingness_threshold`, except datasets with `n <= 88` use a stricter 5% threshold
- Removes metabolites above the missingness threshold
- Imputes remaining missing values with the metabolite-specific median
- Removes zero-variance metabolites
- Runs `prcomp()` with centering and scaling
- Returns PC scores, loadings, variance explained, and the full PCA object

---

### 3. metabolon_pca_plots

Generate PCA plots colored by phenotype variables.

```r
result <- metabolon_pca_plots(pca_result = pca_result,
                              phenotype_df = qc_result$analysis_pheno,
                              sample_id_col = "SampleID",
                              color_vars = c("PlateID", "Batch", "Sex", "Age"),
                              pcs = c(1, 2))
```

**Parameters:**

- `pca_result`: Output from `metabolon_calculate_pca()`
- `phenotype_df`: Phenotype or metadata table containing the sample ID column and requested color variables
- `sample_id_col`: Sample ID column in `phenotype_df`
- `color_vars`: Variables to color plots by (default: `PlateID`, `Batch`, `Sex`, `Age`)
- `pcs`: Which PCs to plot (default: `c(1, 2)`)

**Returns:** List containing:

- `plots`: Named list of ggplot objects, one per available color variable
- `scores`: PCA scores merged with phenotype metadata
- `variance_explained`: Variance explained by each PC

**Behavior:**

- Merges PCA scores with phenotype metadata by sample ID
- Skips requested color variables that are not present in the phenotype data
- Stops if none of the requested color variables are available
- For categorical color variables with more than 50 unique levels, omits the legend to avoid visual clutter

## Running the Pipeline

See `Main.R` for a complete example of running all QC steps in sequence.

Key user-controlled settings in `Main.R`:

```r
output_dir <- "output"
processed_metabolomics_data_dir_name <- "processed_metabolomics_data"
MASK_IDS <- TRUE
SAVE_FULL_PCA_OBJECT <- FALSE
full_pca_output_dir_name <- "restricted_pca_output"
SAMPLE_ID_COL <- "SampleID"
```

**Standard output directory:** `output_dir` controls where shareable CSV files and plots are saved.

**Processed metabolomics data directory:** `processed_metabolomics_data_dir_name` controls the sibling folder where the matched/pruned analysis phenotype table and metabolite matrix are saved.

**Full PCA object directory:** When `SAVE_FULL_PCA_OBJECT = TRUE`, `metabolon_pca_result.rds` is saved to a sibling folder next to `output_dir`. The sibling folder name is controlled by `full_pca_output_dir_name`.

For example, with the defaults above:

```
output/
processed_metabolomics_data/
restricted_pca_output/
```

The full PCA object can contain participant-specific structure through the underlying `prcomp` object. Keep `SAVE_FULL_PCA_OBJECT = FALSE` when collaborators should receive only the standard shareable outputs.

## Output Files

Standard outputs are written to `output_dir`:

```
output/
├── metabolon_sample_qc_summary.csv        # QC summary statistics
├── metabolon_per_sample_masked.csv        # Per-sample QC results when MASK_IDS = TRUE
├── metabolon_variance_explained.csv       # Variance explained per PC
├── metabolon_pca_loadings_masked.csv      # Metabolite loadings
├── metabolon_pca_scores_masked.csv        # PCA scores, masked when MASK_IDS = TRUE
└── plots/
    ├── metabolon_pca_PlateID.png          # PCA colored by plate
    ├── metabolon_pca_Batch.png            # PCA colored by batch
    ├── metabolon_pca_Sex.png              # PCA colored by sex
    └── metabolon_pca_Age.png              # PCA colored by age
```

When `MASK_IDS = FALSE`, filenames that include `_masked` are written without `_masked`.

Processed metabolomics data are written to the sibling folder controlled by `processed_metabolomics_data_dir_name`:

```
processed_metabolomics_data/
├── metabolon_analysis_pheno.csv           # Matched/pruned phenotype table for analysis
└── metabolon_analysis_matrix.rds          # Matched/pruned metabolite matrix for analysis
```

Optional restricted output when `SAVE_FULL_PCA_OBJECT = TRUE`:

```
restricted_pca_output/
└── metabolon_pca_result.rds               # Full PCA result object, including pca_object
```

## Privacy

Set `MASK_IDS = TRUE` in `Main.R` to mask sample IDs in saved per-sample QC and PCA score CSV files. The functions keep real IDs internally so scores can be merged with phenotype metadata.

The full PCA object is more sensitive than the standard CSV outputs because it includes the underlying PCA object and sample-level structure. Use `SAVE_FULL_PCA_OBJECT = FALSE` unless the collaborator specifically needs the full object. If it is needed, set `SAVE_FULL_PCA_OBJECT = TRUE` and choose an appropriate restricted folder name with `full_pca_output_dir_name`.

## Contact

For questions or issues, contact: cpr2139@cumc.columbia.edu
