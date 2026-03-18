# Olink QC Pipeline

Quality control pipeline for Olink proteomics data.

## Installation

1. Install R (version 4.0 or higher recommended)
2. Install required packages:
   ```r
   install.packages(c("ggplot2", "tidyr"))
   ```
3. Load the pipeline functions:
   ```r
   source("olink_sample_qc.R")
   source("olink_calculate_pca.R")
   source("olink_pca_plots.R")
   ```

## Data Format

The pipeline expects Olink long-format data with these columns:

| Column | Description | Examples |
|--------|-------------|----------|
| `SAMPLE_ID` | Sample identifier | TOP140499 |
| `NPX` | Protein expression values | Numeric |
| `AssayType` | Type of assay | assay, ext_ctrl, inc_ctrl |
| `SampleType` | Type of sample | SAMPLE, CONTROL, etc. |
| `PlateID` | Plate identifier | 1, 2, 3 |
| `SampleQC` | Sample quality flag | PASS, WARN |
| `OlinkID` | Protein identifier | OID45518 |
| `Assay` | Protein name | A1BG, APOA1 |

**Note:** By default, only rows where `AssayType == "assay"` AND `SampleType == "SAMPLE"` are included. This excludes all control samples (PLATE_CONTROL, NEGATIVE_CONTROL, SAMPLE_CONTROL, etc.).

### Example Structure
```
SAMPLE_ID  | AssayType | NPX     | PlateID | SampleQC
-----------|-----------|---------|---------|---------
S001       | assay     | 5.432   | 1       | PASS
S001       | assay     | 3.221   | 1       | PASS
S002       | assay     | 4.891   | 1       | PASS
```

## Main Functions

### 1. olink_sample_qc

Sample-level quality control based on median and IQR of NPX values.

```r
result <- olink_sample_qc(df,
                          sample_col = "SAMPLE_ID",
                          value_col = "NPX",
                          assay_type_col = "AssayType",
                          assay_keep = "assay",
                          cutoff = 3,
                          max_prune = 0,
                          output_dir = "output",
                          mask_sample_ids = FALSE)
```

**Parameters:**
- `cutoff`: Number of standard deviations for outlier detection (default: 3)
- `max_prune`: Number of samples to remove if >15% fail (default: 0)
- `output_dir`: Directory for saving summary CSV (default: "output")
- `mask_sample_ids`: If TRUE, mask sample IDs in saved outputs

**Returns:** List containing:
- `summary`: Summary statistics (total, pass/fail counts)
- `per_sample`: Per-sample statistics with pass/fail flags, including:
  - `sample_median`, `sample_iqr`: Sample-specific values
  - `global_median`, `global_iqr`: Population-level values
  - `median_z`, `iqr_z`: Z-scores relative to population
  - `median_fail`, `iqr_fail`, `any_fail`: Pass/fail flags
- `passed_samples`: Vector of passing sample IDs
- `failed_samples`: Vector of failing sample IDs
- `removed_samples`: Vector of removed sample IDs
- `summary_file`: Path to saved summary CSV

**Behavior:**
- If >15% of samples fail, a warning is issued with instructions to re-run with `max_prune`
- When `max_prune > 0`, iteratively removes worst samples until fail rate improves
- Automatically saves summary to CSV in `output_dir`
- Only prints sample lists (passed/failed/removed) when pruning is actually used

---

### 2. olink_calculate_pca

Calculate principal components from the protein matrix.

```r
result <- olink_calculate_pca(df,
                               sample_col = "SAMPLE_ID",
                               olink_id_col = "OlinkID",
                               value_col = "NPX",
                               assay_type_col = "AssayType",
                               assay_keep = "assay",
                               sample_type_col = "SampleType",
                               sample_type_keep = "SAMPLE",
                               impute_missing = TRUE)
```

**Parameters:**
- `olink_id_col`: Column containing protein identifiers (default: "OlinkID")
- `sample_type_col`: Column containing sample type (default: "SampleType")
- `sample_type_keep`: Value indicating sample rows to keep (default: "SAMPLE")
- `impute_missing`: Replace remaining NAs with median (default: TRUE)
- `missingness_threshold`: Remove proteins with >10% missing values (default: 0.10, adjusts to 5% for n <= 88)

**Note:** Filters to rows where `AssayType == assay_keep` AND `SampleType == sample_type_keep` to exclude control samples.

**Returns:** List containing:
- `scores`: Sample scores on PC1-PC5
- `variance_explained`: Proportion of variance explained by each PC
- `loadings`: Protein loadings on PC1-PC5
- `n_proteins`: Number of proteins used after filtering
- `n_proteins_total`: Total number of proteins in data
- `n_samples`: Number of samples used
- `proteins_removed`: Number of proteins removed due to missingness
- `pca_object`: The full prcomp object

**Note:** Handles missing data by:
1. Removing proteins with >10% missingness (or >5% for small datasets)
2. Imputing remaining NAs with protein-specific median

---

### 3. olink_pca_plots

Generate PCA plots colored by various variables.

```r
result <- olink_pca_plots(pca_result,
                          metadata,
                          sample_col = "SAMPLE_ID",
                          color_vars = c("PlateID", "SampleQC"),
                          pcs = c(1, 2))
```

**Parameters:**
- `pca_result`: Output from `olink_calculate_pca()`
- `metadata`: Data frame with sample metadata (must contain sample_col)
- `color_vars`: Variables to color plots by (default: PlateID, SampleQC)
- `pcs`: Which PCs to plot (default: c(1, 2))

**Note:** For categorical color variables with >50 unique levels, the legend is omitted to avoid visual clutter, but the legend title is preserved.

**Returns:** List containing:
- `plots`: Named list of ggplot objects (one per color_var)
- `scores`: PCA scores merged with metadata
- `variance_explained`: Variance explained by each PC

## Running the Pipeline

See `Main.R` for a complete example of running all QC steps in sequence.

## Output Files

When using the save functions, outputs are organized as:

```
output/
├── olink_sample_qc_summary.csv     # QC summary statistics
├── sample_qc_per_sample.csv        # Per-sample QC results
├── pca_result.rds                # Full PCA result object
├── variance_explained.csv        # Variance explained per PC
├── pca_loadings.csv              # Protein loadings
├── pca_scores.csv               # Sample scores with metadata
└── plots/
    ├── pca_PlateID.png          # PCA colored by plate
    └── pca_SampleQC.png         # PCA colored by sample quality
```

When `MASK_IDS = TRUE` in Main.R, CSV filenames will include `_masked`.

## Privacy

Set `MASK_IDS = TRUE` in `Main.R` to mask sample IDs when saving CSV files. The functions themselves keep real IDs internally to enable merging with metadata. Masking is applied only when writing output files.

## Contact

For questions or issues, contact: cpr2139@cumc.columbia.edu
