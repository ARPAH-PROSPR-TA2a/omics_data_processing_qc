# SomaScan QC Pipeline

Quality control pipeline for SomaScan proteomics data.

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

## Data Format

The pipeline expects an ADAT-style data frame with:

| Column Type | Description | Examples |
|-------------|-------------|----------|
| Sample metadata | Identifying information | SampleId, PlateId, SlideId |
| SampleType | Type of well | Sample, Buffer, Calibrator, QC |
| Timepoint (optional) | Visit/follow-up time | Followup, Visit, Timepoint |
| Biological covariates | For PCA coloring | Sex, Age, Treatment, SampleGroup |
| Analyte columns | Protein measurements | Columns starting with `seq.` |

### Example Structure
```
SampleId | SampleType | PlateId | SlideId | Followup | Sex | Age | seq.12345 | seq.12346 | ...
---------|------------|---------|---------|----------|-----|-----|-----------|-----------| ...
S001     | Sample     | P1      | A1      | 0        | M   | 45  | 1234.5    | 5678.9    | ...
S001     | Sample     | P1      | A1      | 6        | M   | 45  | 1245.6    | 5700.1    | ...
Buffer1  | Buffer     | P1      | A1      | NA       | NA  | NA  | 100.2     | 50.3      | ...
QC_1     | QC         | P1      | A1      | NA       | NA  | NA  | 5000.0    | 3000.0    | ...
```

## Main Functions

### 1. somascan_lod_qc

Limit of detection (LOD) quality control based on Buffer wells.

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

**Parameters:**
- `snr_thresh`: Signal-to-noise ratio threshold (default: 2)
- `target_protein_prop`: Minimum proportion of proteins passing SNR per sample (default: 0.80)
- `target_sample_prop`: Minimum proportion of samples passing (default: 0.80)
- `max_prune`: Number of worst samples to remove (default: 0)
- `mask_sample_ids`: If TRUE, replace sample IDs with "Sample_1", "Sample_2", etc.

**Returns:** List containing:
- `pct_samples_with_ge80pct_proteins`: Percentage of samples passing
- `n_samples_initial`, `n_samples_final`: Sample counts
- `n_pruned`: Number of samples removed
- `per_sample`: Per-sample pass/fail data

**Note:** If samples are pruned, a message will print showing how to extract passing sample IDs:
```r
passed_sample_ids <- result$per_sample$SampleId[result$per_sample$sample_pass]
```

---

### 2. somascan_norm_qc

Normalization scale factor quality control for the three dilutions.

```r
result <- somascan_norm_qc(dat,
                          lower = 0.4,
                          upper = 2.5,
                          sample_type_col = "SampleType",
                          sample_id_col = "SampleId",
                          sample_label = "Sample",
                          mask_sample_ids = FALSE)
```

**Parameters:**
- `lower`, `upper`: Acceptable range for normalization scales (default: 0.4 - 2.5)
- `mask_sample_ids`: If TRUE, replace sample IDs with "Sample_1", "Sample_2", etc.

**Returns:** List containing:
- `fail_counts`: Number of samples failing per dilution
- `pass_counts`: Number of samples passing per dilution
- `per_sample`: Per-sample pass/fail for each dilution

---

### 3. somascan_techrep_cor

Technical replicate correlations.

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

**Parameters:**
- `use_qc`: Include QC wells in correlation analysis (default: TRUE)
- `replicate_ids`: Optional vector of SampleIds to include (e.g., biological duplicates)
- `include_only_baseline`: Filter to baseline timepoint only (default: TRUE)
- `time_col`: Column containing timepoint information
- `baseline_values`: Values in time_col that indicate baseline
- `mask_sample_ids`: If TRUE, replace sample IDs with "Sample_1", "Sample_2", etc.

**Returns:** List containing:
- `n_rows_used`: Number of rows analyzed
- `n_ids_with_reps`: Number of unique SampleIds with replicates
- `results`: Pairwise correlations with plate information

---

### 4. somascan_pca_plots

PCA visualization for batch effect assessment.

```r
result <- somascan_pca_plots(dat,
                             sample_type_col = "SampleType",
                             sample_label = "Sample",
                             color_vars = c("PlateId", "SlideId", "Sex", "Age"),
                             pcs = c(1, 2),
                             mask_sample_ids = FALSE)
```

**Parameters:**
- `color_vars`: Variables to color PCA plots by
- `pcs`: Which principal components to plot (default: c(1, 2))
- `log2_transform`: Log2 transform data before PCA (default: TRUE)
- `mask_sample_ids`: If TRUE, replace sample IDs with "Sample_1", "Sample_2", etc.

**Note:** For categorical color variables with >50 unique levels, the legend is omitted to avoid visual clutter, but the legend title is preserved.

**Returns:** List containing:
- `pca`: The prcomp object
- `variance_explained`: Variance explained by each PC
- `scores`: Sample scores with metadata
- `plots`: List of ggplot objects

## Running the Pipeline

See `Main.R` for a complete example of running all QC steps in sequence.

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

When `mask_sample_ids = TRUE`, filenames will include `_masked`.

## Privacy

Set `mask_sample_ids = TRUE` in any function to prevent individual-level sample identifiers from appearing in outputs. This replaces all SampleIds with generic labels ("Sample_1", "Sample_2", etc.).

## Contact

For questions or issues, contact: cpr2139@cumc.columbia.edu
