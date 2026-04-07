# Olink QC Pipeline

Quality control pipeline for Olink proteomics data.

## Installation

1. Install R (version 4.0 or higher recommended)
2. Install required packages:
   ```r
   install.packages(c("dplyr", "tidyr"))
   ```
3. Load the pipeline functions:
   ```r
   source("olink_sample_qc.R")
   source("plate_control_qc.R")
   source("save_sample_qc_outputs.R")
   source("save_plate_control_qc_outputs.R")
   ```

## Data Format

The pipeline expects Olink long-format data (CSV) with these columns:

| Column | Description | Required |
|--------|-------------|----------|
| `SAMPLE_ID` | Sample identifier | Yes |
| `SampleType` | Type of sample (e.g., SAMPLE, PLATE_CONTROL) | Yes |
| `AssayType` | Type of assay (assay, ext_ctrl, inc_ctrl, amp_ctrl) | Yes |
| `Count` | Read count | Yes |
| `PlateID` | Plate identifier | Yes |
| `Block` | Block number | For plate control QC |
| `WellID` | Replicate identifier | For plate control QC |
| `NPX` | Normalized protein expression | For plate control QC |

## Main Functions

### 1. olink_sample_qc

Sample-level quality control based on read counts for assay and control measurements.

```r
result <- olink_sample_qc(df,
                          sample_id_col = "SAMPLE_ID",
                          sample_type_col = "SampleType",
                          sample_type_keep = "SAMPLE",
                          assay_type_col = "AssayType",
                          count_col = "Count",
                          assay_label = "assay",
                          ext_label = "ext_ctrl",
                          inc_label = "inc_ctrl",
                          amp_label = "amp_ctrl",
                          mask_sample_id = FALSE)
```

**Parameters:**
- `sample_type_keep`: Sample type to keep (default: "SAMPLE")
- `assay_label`, `ext_label`, `inc_label`, `amp_label`: Names for assay types in data
- `mask_sample_id`: If TRUE, replace sample IDs with masked IDs

**QC Thresholds:**
- Assay counts > 10,000
- Extension control counts > 500
- Incubation control counts > 500
- Amplification control counts > 500

**Returns:** List containing:
- `qc_counts`: Per-sample counts and pass/fail for each assay type
- `qc_summary`: Summary statistics (total samples, failures by type)
- `passed_samples`: Vector of sample IDs passing all QC checks
- `id_map`: Mapping between original and masked IDs (if mask_sample_id = TRUE)

---

### 2. plate_control_qc

Plate-level quality control using PLATE_CONTROL samples.

```r
result <- plate_control_qc(df,
                           sample_type_col = "SampleType",
                           sample_type_keep = "PLATE_CONTROL",
                           block_col = "Block",
                           plate_col = "PlateID",
                           replicate_col = "WellID",
                           assay_type_col = "AssayType",
                           count_col = "Count",
                           npx_col = "NPX")
```

**Parameters:**
- `sample_type_keep`: Sample type for plate controls (default: "PLATE_CONTROL")
- `block_col`, `plate_col`, `replicate_col`: Column names for grouping

**QC Thresholds:**
- Assay counts > 10,000
- Extension control counts > 1,000
- Incubation control counts > 1,000
- Amplification control counts > 500

**Plate Pass Criteria:** At least 3 of 4 replicates must pass

**Returns:** List containing:
- `qc_counts`: Per-block/per-plate counts and pass/fail
- `qc_summary`: Plate-level summary (pass if >= 3 replicates pass)
- `failed_plates`: List of failed plates (user-friendly labels)
- `failed_plate_npx_check`: Whether failed plates have reported NPX values
- `failed_plate_npx_message`: Summary message

---

### 3. save_sample_qc_outputs

Save sample QC results to CSV files.

```r
paths <- save_sample_qc_outputs(qc_res, study_name, base_dir = ".")
```

**Parameters:**
- `qc_res`: Output from `olink_sample_qc()`
- `study_name`: Name for output files
- `base_dir`: Base directory for output (default: ".")

**Returns:** List with paths to created files.

---

### 4. save_plate_control_qc_outputs

Save plate control QC results to CSV files.

```r
paths <- save_plate_control_qc_outputs(qc_res, study_name, base_dir = ".")
```

**Parameters:**
- `qc_res`: Output from `plate_control_qc()`
- `study_name`: Name for output files
- `base_dir`: Base directory for output (default: ".")

**Returns:** List with paths to created files.

---

## Example Usage

```r
# Load functions
source("olink_sample_qc.R")
source("plate_control_qc.R")
source("save_sample_qc_outputs.R")
source("save_plate_control_qc_outputs.R")

# Load data
df <- read.csv("your_olink_data.csv", stringsAsFactors = FALSE)

# Run sample QC
sample_qc_result <- olink_sample_qc(df, mask_sample_id = FALSE)
save_sample_qc_outputs(sample_qc_result, "MyStudy")

# Run plate control QC
plate_qc_result <- plate_control_qc(df)
save_plate_control_qc_outputs(plate_qc_result, "MyStudy")
```

## Output Files

When running the save functions, outputs are organized as:

```
QC_StudyName/
├── StudyName_qc_counts.csv
├── StudyName_qc_summary.csv
├── StudyName_plate_control_qc_counts.csv
├── StudyName_plate_control_qc_summary.csv
├── StudyName_failed_plate_npx_check.csv
└── StudyName_failed_plates.txt
```

## Privacy

Set `mask_sample_id = TRUE` in `olink_sample_qc()` to replace sample IDs with generic labels (MASKED_SAMPLE_0001, etc.). An ID mapping table is returned in `$id_map` for reference if needed later.

## Contact

For questions or issues, contact: cpr2139@cumc.columbia.edu
