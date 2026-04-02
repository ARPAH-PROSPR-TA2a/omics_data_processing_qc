# Code Walkthrough

This document provides technical details on how each function works.

## Overview

The Olink QC pipeline consists of two main QC functions and two helper functions for saving outputs:

1. **olink_sample_qc** - Sample-level quality control
2. **plate_control_qc** - Plate-level quality control using plate controls
3. **save_sample_qc_outputs** - Save sample QC results
4. **save_plate_control_qc_outputs** - Save plate control QC results

---

## Function Details

### olink_sample_qc

**Purpose:** Assess individual sample quality using read counts for assay and control measurements.

**Method:**
1. Filter data to requested sample type (default: SAMPLE)
2. Group by sample ID and assay type
3. Sum total counts per sample per assay type
4. Pivot to wide format (one row per sample, columns for each assay type)
5. Apply pass/fail thresholds to each assay type
6. Calculate overall pass/fail (all must pass)

**Thresholds:**
| Assay Type | Threshold | Rationale |
|------------|-----------|-----------|
| assay | > 10,000 | Sufficient sequencing depth |
| ext_ctrl | > 500 | Extension control check |
| inc_ctrl | > 500 | Incubation control check |
| amp_ctrl | > 500 | Amplification control check |

**Missing Assay Types:** If a sample has no data for a particular assay type (e.g., no extension controls run), that control is treated as passing (set to 0).

**Key Decisions:**
- Uses Count column (not NPX) for QC - more direct measure of sample quality
- Requires ALL four assay types to pass (strict)
- Returns both detailed per-sample results and summary statistics

---

### plate_control_qc

**Purpose:** Assess plate-level quality using PLATE_CONTROL samples.

**Method:**
1. Filter data to PLATE_CONTROL sample type
2. Group by Block + Plate + Replicate (Investigator_ID) + AssayType
3. Sum counts for each group
4. Pivot to wide format
5. Apply pass/fail thresholds to each assay type
6. For each plate: require at least 3 of 4 replicates to pass
7. Check if failed plates have any reported NPX values

**Thresholds:**
| Assay Type | Threshold | Rationale |
|------------|-----------|-----------|
| assay | > 10,000 | Higher for plate controls |
| ext_ctrl | > 1,000 | Higher for plate controls |
| inc_ctrl | > 1,000 | Higher for plate controls |
| amp_ctrl | > 500 | Same as sample QC |

**Plate Pass Criteria:**
- At least 3 of 4 replicates must pass all four assay types
- This accounts for potential random failures while ensuring overall plate quality

**Failed Plate NPX Check:**
- After identifying failed plates, checks whether those plates have any reported NPX values
- If failed plates have NPX data, those samples may need review

**Key Decisions:**
- Uses Block + Plate as the grouping unit (matches Olink reporting)
- Higher thresholds for plate controls vs. samples (more stringent)
- Returns user-friendly "Block_X_Plate_Y" labels for failed plates

---

### save_sample_qc_outputs

**Purpose:** Save sample QC results to CSV files in a standardized location.

**Output Files:**
- `{study_name}_qc_counts.csv` - Detailed per-sample results
- `{study_name}_qc_summary.csv` - Summary statistics

---

### save_plate_control_qc_outputs

**Purpose:** Save plate control QC results to CSV files.

**Output Files:**
- `{study_name}_plate_control_qc_counts.csv` - Detailed per-replicate results
- `{study_name}_plate_control_qc_summary.csv` - Plate-level summary
- `{study_name}_failed_plate_npx_check.csv` - NPX check results
- `{study_name}_failed_plates.txt` - Simple list of failed plates

---

## Pipeline Integration

To run both QC checks:

```r
# Load data
df <- read.csv("your_data.csv", stringsAsFactors = FALSE)

# Sample QC
sample_qc <- olink_sample_qc(df)
save_sample_qc_outputs(sample_qc, "MyStudy")

# Plate Control QC
plate_qc <- plate_control_qc(df)
save_plate_control_qc_outputs(plate_qc, "MyStudy")

# Check results
print(sample_qc$qc_summary)
print(plate_qc$failed_plates)
```

## Dependencies

- R (>= 4.0)
- dplyr (for data manipulation)
- tidyr (for pivoting)

No Olink-specific packages required.
