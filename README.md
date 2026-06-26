# data/

## Overview

Place your raw extraction spreadsheet here as:

```
data/raw/MEDIATORS_DATA_EXTRACTION.xlsx
data/raw/NITS_AS_SURROGATE_DATA_EXTRACTION.xlsx   # MRI-PDFF sheet
```

A synthetic template (`MEDIATORS_SYNTHETIC.xlsx`) is included in this folder for transparency and reproducibility testing. It mirrors the required column structure but contains no real patient data.

---

## Required sheets and columns

### Sheet: `Weight` (and `HOMA-IR`)

| Column | Description |
|--------|-------------|
| `STUDY` | Trial identifier (must match `study_id` in `R/02_study_metadata.R`) |
| `ARM` | Treatment arm label (used for drug-class classification) |
| `n(PO1)` | Sample size for primary outcome 1 |
| `n(PO2)` | Sample size for primary outcome 2 |
| `PO1` | Event rate (%) for outcome 1 |
| `PO2` | Event rate (%) for outcome 2 |
| `BASELINE WEIGHT (SD)` | Baseline value in format `mean (SD)`, e.g. `"106 (21.5)"` |
| `CHANGE AS REPORTED` | Reported change, format varies by type column |
| `HOW IS CHANGE DATA REPORTED?` | One of: `CHANGE FROM BASELINE WITH SD`, `CHANGE FROM BASELINE WITH SE`, `CHANGE FROM BASELINE WITHOUT SD`, `MEAN DIFFERENCE WITH SD`, `LEAST SQUARE MEAN WITH 95% CI`, etc. |

### Sheet: `NIT = MRI-PDFF`

| Column | Description |
|--------|-------------|
| `STUDY`, `ARM` | As above |
| `n(PO1)`, `n(PO2)` | Sample sizes |
| `PO1 (%)`, `PO2 (%)` | Event rates |
| `n(NIT)` | Sample size for the NIT sub-study |
| `BASELINE NIT (MEAN, SD)` | Format `mean (SD)` |
| `FINAL NIT (MEAN, SD)` | Format `mean (SD)` |
| `CHANGE AS REPORTED` | Reported change (SD or CI format) |
| `NIT CFB (%; MEAN)` | Directly reported % change from baseline (if available) |
| `SD OF CFB` | Directly reported SD of CFB (if available) |

---

## SD harmonisation hierarchy

For Weight and HOMA-IR sheets:
1. Directly reported SD (`WITH SD`)
2. SE × √n (`WITH SE`)
3. (CI_upper − CI_lower) / (2 × 1.96) × √n (`WITH 95% CI`)
4. Mean imputation of all available SDs (flagged as `imputed` in `SD_CFB_source`)

For the MRI-PDFF sheet (Cochrane paired-measurement formula):
1. √(SD_baseline² + SD_final² − 2ρ·SD_baseline·SD_final), default ρ = 0.5
2. Reported change SD
3. CI → SE → SD
4. Mean imputation

Absolute changes are converted to % CFB via `(change / baseline_mean) × 100`.

---

## Output files (git-ignored)

After running the analysis scripts, `output/` will contain:

| File | Description |
|------|-------------|
| `mediation_all_results.csv` | Full mediation table across all datasets × endpoints |
| `bpath_by_class_all.csv` | Per-drug-class b-path slopes (interaction model) |
| `FOREST_*.pdf` | Forest + proportion-mediated panel plots |
| `BUBBLE_*.pdf` | Trial-level scatter plots (weight loss vs. response) |
| `bfit_*_cache.*` | Cached `brms` model objects (Stan compiled fits) |
