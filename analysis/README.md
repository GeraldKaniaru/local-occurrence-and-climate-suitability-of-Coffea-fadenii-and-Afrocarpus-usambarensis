# Taita Hills forest and climatic suitability analyses

This folder supports the statistical analyses in 'Occurrence of Afrocarpus usambarensis and Coffea fadenii in the Taita Hills and predicted climatic suitability across Kenya to guide restoration and surveys' manuscript. The analyses compare woody vegetation and focal-species records in Ngangao and Mbololo, model climatic suitability across Kenya, and estimate its overlap with mapped protected areas.

## Run order

Keep this directory named `analysis` inside the project folder. Open `analysis.Rproj` or run the commands below from this directory. Each command runs in its own R session.

| Script | Analysis |
| --- | --- |
| `scripts/01_woody_composition.R` | Woody assemblage composition, dissimilarity, ordination and family representation |
| `scripts/02_vegetation_structure.R` | Stem structure, cover, growth forms and taxon dominance |
| `scripts/03_focal_abundance.R` | Focal-species counts, detection frequencies and concentration |
| `scripts/04_size_class_overlap.R` | Co-occurrence of juveniles and listed larger stems |
| `scripts/05_climatic_suitability.R` | Maxent calibration and thinning scenarios, evaluation and candidate-area mapping |
| `scripts/06_protected_area_overlap.R` | Protected-area overlap for the baseline models and shortlist from step 05 |

```powershell
Rscript --vanilla scripts/01_woody_composition.R
Rscript --vanilla scripts/02_vegetation_structure.R
Rscript --vanilla scripts/03_focal_abundance.R
Rscript --vanilla scripts/04_size_class_overlap.R
$env:TAITA_HABITAT_MODE = 'taita_screened'
Rscript --vanilla scripts/05_climatic_suitability.R
Rscript --vanilla scripts/06_protected_area_overlap.R 2026-10
```

The numbered filenames follow the six analyses in the manuscript. The descriptive filenames `data_collection.R`, `study_area.R`, `sampling_design.R` and `focal_species.R` identify the separate workbook and figure scripts.

Each numbered script saves figures in `outputs/research question N/figures` and tables in `outputs/research question N/tables`, where N is the question number.

Run `Rscript --vanilla scripts/focal_species.R` to assemble the photo plate from the six JPEGs in `focal_species_photos/`. This script requires the `jpeg` package.

## Inputs, outputs and software

- See [Datasets/README.md](Datasets/README.md) for input locations and known gaps. Retain original field workbooks as source files. Scripts use `Datasets/derived` for intermediate data and `outputs` for reporting copies. Existing output folder names are used by downstream scripts.
- The manuscript uses R 4.6.1. The scripts use dplyr, ggh4x, ggplot2, jsonlite, maxnet, openxlsx, patchwork, png, readxl, rWCVP, rWCVPdata, sf, terra and vegan. Some scripts use an existing local library at `scripts/.r-library`. Package versions are not pinned in this folder.
- Check each run's status before using its outputs. The forest scripts can use provisional worksheet-based analyses when reconciled audit inputs are unavailable. Step 06 requires the six saved suitability rasters from step 05, the October 2026 Protected Planet extract and the Kenya boundary.

The organisation follows the British Ecological Society's *Reproducible Code* guide (2025).
