# Taita analysis datasets

Current R analyses read project-local data from this folder.

- field contains the original forest survey workbooks, the local Species.list reference, the field protocol and associated documentation.
- mapping contains the GBIF exports, study-site reference points and the boundary shapefiles with their required sidecars.
- habitat_suitability contains the prepared WorldClim predictor grid, aligned masks, screened occurrences and scenario design. The selected-variable WorldClim archive used to prepare the inputs is in source.
- historical contains existing datasets queried by the older scripts and the historical IVI workbook.
- derived contains working datasets shared between scripts. Question 1 refreshes the juvenile matrix, WCVP taxonomic lookup and family assignment audit, Question 3 refreshes its worksheet count ledger, and the study-area script refreshes the worksheet coordinates.
- audits is the input location for independently reconciled survey audits, if supplied. These were unavailable during migration. Audited workflows can also write their downstream handoffs here.

Reporting copies are written to analysis/outputs. The source inventory records where the initial copies came from, with SHA-256 checksums at migration. The checksums for derived files may change after regeneration. WCVP reference data are loaded by the installed rWCVPdata package; the extracted tables are under derived.

## Inputs not found during migration

- Study-area raster regeneration still requires mapping/environment/elevation.tif, mapping/environment/tree_cover.tif and derived/worksheet_locations.csv. The existing study-area figure is retained.
- The archived publication ready plots.R script references historical/parameter of individuals.csv, which was not found.
- Archived cached habitat modelling inputs are unavailable. Its legacy code expects them under legacy_habitat. The active Question 5 workflow uses habitat_suitability/exploratory.
