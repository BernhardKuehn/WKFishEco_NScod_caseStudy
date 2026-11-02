# Northwest-shelf cod EMSRR case study

This repository is a case study for developing environmentally modified stock-recruitment relationships (EMSRRs) for the three Northwest-shelf cod substocks and using those relationships in a reference-point calculation workflow. The intended downstream workflow uses the `FLasherEMSRR` and `mseEMSRR` R packages.

The analysis is organized around three substocks: **Southern**, **Northwestern**, and **Viking**. It combines stock-assessment recruitment and spawning-stock biomass (SSB) time series with environmental covariates derived from AHOI data, then aims to carry the fitted relationships into reference-point calculations under historical and projected environmental conditions.

> **Project status:** The repository currently contains scripts for processing AHOI data, calculating substock SST time series, joining SST to stock data, and several candidate stock-recruitment fitting helpers. The orchestration script for fitting EMSRRs is only a placeholder, and the climate-projection, data-conditioning, and reference-point steps in `to_do.txt` are not yet implemented here. The README describes the planned workflow as well as the present code; it does not imply that the full analysis is already executable.

## Analysis workflow

The intended sequence, following `to_do.txt`, is:

1. Prepare AHOI environmental fields and derive annual SST for the three substock areas.
2. Join the environmental time series to recruitment and SSB data, then fit candidate EMSRRs.
3. Collate POLCOMS-ERSEM climate runs and average them over the substock areas.
4. Evaluate and, if appropriate, apply bias correction to the climate projections.
5. Condition the stock and environmental data for the reference-point calculation.
6. Run the reference-point routine using the fitted EMSRRs.

Steps 1 and part of step 2 have scripts in the repository. Steps 3-6 remain planned work.

## Repository contents

| Path | Purpose |
| --- | --- |
| `data_00_Handle_AHOI_data.R` | Reads AHOI NetCDF files, extracts environmental variables, combines monthly fields, and writes time-series rasters. |
| `data_01_prepare_SST_by_stubstock.R` | Uses ICES statistical rectangles and the stock-area lookup table to calculate annual SST summaries for the three substocks. |
| `data_02_join_SST_with_substock.R` | Extracts recruitment and age-lagged SSB from FLCore stock objects, lags SST by one year, joins the series, and saves the combined data. |
| `model_01_fit_EMSRR.R` | Intended entry point for fitting EMSRRs; currently contains only a header. |
| `functions/MO.fitness.func.BevertonHolt_MARS.R` | Multiobjective fitness and fitting helpers combining Beverton-Holt stock-recruitment curves with MARS covariate effects. |
| `functions/MO.fitness.func.envHockeyStick_MARS.R` | Multiobjective fitness and fitting helpers combining hockey-stick curves with MARS covariate effects. |
| `functions/MO.fitness.func.envHockeyStick.RANSACK_MARS.R` | Experimental robust RANSAC hockey-stick and MARS helpers. |
| `to_do.txt` | Short list of planned analysis steps. |

The fitting helpers are not currently invoked by `model_01_fit_EMSRR.R`. Integration with FLasherEMSRR and mseEMSRR, including the expected model-object interfaces and reference-point inputs, still needs to be specified and implemented.

## Data requirements

The scripts expect input data to be available locally under `./data/`, relative to the repository root. Those inputs are not included in the repository. Confirm access permissions, provenance, and licensing before obtaining or sharing assessment or environmental data.

### Inputs expected by the existing scripts

| Input | Used by |
| --- | --- |
| AHOI NetCDF files under `data/Env/AHOI/v23.05/` | `data_00_Handle_AHOI_data.R` |
| ICES rectangle shapefile `data/ICES_stat_rect/StatRec_map_Areas_Full_20170124.shp` and its companion files | `data_01_prepare_SST_by_stubstock.R` |
| `data/stock/Area_NSCOD3_table_6a_inshore_offshore.xlsx` | `data_01_prepare_SST_by_stubstock.R` |
| `data/stock/cod46a7d20_FLstockObjects_assessment_2026.Rdata` | `data_02_join_SST_with_substock.R` |
| POLCOMS-ERSEM climate runs | Planned projection and reference-point steps; no processing script is present yet. |

The directory names and filenames above are hard-coded in the R scripts. Either place the input files at those paths or update the scripts to match the local data layout. The AHOI script comments describe historical monthly data for 1948-2024, while the SST preparation script expects a raster named for 1945-2024. Verify the actual date coverage and generated filename before running the later steps.

### Outputs produced by the existing scripts

- `data_00_Handle_AHOI_data.R` writes raster files under `data/Env/processed/`. The output names include the variable and date range found in the NetCDF filenames (for example, `AHOI_SST_raster_1945_2024.tif`).
- `data_01_prepare_SST_by_stubstock.R` saves `data/Env/processed/SST_yr_by_cod_substock.rda`, containing the annual substock SST table `sst_substk`.
- `data_02_join_SST_with_substock.R` saves `data/Cod_Env_data.rda`, containing the joined stock-environment table `env_stk`.

These are local analysis products, not versioned inputs or published results. The date range in a generated raster depends on the files available in the AHOI input directory.

## Running the current data-preparation scripts

Run R from the repository root so the relative `./data/...` paths resolve. In RStudio, open the repository folder and set it as the working directory. Use a fresh R session when running scripts that call `rm(list = ls())`.

Run the scripts in order:

```r
source("data_00_Handle_AHOI_data.R")
source("data_01_prepare_SST_by_stubstock.R")
source("data_02_join_SST_with_substock.R")
```

The scripts assume that required input files already exist. They do not download the data, install packages, or validate every input before processing. The processing scripts also contain plotting and interactive diagnostics.

### SST depth definition

In `data_00_Handle_AHOI_data.R`, the lookup table sets the extraction depth for `Temperature` (the SST field) to **10 m**. The script averages the available vertical levels at or above that depth limit to create the SST field. This 10 m average is the intended SST covariate for this case study.

The data-preparation code currently:

- averages monthly raster values to annual values;
- derives substock regions by joining ICES rectangles to the stock-area lookup table, grouping Northwestern subareas under `Northwestern`; and
- computes a spatial mean for the cells retained in each substock region.

Review the spatial aggregation, coordinate reference systems, missing values, and time coverage for the scientific requirements of the final analysis.

## Stock and environmental data alignment

`data_02_join_SST_with_substock.R` extracts recruitment and SSB from the FLCore stock objects. It uses the stock object's age dimension to lag SSB and labels the result by age lag (for example, `ssb_lag1`). It separately shifts each substock SST series by one year, then joins recruitment, lagged SSB, and `SST_lag1` by year and substock.

The script normalizes substock names containing "Northwest" to `Northwestern`, and names containing "South" to `Southern`. Check the resulting labels, year coverage, and lag conventions before fitting. The first observations may be missing after lagging, and the joins can exclude observations without matching years or substock labels.

## Candidate EMSRR model components

The helper functions implement candidate model structures in which the stock-recruitment component describes the baseline relationship and a MARS (`earth`) model describes environmental effects on log-scale residuals. The prediction is formed by applying the residual effect multiplicatively to the baseline prediction.

The checked-in helper files include:

- Beverton-Holt plus MARS residual effects;
- hockey-stick plus MARS residual effects; and
- an experimental robust/continuous hockey-stick approach using RANSAC and MARS.

The multiobjective fitness functions support repeated cross-validation and optional holdout evaluation. They can return predicted-versus-observed cross-validation results or evaluate a caller-supplied metric. These are research helpers rather than a finished, reproducible model-fitting pipeline: the model script does not select covariates, fit final substock models, save model objects, or pass fitted objects to either downstream package.

Before using any candidate model for management reference points, document and review the model formula, covariate and lag selection, fitting constraints, cross-validation design, performance metric, holdout design, random seeds, and final model-selection criteria.

## Dependencies and reproducibility

There is no package manifest or lockfile in the repository, so exact package versions are not pinned. The scripts call packages including `raster`, `terra`, `RNetCDF`, `sf`, `sp`, `maptools`, `maps`, `rgeos`, `rnaturalearth`, `FLCore`, `ggplot2`, `readxl`, `marmalaid`, `DescTools`, `foreach`, `doParallel`, `stringr`, `collapse`, `earth`, `optimx`, `minqa`, `caret`, `StockRecruitSET`, `pracma`, `data.table`, and `SimDesign`, along with other spatial and model-fitting dependencies. Some scripts also rely on helpers supplied by packages or loaded in the user's R environment.

Install and validate dependencies in a project-specific R library where possible. Record the R version, package versions, input-data versions, and processing parameters used for each analysis run. The scripts use parallel workers in AHOI processing and may require adjustment for the available machine and memory.

## Work still required

- Add input validation and document data provenance and access requirements.
- Finalize and run the stock-environment data join, including checks for substock names, lag definitions, missing observations, and time alignment.
- Complete `model_01_fit_EMSRR.R` to fit and compare EMSRRs for all three substocks and save reproducible model outputs.
- Define and verify how the fitted EMSRRs are represented and consumed by FLasherEMSRR and mseEMSRR.
- Implement climate-run collation, substock-area aggregation, and bias correction checks for POLCOMS-ERSEM projections.
- Implement data conditioning and the reference-point calculation workflow.
- Add validation or tests for model predictions, reference-point outputs, and reproducibility.

## Citation and interpretation

This repository is a case-study analysis workflow, not an official stock assessment or management recommendation. Cite the source data, assessment products, model methods, and software packages used in any resulting report or publication. Reference-point results should be interpreted only after the data, model, projection, and calculation steps have been reviewed and validated.