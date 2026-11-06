# Northwest-shelf cod EMSRR case study

This repository is a case study for developing environmentally modified stock-recruitment relationships (EMSRRs) for the three Northwest-shelf cod substocks and using those relationships in a reference-point calculation workflow. The intended downstream workflow uses the `FLasherEMSRR` and `mseEMSRR` R packages.

The analysis is organized around three substocks: **Southern**, **Northwestern**, and **Viking**. It combines stock-assessment recruitment and spawning-stock biomass (SSB) time series with environmental covariates derived from AHOI data, then aims to carry the fitted relationships into reference-point calculations under historical and projected environmental conditions.

> **Project status:** The repository contains scripts for processing AHOI bottom temperature, deriving substock and cluster time series, fitting candidate stock-recruitment models, and preparing POLCOMS-ERSEM projections. Data conditioning and reference-point calculations remain to be implemented. Input data and generated outputs are local and are not included in the repository.

## Analysis workflow

The intended sequence, following `to_do.txt`, is:

1. Prepare AHOI environmental fields and derive annual bottom-temperature means and clusters for the three substock areas.
2. Join the environmental time series to recruitment and SSB data, then fit candidate EMSRRs.
3. Bias-correct POLCOMS-ERSEM climate runs and generate spatial BVAR realizations.
4. Compare historical and projected annual substock and cluster time series.
5. Condition the stock and environmental data for the reference-point calculation.
6. Run the reference-point routine using the fitted EMSRRs.

Scripts for steps 1-4 are present. Steps 5-6 remain planned work.

## Repository contents

| Path | Purpose |
| --- | --- |
| `data_00.a_Handle_AHOI_data.R` | Reads AHOI NetCDF files and writes monthly environmental rasters. |
| `data_00.b_Bias_correct_projections.R` | Bias-corrects POLCOMS-ERSEM bottom-temperature projections against AHOI. |
| `data_00.c_Generate_projection_realisations_BVAR_500.R` | Generates monthly spatial BVAR realizations for each RCP. |
| `data_00.d_Join_historical_and_projection_pattern.R` | Compares historical and projected annual substock and cluster means on a common footprint. |
| `data_01_prepare_bottomT_by_stubstock.R` | Derives historical annual bottom-temperature means, cluster maps, and cluster time series. |
| `data_02_join_bottomT_with_substock.R` | Joins lagged environmental covariates to recruitment and age-lagged SSB. |
| `model_01_fit_EMSRR_BH.R`, `model_01_fit_EMSRR_Segreg.R` | Fit candidate EMSRRs and evaluate them against stock-recruitment baselines. |
| `functions/MO.fitness.func.BevertonHolt_MARS.R` | Multiobjective fitness and fitting helpers combining Beverton-Holt stock-recruitment curves with MARS covariate effects. |
| `functions/MO.fitness.func.envHockeyStick_MARS.R` | Multiobjective fitness and fitting helpers combining hockey-stick curves with MARS covariate effects. |
| `functions/MO.fitness.func.envHockeyStick.RANSACK_MARS.R` | Experimental robust RANSAC hockey-stick and MARS helpers. |
| `to_do.txt` | Short list of planned analysis steps. |

The helper functions remain research candidates. Integration with FLasherEMSRR and mseEMSRR, including the expected model-object interfaces and reference-point inputs, still needs to be specified and implemented.

## Data requirements

The scripts expect input data to be available locally under `./data/`, relative to the repository root. Those inputs are not included in the repository. Confirm access permissions, provenance, and licensing before obtaining or sharing assessment or environmental data.

### Inputs expected by the existing scripts

| Input | Used by |
| --- | --- |
| AHOI NetCDF files under `data/Env/AHOI/v23.05/` | `data_00.a_Handle_AHOI_data.R` |
| ICES rectangle shapefile `data/ICES_stat_rect/StatRec_map_Areas_Full_20170124.shp` and its companion files | `data_01_prepare_bottomT_by_stubstock.R` |
| `data/stock/Area_NSCOD3_table_6a_inshore_offshore.xlsx` | `data_01_prepare_bottomT_by_stubstock.R` |
| `data/stock/cod46a7d20_FLstockObjects_assessment_2026.Rdata` | `data_02_join_bottomT_with_substock.R` |
| POLCOMS-ERSEM bottom-temperature climate runs in `data/Env/POLCOMS-ERSEM/` | `data_00.b_Bias_correct_projections.R` |

The directory names and filenames above are hard-coded in the R scripts. Either place the input files at those paths or update the scripts to match the local data layout. The AHOI script comments describe historical monthly data for 1948-2024, while the bottom-temperature preparation script expects a raster named for 1945-2024. Verify the actual date coverage and generated filename before running the later steps.

### Outputs produced by the existing scripts

- `data_00.a_Handle_AHOI_data.R` writes monthly rasters under `data/Env/processed/`.
- `data_01_prepare_bottomT_by_stubstock.R` saves the annual mean and cluster tables and their corresponding spatial fields and cluster maps under `data/Env/processed/`.
- `data_02_join_bottomT_with_substock.R` saves `data/Cod_Env_data.rda`, containing the joined stock-environment table `env_stk`.
- `data_00.d_Join_historical_and_projection_pattern.R` saves annual projections and historical-comparison diagnostics under `data/Env/processed/bottomT_projections/`; see below.

These are local analysis products, not versioned inputs or published results. The date range in a generated raster depends on the files available in the AHOI input directory.

## Running the current data-preparation scripts

Run R from the repository root so the relative `./data/...` paths resolve. In RStudio, open the repository folder and set it as the working directory. Use a fresh R session when running scripts that call `rm(list = ls())`.

Run the scripts in order:

```r
source("data_00.a_Handle_AHOI_data.R")
source("data_01_prepare_bottomT_by_stubstock.R")
source("data_02_join_bottomT_with_substock.R")
```

The scripts assume that required input files already exist. They do not download the data, install packages, or validate every input before processing. The processing scripts also contain plotting and interactive diagnostics.

### Environmental depth definitions

In `data_00.a_Handle_AHOI_data.R`, the lookup table sets the extraction depth for `Temperature` (the SST field) to **10 m**. Bottom-temperature summaries instead use the `Temperature_Bottom` field.

The data-preparation code currently:

- averages monthly raster values to annual values;
- derives substock regions by joining ICES rectangles to the stock-area lookup table, grouping Northwestern subareas under `Northwestern`; and
- computes a spatial mean for the cells retained in each substock region.

Review the spatial aggregation, coordinate reference systems, missing values, and time coverage for the scientific requirements of the final analysis.

## Stock and environmental data alignment

`data_02_join_bottomT_with_substock.R` extracts recruitment and SSB from the FLCore stock objects. It uses the stock object's age dimension to lag SSB and labels the result by age lag (for example, `ssb_lag1`). It shifts each bottom-temperature mean and cluster series by one year, then joins recruitment, lagged SSB, and the environmental variables by year and substock.

The script normalizes substock names containing "Northwest" to `Northwestern`, and names containing "South" to `Southern`. Check the resulting labels, year coverage, and lag conventions before fitting. The first observations may be missing after lagging, and the joins can exclude observations without matching years or substock labels.

## Annual historical and projected bottom temperature

After creating the historical products with `data_01_prepare_bottomT_by_stubstock.R` and the bias-corrected BVAR TIFFs with `data_00.b` and `data_00.c`, run from the repository root:

```r
source("data_00.d_Join_historical_and_projection_pattern.R")
```

Re-run the updated `data_01` first if the historical spatial file was created with an older version: the old script saved `SpatRaster` external pointers directly, which may not survive loading. The updated script saves portable in-memory raster fields and cluster maps.

The projection script reads annual historical substock fields, fixed historical cluster maps and saved historical time series, plus files named `POLCOM_ERSEM_rcp45_bottomT_<realization>.tiff` and `POLCOM_ERSEM_rcp85_bottomT_<realization>.tiff` under `data/Env/POLCOMS-ERSEM_biasRM_BVARsims/bottomT/`. It requires `terra` and `raster`, matching longitude/latitude CRS, resolution and grid origin, and one ordered observation per month for every complete projection year. It does not resample or bias-correct the inputs. Every realization is scanned twice, so a full ensemble can take considerable time and disk I/O.

The common footprint contains cells with complete historical annual values and finite monthly values in **every** realization of both scenarios. The same cells are used for history and projections; the original cluster membership is never refitted. Within each substock, annual fields are averaged without spatial weights, matching the saved substock means. Cluster means use fixed `raster::area()` weights, matching the original clustering routine. Cells outside the common footprint are excluded and weights are renormalized within each cluster. The script stops if any mean or cluster loses all its cells. Since cropping can change the historical covariates used to fit the models, inspect the coverage and historical comparisons before using these projections with fitted models.

Outputs in `data/Env/processed/bottomT_projections/` are:

- `bottomT_historical_projection_summaries.rds`: a list containing `historical_full`, `historical_common`, `historical_comparison`, `coverage`, `common_masks`, `overlap`, `projections`, and `projection_files`. `projections$rcp45` and `$rcp85` are named by realization number; each contains an unlagged wide table with `year`, `bottomT_<substock>`, and `bottomT_<substock>_Cluster<n>` columns, matching the original saved historical table names.
- `historical_comparison.csv`: saved versus recomputed full-footprint and common-footprint historical values by variable and year.
- `common_footprint_coverage.csv`: historical and retained cell counts for each mean and cluster.
- `historical_projection_overlap.csv`: historical values and the ensemble 10th percentile, median and 90th percentile for each overlapping year, scenario and variable.
- `historical_projection_comparison.pdf`: saved versus recalculated history, plus projection medians and 10-90% ensemble bands, including the overlap.

Historical and projected records are **not spliced together**. The projections retain their overlap and are not lagged. Before using them with models, account for the one-year lag and `bottomT_mean_<substock>` renaming performed by `data_02`; no projection-to-model join or reference-point interface is implemented here. Bias correction was applied before BVAR simulations, so overlap diagnostics assess the simulated ensemble, not a guarantee that every realization matches historical observations.

## Candidate EMSRR model components

The helper functions implement candidate model structures in which the stock-recruitment component describes the baseline relationship and a MARS (`earth`) model describes environmental effects on log-scale residuals. The prediction is formed by applying the residual effect multiplicatively to the baseline prediction.

The checked-in helper files include:

- Beverton-Holt plus MARS residual effects;
- hockey-stick plus MARS residual effects; and
- an experimental robust/continuous hockey-stick approach using RANSAC and MARS.

The multiobjective fitness functions support repeated cross-validation and optional holdout evaluation. They can return predicted-versus-observed cross-validation results or evaluate a caller-supplied metric. The fitting scripts evaluate lagged bottom-temperature means and clusters as candidate predictors for Southern, Northwestern, and Viking. Their MARS residual effect is fitted on log recruitment and applied multiplicatively to the baseline prediction.

The script selects the latest 10 distinct calendar years in `env_stk` as a chronological holdout and uses every earlier complete, positive observation for training. Therefore the realized training fraction depends on the joined data coverage and may differ from 75%; no years between the training period and holdout are discarded. All compared models use the same matched observations. The script reports RMSE, MAE, log-RMSE, and log bias; positive percentage improvement in log-RMSE means the EMSRR outperformed that baseline. This is a holdout comparison, not a separately tuned cross-validation estimate. The model objects and run metadata are saved for reproducibility.

Run the fitting scripts from the repository root after preparing `data/Cod_Env_data.rda`:

```r
source("model_01_fit_EMSRR_BH.R")
source("model_01_fit_EMSRR_Segreg.R")
```

The scripts require the `earth` package and expect `env_stk` to contain `year`, `rec`, `ssb_lag1`, `variable`, `value`, and `substock`. They write results under `results/model_01_fit_EMSRR/`. Review the data alignment, holdout year coverage, fitted curves, and performance before using any candidate model for management reference points. The fitted objects are not yet adapted to the interfaces required by FLasherEMSRR or mseEMSRR.

## Dependencies and reproducibility

There is no package manifest or lockfile in the repository, so exact package versions are not pinned. The scripts call packages including `raster`, `terra`, `RNetCDF`, `sf`, `sp`, `maptools`, `maps`, `rgeos`, `rnaturalearth`, `FLCore`, `ggplot2`, `readxl`, `marmalaid`, `DescTools`, `foreach`, `doParallel`, `stringr`, `collapse`, `earth`, `optimx`, `minqa`, `caret`, `StockRecruitSET`, `pracma`, `data.table`, and `SimDesign`, along with other spatial and model-fitting dependencies. Some scripts also rely on helpers supplied by packages or loaded in the user's R environment.

Install and validate dependencies in a project-specific R library where possible. Record the R version, package versions, input-data versions, and processing parameters used for each analysis run. The scripts use parallel workers in AHOI processing and may require adjustment for the available machine and memory.

## Work still required

- Add input validation and document data provenance and access requirements.
- Finalize and run the stock-environment data join, including checks for substock names, lag definitions, missing observations, and time alignment.
- Review and validate the model fits for all three substocks.
- Define and verify how the fitted EMSRRs are represented and consumed by FLasherEMSRR and mseEMSRR.
- Review the POLCOMS-ERSEM bias correction and common-footprint comparison on the actual inputs.
- Implement data conditioning and the reference-point calculation workflow.
- Add validation or tests for model predictions, reference-point outputs, and reproducibility.

## Citation and interpretation

This repository is a case-study analysis workflow, not an official stock assessment or management recommendation. Cite the source data, assessment products, model methods, and software packages used in any resulting report or publication. Reference-point results should be interpreted only after the data, model, projection, and calculation steps have been reviewed and validated.