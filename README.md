# Northwest-shelf cod EMSRR case study

This repository develops environmentally modified stock-recruitment
relationships (EMSRRs) for the **Southern**, **Northwestern**, and
**Viking** cod substocks. It combines stock-assessment recruitment and
spawning-stock biomass (SSB) with historical and projected bottom
temperature, then uses the fitted relationships in a climate-dependent
reference-point workflow.

> **Project status:** Scripts are available for the workflow from
> environmental data preparation through operating-model conditioning
> and reference-point simulations. This is a research case study:
> inputs, model fits, assumptions, and reference-point calculations
> require review and validation before the results are used for
> assessment or management.

## Workflow overview

| Step | Script | Purpose |
|------------------------|------------------------|------------------------|
| Historical environment | `data_00.a_Handle_AHOI_data.R` | Process AHOI NetCDF data into monthly environmental rasters. |
| Historical substock summaries | `data_01_prepare_bottomT_by_stubstock.R` | Calculate annual bottom-temperature means and spatial clusters for each substock. |
| Stock-environment join | `data_02_join_bottomT_with_substock.R` | Join recruitment, lagged SSB, and lagged historical environmental covariates. |
| EMSRR fitting | `model_01_fit_EMSRR_BH.R` and `model_01_fit_EMSRR_Segreg.R` | Select and fit Beverton-Holt and segmented EMSRRs. |
| Projection bias correction | `data_00.b_Bias_correct_projections.R` | Bias-correct POLCOMS-ERSEM bottom-temperature projections against AHOI. |
| Projection realizations | `data_00.c_Generate_projection_realisations_BVAR_500.R` | Fit spatial BVARs and generate projected bottom-temperature realizations. |
| Historical-projection comparison | `data_00.d_Join_historical_and_projection_pattern.R` | Summarize history and projections on a common spatial footprint. |
| Reference-point lookup tables | `data_03_ccrefpts_make_lookup_tables.R` | Create EqSim settings and current biological reference-point tables. |
| Model and data preparation | `model_02_Prepare_EMSRR_model_and_data.R` | Refit selected EMSRRs on common-footprint history and prepare stock and environmental inputs. |
| Operating-model conditioning | `model_03_ccrefpts_conditionOM.R` | Create operating models for each substock, SR type, and climate setting. |
| Reference-point simulations | `model_04_ccrefpts_runRefptSim.R` | Simulate fishing-mortality targets and calculate candidate reference points. |

The `functions/` directory contains helper functions for spatial
clustering, spatial BVAR, partial-effect, and simulation utility
functions.

## Requirements

### Local input data

Input data are not included in the repository. The scripts use
hard-coded paths relative to the repository root:

| Input | Expected location |
|------------------------------------|------------------------------------|
| AHOI NetCDF files | `data/Env/AHOI/v23.05/` |
| ICES rectangle shapefile and companion files | `data/ICES_stat_rect/StatRec_map_Areas_Full_20170124.shp` |
| Substock-area lookup | `data/stock/Area_NSCOD3_table_6a_inshore_offshore.xlsx` |
| Assessment stock objects | `data/stock/cod46a7d20_FLstockObjects_assessment_2026.Rdata` |
| RCP4.5 POLCOMS-ERSEM bottom temperature | `data/Env/POLCOMS-ERSEM/POLCOM_ERSEM_GreaterNorthSea_rcp45_BottomTemp_raster_2006_2099.RData` |
| RCP8.5 POLCOMS-ERSEM bottom temperature | `data/Env/POLCOMS-ERSEM/POLCOM_ERSEM_GreaterNorthSea_rcp85_BottomTemp_raster_2006_2099.RData` |

Confirm data access permissions, provenance, and licensing before
obtaining or sharing the inputs.

The AHOI processing script describes historical data for 1948-2024,
whereas the later scripts expect a raster named
`AHOI_BottomTemp_raster_1945_2024.tif`. The generated filename reflects
the actual input dates. Check the coverage and filename, and update the
paths in the later scripts if necessary.

### R packages and environment

There is no package manifest or lockfile. Package versions are therefore
not pinned. Dependencies include:

-   Spatial and environmental processing: `raster`, `terra`, `RNetCDF`,
    `sf`, `readxl`, `marmalaid`, `maps`, `maptools`, `rgeos`,
    `rnaturalearth`, and `mclust`.
-   Projection and parallel processing: `BVAR`, `coda`, `corrplot`,
    `tsDyn`, `Rcpp`, `Rcpp2doParallel`, `foreach`, `doParallel`, and
    `doSNOW`.
-   Model fitting and simulation: `earth`, `minqa`, `collapse`,
    `ggplot2`, `DescTools`, `stringr`, `dplyr`, `data.table`, `FLCore`,
    `FLasherEMSRR`, and `mseEMSRR`.

Other dependencies are called by individual scripts and helper
functions. Some older spatial packages may need a compatible R
environment. Install and validate the required packages in a
project-specific library, and record the R and package versions used.

**For `model_03` and `model_04`, start a fresh R session without the
original `FLasher` or `mse` packages attached.** The scripts warn that
those packages can conflict with classes provided by `FLasherEMSRR` and
`mseEMSRR`.

## Running the workflow

Run scripts from the repository root so relative paths resolve
correctly. Several scripts begin with `rm(list = ls())`, which clears
objects from the current workspace. Use a fresh R session and save any
work beforehand.

The following order shows the dependencies. Existing intermediate
products can be reused when they are compatible with the current scripts
and inputs.

``` r
# Historical environmental data and stock-environment join
source("data_00.a_Handle_AHOI_data.R")
source("data_01_prepare_bottomT_by_stubstock.R")
source("data_02_join_bottomT_with_substock.R")

# Select and fit EMSRRs
source("model_01_fit_EMSRR_BH.R")
source("model_01_fit_EMSRR_Segreg.R")

# Climate projections and common-footprint summaries
source("data_00.b_Bias_correct_projections.R")
source("data_00.c_Generate_projection_realisations_BVAR_500.R")
source("data_00.d_Join_historical_and_projection_pattern.R")

# Lookup tables and reference-point inputs
source("data_03_ccrefpts_make_lookup_tables.R")
source("model_02_Prepare_EMSRR_model_and_data.R")

# In a fresh session without FLasher or mse attached:
source("model_03_ccrefpts_conditionOM.R")

# Attach data.table for the unqualified rbindlist() calls in model_04.
library(data.table)
source("model_04_ccrefpts_runRefptSim.R")
```

The scripts do not download data or install packages. Some contain plots
and interactive diagnostics, and not every input is validated in
advance. Review configuration and intermediate outputs before running
the full workflow.

Processing can require substantial time, memory, and disk space. In
particular:

-   `data_00.a` creates a cluster with 10 workers.
-   `data_00.c` is configured for 500 realizations per RCP, in batches
    of 25, using four cores for simulation.
-   `data_00.d` uses two workers and writes worker messages to
    `bottomT_workers.log`.
-   `model_03` currently uses 100 iterations and a final projection year
    of 2099.
-   `model_04` currently requests 25 workers. Adjust `ncores` in the
    script to suit the available machine.

## Historical environmental data

`data_00.a` reads the AHOI fields and writes monthly rasters under
`data/Env/processed/`. Its lookup table uses a 10 m depth limit for the
`Temperature` surface field. Bottom-temperature summaries use the
`Temperature_Bottom` field.

`data_01` averages monthly bottom temperature to annual fields,
constructs the substock regions from ICES rectangles and the area
lookup, and calculates substock means and spatial clusters. Clustering
uses temporal patterns in complete historical grid cells; it does not
enforce spatial connectivity. Inspect the selected cluster maps and
diagnostics.

The script saves these products under `data/Env/processed/`:

-   `bottomT_yr_by_cod_substock.rda`
-   `Cluster_bottomT_yr_by_cod_substock.rda`
-   `Spatial_bottomT_yr_by_cod_substock.rda`
-   `SpatialCluster_bottomT_yr_by_cod_substock.rda`

If spatial products were created with an older script that saved
`SpatRaster` external pointers directly, rerun the current `data_01`
script. The current script saves portable in-memory raster fields and
cluster maps.

## Stock-environment alignment and EMSRR fitting

`data_02` extracts recruitment and SSB from the FLCore stock objects. It
uses the recruitment age to lag SSB and shifts the environmental series
by one year before joining by year and substock. Check substock labels,
lag conventions, missing values, and retained years before fitting.

The joined table is saved as `data/Cod_Env_data.rda`, containing
`env_stk`.

The two `model_01` scripts fit Beverton-Holt (`BH`) and segmented
(`SEG`) stock-recruitment baselines with MARS environmental effects. For
each substock they use bootstrap out-of-bag performance to select a
candidate environmental variable and MARS form. The latest 10 distinct
calendar years are reserved for holdout evaluation; earlier complete
observations are used for training.

Outputs include selection summaries, holdout metrics and predictions,
partial effects, plots, and selected model objects under
`results/model_01_fit_EMSRR/`. Each script also refits its selected
model on all available fitting years and saves it under:

-   `data/EMSRRs/Cod.EMSRR_BH.RData`
-   `data/EMSRRs/Cod.EMSRR_SEG.RData`

The holdout scores refer to the training-only evaluation, not an
independent evaluation of the full-data refits.

## Climate projections and the common footprint

`data_00.b` bias-corrects the POLCOMS-ERSEM projections. `data_00.c`
fits a spatial BVAR and generates monthly realizations for RCP4.5 and
RCP8.5. The realization files are written under
`data/Env/POLCOMS-ERSEM_biasRM_BVARsims/bottomT/`, with names such as
`POLCOM_ERSEM_rcp45_bottomT_1.tiff`.

`data_00.d` uses historical fields and fixed historical cluster maps to
calculate comparable annual historical and projected series. It requires
matching longitude/latitude CRS, resolution, and grid origin, plus
ordered, complete monthly projection years. It does not resample or
bias-correct the inputs.

The common footprint retains cells with complete historical annual
values and finite monthly values in **every** realization of both
scenarios. Substock means are unweighted; cluster means use fixed
`raster::area()` weights. Inspect coverage and historical comparisons:
restricting the footprint can change the covariates used for model
fitting.

Outputs under `data/Env/processed/bottomT_projections/` include:

-   `bottomT_historical_projection_summaries.rds`: historical summaries,
    common masks, coverage, overlap diagnostics, projection tables, and
    file information.
-   `historical_comparison.csv` and `common_footprint_coverage.csv`.
-   `historical_projection_overlap.csv` and
    `historical_projection_comparison.pdf`.

History and projections remain separate in this output; they are not
spliced or lagged. Bias correction was performed before BVAR simulation,
so overlap diagnostics do not guarantee that each stochastic realization
matches history.

## Reference-point inputs and operating models

`data_03` writes:

-   `data/eqsim_lut.Rdata`, containing substock-specific EqSim settings.
-   `data/BRPs.Rdata`, containing the current biological reference-point
    table.

The script identifies ICES 2026 advice and reference-point scripts as
its sources. Review those values and assumptions against the intended
assessment.

`model_02` reads the selected BH and SEG models, stock data, historical
common-footprint summaries, and both projection scenarios. For each
substock and SR type, it retains the selected variable and MARS form,
then refits the baseline and environmental effect using common-footprint
historical values for the original fitting years.

It applies a one-year environmental lag to history and projections,
keeps them separate, and checks that the packaged historical inputs
reproduce the refitted predictions. Model selection and holdout metrics
still refer to the original historical footprint; they are not
recalculated on the common footprint.

Outputs are:

-   `data/EMSRRs/Cod.EMSRR_<BH|SEG>_common.RData`: common-footprint
    refits.
-   `model/Cod.EMSRR_<BH|SEG>_model_and_data.RData`: prepared models,
    stock data, historical covariates, and projections.

`model_03` uses these files to condition operating models for the three
substocks and both SR types. Its current setup uses a single fitted
model rather than bootstrapped model versions, with autocorrelated
lognormal recruitment deviances.

For each substock and SR type, it creates seven climate settings: a
no-climate- change setting and RCP4.5/RCP8.5 settings for the 2030,
2040, and 2050 periods. The no-climate-change series is derived from
RCP4.5 variability with its median trend removed and a historical mean
offset. The RCP period settings similarly remove the annual median trend
and use the respective decadal mean level.

Conditioned inputs are saved as:

``` text
model/Cod_<substock>_EMSRR_<BH|SEG>_ccrefpts_input.rda
```

## Reference-point simulations

`model_04` loads the conditioned operating models and runs simulations
over target fishing mortalities from 0 to 1 in steps of 0.02. It first
runs with a zero biomass trigger, then runs with a biomass-trigger
harvest control rule. It derives candidate Fmsy, Fmsy-range, Btrigger,
Blim, and Fpa values from the simulated outputs.

The script creates `model/refpts/` and saves one file per substock, SR
type, climate scenario, and period:

``` text
model/refpts/<substock>_<BH|SEG>_<scenario>_<period>_EMSRR_model.rda
```

Each file contains `df`, `df2`, and `tt`: the simulation summaries for
the two runs and their recorded elapsed time. The combined
reference-point tables are saved as `newBRPs` and `newBRPs_list` in:

``` text
model/newBRPs_EMSRR_model.rda
```

These are candidate calculations, not validated management reference
points. In particular, review the following implementation choices:

-   `ey <- 20` is used with the inclusive range `(fy - ey):fy`, which
    selects 21 calendar years.
-   The reported Fmsy range is based on predicted landings above 95% of
    the maximum; it is not a statistical confidence interval.
-   Fpa uses a smoothed curve of per-iteration fifth percentiles over
    the evaluation years, compared with `0.3 * Bmsy`. Verify that this
    matches the intended risk definition. Fpa can be missing if no point
    on the grid crosses the threshold.
-   `Bpa`, `Flim`, and `Fmgt` are not newly calculated by this script.
-   Check scenario construction, projection coverage, stochastic
    assumptions, units, and the behavior of parallel workers on the
    target platform.

## Interpretation and citation

This is a research case study, not an official stock assessment or
management recommendation. Cite the source data, assessment products,
model methods, and software used in any resulting report or publication.
Interpret results only after the data, model fits, projections, and
reference-point calculations have been reviewed and validated.

