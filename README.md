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

The diagram summarizes the scientific workflow. Script names and
run-order details are listed below.

![Framework for estimating reference points with environmental
influence](docs/images/workflow-framework.png)

*Conceptual framework for combining stock data, environmental
observations and projections, EMSRRs, and simulations for
reference-point estimation. The repository implements the workflow in
the scripts described below.*

## Repository layout

Run scripts from the repository root so their relative paths resolve
correctly.

``` text
.
├── data_00.*.R                         # Input processing and projections
├── data_01_prepare_bottomT_by_stubstock.R
├── data_02_join_bottomT_with_substock.R
├── data_03_ccrefpts_make_lookup_tables.R
├── model_01_fit_EMSRR_*.R              # EMSRR fitting and selection
├── model_02_Prepare_EMSRR_model_and_data.R
├── model_03_ccrefpts_conditionOM.R
├── model_04_ccrefpts_runRefptSim*.R
├── output_*.R                          # Results and diagnostic plots
├── functions/                          # Shared analysis helpers
├── data/                               # External inputs and generated data
├── model/                              # Prepared models and simulations
└── results/                            # Model-selection outputs and plots
```

Input data and most generated products are not included in the
repository. The scripts create intermediate and result files under
`data/`, `model/`, and `results/`.

## Script guide

### 1. Prepare historical data

| Script | Purpose |
|------------------------------------|------------------------------------|
| `data_00.a_Handle_AHOI_data.R` | Process monthly AHOI NetCDF files into environmental rasters cutted to the region of interest and extract bottom temperature (and other variables as well if wanted) |
| `data_01_prepare_bottomT_by_stubstock.R` | Calculate annual bottom-temperature fields, substock time series, and spatial clusters. |
| `data_02_join_bottomT_with_substock.R` | Join recruitment and lagged SSB with lagged historical environmental covariates. |

### 2. Fit and select EMSRRs

| Script | Purpose |
|------------------------------------|------------------------------------|
| `model_01_fit_EMSRR_Segreg.R` | Fit and select segmented EMSRRs. |
| `model_01_fit_EMSRR_Segreg_trunc.ICES.SRR.R` | Fit the truncated-ICES segmented EMSRR variant. |
| `output_05_plot_spatial_clusters.R` | Plot spatial cluster maps and, when available, the spatial representation of the covariates selected by each model. |

The fitting scripts write model-selection summaries and other outputs
under `results/model_01_fit_EMSRR/` and
`results/model_01_fit_EMSRR_SEG_ICES/`. `output_05` uses the generated
spatial cluster raster and selection summaries; run the relevant model
scripts before plotting selected covariates.

### 3. Prepare climate projections and reference-point inputs

| Script | Purpose |
|------------------------------------|------------------------------------|
| `data_00.b_Bias_correct_projections.R` | Bias-correct POLCOMS-ERSEM bottom-temperature projections against AHOI. |
| `data_00.c_Generate_projection_realisations_BVAR_500.R` | Fit spatial BVARs and generate projected bottom-temperature realizations. |
| `data_00.d_Join_historical_and_projection_pattern.R` | Summarize historical and projected bottom temperature on a common spatial footprint. |
| `data_03_ccrefpts_make_lookup_tables.R` | Create EqSim settings and current biological reference-point tables. |
| `model_02_Prepare_EMSRR_model_and_data.R` | Refit selected EMSRRs on common-footprint history and prepare stock and environmental inputs. |

### 4. Condition operating models and simulate reference points

| Script | Purpose |
|------------------------------------|------------------------------------|
| `model_03_ccrefpts_conditionOM.R` | Create operating models for each substock, SR type, and climate setting. |
| `model_04_ccrefpts_runRefptSim.R` | Simulate fishing-mortality targets and calculate candidate reference points. |
| `model_04_ccrefpts_runRefptSim_truncICES.R` | Run the corresponding truncated-ICES reference-point simulations. |

The `functions/` directory contains shared utilities, including spatial
clustering, spatial BVAR, partial-effect, and simulation functions.

## Spatial bottom-temperature patterns

`data_01_prepare_bottomT_by_stubstock.R` calculates annual
bottom-temperature fields by substock and clusters cells by their
temporal patterns. Clustering does not enforce spatial connectivity;
inspect the resulting maps and diagnostics.

The figure below shows the selected spatial cluster map for each
substock. Each color represents a cluster; the number of clusters may
differ among substocks.

![Selected spatial bottom-temperature clusters for the Northwestern,
Southern, and Viking
substocks](docs/images/bottomT-clusters-per-substock.png)

*Spatial clustering of annual bottom-temperature patterns by substock.*

The script writes these products under `data/Env/processed/`:

-   `bottomT_yr_by_cod_substock.rda`
-   `Cluster_bottomT_yr_by_cod_substock.rda`
-   `Spatial_bottomT_yr_by_cod_substock.rda`
-   `SpatialCluster_bottomT_yr_by_cod_substock.rda`

`output_05_plot_spatial_clusters.R` plots the cluster maps and can
compare them with the covariates selected by the `Segreg` and
`Segreg_trunc_ICES` models. For best results, run it from the repository
root after generating the spatial products and model-selection
summaries. The plots can be saved under `output/figs` and the script can
also capture base-R plots for replay in an R session.

## Data and generated products

### External inputs

Input datasets are not included in this repository. A Figshare archive
for the case study has been configured here:
<https://doi.org/10.6084/m9.figshare.34961302> For the automatic
download run the downloader-Script:

``` r
source("data_00_download_inputs.R")
```

The ZIP should contain the following paths directly at its root, without
a wrapping folder:

| Input | Expected location |
|------------------------------------|------------------------------------|
| Monthly bottom T AHOI v23.05 files | `data/Env/processed/AHOI_BottomTemp_raster_1945_2024.tif` |
| ICES rectangle shapefile and companion files | `data/ICES_stat_rect/StatRec_map_Areas_Full_20170124.{shp,shx,dbf,prj}` |
| Substock-area lookup | `data/stock/Area_NSCOD3_table_6a_inshore_offshore.xlsx` |
| Assessment stock objects | `data/stock/cod46a7d20_FLstockObjects_assessment_2026.Rdata` |
| RCP4.5 POLCOMS-ERSEM bottom temperature | `data/Env/POLCOMS-ERSEM/POLCOM_ERSEM_GreaterNorthSea_rcp45_BottomTemp_raster_2006_2099.RData` |
| RCP8.5 POLCOMS-ERSEM bottom temperature | `data/Env/POLCOMS-ERSEM/POLCOM_ERSEM_GreaterNorthSea_rcp85_BottomTemp_raster_2006_2099.RData` |

The original monthly AHOI files available for the North East Atlantic
are available on request, as those were too big to store directly. Those
can be processed with `data_00_Handle_AHOI_data.R`. If you already
downloaded `data/Env/processed/AHOI_BottomTemp_raster_1945_2024.tif`
than there is no need to run `data_00_Handle_AHOI_data.R`.

### Generated data and model outputs

The analysis scripts generate intermediate products under `data/` and
prepared models and simulations under `model/`. Important outputs
include:

-   `data/Cod_Env_data.rda`: joined stock and environmental data
    (`env_stk`).
-   `results/model_01_fit_EMSRR/`:segmented regression model-selection
    summaries, holdout metrics, predictions, partial effects, and plots.
-   `results/model_01_fit_EMSRR_SEG_ICES/`: truncated-ICES segmented
    model-selection outputs.
-   `data/EMSRRs/`: selected EMSRR model objects.
-   `data/Env/processed/bottomT_projections/`: historical and projection
    summaries and common-footprint diagnostics.
-   `model/`: prepared models, conditioned operating models, and
    simulation outputs.

Review intermediate products and their coverage before continuing
through the workflow.

## Running the workflow

Run scripts from the repository root. Several scripts begin with
`rm(list = ls())`, which clears objects in the current R workspace. Use
a fresh R session and save any work beforehand.

``` r
# Historical environmental processing and stock-environment join
source("data_00.a_Handle_AHOI_data.R")
source("data_01_prepare_bottomT_by_stubstock.R")
source("data_02_join_bottomT_with_substock.R")

# Fit/select EMSRRs
source("model_01_fit_EMSRR_BH.R")
source("model_01_fit_EMSRR_Segreg.R")
source("model_01_fit_EMSRR_Segreg_trunc.ICES.SRR.R")

# Optional spatial-cluster plots; requires spatial products and
# model-selection summaries for the selected-covariate comparison
source("output_05_plot_spatial_clusters.R")

# Climate projections and common-footprint summaries
source("data_00.b_Bias_correct_projections.R")
source("data_00.c_Generate_projection_realisations_BVAR_500.R")
source("data_00.d_Join_historical_and_projection_pattern.R")

# Reference-point inputs and operating models
source("data_03_ccrefpts_make_lookup_tables.R")
source("model_02_Prepare_EMSRR_model_and_data.R")

# In a fresh session without FLasher or mse attached:
source("model_03_ccrefpts_conditionOM.R")

# Attach data.table for unqualified rbindlist() calls in model_04.
library(data.table)
source("model_04_ccrefpts_runRefptSim.R")
```

The exact order and inputs for the truncated-ICES workflow may differ
from the standard workflow; review that script's configuration before
running it. Some scripts contain interactive diagnostics and plots.

## Method notes and interpretation

### Historical environment and stock alignment

`data_01` calculates annual bottom-temperature fields and constructs the
substock regions from ICES rectangles and the area lookup. It calculates
substock averages and spatial clusters from complete historical grid
cells. The selected cluster count is based on temporal-pattern
stability; the clustering does not enforce spatial connectivity.

`data_02` extracts recruitment and SSB from the FLCore stock objects. It
lags SSB by recruitment age and shifts environmental series by one year
before joining by year and substock. Check substock labels, lag
conventions, missing values, and retained years before fitting.

### EMSRR model selection

The BH and segmented model scripts use bootstrap out-of-bag performance
to select a candidate environmental variable and MARS form for each
substock. The latest 10 distinct calendar years are reserved for holdout
evaluation; earlier complete observations are used for training. The
holdout scores are for training-only model selection and do not
constitute an independent evaluation of the final full-data refits.

`model_02` uses the selected variable and MARS form to refit models
using common-footprint historical environmental values for the original
fitting years. Model selection and holdout metrics are not recalculated
on that common footprint.

### Projections and common footprint

The projection workflow bias-corrects POLCOMS-ERSEM projections,
generates stochastic realizations with spatial BVARs, and summarizes
historical and projected annual bottom temperature on a common spatial
footprint. The common footprint retains cells with complete historical
annual values and finite monthly values in every realization of both
scenarios. Substock means are unweighted; cluster means use fixed
`raster::area()` weights.

History and projections remain separate in the summary output; they are
not spliced or lagged there. Inspect coverage, historical comparisons,
and overlap diagnostics. Restricting the footprint can change
environmental covariates.

### Reference-point simulations

`model_03` conditions operating models for the three substocks and
configured SR types and climate settings. `model_04` simulates target
fishing mortalities and derives candidate reference-point values from
the simulation outputs.

These are candidate calculations, not validated management reference
points. Review scenario construction, projection coverage, stochastic
assumptions, units, and implementation choices.

## Requirements

There is no package manifest or lockfile, so package versions are not
pinned. Dependencies include:

-   Spatial and environmental processing: `raster`, `terra`, `RNetCDF`,
    `sf`, `readxl`, `marmalaid`, `maps`, `maptools`, `rgeos`,
    `rnaturalearth`, and `mclust`.
-   Projection and parallel processing: `BVAR`, `coda`, `corrplot`,
    `tsDyn`, `Rcpp`, `Rcpp2doParallel`, `foreach`, `doParallel`, and
    `doSNOW`.
-   Model fitting and simulation: `earth`, `minqa`, `collapse`,
    `ggplot2`, `DescTools`, `stringr`, `dplyr`, `data.table`, `FLCore`,
    `FLasherEMSRR`, and `mseEMSRR`.

Other dependencies may be used by individual scripts and helper
functions. Some older spatial packages may need a compatible R
environment. Install and validate dependencies in a project-specific
library, and record the R and package versions used.

For `model_03` and `model_04`, start a fresh R session without the
original `FLasher` or `mse` packages attached; they may conflict with
classes from `FLasherEMSRR` and `mseEMSRR`.

Processing can require substantial time, memory, and disk space. In
particular, `data_00.a` creates a 10-worker cluster, `data_00.c` is
configured for 500 realizations per RCP, and `model_04` currently
requests 25 workers. Adjust worker settings to match available
resources.

## Interpretation and citation

This is a research case study, not an official stock assessment or
management recommendation. Cite the source data, assessment products,
model methods, and software used in any resulting report or publication.
Interpret results only after the data, model fits, projections, and
reference-point calculations have been reviewed and validated.
