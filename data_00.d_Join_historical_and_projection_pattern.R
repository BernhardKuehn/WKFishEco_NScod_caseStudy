# Annual bottom-temperature summaries for the historical and projected fields.
# Run data_01_prepare_bottomT_by_stubstock.R and data_00.c first.
# Run from the repository root; the inputs are local and are not versioned.

source("functions/prepare_bottomT_projections.R")

bottomT_projection_results <- prepare_bottomT_projections()
