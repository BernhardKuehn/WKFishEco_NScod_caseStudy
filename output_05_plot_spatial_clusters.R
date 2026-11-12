# Plot all spatial clusters, then highlight clusters selected by each model
# Run from the repository root after the model scripts have produced their
# variable-selection summary CSV files.

rm(list = ls())

library(terra)
library(raster)
library(wesanderson)

load(file.path(
  "./data/Env/processed",
  "SpatialCluster_bottomT_yr_by_cod_substock.rda"
))

substocks <- c("Northwestern", "Southern", "Viking")
cluster_maps <- cluster_bottomT_substk_spatial

if (!all(substocks %in% names(cluster_maps))) {
  stop("Expected cluster maps for: ", paste(substocks, collapse = ", "))
}

# Use one palette across both graphics so cluster colors remain consistent.
base_palette <- c(
  wesanderson::wes_palette("GrandBudapest1", n = 4),
  wesanderson::wes_palette("Darjeeling1", n = 5)
)

n_clusters <- vapply(substocks, function(substock) {
  max_cluster <- raster::cellStats(
    cluster_maps[[substock]],
    stat = "max",
    na.rm = TRUE
  )
  
  if (!is.finite(max_cluster)) {
    stop("Could not determine the number of clusters for ", substock, ".")
  }
  
  as.integer(max_cluster)
}, integer(1))

pal <- grDevices::colorRampPalette(base_palette)(max(n_clusters))

# store plots
plots_out <- list()

# Graphic 1: all clusters for each substock.
old_par <- par(no.readonly = TRUE)
par(mfrow = c(1, 3), mar = c(3, 3, 3, 4))

for (substock in substocks) {
  k <- n_clusters[[substock]]
  
  xlim = extent(cluster_maps[[substock]])[1:2] + c(-2,2)
  ylim = extent(cluster_maps[[substock]])[3:4] + c(-2,2)
  
  raster::image(
    cluster_maps[[substock]],
    xlim = xlim,
    ylim = ylim,
    main = substock,
    col = pal[seq_len(k)],
    breaks = seq(0.5, k + 0.5, by = 1)
  )
  maps::map(add = TRUE)
  legend(
    "topright",
    legend = seq_len(k),
    fill = pal[seq_len(k)],
    title = "Cluster",box.col = "white",
    cex = 0.8
  )
  box()
}

# record plot
plots_out$all_clusters <- grDevices::recordPlot()
par(old_par)

# Each CSV contains one selected covariate per substock.
selection_files <- c(
  Segreg = "./results/model_01_fit_EMSRR/SEG_variable_selection_summary.csv",
  BH = "./results/model_01_fit_EMSRR/BH_variable_selection_summary.csv",
  Segreg_trunc_ICES =
    "./results/model_01_fit_EMSRR_SEG_ICES/SEG_variable_selection_summary.csv"
)

selection_tables <- lapply(selection_files, function(path) {
  if (!file.exists(path)) {
    stop("Selection summary not found: ", path)
  }
  
  selection <- read.csv(path, stringsAsFactors = FALSE)
  required_columns <- c("substock", "variable", "selected")
  
  if (!all(required_columns %in% names(selection))) {
    stop(
      "Expected columns substock, variable, and selected in ",
      path
    )
  }
  
  is_selected <- tolower(trimws(as.character(selection$selected))) %in%
    c("true", "t", "1")
  selected <- selection[is_selected, , drop = FALSE]
  
  if (anyDuplicated(selected$substock) ||
      !setequal(selected$substock, substocks)) {
    stop(
      "Expected exactly one selected covariate per substock in ",
      path
    )
  }
  
  selected[match(substocks, selected$substock), , drop = FALSE]
})

# Graphic 2: show only the selected cluster cells, for each model and substock.
par(
  mfrow = c(length(selection_files), length(substocks)),
  mar = c(2, 2, 4, 1)
)

for (model_name in names(selection_files)) {
  selection <- selection_tables[[model_name]]
  
  for (substock in substocks) {
    variable <- selection$variable[selection$substock == substock]
    
    # Model covariates such as bottomT_mean_Southern have no cluster map.
    if (!grepl("Cluster[0-9]+$", variable)) {
      k <- n_clusters[[substock]]
      all = cluster_maps[[substock]] != 0
      
      raster::image(
        all,
        main = paste(model_name, substock, "\n", variable, "(whole area shown)"),
        col = "red",
        breaks = seq(0,1, by = 1),
        axes = FALSE,
        box = FALSE
      )
      maps::map(add = TRUE)
      next
    }
    
    cluster_id <- as.integer(sub(".*Cluster([0-9]+)$", "\\1", variable))
    
    if (cluster_id < 1L || cluster_id > n_clusters[[substock]]) {
      stop(
        "Selected cluster ", cluster_id, " for ", substock,
        " is not present in its spatial cluster map."
      )
    }
    
    selected_map <- cluster_maps[[substock]]
    map_values <- raster::getValues(selected_map)
    map_values[!is.na(map_values) & map_values != cluster_id] <- NA
    selected_map <- raster::setValues(selected_map, map_values)
    
    raster::image(
      selected_map,
      main = paste(model_name, substock, "\n", variable),
      col = pal,
      breaks = seq(0.5, max(n_clusters) + 0.5, by = 1),
      legend = FALSE,
      axes = FALSE,
      box = FALSE
    )
    maps::map(add = TRUE)
  }
}

# record plot
plots_out$selected_by_model <- grDevices::recordPlot()
par(old_par)

### save plots ###

png(file.path("./output/figs","BottomT_clusters_per_substock.png"),res = 400,units = "in",width = 9,height = 6)
grDevices::replayPlot(plots_out$all_clusters)
dev.off()

png(file.path("./output/figs","BottomT_clusters_per_substock_selectedModel.png"),res = 400,units = "in",width = 7,height = 6)
grDevices::replayPlot(plots_out$selected_by_model)
dev.off()


