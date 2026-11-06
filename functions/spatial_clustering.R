# ----------------------------------------------- # 
# functions to perform a quick spatial clustering
# ----------------------------------------------- #

# Cluster a longitude/latitude environmental raster by temporal evolution.
#
# Rows are grid cells and layers are ordered time steps. Clustering does not
# enforce spatial connectivity; inspect the selected cluster map.
spatial_kmeans <- function(
    sst,
    k_values = 2:10,
    time = names(sst),
    nstart = 30,
    seed = 123
) {
  if (!requireNamespace("raster", quietly = TRUE)) {
    stop("Install the 'raster' package.")
  }
  if (!requireNamespace("mclust", quietly = TRUE)) {
    stop("Install the 'mclust' package.")
  }
  
  sst <- raster::brick(sst)
  
  if (!isTRUE(raster::isLonLat(sst))) {
    stop("sst must use longitude/latitude coordinates.")
  }
  if (length(time) != raster::nlayers(sst) ||
      anyNA(time) || anyDuplicated(time)) {
    stop("time must contain one unique label per raster layer, in order.")
  }
  if (raster::nlayers(sst) < 6L) {
    stop("At least six time steps are needed for the stability diagnostic.")
  }
  if (!length(k_values) || any(!is.finite(k_values)) ||
      any(k_values < 2) || any(k_values != floor(k_values))) {
    stop("k_values must contain integers of at least 2.")
  }
  if (length(nstart) != 1L || !is.finite(nstart) ||
      nstart < 1 || nstart != floor(nstart)) {
    stop("nstart must be a positive integer.")
  }
  
  k_values <- sort(unique(k_values))
  
  # Keep a fixed set of ocean cells with complete time series.
  values <- as.matrix(raster::getValues(sst))
  area <- raster::getValues(raster::area(sst[[1]]))
  
  complete <- rowSums(is.finite(values)) == ncol(values) &
    is.finite(area) & area > 0
  
  cell_sd <- rep(NA_real_, nrow(values))
  cell_sd[complete] <- apply(
    values[complete, , drop = FALSE], 1, stats::sd
  )
  cells <- which(complete & is.finite(cell_sd) & cell_sd > 0)
  
  if (length(cells) <= max(k_values)) {
    stop("Too few complete, variable cells for the requested k_values.")
  }
  
  X <- values[cells, , drop = FALSE]
  cell_area <- area[cells]
  
  standardize_cells <- function(x) {
    cell_mean <- rowMeans(x)
    cell_sd <- apply(x, 1, stats::sd)
    
    if (any(!is.finite(cell_sd) | cell_sd == 0)) {
      stop("A cell has zero or invalid variability in a time block.")
    }
    
    sweep(sweep(x, 1, cell_mean, "-"), 1, cell_sd, "/")
  }
  
  Z <- standardize_cells(X)
  
  # Fit to contiguous halves of the record for a stability comparison.
  n_time <- ncol(X)
  split_at <- floor(n_time / 2)
  Z_first <- standardize_cells(X[, seq_len(split_at), drop = FALSE])
  Z_second <- standardize_cells(
    X[, seq.int(split_at + 1, n_time), drop = FALSE]
  )
  
  results <- setNames(
    vector("list", length(k_values)),
    paste0("k", k_values)
  )
  
  for (k in k_values) {
    set.seed(seed + k)
    full <- stats::kmeans(Z, centers = k, nstart = nstart)
    
    set.seed(seed + k)
    first <- stats::kmeans(Z_first, centers = k, nstart = nstart)
    
    set.seed(seed + k)
    second <- stats::kmeans(Z_second, centers = k, nstart = nstart)
    
    stability <- mclust::adjustedRandIndex(
      first$cluster, second$cluster
    )
    
    # Raster cells excluded from clustering remain NA.
    map_values <- rep(NA_integer_, raster::ncell(sst))
    map_values[cells] <- full$cluster
    cluster_map <- raster::setValues(
      raster::raster(sst[[1]]), map_values
    )
    names(cluster_map) <- "cluster"
    
    # Fixed area weights applied to the original, unstandardized SST.
    cluster_means <- sapply(seq_len(k), function(j) {
      idx <- which(full$cluster == j)
      weights <- cell_area[idx] / sum(cell_area[idx])
      drop(crossprod(weights, X[idx, , drop = FALSE]))
    })
    cluster_means <- matrix(cluster_means, nrow = n_time, ncol = k)
    colnames(cluster_means) <- paste0("Cluster", seq_len(k))
    
    cluster_area <- vapply(seq_len(k), function(j) {
      sum(cell_area[full$cluster == j])
    }, numeric(1))
    
    results[[paste0("k", k)]] <- list(
      cluster_map = cluster_map,
      timeseries = data.frame(time = time, cluster_means),
      cell = cells,
      membership = full$cluster,
      cluster_area_km2 = cluster_area,
      withinss = full$tot.withinss,
      between_fraction = full$betweenss / full$totss,
      block_stability_ari = stability
    )
  }
  
  diagnostics <- data.frame(
    k = k_values,
    withinss = vapply(results, `[[`, numeric(1), "withinss"),
    between_fraction = vapply(
      results, `[[`, numeric(1), "between_fraction"
    ),
    block_stability_ari = vapply(
      results, `[[`, numeric(1), "block_stability_ari"
    )
  )
  
  if (any(!is.finite(diagnostics$block_stability_ari))) {
    stop("A cluster stability score is missing or non-finite.")
  }
  
  # Highest block stability; prefer the smaller k in an exact tie.
  best <- order(-diagnostics$block_stability_ari, diagnostics$k)[1]
  selected_k <- diagnostics$k[best]
  
  list(
    selected_k = selected_k,
    selected = results[[paste0("k", selected_k)]],
    diagnostics = diagnostics,
    results = results
  )
}

choose_cluster_count <- function(cluster_output) {
  diagnostics <- cluster_output$diagnostics
  
  if (any(!is.finite(diagnostics$block_stability_ari))) {
    stop("Cluster stability contains missing or non-finite values.")
  }
  
  # Prefer the smaller k if stability scores are identical.
  best <- order(
    -diagnostics$block_stability_ari,
    diagnostics$k
  )[1]
  
  k <- diagnostics$k[best]
  
  list(
    k = k,
    stability = diagnostics$block_stability_ari[best],
    result = cluster_output$results[[paste0("k", k)]]
  )
}