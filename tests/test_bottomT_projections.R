stopifnot(requireNamespace("terra"), requireNamespace("raster"))

test_bottomT_projections <- function() {
  script <- normalizePath("data_00.d_Join_historical_and_projection_pattern.R")
  original_dir <- getwd()
  root <- tempfile("bottomT-test-")
  dir.create(root)
  on.exit({setwd(original_dir); unlink(root, recursive = TRUE)}, add = TRUE)
  processed <- file.path(root, "data/Env/processed")
  projections <- file.path(root, "data/Env/POLCOMS-ERSEM_biasRM_BVARsims/bottomT")
  output <- file.path(processed, "bottomT_projections")
  dir.create(processed, recursive = TRUE)
  dir.create(projections, recursive = TRUE)
  setwd(root)
  run_script <- function() {
    env <- new.env(parent = globalenv())
    source(script, local = env)
    env
  }

  h <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2,
                   ymin = 50, ymax = 52, crs = "EPSG:4326", nlyrs = 3)
  terra::values(h) <- cbind(c(1, 2, 10, 20), c(2, 3, 11, 21),
                           c(3, 4, 12, 22))
  names(h) <- as.character(2006:2008)
  m <- terra::rast(h[[1]])
  terra::values(m) <- c(1, 1, 2, 2)
  hn <- terra::crop(h, terra::ext(1, 2, 50, 52))
  mn <- terra::crop(m, terra::ext(1, 2, 50, 52))
  bottomT_yr_by_substock_spatial <- list(Southern = raster::brick(h),
                                         Northwestern = raster::brick(hn))
  cluster_bottomT_substk_spatial <- list(Southern = raster::raster(m),
                                          Northwestern = raster::raster(mn))
  save(bottomT_yr_by_substock_spatial,
       file = file.path(processed, "Spatial_bottomT_yr_by_cod_substock.rda"))
  save(cluster_bottomT_substk_spatial,
       file = file.path(processed, "SpatialCluster_bottomT_yr_by_cod_substock.rda"))

  area <- raster::getValues(raster::area(raster::raster(h[[1]])))
  vals <- terra::values(h)
  bottomT_substk <- data.frame(year = 2006:2008,
                               bottomT_Southern = colMeans(vals),
                               bottomT_Northwestern = colMeans(vals[c(2, 4), ]))
  cluster_bottomT_substk <- data.frame(
    year = 2006:2008,
    bottomT_Southern_Cluster1 = drop(crossprod(area[1:2] / sum(area[1:2]), vals[1:2, ])),
    bottomT_Southern_Cluster2 = drop(crossprod(area[3:4] / sum(area[3:4]), vals[3:4, ])),
    bottomT_Northwestern_Cluster1 = vals[2, ],
    bottomT_Northwestern_Cluster2 = vals[4, ]
  )
  save(bottomT_substk, file = file.path(processed, "bottomT_yr_by_cod_substock.rda"))
  save(cluster_bottomT_substk, file = file.path(processed, "Cluster_bottomT_yr_by_cod_substock.rda"))

  dates <- seq(as.Date("2007-01-15"), as.Date("2009-12-15"), by = "month")
  for (scenario in c("rcp45", "rcp85")) {
    for (id in c(1, 10)) {
      p <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2,
                       ymin = 50, ymax = 52, crs = "EPSG:4326", nlyrs = length(dates))
      base <- c(10, 20, 30, 50) + id + if (scenario == "rcp85") 5 else 0
      projection_data <- sapply(seq_along(dates), function(i) base + (i - 1) %/% 12)
      if (scenario == "rcp85" && id == 10) projection_data[1, 5] <- NA
      terra::values(p) <- projection_data
      names(p) <- as.character(dates)
      terra::time(p) <- dates
      terra::writeRaster(p, file.path(projections, paste0("POLCOM_ERSEM_", scenario,
                                                         "_bottomT_", id, ".tiff")))
    }
  }
  run <- run_script()
  result <- run$bottomT_projection_results
  stopifnot(identical(names(result$projections$rcp45), c("1", "10")),
            identical(result$historical_common$year, 2006:2008),
            identical(result$projections$rcp45[["1"]]$year, 2007:2009),
            max(abs(result$historical_comparison$full_minus_saved)) < 1e-10,
            any(abs(result$historical_comparison$common_minus_saved) > 0),
            result$coverage$common_cells[result$coverage$substock == "Southern" &
                                           result$coverage$group == "mean"] == 3,
            is.na(raster::getValues(result$common_masks$Southern)[1]),
            all.equal(result$historical_common$bottomT_Southern,
                      unname(colMeans(vals[2:4, ]))),
            all.equal(result$projections$rcp45[["1"]]$bottomT_Southern,
                      c(mean(c(21, 31, 51)), mean(c(22, 32, 52)), mean(c(23, 33, 53)))),
            all.equal(result$projections$rcp45[["1"]]$bottomT_Southern_Cluster2,
                      c(41, 42, 43)),
            all.equal(result$projections$rcp45[["1"]]$bottomT_Northwestern,
                      c(36, 37, 38)),
            all(result$overlap$year %in% 2007:2008),
            nrow(result$overlap) == 2 * 6 * 2,
            file.exists(file.path(output, "historical_projection_comparison.pdf")),
            file.exists(file.path(output, "historical_projection_overlap.csv")),
            identical(readRDS(file.path(output, "bottomT_historical_projection_summaries.rds"))$projections,
                      result$projections))

  # Cluster 1 spans latitudes in this variant; its area weights are not equal.
  m2 <- m
  terra::values(m2) <- c(1, 2, 1, 2)
  cluster_bottomT_substk_spatial <- list(Southern = raster::raster(m2),
                                          Northwestern = raster::raster(mn))
  save(cluster_bottomT_substk_spatial,
       file = file.path(processed, "SpatialCluster_bottomT_yr_by_cod_substock.rda"))
  cluster_bottomT_substk$bottomT_Southern_Cluster1 <-
    drop(crossprod(area[c(1, 3)] / sum(area[c(1, 3)]), vals[c(1, 3), ]))
  cluster_bottomT_substk$bottomT_Southern_Cluster2 <-
    drop(crossprod(area[c(2, 4)] / sum(area[c(2, 4)]), vals[c(2, 4), ]))
  save(cluster_bottomT_substk, file = file.path(processed, "Cluster_bottomT_yr_by_cod_substock.rda"))
  run2 <- run_script()
  result2 <- run2$bottomT_projection_results
  expected <- (area[2] * 21 + area[4] * 51) / sum(area[c(2, 4)])
  stopifnot(abs(result2$projections$rcp45[["1"]]$bottomT_Southern_Cluster2[1] -
                  expected) < 1e-5,
            max(abs(result2$historical_comparison$full_minus_saved)) < 1e-10)

  p <- terra::rast(file.path(projections, "POLCOM_ERSEM_rcp45_bottomT_1.tiff"))
  stopifnot(inherits(try(run$get.proj.time(p[[-1]], "incomplete"), silent = TRUE), "try-error"))
  stopifnot(inherits(try(run$get.proj.time(p[[-(13:24)]], "missing year"), silent = TRUE),
                     "try-error"))
  names(p) <- gsub("-", ".", paste0("X", dates))
  terra::time(p) <- NULL
  stopifnot(identical(format(run$get.proj.time(p, "layer names"), "%Y-%m"),
                      format(dates, "%Y-%m")))

  # A cluster with no projection coverage must not produce a time series.
  no_cluster <- terra::rast(file.path(projections, "POLCOM_ERSEM_rcp45_bottomT_1.tiff"))
  no_cluster_values <- terra::values(no_cluster)
  no_cluster_values[3:4, ] <- NA
  terra::values(no_cluster) <- no_cluster_values
  missing_file <- file.path(projections, "POLCOM_ERSEM_rcp45_bottomT_2.tiff")
  terra::writeRaster(no_cluster, missing_file)
  err <- try(run_script(), silent = TRUE)
  stopifnot(inherits(err, "try-error"), grepl("no cells on the common footprint", err))
  unlink(missing_file)

  # A shifted projection grid must be rejected rather than silently resampled.
  bad <- terra::rast(nrows = 2, ncols = 2, xmin = 0.1, xmax = 2.1,
                     ymin = 50, ymax = 52, crs = "EPSG:4326", nlyrs = length(dates))
  terra::values(bad) <- 1
  terra::time(bad) <- dates
  bad_file <- file.path(projections, "POLCOM_ERSEM_rcp45_bottomT_2.tiff")
  terra::writeRaster(bad, bad_file)
  err <- try(run_script(), silent = TRUE)
  stopifnot(inherits(err, "try-error"), grepl("grid origin differs", err))
  unlink(bad_file)

  # Old direct SpatRaster saves cannot be restored; report how to repair them.
  bottomT_yr_by_substock_spatial$Southern <- h
  save(bottomT_yr_by_substock_spatial,
       file = file.path(processed, "Spatial_bottomT_yr_by_cod_substock.rda"))
  err <- try(run_script(), silent = TRUE)
  stopifnot(inherits(err, "try-error"), grepl("Regenerate it with the updated data_01", err))
  cat("Bottom-temperature projection tests passed.\n")
}

test_bottomT_projections()
