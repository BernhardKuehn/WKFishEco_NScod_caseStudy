source("functions/prepare_bottomT_projections.R")
stopifnot(requireNamespace("terra"), requireNamespace("raster"))

test_bottomT_projections <- function() {
  root <- tempfile("bottomT-test-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  processed <- file.path(root, "processed")
  projections <- file.path(root, "projections")
  output <- file.path(root, "output")
  dir.create(processed)
  dir.create(projections)

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
  result <- prepare_bottomT_projections(processed, projections, output)
  stopifnot(identical(names(result$projections$rcp45), c("1", "10")),
            identical(result$historical_common$year, 2006:2008),
            identical(result$projections$rcp45[["1"]]$year, 2007:2009),
            max(abs(result$historical_comparison$full_minus_saved)) < 1e-10,
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
  result2 <- prepare_bottomT_projections(processed, projections, output)
  expected <- (area[2] * 21 + area[4] * 51) / sum(area[c(2, 4)])
  stopifnot(abs(result2$projections$rcp45[["1"]]$bottomT_Southern_Cluster2[1] -
                  expected) < 1e-5,
            max(abs(result2$historical_comparison$full_minus_saved)) < 1e-10)

  bad <- terra::rast(nrows = 2, ncols = 2, xmin = 0.1, xmax = 2.1,
                     ymin = 50, ymax = 52, crs = "EPSG:4326")
  stopifnot(inherits(try(bottomT_same_grid(h, bad, "shifted"), silent = TRUE), "try-error"))
  p <- terra::rast(file.path(projections, "POLCOM_ERSEM_rcp45_bottomT_1.tiff"))
  stopifnot(inherits(try(bottomT_months(p[[-1]], "incomplete"), silent = TRUE), "try-error"))
  stopifnot(inherits(try(bottomT_months(p[[-(13:24)]], "missing year"), silent = TRUE),
                     "try-error"))
  names(p) <- gsub("-", ".", paste0("X", dates))
  terra::time(p) <- NULL
  stopifnot(identical(bottomT_months(p, "layer names")$months, format(dates, "%Y-%m")))
  overlap <- bottomT_overlap(data.frame(year = 2007L, bottomT_Southern = 3),
                             list(rcp45 = list("1" = data.frame(year = 2007L,
                                                                bottomT_Southern = 4))))
  stopifnot(nrow(overlap) == 1, overlap$median_minus_historical == 1)
  bottomT_plot_diagnostics(file.path(root, "single_year.pdf"),
                           data.frame(year = 2007L, bottomT_Southern = 3),
                           data.frame(year = 2007L, variable = "bottomT_Southern", saved = 3),
                           list(rcp45 = list("1" = data.frame(year = 2007L,
                                                              bottomT_Southern = 4))))

  # Old direct SpatRaster saves cannot be restored; report how to repair them.
  bottomT_yr_by_substock_spatial$Southern <- h
  save(bottomT_yr_by_substock_spatial,
       file = file.path(processed, "Spatial_bottomT_yr_by_cod_substock.rda"))
  err <- try(prepare_bottomT_projections(processed, projections, output), silent = TRUE)
  stopifnot(inherits(err, "try-error"), grepl("Regenerate the spatial product", err))
  cat("Bottom-temperature projection tests passed.\n")
}

test_bottomT_projections()
