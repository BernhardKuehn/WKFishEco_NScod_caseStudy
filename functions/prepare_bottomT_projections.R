# Prepare annual substock and fixed-cluster means on a common spatial footprint.

read_bottomT_object <- function(path, object) {
  if (!file.exists(path)) stop("Missing input: ", path)
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (!exists(object, envir = env, inherits = FALSE)) {
    stop("Expected object '", object, "' in ", path)
  }
  get(object, envir = env, inherits = FALSE)
}

bottomT_years <- function(x, context) {
  labels <- sub("^X", "", names(x))
  if (any(!grepl("^[0-9]{4}$", labels))) {
    stop(context, ": historical raster layers must be named by year.")
  }
  years <- as.integer(labels)
  if (anyDuplicated(years) || is.unsorted(years, strictly = TRUE)) {
    stop(context, ": historical years must be unique and increasing.")
  }
  years
}

bottomT_months <- function(x, context) {
  dates <- terra::time(x)
  if (length(dates) != terra::nlyr(x) || anyNA(dates)) {
    labels <- gsub("\\.", "-", sub("^X", "", names(x)))
    if (any(!grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", labels))) {
      stop(context, ": missing monthly dates in raster time and layer names.")
    }
    dates <- as.Date(labels, format = "%Y-%m-%d")
  } else {
    dates <- as.Date(dates)
  }
  if (anyNA(dates)) stop(context, ": invalid monthly dates.")
  months <- format(dates, "%Y-%m")
  if (anyDuplicated(months) || is.unsorted(months, strictly = TRUE)) {
    stop(context, ": months must be unique and increasing.")
  }
  years <- as.integer(format(dates, "%Y"))
  if (any(table(years) != 12L)) {
    stop(context, ": every projection year must contain all 12 months.")
  }
  month_number <- years * 12L + as.integer(format(dates, "%m"))
  if (any(diff(month_number) != 1L)) {
    stop(context, ": projection months must be consecutive, without missing years.")
  }
  list(dates = dates, months = months, years = years)
}

bottomT_same_grid <- function(x, y, context) {
  same <- terra::compareGeom(x[[1]], y[[1]], crs = TRUE, ext = FALSE,
                             rowcol = FALSE, res = TRUE, stopOnError = FALSE)
  # Equal resolution alone does not ensure the grid cells have the same origin.
  offset <- (terra::origin(x) - terra::origin(y)) / terra::res(x)
  if (!isTRUE(same) || any(abs(offset - round(offset)) > 1e-6)) {
    stop(context, ": CRS, resolution or grid origin differs from history; ",
         "resample explicitly before running this script.")
  }
}

bottomT_projection_files <- function(directory, scenarios) {
  if (!dir.exists(directory)) stop("Missing projection directory: ", directory)
  files <- list.files(directory, pattern = "\\.tiff?$", full.names = TRUE,
                      ignore.case = TRUE)
  pattern <- "^POLCOM_ERSEM_(rcp45|rcp85)_bottomT_([0-9]+)\\.tiff?$"
  if (any(!grepl(pattern, basename(files), ignore.case = TRUE))) {
    stop("Unexpected projection TIFF filename(s): ",
         paste(basename(files[!grepl(pattern, basename(files),
                                      ignore.case = TRUE)]), collapse = ", "))
  }
  result <- setNames(vector("list", length(scenarios)), scenarios)
  for (scenario in scenarios) {
    selected <- files[grepl(paste0("^POLCOM_ERSEM_", scenario,
                                   "_bottomT_[0-9]+\\.tiff?$"),
                            basename(files), ignore.case = TRUE)]
    if (!length(selected)) stop("No projection files found for ", scenario)
    ids <- as.integer(sub(pattern, "\\2", basename(selected), ignore.case = TRUE))
    if (anyNA(ids) || anyDuplicated(ids)) {
      stop("Invalid or duplicate realization numbers for ", scenario)
    }
    selected <- selected[order(ids)]
    ids <- sort(ids)
    result[[scenario]] <- setNames(selected, as.character(ids))
  }
  result
}

bottomT_summarize <- function(values, years, regions, common = TRUE) {
  out <- data.frame(year = years)
  for (region in names(regions)) {
    r <- regions[[region]]
    cells <- if (common) r$common else r$historical
    data <- values[[region]]
    if (nrow(data) != length(r$historical)) {
      stop(region, ": unexpected number of raster cells.")
    }
    selected <- which(cells)
    if (!length(selected)) stop(region, ": no cells available for the mean.")
    out[[paste0("bottomT_", region)]] <-
      colMeans(data[selected, , drop = FALSE], na.rm = !common)
    for (id in r$cluster_ids) {
      selected <- which(cells & !is.na(r$cluster) & r$cluster == id)
      if (!length(selected)) stop(region, ": no cells remain in Cluster", id)
      weights <- r$area[selected] / sum(r$area[selected])
      out[[paste0("bottomT_", region, "_Cluster", id)]] <-
        drop(crossprod(weights, data[selected, , drop = FALSE]))
    }
  }
  out
}

bottomT_compare <- function(original, full, common) {
  columns <- setdiff(names(full), "year")
  if (!all(columns %in% names(original))) {
    stop("Saved historical tables lack columns: ",
         paste(setdiff(columns, names(original)), collapse = ", "))
  }
  if (anyDuplicated(original$year)) stop("Saved historical years are duplicated.")
  do.call(rbind, lapply(columns, function(column) {
    years <- sort(unique(c(original$year, full$year, common$year)))
    saved <- original[[column]][match(years, original$year)]
    uncropped <- full[[column]][match(years, full$year)]
    cropped <- common[[column]][match(years, common$year)]
    data.frame(year = years, variable = column, saved = saved,
               recomputed_full = uncropped, recomputed_common = cropped,
               full_minus_saved = uncropped - saved,
               common_minus_saved = cropped - saved)
  }))
}

bottomT_overlap <- function(historical, projections) {
  rows <- list()
  for (scenario in names(projections)) {
    runs <- projections[[scenario]]
    common_years <- intersect(historical$year, runs[[1]]$year)
    if (!length(common_years)) {
      warning("Historical and ", scenario, " time series have no overlapping years.")
      next
    }
    for (variable in setdiff(names(historical), "year")) {
      indices <- match(common_years, runs[[1]]$year)
      values <- vapply(runs, function(x) x[[variable]][indices],
                       numeric(length(common_years)))
      values <- matrix(values, nrow = length(common_years))
      quantiles <- t(apply(values, 1, stats::quantile, probs = c(0.1, 0.5, 0.9),
                           names = FALSE))
      observed <- historical[[variable]][match(common_years, historical$year)]
      rows[[length(rows) + 1L]] <- data.frame(
        scenario = scenario, variable = variable, year = common_years,
        historical = observed, projection_p10 = quantiles[, 1],
        projection_median = quantiles[, 2], projection_p90 = quantiles[, 3],
        median_minus_historical = quantiles[, 2] - observed
      )
    }
  }
  if (!length(rows)) {
    return(data.frame(scenario = character(), variable = character(), year = integer(),
                      historical = numeric(), projection_p10 = numeric(),
                      projection_median = numeric(), projection_p90 = numeric(),
                      median_minus_historical = numeric()))
  }
  do.call(rbind, rows)
}

bottomT_plot_diagnostics <- function(path, historical, comparison, projections) {
  grDevices::pdf(path, width = 10, height = 6)
  on.exit(grDevices::dev.off(), add = TRUE)
  colors <- c(rcp45 = "#1B7837", rcp85 = "#B2182B")
  for (variable in setdiff(names(historical), "year")) {
    ranges <- lapply(projections, function(runs) {
      years <- runs[[1]]$year
      values <- matrix(vapply(runs, function(x) x[[variable]], numeric(length(years))),
                       nrow = length(years))
      qs <- t(apply(values, 1, stats::quantile, probs = c(0.1, 0.5, 0.9),
                    names = FALSE))
      data.frame(year = years, low = qs[, 1], median = qs[, 2], high = qs[, 3])
    })
    saved <- comparison[comparison$variable == variable, ]
    xlim <- range(c(historical$year, unlist(lapply(ranges, `[[`, "year"))))
    ylim <- range(c(historical[[variable]], saved$saved,
                    unlist(lapply(ranges, function(x) c(x$low, x$high)))),
                  finite = TRUE)
    graphics::plot(historical$year, historical[[variable]], type = "n",
                   xlim = xlim, ylim = ylim, xlab = "Year",
                   ylab = "Annual bottom temperature", main = variable)
    for (scenario in names(ranges)) {
      x <- ranges[[scenario]]
      graphics::polygon(c(x$year, rev(x$year)), c(x$low, rev(x$high)),
                        border = NA, col = grDevices::adjustcolor(colors[scenario],
                                                                 alpha.f = 0.16))
      graphics::lines(x$year, x$median, col = colors[scenario], lwd = 2)
    }
    graphics::lines(saved$year, saved$saved, col = "grey40", lty = 2)
    graphics::lines(historical$year, historical[[variable]], lwd = 2)
    graphics::legend("topleft", bty = "n", lwd = 2, lty = c(1, 2, 1, 1),
                     col = c("black", "grey40", colors),
                     legend = c("Historical, common footprint", "Historical, saved",
                                "RCP4.5 median (10-90%)", "RCP8.5 median (10-90%)"))
  }
}

prepare_bottomT_projections <- function(
    processed_dir = "./data/Env/processed",
    projection_dir = "./data/Env/POLCOMS-ERSEM_biasRM_BVARsims/bottomT",
    output_dir = "./data/Env/processed/bottomT_projections") {
  for (package in c("terra", "raster")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Install the '", package, "' package.")
    }
  }
  history <- read_bottomT_object(file.path(processed_dir, "Spatial_bottomT_yr_by_cod_substock.rda"),
                                 "bottomT_yr_by_substock_spatial")
  maps <- read_bottomT_object(file.path(processed_dir, "SpatialCluster_bottomT_yr_by_cod_substock.rda"),
                              "cluster_bottomT_substk_spatial")
  saved_mean <- read_bottomT_object(file.path(processed_dir, "bottomT_yr_by_cod_substock.rda"),
                                    "bottomT_substk")
  saved_clusters <- read_bottomT_object(file.path(processed_dir, "Cluster_bottomT_yr_by_cod_substock.rda"),
                                        "cluster_bottomT_substk")
  if (!is.list(history) || !length(history) || is.null(names(history)) ||
      any(!nzchar(names(history))) || anyDuplicated(names(history)) ||
      !is.list(maps) || is.null(names(maps)) || anyDuplicated(names(maps)) ||
      !setequal(names(history), names(maps))) {
    stop("Historical fields and cluster maps must have the same unique substock names.")
  }
  if (!is.data.frame(saved_mean) || !is.data.frame(saved_clusters) ||
      !"year" %in% names(saved_mean) || !"year" %in% names(saved_clusters) ||
      anyNA(saved_mean$year) || anyNA(saved_clusters$year) ||
      anyDuplicated(saved_mean$year) || anyDuplicated(saved_clusters$year)) {
    stop("Saved historical mean and cluster tables must have unique year columns.")
  }
  original <- merge(saved_mean, saved_clusters, by = "year", all = TRUE)
  files <- bottomT_projection_files(projection_dir, c("rcp45", "rcp85"))

  regions <- list()
  history_values <- list()
  all_extents <- list()
  years <- NULL
  for (region in names(history)) {
    h <- tryCatch({
      field <- terra::rast(history[[region]])
      terra::nlyr(field)
      field
    }, error = function(e) {
      stop(region, ": could not restore the historical raster (", conditionMessage(e),
           "). Regenerate the spatial product with the updated ",
           "data_01_prepare_bottomT_by_stubstock.R.", call. = FALSE)
    })
    m <- terra::rast(maps[[region]])
    if (terra::nlyr(m) != 1L ||
        !isTRUE(terra::compareGeom(h[[1]], m, crs = TRUE, ext = TRUE,
                                   rowcol = TRUE, res = TRUE, stopOnError = FALSE))) {
      stop(region, ": cluster map must match the historical raster grid and extent.")
    }
    if (!isTRUE(terra::is.lonlat(h))) stop(region, ": historical grid must be lon/lat.")
    if (is.null(years)) years <- bottomT_years(h, region)
    if (!identical(years, bottomT_years(h, region))) {
      stop(region, ": historical years differ among substocks.")
    }
    if (length(regions)) bottomT_same_grid(h, terra::rast(regions[[1]]$template), region)
    values <- terra::values(h, mat = TRUE)
    cluster <- as.vector(terra::values(m, mat = FALSE))
    ids <- sort(unique(cluster[!is.na(cluster)]))
    if (!length(ids) || any(!is.finite(ids) | ids < 1 | ids != floor(ids))) {
      stop(region, ": cluster map must contain positive integer cluster IDs.")
    }
    area <- raster::getValues(raster::area(raster::raster(h[[1]])))
    complete <- rowSums(is.finite(values)) == ncol(values) & is.finite(area) & area > 0
    if (any(!is.na(cluster) & !complete)) {
      stop(region, ": cluster map includes cells without complete historical data.")
    }
    regions[[region]] <- list(historical = rowSums(is.finite(values)) > 0,
                              common = complete, cluster = cluster,
                              cluster_ids = ids, area = area,
                              xy = terra::xyFromCell(h, seq_len(terra::ncell(h))),
                              template = raster::raster(h[[1]]))
    history_values[[region]] <- values
    all_extents[[region]] <- as.vector(terra::ext(h))
  }
  bounds <- do.call(rbind, all_extents)
  bounds <- terra::ext(min(bounds[, 1]), max(bounds[, 2]),
                       min(bounds[, 3]), max(bounds[, 4]))

  # Scan every realization and retain cells finite for every month of both RCPs.
  months_by_scenario <- list()
  for (scenario in names(files)) {
    for (id in names(files[[scenario]])) {
      path <- files[[scenario]][[id]]
      message("Checking footprint: ", basename(path))
      p <- terra::rast(path)
      bottomT_same_grid(p, terra::rast(regions[[1]]$template), path)
      t <- bottomT_months(p, path)
      if (is.null(months_by_scenario[[scenario]])) {
        months_by_scenario[[scenario]] <- t$months
      } else if (!identical(months_by_scenario[[scenario]], t$months)) {
        stop(path, ": months differ among realizations of ", scenario)
      }
      p <- terra::crop(p, bounds)
      valid <- terra::app(is.finite(p), "sum") == terra::nlyr(p)
      valid_values <- terra::values(valid, mat = FALSE)
      for (region in names(regions)) {
        index <- terra::cellFromXY(valid, regions[[region]]$xy)
        present <- valid_values[index]
        regions[[region]]$common <- regions[[region]]$common &
          !is.na(present) & present == 1
      }
    }
  }

  coverage <- do.call(rbind, lapply(names(regions), function(region) {
    r <- regions[[region]]
    groups <- c(list(mean = r$historical),
                setNames(lapply(r$cluster_ids, function(id) !is.na(r$cluster) &
                                  r$cluster == id), paste0("Cluster", r$cluster_ids)))
    do.call(rbind, lapply(names(groups), function(group) {
      original_cells <- sum(groups[[group]])
      common_cells <- sum(groups[[group]] & r$common)
      if (!common_cells) stop(region, " ", group, ": no cells on the common footprint.")
      data.frame(substock = region, group = group, historical_cells = original_cells,
                 common_cells = common_cells, retained_fraction = common_cells / original_cells)
    }))
  }))
  if (any(coverage$retained_fraction < 1)) {
    warning("The common footprint excludes historical cells; inspect ",
            "common_footprint_coverage.csv and historical_comparison.csv.")
  }
  historical_full <- bottomT_summarize(history_values, years, regions, common = FALSE)
  historical_common <- bottomT_summarize(history_values, years, regions)
  comparison <- bottomT_compare(original, historical_full, historical_common)
  missing_comparison <- !is.finite(comparison$saved) | !is.finite(comparison$recomputed_full)
  if (any(missing_comparison)) {
    warning(sum(missing_comparison), " saved/full historical observations do not match by year; ",
            "inspect historical_comparison.csv.")
  }
  differences <- abs(comparison$full_minus_saved)
  if (any(differences > 1e-6, na.rm = TRUE)) {
    warning("Recomputed full-footprint historical values differ from saved values; ",
            "inspect historical_comparison.csv.")
  }

  projections <- setNames(vector("list", length(files)), names(files))
  for (scenario in names(files)) {
    projections[[scenario]] <- setNames(vector("list", length(files[[scenario]])),
                                        names(files[[scenario]]))
    for (id in names(files[[scenario]])) {
      path <- files[[scenario]][[id]]
      message("Aggregating: ", basename(path))
      p <- terra::crop(terra::rast(path), bounds)
      t <- bottomT_months(p, path)
      annual <- terra::tapp(p, t$years, "mean")
      projection_values <- terra::values(annual, mat = TRUE)
      values <- lapply(regions, function(r) {
        index <- terra::cellFromXY(annual, r$xy)
        projection_values[index, , drop = FALSE]
      })
      if (any(vapply(names(regions), function(region) {
        any(!is.finite(values[[region]][regions[[region]]$common, , drop = FALSE]))
      }, logical(1)))) stop(path, ": missing data on the common footprint.")
      projections[[scenario]][[id]] <-
        bottomT_summarize(values, sort(unique(t$years)), regions)
    }
  }

  common_masks <- lapply(regions, function(r) {
    raster::setValues(r$template, ifelse(r$common, 1, NA_real_))
  })
  overlap <- bottomT_overlap(historical_common, projections)
  result <- list(historical_full = historical_full, historical_common = historical_common,
                 historical_comparison = comparison, coverage = coverage,
                 common_masks = common_masks, projections = projections,
                 overlap = overlap, projection_files = files)
  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) {
    stop("Could not create output directory: ", output_dir)
  }
  saveRDS(result, file.path(output_dir, "bottomT_historical_projection_summaries.rds"))
  utils::write.csv(comparison, file.path(output_dir, "historical_comparison.csv"), row.names = FALSE)
  utils::write.csv(coverage, file.path(output_dir, "common_footprint_coverage.csv"), row.names = FALSE)
  utils::write.csv(overlap, file.path(output_dir, "historical_projection_overlap.csv"), row.names = FALSE)
  bottomT_plot_diagnostics(file.path(output_dir, "historical_projection_comparison.pdf"),
                           historical_common, comparison, projections)
  invisible(result)
}
