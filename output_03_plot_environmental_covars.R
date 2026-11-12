# Historical and projected environmental inputs used by the EMSRR models.
# Run from the repository root in a fresh R session.
library(FLCore)
library(earth)
library(ggplot2)

fig_dir <- file.path("output", "figs")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# plot only for one model type as the selected covars are identical
sr_types <- "SEG" #c("BH", "SEG")
scenarios <- c("rcp45", "rcp85")
periods <- c(2030L, 2040L, 2050L)

load_object <- function(path, object) {
  if (!file.exists(path)) stop("Missing input file: ", path)
  
  e <- new.env(parent = baseenv())
  load(path, envir = e)
  
  if (!exists(object, envir = e, inherits = FALSE))
    stop("Missing object '", object, "' in ", path)
  
  get(object, envir = e, inherits = FALSE)
}

check_years <- function(years, label) {
  if (!is.numeric(years) || !length(years) ||
      any(!is.finite(years)) || any(years != floor(years)) ||
      anyDuplicated(years))
    stop(label, ": years must be unique finite integers.")
}

projection_matrix <- function(x, variable, label) {
  if (!is.data.frame(x) || !"year" %in% names(x))
    stop(label, ": missing environmental projection table or year column.")
  
  check_years(x$year, label)
  x <- x[order(x$year), , drop = FALSE]
  
  if (any(diff(x$year) != 1L))
    stop(label, ": projection years must be consecutive.")
  
  cols <- setdiff(names(x), "year")
  prefix <- paste0(variable, "_it")
  ids <- substring(cols, nchar(prefix) + 1L)
  
  if (!length(cols) || any(!startsWith(cols, prefix)) ||
      any(!grepl("^[0-9]+$", ids)) || anyDuplicated(ids))
    stop(label, ": invalid environmental realization columns.")
  
  values <- as.matrix(x[, cols, drop = FALSE])
  if (!is.numeric(values) || any(!is.finite(values)))
    stop(label, ": environmental realizations must be finite numeric values.")
  
  list(year = x$year, values = values)
}

summarise_realizations <- function(year, values, stock, type, variable, scenario) {
  q <- t(apply(
    values, 1, quantile,
    probs = c(0.05, 0.5, 0.95),
    names = FALSE
  ))
  
  data.frame(
    stock = stock,
    type = type,
    variable = variable,
    scenario = scenario,
    year = year,
    low = q[, 1],
    median = q[, 2],
    high = q[, 3]
  )
}

models <- setNames(lapply(sr_types, function(type) {
  load_object(
    file.path("model", paste0("Cod.EMSRR_", type, "_model_and_data.RData")),
    "Cod.EMSRR.obj"
  )
}), sr_types)

stocks <- names(models[["SEG"]])
# if (is.null(stocks) || anyDuplicated(stocks) ||
#     !setequal(stocks, names(models[["SEG"]])))
#   stop("The BH and SEG files must contain the same named substocks.")

history_parts <- list()
projection_parts <- list()
period_parts <- list()
offset_parts <- list()
k <- 0L

for (stock in stocks) {
  for (type in sr_types) {
    k <- k + 1L
    obj <- models[[type]][[stock]]
    variable <- obj$model$env.vars
    lag <- obj$model$env.lag
    label <- paste(stock, type)
    
    if (length(variable) != 1L || is.na(variable) || !nzchar(variable) ||
        length(lag) != 1L || is.na(lag) || lag != 1L)
      stop(label, ": expected one environmental variable with a one-year lag.")
    
    history <- obj$env.variables$historical
    if (!is.data.frame(history) ||
        !all(c("year", variable) %in% names(history)))
      stop(label, ": missing historical environmental covariate.")
    
    check_years(history$year, paste(label, "history"))
    history <- history[order(history$year), , drop = FALSE]
    
    if (!is.numeric(history[[variable]]) ||
        any(!is.finite(history[[variable]])))
      stop(label, ": historical covariate must contain finite numeric values.")
    
    history_parts[[k]] <- data.frame(
      stock = stock,
      type = type,
      variable = variable,
      year = history$year,
      value = history[[variable]]
    )
    
    # model_03 offsets noCC by the 2006-2020 historical environmental mean.
    # Values in model_02 have already been shifted to recruitment year.
    offset_years <- (2006:2020) + lag
    offset_values <- history[[variable]][match(offset_years, history$year)]
    
    if (any(!is.finite(offset_values)))
      stop(label, ": missing historical environmental values for noCC offset.")
    
    offset <- mean(offset_values)
    offset_parts[[k]] <- data.frame(
      stock = stock,
      type = type,
      variable = variable,
      noCC_offset = offset
    )
    
    projections <- setNames(lapply(scenarios, function(scenario) {
      projection_matrix(
        obj$env.variables$climate.projections[[scenario]][[variable]],
        variable,
        paste(label, scenario)
      )
    }), scenarios)
    
    if (!identical(projections$rcp45$year, projections$rcp85$year))
      stop(label, ": the two RCP scenarios have different projection years.")
    
    # Match model_03: remove the RCP4.5 annual realization median and
    # replace it with the fixed historical offset, retaining variability.
    noCC <- projections$rcp45
    noCC$values <- sweep(
      noCC$values, 1, apply(noCC$values, 1, median), "-"
    ) + offset
    
    all_projections <- c(list(noCC = noCC), projections)
    projection_parts[[k]] <- do.call(
      rbind,
      lapply(names(all_projections), function(sc) {
        x <- all_projections[[sc]]
        summarise_realizations(
          x$year, x$values, stock, type, variable, sc
        )
      })
    )
    
    # These are means, not medians: model_03 averages each realization
    # over the decade and then averages across realizations.
    period_parts[[k]] <- do.call(
      rbind,
      lapply(scenarios, function(sc) {
        x <- projections[[sc]]
        
        do.call(rbind, lapply(periods, function(start) {
          years <- seq.int(start, start + 9L)
          i <- match(years, x$year)
          
          if (anyNA(i))
            stop(
              label, " ", sc, ": missing years in period ",
              start, "-", start + 9L
            )
          
          level <- mean(colMeans(x$values[i, , drop = FALSE]))
          
          data.frame(
            stock = stock,
            type = type,
            variable = variable,
            scenario = sc,
            start = start,
            end = start + 9L,
            level = level
          )
        }))
      })
    )
  }
}

historical <- do.call(rbind, history_parts)
projected <- do.call(rbind, projection_parts)
period_levels <- do.call(rbind, period_parts)
offsets <- do.call(rbind, offset_parts)

panel_label <- function(stock, type, variable) {
  paste0(stock, " | EMSRR ", type, "\n", variable)
}

panel_levels <- unique(panel_label(
  historical$stock, historical$type, historical$variable
))

historical$panel <- factor(
  panel_label(historical$stock, historical$type, historical$variable),
  levels = panel_levels
)
projected$panel <- factor(
  panel_label(projected$stock, projected$type, projected$variable),
  levels = panel_levels
)
period_levels$panel <- factor(
  panel_label(period_levels$stock, period_levels$type, period_levels$variable),
  levels = panel_levels
)

series_colours <- c(
  "Historical" = "#333333",
  noCC = "#009E73",
  rcp45 = "#0072B2",
  rcp85 = "#D55E00"
)
series_labels <- c(
  "Historical" = "Historical",
  noCC = "No climate change",
  rcp45 = "RCP 4.5",
  rcp85 = "RCP 8.5"
)
series_linetypes <- c(
  "Historical" = "solid",
  noCC = "dashed",
  rcp45 = "solid",
  rcp85 = "solid"
)

p <- ggplot() +
  geom_ribbon(
    data = projected,
    aes(x = year, ymin = low, ymax = high, fill = scenario,
        group = scenario),
    alpha = 0.14, colour = NA
  ) +
  geom_line(
    data = projected,
    aes(x = year, y = median, colour = scenario,
        linetype = scenario, group = scenario),
    linewidth = 0.55
  ) +
  geom_line(
    data = historical,
    aes(x = year, y = value, colour = "Historical",
        linetype = "Historical"),
    linewidth = 0.65
  ) +
  geom_segment(
    data = period_levels,
    aes(x = start, xend = end, y = level, yend = level),
    colour = "white", linewidth = 2.8, show.legend = FALSE
  ) +
  geom_segment(
    data = period_levels,
    aes(x = start, xend = end, y = level, yend = level,
        colour = scenario),
    linewidth = 1.5, show.legend = FALSE
  ) +
  facet_wrap(vars(panel), scales = "free_y") +
  scale_colour_manual(
    values = series_colours, labels = series_labels,
    breaks = names(series_colours), name = "Series"
  ) +
  scale_linetype_manual(
    values = series_linetypes, labels = series_labels,
    breaks = names(series_colours), name = "Series"
  ) +
  scale_fill_manual(
    values = series_colours, labels = series_labels,
    breaks = names(series_colours), name = "Series"
  ) +
  labs(
    x = "Recruitment year (environment from the preceding year)",
    y = "Environmental covariate (model units)",
    title = "Historical and projected EMSRR environmental covariates",
    subtitle = paste(
      "Lines: annual realization medians; ribbons: 5th-95th percentiles",
      "of all model_02 realizations.",
      "Thick bars: model_03 decade-run mean levels."
    ),
    caption = paste(
      "noCC retains RCP4.5 variability around its historical 2006-2020 mean.",
      "Decade bars cover recruitment years 2030-2039, 2040-2049 and 2050-2059."
    )
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(
  file.path(fig_dir, "environmental_covariates.png"),
  p, width = 13, height = 10, dpi = 300
)

write.csv(
  period_levels[, c("stock", "type", "variable", "scenario",
                    "start", "end", "level")],
  file.path("output", "environmental_period_levels.csv"),
  row.names = FALSE
)
write.csv(
  offsets,
  file.path("output", "environmental_noCC_offsets.csv"),
  row.names = FALSE
)

cat("Saved environmental figure and period levels under output/.\n")