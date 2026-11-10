# Compare the model_04 simulations, ICES advice, and fitted SR relationships.
# Run from the repository root in a fresh R session.
library(FLCore)
library(ggplot2)
library(dplyr)
library(tidyr)
library(earth)
library(data.table)

fig_dir <- file.path("output", "figs")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# model_04 uses the inclusive range (fy - ey):fy: 21 years when ey = 20.
ey <- 20L
sr_types <- c("BH", "SEG")
scenario_order <- c("noCC", "rcp45", "rcp85")
scenario_labels <- c(noCC = "No climate change",
                     rcp45 = "RCP 4.5", rcp85 = "RCP 8.5")
scenario_colours <- c(noCC = "#333333", rcp45 = "#0072B2",
                      rcp85 = "#D55E00", Benchmark = "#009E73")
period_labels <- c("None", "2030", "2040", "2050")
set.seed(2026)

load_object <- function(path, object) {
  if (!file.exists(path)) stop("Missing input file: ", path)
  e <- new.env(parent = baseenv())
  load(path, envir = e)
  if (!exists(object, envir = e, inherits = FALSE))
    stop("Missing object '", object, "' in ", path)
  get(object, envir = e, inherits = FALSE)
}

require_columns <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  if (length(missing))
    stop(label, " is missing columns: ", paste(missing, collapse = ", "))
}

as_year <- function(x, label) {
  y <- suppressWarnings(as.numeric(as.character(x)))
  if (anyNA(y) || any(!is.finite(y)) || any(y != floor(y)))
    stop("Invalid year in ", label)
  y
}

# ---------------------------------------------------------------------------
# Reference points. BRPs contains published ICES advice values, not the
# output of a separately rerun original ICES EqSim routine.
# ---------------------------------------------------------------------------
ices <- load_object(file.path("data", "BRPs.Rdata"), "BRPs")
bench <- load_object(file.path("model", "newBRPs_SRR_benchmark_model.rda"),
                     "newBRPs_bench")
climate <- load_object(file.path("model", "newBRPs_EMSRR_model.rda"),
                       "newBRPs")
stocks <- rownames(ices)
if (is.null(stocks) || anyDuplicated(stocks) ||
    !setequal(stocks, rownames(bench)))
  stop("ICES and benchmark reference-point tables have different stocks.")
require_columns(ices, c("Fmsy", "Btrigger"), "ICES BRPs")
require_columns(bench, c("Fmsy", "Btrigger"), "Benchmark BRPs")
require_columns(climate, c("stock", "srr.type", "climScen", "period",
                           "Fmsy", "FmsyL", "FmsyH", "Btrigger"),
                "Climate BRPs")

climate <- as.data.frame(climate) |>
  mutate(stock = as.character(stock), srr.type = as.character(srr.type),
         climScen = as.character(climScen), period = as.character(period))
expected <- expand.grid(stock = stocks, srr.type = sr_types,
                        setting = c("noCC_None", paste0("rcp45_", c(2030, 2040, 2050)),
                                    paste0("rcp85_", c(2030, 2040, 2050))),
                        stringsAsFactors = FALSE)
actual <- paste(climate$stock, climate$srr.type,
                paste(climate$climScen, climate$period, sep = "_"))
wanted <- paste(expected$stock, expected$srr.type, expected$setting)
if (anyDuplicated(actual) || !setequal(actual, wanted))
  stop("Climate reference points do not contain exactly the expected stock/model/scenario combinations.")

ice_long <- data.frame(stock = stocks, ices, row.names = NULL) |>
  pivot_longer(c(Fmsy, Btrigger), names_to = "metric", values_to = "value") |>
  mutate(reference = "ICES advice")
bench_long <- data.frame(stock = stocks, bench[stocks, , drop = FALSE],
                         row.names = NULL) |>
  pivot_longer(c(Fmsy, Btrigger), names_to = "metric", values_to = "value") |>
  mutate(reference = "SRR benchmark")
if (any(!is.finite(ice_long$value)) || any(!is.finite(bench_long$value)) ||
    any(!is.finite(climate$Fmsy)) || any(!is.finite(climate$Btrigger)))
  stop("Non-finite Fmsy or Btrigger in reference-point tables.")
baseline_long <- bind_rows(ice_long, bench_long)

ref_long <- climate |>
  pivot_longer(c(Fmsy, Btrigger), names_to = "metric", values_to = "value") |>
  mutate(period = factor(period, levels = period_labels),
         climScen = factor(climScen, levels = scenario_order),
         srr.type = factor(srr.type, levels = sr_types))

plot_refpts <- function(metric, y_label, filename) {
  p <- ggplot(filter(ref_long, .data$metric == .env$metric),
              aes(x = period, y = value, colour = climScen,
                  group = climScen)) +
    geom_hline(data = filter(baseline_long, .data$metric == .env$metric),
               aes(yintercept = value, linetype = reference),
               colour = "grey35") +
    geom_line(linewidth = 0.65) +
    geom_point(size = 2.3) +
    facet_grid(stock ~ srr.type, scales = "free_y") +
    scale_colour_manual(values = scenario_colours[scenario_order],
                        labels = scenario_labels, drop = FALSE,
                        name = "Climate") +
    scale_linetype_manual(values = c("ICES advice" = "dashed",
                                     "SRR benchmark" = "dotdash"),
                          name = "Reference") +
    labs(x = "Climate period (starting year)", y = y_label,
         title = paste(metric, "across climate settings"),
         subtitle = "Horizontal lines: published ICES advice and no-environment SRR benchmark",
         caption = "No climate change has no future period; climate periods represent decade-long environmental levels.") +
    theme_bw() +
    theme(legend.position = "bottom")
  if (metric == "Fmsy") {
    p <- p +
      geom_errorbar(aes(ymin = FmsyL, ymax = FmsyH), width = 0.12) +
      labs(caption = paste("Bars: target-F range yielding >95% of maximum landings,",
                           "not confidence intervals. NoCC has no future period."))
  }
  ggsave(file.path(fig_dir, filename), p, width = 11, height = 8, dpi = 300)
}
plot_refpts("Fmsy", "Fmsy (target F)", "refpts_Fmsy.png")
plot_refpts("Btrigger", "Btrigger (tonnes)", "refpts_Btrigger.png")

# ---------------------------------------------------------------------------
# No-HCR yield and SSB curves. Summarise the same years as model_04:
# first the median across years per iteration, then across iterations.
# Bands are the 5th–95th percentiles of those iteration medians.
# ---------------------------------------------------------------------------
summarise_run <- function(path, stock, srr_type, scenario, period) {
  x <- load_object(path, "df")
  require_columns(x, c("year", "iter", "Ftarg", "slot", "data",
                       "stock", "climScen", "period"),
                  basename(path))
  expected_scenario <- if (scenario == "Benchmark") "benchmark" else scenario
  if (anyNA(x$stock) || anyNA(x$climScen) || anyNA(x$period) ||
      any(as.character(x$stock) != stock) ||
      any(as.character(x$climScen) != expected_scenario) ||
      any(as.character(x$period) != period))
    stop("Simulation metadata does not match the filename: ", path)
  x <- as.data.frame(x[, c("year", "iter", "Ftarg", "slot", "data")])
  x$year <- as_year(x$year, basename(path))
  x <- x[x$slot %in% c("landings", "ssb"), , drop = FALSE]
  if (!nrow(x) || !all(c("landings", "ssb") %in% x$slot))
    stop("Missing landings or SSB in ", path)
  fy <- max(x$year)
  years <- seq.int(fy - ey, fy)
  if (!all(years %in% x$year))
    stop("Incomplete evaluation-year range in ", path)
  x <- x[x$year %in% years, , drop = FALSE]
  if (any(!is.finite(x$data)) || any(!is.finite(x$Ftarg)) || anyNA(x$iter))
    stop("Invalid simulation values in ", path)
  if (anyDuplicated(x[, c("year", "iter", "Ftarg", "slot")]))
    stop("Duplicate simulation year/iteration/F/slot values in ", path)
  per_iter <- x |>
    group_by(slot, Ftarg, iter) |>
    summarise(value = median(data), n_years = n_distinct(year),
              .groups = "drop")
  if (any(per_iter$n_years != length(years)))
    stop("Missing evaluation years for an iteration/F/slot in ", path)
  per_iter |>
    group_by(slot, Ftarg) |>
    summarise(median = median(value), low = quantile(value, 0.05),
              high = quantile(value, 0.95), .groups = "drop") |>
    mutate(stock = stock, srr.type = srr_type, scenario = scenario,
           period = period)
}

curve_parts <- list()
k <- 0L
for (stock in stocks) {
  cat("Calc. summary for",stock,"...","\n")
  bench_file <- file.path("model", "refpts",
                          paste0(stock, "_SRR_benchmark_model.rda"))
  benchmark_curves <- summarise_run(bench_file, stock, "Benchmark",
                                    "Benchmark", "None")
  for (type in sr_types) {
    k <- k + 1L
    curve_parts[[k]] <- mutate(benchmark_curves, srr.type = type)
    rows <- climate[climate$stock == stock & climate$srr.type == type, ]
    for (i in seq_len(nrow(rows))) {
      sc <- rows$climScen[i]
      period <- rows$period[i]
      path <- file.path("model", "refpts",
                        paste(stock, type, sc, period, "EMSRR", "model.rda",
                              sep = "_"))
      k <- k + 1L
      curve_parts[[k]] <- summarise_run(path, stock, type, sc, period)
    }
  }
}
curves <- bind_rows(curve_parts) |>
  mutate(scenario = factor(scenario, levels = c("Benchmark", scenario_order)),
         period = factor(period, levels = period_labels),
         srr.type = factor(srr.type, levels = sr_types))

# Separate period facets avoid overlapping seven ribbons. The benchmark and
# noCC curves are repeated in future-period panels as references.
curve_plot_data <- bind_rows(
  lapply(period_labels, function(p) {
    curves |>
      filter((scenario %in% c("Benchmark", "noCC") & period == "None") |
               (p != "None" & period == p)) |>
      mutate(panel_period = if (p == "None") "Baseline" else p)
  })
) |>
  mutate(panel_period = factor(panel_period,
                               levels = c("Baseline", "2030", "2040", "2050")))

plot_curves <- function(slot_name, y_label, filename) {
  for (type in sr_types) {
    p <- ggplot(filter(curve_plot_data, slot == slot_name, srr.type == type),
                aes(x = Ftarg, y = median, colour = scenario, fill = scenario,
                    group = scenario)) +
      geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.12, colour = NA) +
      geom_line(linewidth = 0.7) +
      facet_grid(stock ~ panel_period, scales = "free_y") +
      scale_colour_manual(values = scenario_colours, name = "Scenario",
                          labels = c(Benchmark = "SRR benchmark",
                                     scenario_labels), drop = FALSE) +
      scale_fill_manual(values = scenario_colours, name = "Scenario",
                        labels = c(Benchmark = "SRR benchmark",
                                   scenario_labels), drop = FALSE) +
      labs(x = "Target fishing mortality (Ftarg)", y = y_label,
           title = paste(type, "no-HCR simulations"),
           subtitle = "Lines: median of iteration medians; ribbons: 5th–95th percentiles",
           caption = paste0("Each iteration is summarised over (fy - ", ey,
                            "):fy, inclusive. Yield is landings, not catch.")) +
      theme_bw() + theme(legend.position = "bottom")
    ggsave(file.path(fig_dir, paste0(filename, "_", type, ".png")),
           p, width = 14, height = 8, dpi = 300)
  }
}
plot_curves("landings", "Landings (assessment units)", "yield_F")
plot_curves("ssb", "SSB (assessment units)", "SSB_F")

# ---------------------------------------------------------------------------
# Historical recruitment against lagged SSB and fitted relationships.
# model_03's noCC offset is the mean of lagged historical environmental
# values for years 2007–2021 (2006–2020 plus the one-year lag).
# ---------------------------------------------------------------------------
models <- setNames(lapply(sr_types, function(type) {
  load_object(file.path("model", paste0("Cod.EMSRR_", type,
                                        "_model_and_data.RData")),
              "Cod.EMSRR.obj")
}), sr_types)

sr_history <- list()
sr_predictions <- list()
offsets <- list()
for (stock in stocks) {
  hist <- models[["BH"]][[stock]]$stock.data
  require_columns(hist, c("year", "ssb_lag1", "rec"), paste(stock, "history"))
  hist <- hist[is.finite(hist$ssb_lag1) & is.finite(hist$rec) &
                 hist$ssb_lag1 > 0 & hist$rec > 0, ]
  if (nrow(hist) < 5L) stop("Insufficient historical stock data for ", stock)
  hist$stock <- stock
  sr_history[[stock]] <- hist
  ssb_grid <- seq(0, max(hist$ssb_lag1), length.out = 201L)
  
  for (type in sr_types) {
    obj <- models[[type]][[stock]]
    model <- obj$model
    var <- model$env.vars
    if (length(var) != 1L || !identical(model$env.lag, 1L) ||
        length(model$SR.model$vars) != 2L ||
        !identical(model$SR.model$vars[2], var))
      stop("Unexpected EMSRR environmental inputs for ", stock, " ", type)
    env <- obj$env.variables$historical
    require_columns(env, c("year", var), paste(stock, type, "environment"))
    env$year <- as_year(env$year, paste(stock, type, "environment"))
    if (anyDuplicated(env$year)) stop("Duplicate environmental years for ", stock)
    noCC_years <- 2006:2020 + model$env.lag
    env_values <- env[[var]][match(noCC_years, env$year)]
    if (length(env_values) != 15L || any(!is.finite(env_values)))
      stop("Missing historical noCC offset values for ", stock, " ", type)
    env_mean <- mean(env_values)
    offsets[[paste(stock, type)]] <- data.frame(stock, type, variable = var,
                                                noCC_offset = env_mean)
    newdata <- data.frame(ssb_grid, rep(env_mean, length(ssb_grid)))
    names(newdata) <- model$SR.model$vars
    pred <- as.numeric(model$predict.func(model$SR.model, newdata))
    sigma <- model$sd.residuals
    if (length(pred) != length(ssb_grid) || any(!is.finite(pred)) ||
        any(pred < 0) || length(sigma) != 1L || !is.finite(sigma) || sigma < 0)
      stop("Invalid EMSRR prediction or residual SD for ", stock, " ", type)
    mult <- exp(qnorm(c(0.025, 0.975)) * sigma)
    sr_predictions[[paste(stock, type)]] <-
      data.frame(stock, relationship = paste("EMSRR", type), ssb = ssb_grid,
                 fit = pred, low = pred * mult[1], high = pred * mult[2])
  }
  
  srpars <- load_object(file.path("data", "SRRs",
                                  paste0("Cod_", stock, "_SRR_bootstrap.rda")),
                        "srpars")
  pars_long <- as.data.frame(srpars) |>
    select(iter, params, data) |>
    filter(params %in% c("a", "b", "m", "sigmaR"))
  if (anyDuplicated(pars_long[, c("iter", "params")]))
    stop("Duplicate SRR bootstrap parameters for ", stock)
  pars <- pivot_wider(pars_long, names_from = params, values_from = data)
  require_columns(pars, c("iter", "a", "b", "m", "sigmaR"),
                  paste(stock, "SRR bootstrap"))
  if (!nrow(pars) || any(!is.finite(as.matrix(pars[, c("a", "b", "m", "sigmaR")]))) ||
      any(!pars$m %in% 1:3) || any(pars$a <= 0) || any(pars$b <= 0) ||
      any(pars$sigmaR < 0))
    stop("Invalid SRR bootstrap parameters for ", stock)
  
  # Include both saved parameter variation and lognormal recruitment noise.
  # As in model_03, the lognormal deviations are not bias-corrected.
  draws_per_iteration <- 400L
  idx <- rep(seq_len(nrow(pars)), each = draws_per_iteration)
  noise <- exp(rnorm(length(idx), sd = pars$sigmaR[idx]))
  srr_at_ssb <- lapply(ssb_grid, function(s) {
    fitted <- numeric(nrow(pars))
    bh <- pars$m == 1
    ricker <- pars$m == 2
    seg <- pars$m == 3
    fitted[bh] <- pars$a[bh] * s / (pars$b[bh] + s)
    fitted[ricker] <- pars$a[ricker] * s * exp(-pars$b[ricker] * s)
    fitted[seg] <- pars$a[seg] * pmin(s, pars$b[seg])
    interval <- quantile(fitted[idx] * noise, c(0.025, 0.975))
    c(fit = median(fitted), low = interval[1], high = interval[2])
  })
  vals <- do.call(rbind, srr_at_ssb)
  sr_predictions[[paste(stock, "SRR")]] <-
    data.frame(stock, relationship = "SRR benchmark", ssb = ssb_grid,
               fit = vals[, 1], low = vals[, 2], high = vals[, 3])
}

history <- bind_rows(sr_history)
sr_plot_data <- bind_rows(sr_predictions)
offset_table <- bind_rows(offsets)
print(offset_table)

# One panel per relationship makes the prediction intervals legible.
history_panels <- tidyr::crossing(
  history, relationship = c("SRR benchmark", "EMSRR BH", "EMSRR SEG")
)
p_sr <- ggplot() +
  geom_ribbon(data = sr_plot_data,
              aes(x = ssb, ymin = low, ymax = high, fill = relationship),
              alpha = 0.2) +
  geom_line(data = sr_plot_data,
            aes(x = ssb, y = fit, colour = relationship), linewidth = 0.8) +
  geom_point(data = history_panels, aes(x = ssb_lag1, y = rec),
             alpha = 0.55, size = 1.4) +
  facet_wrap(~stock, ncol = 3, scales = "free") +
  scale_colour_manual(values = c("SRR benchmark" = "#009E73",
                                 "EMSRR BH" = "#0072B2",
                                 "EMSRR SEG" = "#D55E00"), guide = "none") +
  scale_fill_manual(values = c("SRR benchmark" = "#009E73",
                               "EMSRR BH" = "#0072B2",
                               "EMSRR SEG" = "#D55E00"), guide = "none") +
  labs(x = "SSB",
       y = "Recruitment (in thousands)",
       title = "Historical stock–recruitment data and fitted relationships",
       subtitle = "EMSRR environment fixed at the noCC offset",
       caption = paste("Ribbons: 95% recruitment prediction intervals.",
                       "SRR line: median bootstrap fit; its ribbon also includes bootstrap variation.")) + 
  theme_bw()


ggsave(file.path(fig_dir, "stock_recruitment.png"), p_sr + theme_bw(),
       width = 12, height = 9, dpi = 300)

cat("Figures saved to ", fig_dir, "\n", sep = "")