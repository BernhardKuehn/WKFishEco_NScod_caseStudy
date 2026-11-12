# Recruitment projections at fixed recent SSB, using model_02 inputs.
# Run from the repository root in a fresh R session.
library(FLCore)
library(earth)
library(ggplot2)
library(dplyr)
library(tidyr)

fig_dir <- file.path("output", "figs")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# The 2006 projected environment predicts 2007 recruitment after the lag.
first_projection_year <- 2007L
last_projection_year <- 2099L
sr_types <- c("BH", "SEG","SEG_trunc_ICES")
scenarios <- c("rcp45", "rcp85")

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

check_years <- function(x, label) {
  if (!is.numeric(x) || any(!is.finite(x)) ||
      any(x != floor(x)) || anyDuplicated(x))
    stop(label, ": years must be unique finite integers.")
}

# Pointwise quantiles of an equally weighted mixture of lognormal predictions.
# log_sd can differ between SRR bootstrap fits; zero SD is a point mass.
mixture_quantiles <- function(fit, log_sd, probs = c(0.025, 0.5, 0.975)) {
  if (!length(fit) || length(log_sd) != length(fit) ||
      any(!is.finite(fit)) || any(fit <= 0) ||
      any(!is.finite(log_sd)) || any(log_sd < 0))
    stop("Invalid recruitment predictions or log-scale residual SDs.")
  
  mu <- log(fit)
  lower <- min(mu - 10 * log_sd) - 1
  upper <- max(mu + 10 * log_sd) + 1
  
  result <- vapply(probs, function(p) {
    lo <- lower
    hi <- upper
    
    for (i in seq_len(60L)) {
      mid <- (lo + hi) / 2
      cdf <- numeric(length(mu))
      stochastic <- log_sd > 0
      
      cdf[stochastic] <- pnorm(
        (mid - mu[stochastic]) / log_sd[stochastic]
      )
      cdf[!stochastic] <- as.numeric(mid >= mu[!stochastic])
      
      if (mean(cdf) >= p) hi <- mid else lo <- mid
    }
    
    exp(hi)
  }, numeric(1))
  
  if (any(!is.finite(result))) stop("Recruitment quantiles overflowed.")
  result
}

models <- setNames(lapply(sr_types, function(type) {
  load_object(
    file.path("model", paste0("Cod.EMSRR_", type, "_model_and_data.RData")),
    "Cod.EMSRR.obj"
  )
}), sr_types)

stocks <- names(models[["BH"]])
if (is.null(stocks) || anyDuplicated(stocks) ||
    !setequal(stocks, names(models[["SEG"]])))
  stop("The BH and SEG model files must contain the same named substocks.")

stock_file <- file.path(
  "data", "stock",
  "cod46a7d20_FLstockObjects_assessment_2026.Rdata"
)
assessment <- load_object(stock_file, "stks")

stock_names <- vapply(assessment, function(x) {
  name <- x@name
  if (grepl("Southern|South", name)) "Southern" else
    if (grepl("Northwestern|Northwest", name)) "Northwestern" else
      if (grepl("Viking", name)) "Viking" else
        stop("Unrecognized assessment stock name: ", name)
}, character(1))

if (anyDuplicated(stock_names) || !setequal(stocks, stock_names))
  stop("Assessment stocks do not match the prepared EMSRR models.")
names(assessment) <- stock_names

histories <- list()
ssb_settings <- list()
forecast_parts <- list()
k <- 0L

for (stock in stocks) {
  # model_02's stock.data contains lagged SSB and omits the final SSB year.
  # Obtain the actual last five assessment SSB years from the FLStock.
  s <- FLCore::ssb(assessment[[stock]])
  stock_years <- suppressWarnings(as.numeric(dimnames(s)$year))
  check_years(stock_years, paste(stock, "assessment SSB"))
  
  ssb_values <- as.numeric(s)
  if (length(ssb_values) != length(stock_years))
    stop(stock, ": expected one SSB value per assessment year.")
  
  origin <- max(stock_years)
  recent_years <- seq.int(origin - 4L, origin)
  recent_ssb <- ssb_values[match(recent_years, stock_years)]
  
  if (any(!is.finite(recent_ssb)) || any(recent_ssb <= 0))
    stop(stock, ": the last five assessment SSB years must be complete and positive.")
  
  fixed_ssb <- mean(recent_ssb)
  ssb_settings[[stock]] <- data.frame(
    stock = stock,
    assessment_year = origin,
    first_ssb_year = origin - 4L,
    last_ssb_year = origin,
    mean_ssb = fixed_ssb
  )
  
  history <- models[["BH"]][[stock]]$stock.data
  require_columns(history, c("year", "rec"), paste(stock, "history"))
  check_years(history$year, paste(stock, "historical recruitment"))
  
  if (any(history$year > origin))
    stop(stock, ": recruitment history extends beyond the assessment SSB.")
  
  history <- history[is.finite(history$rec) & history$rec >= 0, ]
  if (!nrow(history)) stop(stock, ": no historical recruitment observations.")
  history$stock <- stock
  histories[[stock]] <- history
  
  # SRR benchmark: each saved iteration has a model and residual SD.
  srpars <- load_object(
    file.path("data", "SRRs", paste0("Cod_", stock, "_SRR_bootstrap.rda")),
    "srpars"
  )
  
  pars_long <- as.data.frame(srpars) |>
    select(iter, params, data) |>
    filter(params %in% c("a", "b", "m", "sigmaR"))
  
  if (anyDuplicated(pars_long[, c("iter", "params")]))
    stop(stock, ": duplicate SRR bootstrap parameters.")
  
  pars <- pivot_wider(pars_long, names_from = params, values_from = data)
  require_columns(
    pars, c("iter", "a", "b", "m", "sigmaR"),
    paste(stock, "SRR bootstrap")
  )
  
  if (!nrow(pars) ||
      any(!is.finite(as.matrix(pars[, c("a", "b", "m", "sigmaR")]))) ||
      any(!pars$m %in% 1:3) || any(pars$a <= 0) || any(pars$b <= 0) ||
      any(pars$sigmaR < 0))
    stop(stock, ": invalid SRR bootstrap parameters.")
  
  fit <- numeric(nrow(pars))
  bh <- pars$m == 1
  ricker <- pars$m == 2
  seg <- pars$m == 3
  
  fit[bh] <- pars$a[bh] * fixed_ssb / (pars$b[bh] + fixed_ssb)
  fit[ricker] <- pars$a[ricker] * fixed_ssb *
    exp(-pars$b[ricker] * fixed_ssb)
  fit[seg] <- pars$a[seg] * pmin(fixed_ssb, pars$b[seg])
  
  srr_quantiles <- mixture_quantiles(fit, pars$sigmaR)
  
  expected_years <- NULL
  
  for (type in sr_types) {
    obj <- models[[type]][[stock]]
    model <- obj$model
    var <- model$env.vars
    
    if (length(var) != 1L || !identical(model$env.lag, 1L) ||
        length(model$SR.model$vars) != 2L ||
        !identical(model$SR.model$vars, c("ssb_lag1", var)))
      stop(stock, " ", type, ": unexpected EMSRR predictor names or lag.")
    
    sigma <- model$sd.residuals
    if (length(sigma) != 1L || !is.finite(sigma) || sigma < 0)
      stop(stock, " ", type, ": invalid EMSRR residual SD.")
    
    for (scenario in scenarios) {
      projection <- obj$env.variables$climate.projections[[scenario]][[var]]
      
      if (!is.data.frame(projection))
        stop(stock, " ", type, " ", scenario, ": missing projection table.")
      
      require_columns(projection, "year", paste(stock, type, scenario))
      check_years(projection$year, paste(stock, type, scenario))
      
      cols <- setdiff(names(projection), "year")
      prefix <- paste0(var, "_it")
      suffix <- substring(cols, nchar(prefix) + 1L)
      
      if (!length(cols) || any(!startsWith(cols, prefix)) ||
          any(!grepl("^[0-9]+$", suffix)) || anyDuplicated(suffix))
        stop(stock, " ", type, " ", scenario, ": invalid realization columns.")
      
      projection <- projection[
        projection$year >= first_projection_year &
          projection$year <= last_projection_year, ,
        drop = FALSE
      ]
      projection <- projection[order(projection$year), , drop = FALSE]
      years <- projection$year
      
      if (!length(years) || years[1] != first_projection_year)
        stop(
          stock, " ", type, " ", scenario,
          ": no lagged environmental projection for ", first_projection_year,
          ". model_02 requires the preceding year's projected environment; ",
          "adjust first_projection_year only if a later start is intended."
        )
      
      if (any(diff(years) != 1L))
        stop(stock, " ", type, " ", scenario,
             ": projection years must be consecutive.")
      
      if (tail(years, 1) != last_projection_year)
        stop(
          stock, " ", type, " ", scenario,
          ": projection does not reach ", last_projection_year,
          ". Adjust last_projection_year if a shorter forecast is intended."
        )
      
      if (is.null(expected_years)) expected_years <- years
      if (!identical(years, expected_years))
        stop(stock, ": forecast years differ across models or climate scenarios.")
      
      env <- as.matrix(projection[, cols, drop = FALSE])
      if (!is.numeric(env) || any(!is.finite(env)))
        stop(stock, " ", type, " ", scenario,
             ": invalid projected environment.")
      
      # Each column is a realization; all years are already lagged by model_02.
      newdata <- data.frame(rep(fixed_ssb, length(env)), as.vector(env))
      names(newdata) <- model$SR.model$vars
      
      predictions <- as.numeric(model$predict.func(model$SR.model, newdata))
      if (length(predictions) != length(env))
        stop(stock, " ", type, " ", scenario, ": unexpected prediction length.")
      
      predictions <- matrix(predictions, nrow = length(years))
      quants <- t(vapply(seq_along(years), function(i) {
        mixture_quantiles(predictions[i, ], rep(sigma, ncol(predictions)))
      }, numeric(3)))
      
      k <- k + 1L
      forecast_parts[[k]] <- data.frame(
        stock = stock,
        model = type,
        scenario = scenario,
        year = years,
        period = ifelse(years <= origin, "historical overlap", "forecast"),
        low = quants[, 1],
        median = quants[, 2],
        high = quants[, 3]
      )
    }
    
    k <- k + 1L
    forecast_parts[[k]] <- data.frame(
      stock = stock,
      model = type,
      scenario = "SRR benchmark",
      year = expected_years,
      period = ifelse(
        expected_years <= origin, "historical overlap", "forecast"
      ),
      low = srr_quantiles[1],
      median = srr_quantiles[2],
      high = srr_quantiles[3]
    )
  }
}

historical <- bind_rows(histories)
forecast <- bind_rows(forecast_parts) |>
  mutate(
    model = factor(model, levels = sr_types),
    scenario = factor(
      scenario, levels = c("rcp45", "rcp85", "SRR benchmark")
    )
  )

ssb_used <- bind_rows(ssb_settings)
print(ssb_used)

write.csv(
  ssb_used, file.path("output", "recruitment_forecast_ssb.csv"),
  row.names = FALSE
)
write.csv(
  forecast, file.path("output", "recruitment_forecast_summary.csv"),
  row.names = FALSE
)

colours <- c(
  rcp45 = "#0072B2", rcp85 = "#D55E00", "SRR benchmark" = "#555555"
)
labels <- c(
  rcp45 = "EMSRR RCP 4.5",
  rcp85 = "EMSRR RCP 8.5",
  "SRR benchmark" = "SRR benchmark"
)

p <- ggplot() +
  geom_ribbon(
    data = filter(forecast, scenario == "SRR benchmark"),
    aes(x = year, ymin = low, ymax = high, fill = scenario),
    alpha = 0.12
  ) +
  geom_ribbon(
    data = filter(forecast, scenario != "SRR benchmark"),
    aes(x = year, ymin = low, ymax = high, fill = scenario),
    alpha = 0.20
  ) +
  geom_vline(
    data = ssb_used, aes(xintercept = assessment_year + 0.5),
    colour = "grey50", linetype = "dotted"
  ) +
  geom_line(
    data = forecast,
    aes(x = year, y = median, colour = scenario, linetype = scenario),
    linewidth = 0.7
  ) +
  geom_line(
    data = historical, aes(x = year, y = rec),
    colour = "grey25", linewidth = 0.4
  ) +
  geom_point(
    data = historical, aes(x = year, y = rec),
    colour = "grey25", size = 0.65, alpha = 0.7
  ) +
  facet_grid(stock ~ model, scales = "free_y") +
  scale_colour_manual(values = colours, labels = labels, name = "Projection") +
  scale_fill_manual(values = colours, labels = labels, name = "Projection") +
  scale_linetype_manual(
    values = c(rcp45 = "solid", rcp85 = "solid", "SRR benchmark" = "dashed"),
    labels = labels, name = "Projection"
  ) +
  labs(
    x = "Recruitment year",
    y = "Recruitment (assessment units)",
    title = "Historical and projected recruitment at fixed recent SSB",
    subtitle = paste(
      "SSB fixed at the mean of the last five assessment years;",
      "lines are predictive medians"
    ),
    caption = paste(
      "Shading: pointwise 95% prediction intervals.",
      "EMSRR: climate realizations and recruitment variability;",
      "SRR: bootstrap fits and recruitment variability.",
      "Dotted line: end of assessment.",
      "These are conditional projections, not stock forecasts."
    )
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(
  file.path(fig_dir, "recruitment_forecast.png"),
  p, width = 13, height = 9, dpi = 300
)
cat("Saved recruitment forecast figure and summary tables under output/.\n")