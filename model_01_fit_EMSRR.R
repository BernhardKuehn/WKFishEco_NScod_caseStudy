# ----------------------------------------------------------------------------- #
# Fit cod substock stock-recruitment models and evaluate the final 10 years.
#
# Models fitted independently for Southern, Northwestern, and Viking:
#   1. Beverton-Holt with a MARS effect of lagged SST on log-scale residuals
#      (the EMSRR);
#   2. Beverton-Holt using SSB only; and
#   3. a hockey-stick stock-recruitment model using SSB only.
#
# Run from the repository root after data_02_join_SST_with_substock.R.
# Input:  ./data/Cod_Env_data.rda (object env_stk)
# Output: ./results/model_01_fit_EMSRR/
# ----------------------------------------------------------------------------- #

input_file <- "./data/Cod_Env_data.rda"
output_dir <- "./results/model_01_fit_EMSRR"
ssb_column <- "ssb_lag1"
sst_column <- "SST_lag1"
expected_substocks <- c("Southern", "Northwestern", "Viking")
holdout_n_years <- 10L

if (!file.exists(input_file)) {
  stop("Input data not found: ", input_file,
       ". Run data_02_join_SST_with_substock.R first.")
}
if (!requireNamespace("earth", quietly = TRUE)) {
  stop("Package 'earth' is required. Install it with install.packages('earth').")
}

loaded_data <- new.env(parent = emptyenv())
loaded_objects <- load(input_file, envir = loaded_data)
if (!"env_stk" %in% loaded_objects) {
  stop("Expected object 'env_stk' was not found in ", input_file, ".")
}
if (!is.data.frame(loaded_data$env_stk)) {
  stop("Object 'env_stk' must be a data frame.")
}

required_columns <- c("year", "rec", ssb_column, sst_column, "substock")
missing_columns <- setdiff(required_columns, names(loaded_data$env_stk))
if (length(missing_columns) > 0L) {
  stop("Missing required column(s) in env_stk: ",
       paste(missing_columns, collapse = ", "))
}

dat <- loaded_data$env_stk[, required_columns]
dat$year <- suppressWarnings(as.integer(as.character(dat$year)))
dat$rec <- suppressWarnings(as.numeric(as.character(dat$rec)))
dat[[ssb_column]] <- suppressWarnings(as.numeric(as.character(dat[[ssb_column]])))
dat[[sst_column]] <- suppressWarnings(as.numeric(as.character(dat[[sst_column]])))
dat$substock <- as.character(dat$substock)

if (anyNA(dat$year) || anyNA(dat$substock) || any(!nzchar(dat$substock))) {
  stop("env_stk contains missing/invalid years or substock labels.")
}

# Normalize the variants that the data-preparation scripts may produce.
dat$substock[grepl("northwest", dat$substock, ignore.case = TRUE)] <- "Northwestern"
dat$substock[grepl("south", dat$substock, ignore.case = TRUE)] <- "Southern"
dat$substock[grepl("viking", dat$substock, ignore.case = TRUE)] <- "Viking"
unexpected_substocks <- setdiff(unique(dat$substock), expected_substocks)
if (length(unexpected_substocks) > 0L) {
  stop("Unexpected substock label(s): ",
       paste(unexpected_substocks, collapse = ", "))
}
if (!setequal(unique(dat$substock), expected_substocks)) {
  stop("env_stk must contain all three substocks: ",
       paste(expected_substocks, collapse = ", "))
}
if (anyDuplicated(dat[c("year", "substock")])) {
  stop("env_stk must contain at most one row per year and substock.")
}

# Define the holdout from calendar years before removing incomplete rows, so a
# missing SST/recruitment value cannot silently move the evaluation period.
available_years <- sort(unique(dat$year))
if (length(available_years) <= holdout_n_years) {
  stop("Need more than ", holdout_n_years, " distinct years to fit and evaluate.")
}
holdout_years <- tail(available_years, holdout_n_years)
first_holdout_year <- min(holdout_years)

# Use matched complete, positive observations for all three models. Positive
# recruitment and SSB are required for the log-scale fitting objective.
complete <- stats::complete.cases(dat[, required_columns]) &
  is.finite(dat$rec) & is.finite(dat[[ssb_column]]) &
  is.finite(dat[[sst_column]]) & dat$rec > 0 & dat[[ssb_column]] > 0
dat <- dat[complete, , drop = FALSE]
names(dat)[names(dat) == ssb_column] <- "ssb"
names(dat)[names(dat) == sst_column] <- "SST_lag1"
dat <- dat[order(dat$substock, dat$year), , drop = FALSE]

fit_beverton_holt <- function(training_data) {
  objective <- function(log_parameters) {
    a <- exp(log_parameters[1L])
    b <- exp(log_parameters[2L])
    prediction <- a * training_data$ssb / (training_data$ssb + b)
    if (any(!is.finite(prediction)) || any(prediction <= 0)) {
      return(.Machine$double.xmax)
    }
    sum((log(training_data$rec) - log(prediction))^2)
  }

  a_starts <- max(training_data$rec) * c(1.05, 2, 5)
  b_starts <- c(stats::median(training_data$ssb),
                stats::quantile(training_data$ssb, 0.75, names = FALSE),
                max(training_data$ssb))
  starts <- expand.grid(a = a_starts, b = b_starts)
  lower <- log(c(a = max(.Machine$double.eps, max(training_data$rec) * 1e-8),
                 b = max(.Machine$double.eps, min(training_data$ssb) * 1e-6)))
  upper <- log(c(a = max(training_data$rec) * 1e4,
                 b = max(training_data$ssb) * 1e4))

  fits <- lapply(seq_len(nrow(starts)), function(i) {
    stats::optim(
      par = log(c(starts$a[i], starts$b[i])),
      fn = objective,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      control = list(maxit = 10000)
    )
  })
  converged <- which(vapply(
    fits,
    function(fit) fit$convergence == 0L && is.finite(fit$value),
    logical(1)
  ))
  if (length(converged) == 0L) {
    stop("Beverton-Holt optimization did not converge.")
  }
  best <- fits[[converged[which.min(vapply(
    fits[converged], `[[`, numeric(1), "value"
  ))]]]
  parameters <- exp(best$par)
  names(parameters) <- c("a", "b")
  parameters
}

predict_beverton_holt <- function(parameters, ssb) {
  parameters[["a"]] * ssb / (ssb + parameters[["b"]])
}

fit_hockey_stick <- function(training_data) {
  # For each candidate breakpoint, the optimal log-scale slope is available
  # analytically; evaluating observed and log-spaced breakpoints is stable and
  # avoids relying on a derivative-based optimizer at the hockey-stick kink.
  ssb_range <- range(training_data$ssb)
  candidates <- sort(unique(c(
    training_data$ssb,
    exp(seq(log(ssb_range[1L]), log(ssb_range[2L]), length.out = 501L))
  )))
  fits <- lapply(candidates, function(breakpoint) {
    shape <- pmin(training_data$ssb, breakpoint)
    a <- exp(mean(log(training_data$rec) - log(shape)))
    prediction <- a * shape
    list(
      a = a,
      b = breakpoint,
      objective = sum((log(training_data$rec) - log(prediction))^2)
    )
  })
  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "objective"))]]
  c(a = best$a, b = best$b)
}

predict_hockey_stick <- function(parameters, ssb) {
  parameters[["a"]] * pmin(ssb, parameters[["b"]])
}

score_predictions <- function(observed, predicted) {
  data.frame(
    n = length(observed),
    rmse = sqrt(mean((observed - predicted)^2)),
    mae = mean(abs(observed - predicted)),
    log_rmse = sqrt(mean((log(observed) - log(predicted))^2)),
    log_bias = mean(log(predicted) - log(observed))
  )
}

model_objects <- setNames(vector("list", length(expected_substocks)),
                          expected_substocks)
prediction_rows <- list()
metric_rows <- list()
parameter_rows <- list()

for (substock in expected_substocks) {
  stock_data <- dat[dat$substock == substock, , drop = FALSE]
  training_data <- stock_data[stock_data$year < first_holdout_year, , drop = FALSE]
  holdout_data <- stock_data[stock_data$year %in% holdout_years, , drop = FALSE]

  if (nrow(training_data) < 10L) {
    stop(substock, " has fewer than 10 complete training observations.")
  }
  if (nrow(holdout_data) == 0L) {
    stop(substock, " has no complete observations in the final-10-year holdout.")
  }

  bh_parameters <- fit_beverton_holt(training_data)
  hs_parameters <- fit_hockey_stick(training_data)
  bh_training_prediction <- predict_beverton_holt(bh_parameters, training_data$ssb)
  log_residuals <- log(training_data$rec) - log(bh_training_prediction)

  mars_data <- data.frame(
    log_resid = log_residuals,
    SST_lag1 = training_data$SST_lag1
  )
  mars_model <- earth::earth(
    log_resid ~ SST_lag1,
    data = mars_data,
    degree = 1L,
    pmethod = "backward",
    nk = min(21L, nrow(training_data) - 1L),
    nprune = min(5L, nrow(training_data) - 1L)
  )

  prediction_table <- data.frame(
    year = holdout_data$year,
    substock = substock,
    observed = holdout_data$rec,
    BH_MARS = exp(
      log(predict_beverton_holt(bh_parameters, holdout_data$ssb)) +
        as.numeric(stats::predict(
          mars_model,
          newdata = data.frame(SST_lag1 = holdout_data$SST_lag1)
        ))
    ),
    BH_SSB_only = predict_beverton_holt(bh_parameters, holdout_data$ssb),
    HS_SSB_only = predict_hockey_stick(hs_parameters, holdout_data$ssb)
  )
  if (any(!is.finite(as.matrix(prediction_table[, c(
    "BH_MARS", "BH_SSB_only", "HS_SSB_only"
  )])))) {
    stop("Non-finite holdout prediction for ", substock, ".")
  }

  long_predictions <- do.call(rbind, lapply(
    c("BH_MARS", "BH_SSB_only", "HS_SSB_only"),
    function(model_name) {
      data.frame(
        year = prediction_table$year,
        substock = substock,
        model = model_name,
        observed = prediction_table$observed,
        predicted = prediction_table[[model_name]]
      )
    }
  ))
  prediction_rows[[substock]] <- long_predictions

  for (model_name in c("BH_MARS", "BH_SSB_only", "HS_SSB_only")) {
    scores <- score_predictions(
      prediction_table$observed,
      prediction_table[[model_name]]
    )
    metric_rows[[length(metric_rows) + 1L]] <- data.frame(
      substock = substock,
      model = model_name,
      scores
    )
  }

  parameter_rows[[substock]] <- data.frame(
    substock = substock,
    model = c("BH_MARS", "BH_SSB_only", "HS_SSB_only"),
    a = c(bh_parameters[["a"]], bh_parameters[["a"]], hs_parameters[["a"]]),
    b = c(bh_parameters[["b"]], bh_parameters[["b"]], hs_parameters[["b"]]),
    n_training = nrow(training_data),
    n_holdout = nrow(holdout_data)
  )
  model_objects[[substock]] <- list(
    BevertonHolt_parameters = bh_parameters,
    MARS_environment_model = mars_model,
    HockeyStick_parameters = hs_parameters,
    training_years = range(training_data$year),
    holdout_years = sort(unique(holdout_data$year))
  )
}

holdout_predictions <- do.call(rbind, prediction_rows)
holdout_metrics <- do.call(rbind, metric_rows)
fitted_parameters <- do.call(rbind, parameter_rows)

# Positive improvement means lower holdout log-RMSE for the EMSRR.
comparison_rows <- lapply(expected_substocks, function(substock) {
  scores <- holdout_metrics[holdout_metrics$substock == substock, ]
  emsrr_error <- scores$log_rmse[scores$model == "BH_MARS"]
  bh_error <- scores$log_rmse[scores$model == "BH_SSB_only"]
  hs_error <- scores$log_rmse[scores$model == "HS_SSB_only"]
  data.frame(
    substock = substock,
    BH_MARS_vs_BH_SSB_only_improvement_pct =
      if (bh_error > 0) 100 * (1 - emsrr_error / bh_error) else NA_real_,
    BH_MARS_vs_HS_SSB_only_improvement_pct =
      if (hs_error > 0) 100 * (1 - emsrr_error / hs_error) else NA_real_
  )
})
holdout_comparison <- do.call(rbind, comparison_rows)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  holdout_predictions,
  file.path(output_dir, "holdout_predictions.csv"),
  row.names = FALSE
)
utils::write.csv(
  holdout_metrics,
  file.path(output_dir, "holdout_metrics.csv"),
  row.names = FALSE
)
utils::write.csv(
  holdout_comparison,
  file.path(output_dir, "emsrr_vs_native_comparison.csv"),
  row.names = FALSE
)
utils::write.csv(
  fitted_parameters,
  file.path(output_dir, "fitted_parameters.csv"),
  row.names = FALSE
)
saveRDS(
  list(
    holdout_years = holdout_years,
    first_holdout_year = first_holdout_year,
    training_fraction = mean(dat$year < first_holdout_year),
    model_objects = model_objects,
    holdout_predictions = holdout_predictions,
    holdout_metrics = holdout_metrics,
    holdout_comparison = holdout_comparison,
    fitted_parameters = fitted_parameters,
    session_info = utils::capture.output(utils::sessionInfo())
  ),
  file.path(output_dir, "model_fit_results.rds")
)

message(
  "Fitting complete. Training through ", first_holdout_year - 1L,
  "; holdout years: ", min(holdout_years), "-", max(holdout_years),
  ". Results saved to ", output_dir, "."
)
print(holdout_comparison, row.names = FALSE)