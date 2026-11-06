# ----------------------------------------------------------------------------- #
# Fit cod substock EMSRRs and evaluate selected models on a final-year holdout.
#
# For each substock and environmental covariate:
#   1. Bootstrap training years and calculate OOB predictions.
#   2. Fit both MARS forms (linpreds = FALSE and TRUE) per bootstrap.
#   3. Keep the form with lower OOB log-RMSE for each bootstrap.
#   4. Select the covariate with the lowest median winning OOB log-RMSE.
#   5. Fit the selected model on the training data and compare its holdout
#      performance with a Beverton-Holt model using SSB only.
#
# Run from the repository root after data_02_join_bottomT_with_substock.R.
# Input:  ./data/Cod_Env_data.rda (object env_stk)
# Output: ./results/model_01_fit_EMSRR/
# ----------------------------------------------------------------------------- #

rm(list = ls())
library(earth)
library(collapse)

# -------------------------- #
# Configuration
# -------------------------- #

input_file <- "./data/Cod_Env_data.rda"
output_dir <- "./results/model_01_fit_EMSRR"
ssb_column <- "ssb_lag1"
expected_substocks <- c("Southern", "Northwestern", "Viking")
holdout_n_years <- 10L # approximately a 75/25% split
n_boot <- 1000L
seed <- 2026L
earth_pmethod <- "exhaustive"

if (!file.exists(input_file)) {
  stop("Input data not found: ", input_file,
       ". Run data_02_join_bottomT_with_substock.R first.")
}
if (n_boot < 1L || holdout_n_years < 1L) {
  stop("n_boot and holdout_n_years must be positive.")
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(seed)

# Load without polluting the global environment.
loaded_data <- new.env(parent = emptyenv())
loaded_objects <- load(input_file, envir = loaded_data)
if (!"env_stk" %in% loaded_objects ||
    !is.data.frame(loaded_data$env_stk)) {
  stop("Expected a data frame named 'env_stk' in ", input_file, ".")
}
dat <- loaded_data$env_stk

required_columns <- c("year", "substock", "rec", ssb_column,
                      "variable", "value")
missing_columns <- setdiff(required_columns, names(dat))
if (length(missing_columns)) {
  stop("Missing required columns in env_stk: ",
       paste(missing_columns, collapse = ", "))
}
if (!is.numeric(dat$year) || !is.numeric(dat$rec) ||
    !is.numeric(dat[[ssb_column]]) || !is.numeric(dat$value)) {
  stop("year, rec, ", ssb_column, ", and value must be numeric.")
}

dat <- dat[, required_columns]
dat$variable <- as.character(dat$variable)
dat$substock <- as.character(dat$substock)

# -------------------------- #
# Model and metric helpers
# -------------------------- #

source("./functions/calc.partial.effects.R")

# define bounds for the BH parameter search
get_bh_bounds <- function(substock, rec, ssb) {
  bounds <- bh_bounds_by_substock[[substock]]
  
  # Fallbacks preserve the former lower bound for b and add broad finite
  # upper limits so BOBYQA can work in a bounded parameter space.
  if (is.null(bounds)) {
    bounds <- list(
      a = c(max(rec) * (1 + 1e-6), bh_upper_factor * max(rec)),
      b = c(0.01 * min(ssb), bh_upper_factor * max(ssb))
    )
  }
  
  for (parameter in c("a", "b")) {
    limit <- bounds[[parameter]]
    if (length(limit) != 2L || any(!is.finite(limit)) ||
        any(limit <= 0) || limit[1] >= limit[2]) {
      stop("Invalid positive finite bounds for BH parameter ",
           parameter, " in ", substock, ".")
    }
  }
  
  bounds
}


# Fit R = a * S / (S + b) by minimizing squared log-scale residuals.
# For any fixed b, the minimizing log(a) has a closed-form solution, leaving
# only a one-dimensional optimization over log(b).
fit_bh <- function(rec, ssb, bounds, n_starts = bh_n_starts) {
  keep <- is.finite(rec) & is.finite(ssb) & rec > 0 & ssb > 0
  rec <- rec[keep]
  ssb <- ssb[keep]
  
  if (length(rec) < 3L) {
    stop("At least three positive recruitment/SSB observations are required.")
  }
  if (n_starts < 2L) {
    stop("Use at least two BH starting points.")
  }
  
  lower <- log(c(a = bounds$a[1], b = bounds$b[1]))
  upper <- log(c(a = bounds$a[2], b = bounds$b[2]))
  
  # Optimize L1-loss on log-scale residuals, using log parameters to guarantee
  # positive a and b.
  objective <- function(log_par) {
    log_pred <- log_par[1] + log(ssb) -
      log(ssb + exp(log_par[2]))
    sum(abs(log(rec) - log_pred)) # L1-loss
  }
  
  # Data-informed center, following the FLR starting-value formulation.
  a0 <- max(quantile(rec, 0.75, na.rm = TRUE))
  b0 <- a0 / max(quantile(rec / ssb, 0.90, na.rm = TRUE))
  
  # Keep the informed start and draw the remaining starts around it.
  starts_raw <- cbind(
    a = abs(rnorm(n_starts - 1L, mean = a0, sd = 1.8 * a0)),
    b = pmax(
      bounds$b[1],
      abs(rnorm(n_starts - 1L, mean = b0, sd = 1.8 * b0))
    )
  )
  starts_raw <- rbind(c(a = a0, b = b0), starts_raw)
  
  # Ensure every start lies within the same bounds used by BOBYQA.
  starts_raw[, "a"] <- pmin(
    bounds$a[2], pmax(bounds$a[1], starts_raw[, "a"])
  )
  starts_raw[, "b"] <- pmin(
    bounds$b[2], pmax(bounds$b[1], starts_raw[, "b"])
  )
  
  # BOBYQA still optimizes log(a), log(b); the objective and bounds remain
  # unchanged.
  starts <- log(starts_raw)
  
  fits <- lapply(seq_len(nrow(starts)), function(i) {
    tryCatch(
      minqa::bobyqa(
        par = starts[i, ],
        fn = objective,
        lower = lower,
        upper = upper,
        control = list(maxfun = bh_maxfun)
      ),
      error = function(e) NULL
    )
  })
  
  valid <- which(!vapply(fits, is.null, logical(1)))
  valid <- valid[vapply(fits[valid], function(fit) {
    is.finite(fit$fval) && all(is.finite(fit$par))
  }, logical(1))]
  
  if (!length(valid)) {
    stop("All bounded multistart BOBYQA fits failed.")
  }
  
  converged <- valid[vapply(fits[valid], function(fit) {
    isTRUE(fit$ierr == 0)
  }, logical(1))]
  
  if (!length(converged)) {
    warning("No BOBYQA start reported normal convergence; using the best ",
            "finite solution.")
    converged <- valid
  }
  
  best_index <- converged[which.min(vapply(
    fits[converged],
    function(fit) fit$fval,
    numeric(1)
  ))]
  best <- fits[[best_index]]
  
  list(
    a = unname(exp(best$par[1])),
    b = unname(exp(best$par[2])),
    objective = best$fval,
    convergence = best$ierr,
    bounds = bounds,
    n_starts = nrow(starts),
    best_start = best_index,
    optimizer = best
  )
}

predict_bh <- function(model, ssb) {
  model$a * ssb / (ssb + model$b)
}

score_predictions <- function(observed, predicted) {
  if (!length(observed) || length(observed) != length(predicted)) {
    stop("Observed and predicted values must have the same non-zero length.")
  }
  if (any(!is.finite(observed)) || any(!is.finite(predicted)) ||
      any(observed <= 0) || any(predicted <= 0)) {
    stop("Scoring requires finite, positive observed and predicted values.")
  }
  
  data.frame(
    n = length(observed),
    rmse = sqrt(mean((observed - predicted)^2)),
    mae = mean(abs(observed - predicted)),
    log_rmse = sqrt(mean((log(observed) - log(predicted))^2)),
    log_bias = mean(log(predicted) - log(observed))
  )
}

fit_mars_residual <- function(log_residuals, covariate, linpreds) {
  earth::earth(
    log_residuals ~ x,
    data = data.frame(log_residuals = log_residuals, x = covariate),
    linpreds = linpreds, fast.k = 20, nprune = 20,
    pmethod = earth_pmethod
  )
}

predict_emsrr <- function(bh_model, mars_model, ssb, covariate) {
  baseline <- predict_bh(bh_model, ssb)
  
  # Match the existing helper's behavior: an intercept-only MARS fit adds no
  # environmental adjustment.
  if (length(stats::coef(mars_model)) <= 1L) {
    residual_effect <- rep(0, length(ssb))
  } else {
    residual_effect <- as.numeric(stats::predict(
      mars_model,
      newdata = data.frame(x = covariate)
    ))
  }
  
  exp(log(baseline) + residual_effect)
}

predict_emsrr_model <- function(model.obj, newdata) {
  predict_emsrr(
    bh_model = model.obj$bh_model,
    mars_model = model.obj$mars_model,
    ssb = newdata[[model.obj$ssb_name]],
    covariate = newdata[[model.obj$covariate_name]]
  )
}

# ------------------------------ #
# starting values for the BH-fit
# ------------------------------ #

bh_n_starts <- 10L
bh_maxfun <- 500L
bh_upper_factor <- 100

# Optional: enter biologically justified bounds in the same units as rec and ssb.
# Example structure:
# bh_bounds_by_substock <- list(
#   Southern = list(a = c(lower, upper), b = c(lower, upper)),
#   Northwestern = list(a = c(lower, upper), b = c(lower, upper)),
#   Viking = list(a = c(lower, upper), b = c(lower, upper))
# )
bh_bounds_by_substock <- list(
  Southern = NULL,
  Northwestern = NULL,
  Viking = NULL
)

# -------------------------- #
# Bootstrap selection
# -------------------------- #

candidate_rows <- list()
candidate_row_index <- 1L

for (substock in expected_substocks) {
  message("Processing substock: ", substock)
  
  stock_data <- dat[dat$substock == substock, , drop = FALSE]
  if (!nrow(stock_data)) {
    stop("No rows found for substock: ", substock)
  }
  
  # Recruitment and SSB repeat across the long-format covariate rows. Confirm
  # there is exactly one consistent stock observation per year.
  stock_by_year <- unique(stock_data[, c("year", "rec", ssb_column)])
  if (anyDuplicated(stock_by_year$year)) {
    stop("Conflicting recruitment/SSB values for a year in ", substock, ".")
  }
  
  base_ok <- is.finite(stock_by_year$year) &
    is.finite(stock_by_year$rec) &
    is.finite(stock_by_year[[ssb_column]]) &
    stock_by_year$rec > 0 &
    stock_by_year[[ssb_column]] > 0
  stock_by_year <- stock_by_year[base_ok, , drop = FALSE]
  available_years <- sort(unique(stock_by_year$year))
  
  if (length(available_years) <= holdout_n_years) {
    stop("Substock ", substock, " has only ", length(available_years),
         " usable years; need more than ", holdout_n_years, ".")
  }
  
  holdout_years <- tail(available_years, holdout_n_years)
  first_holdout_year <- min(holdout_years)
  
  training_stock <- stock_by_year[
    stock_by_year$year < first_holdout_year, ,
    drop = FALSE
  ]
  bh_bounds <- get_bh_bounds(
    substock = substock,
    rec = training_stock$rec,
    ssb = training_stock[[ssb_column]]
  )
  
  covariates <- sort(unique(stock_data$variable[
    !is.na(stock_data$variable) & nzchar(stock_data$variable)
  ]))
  
  if (!length(covariates)) {
    stop("No environmental covariates found for ", substock, ".")
  }
  
  for (variable_name in covariates) {
    covariate_data <- stock_data[
      stock_data$variable == variable_name,
      c("year", "rec", ssb_column, "value"),
      drop = FALSE
    ]
    
    valid <- complete.cases(covariate_data) &
      is.finite(covariate_data$year) &
      is.finite(covariate_data$rec) &
      is.finite(covariate_data[[ssb_column]]) &
      is.finite(covariate_data$value) &
      covariate_data$rec > 0 &
      covariate_data[[ssb_column]] > 0
    
    covariate_data <- covariate_data[valid, , drop = FALSE]
    covariate_data <- covariate_data[
      order(covariate_data$year), , drop = FALSE
    ]
    
    if (anyDuplicated(covariate_data$year)) {
      stop("More than one row per year for ", substock, " / ",
           variable_name, ".")
    }
    
    training_data <- covariate_data[
      covariate_data$year < first_holdout_year, , drop = FALSE
    ]
    
    if (nrow(training_data) < 5L) {
      warning("Skipping ", substock, " / ", variable_name,
              ": fewer than five complete training years.")
      next
    }
    
    message("  Covariate: ", variable_name)
    
    for (bootstrap_id in seq_len(n_boot)) {
      bootstrap_index <- sample.int(
        n = nrow(training_data),
        size = nrow(training_data),
        replace = TRUE
      )
      oob_index <- setdiff(seq_len(nrow(training_data)), unique(bootstrap_index))
      
      # Always create a row for each form so failures and the winning form are
      # visible in the metrics table.
      bootstrap_rows <- vector("list", 2L)
      names(bootstrap_rows) <- c("FALSE", "TRUE")
      
      if (!length(oob_index)) {
        for (linear_form in c(FALSE, TRUE)) {
          bootstrap_rows[[as.character(linear_form)]] <- data.frame(
            substock = substock,
            variable = variable_name,
            bootstrap = bootstrap_id,
            linpreds = linear_form,
            model_form = if (linear_form) "linear" else "MARS",
            n_train = length(bootstrap_index),
            n_oob = 0L,
            rmse = NA_real_,
            mae = NA_real_,
            log_rmse = NA_real_,
            log_bias = NA_real_,
            status = "no_oob",
            error = NA_character_,
            selected = FALSE
          )
        }
      } else {
        fit_result <- tryCatch({
          bootstrap_data <- training_data[bootstrap_index, , drop = FALSE]
          oob_data <- training_data[oob_index, , drop = FALSE]
          
          bh_model <- fit_bh(
            rec = bootstrap_data$rec,
            ssb = bootstrap_data[[ssb_column]],
            bounds = bh_bounds
          )
          fitted_baseline <- predict_bh(
            bh_model,
            bootstrap_data[[ssb_column]]
          )
          log_residuals <- log(bootstrap_data$rec) -
            log(fitted_baseline)
          
          list(
            bootstrap_data = bootstrap_data,
            oob_data = oob_data,
            bh_model = bh_model,
            log_residuals = log_residuals
          )
        }, error = function(e) e)
        
        for (linear_form in c(FALSE, TRUE)) {
          model_form <- if (linear_form) "linear" else "MARS"
          result <- tryCatch({
            if (inherits(fit_result, "error")) {
              stop(conditionMessage(fit_result))
            }
            
            mars_model <- fit_mars_residual(
              log_residuals = fit_result$log_residuals,
              covariate = fit_result$bootstrap_data$value,
              linpreds = linear_form
            )
            predicted <- predict_emsrr(
              bh_model = fit_result$bh_model,
              mars_model = mars_model,
              ssb = fit_result$oob_data[[ssb_column]],
              covariate = fit_result$oob_data$value
            )
            score <- score_predictions(
              observed = fit_result$oob_data$rec,
              predicted = predicted
            )
            
            list(score = score, error = NA_character_)
          }, error = function(e) {
            list(score = NULL, error = conditionMessage(e))
          })
          
          score <- result$score
          bootstrap_rows[[as.character(linear_form)]] <- data.frame(
            substock = substock,
            variable = variable_name,
            bootstrap = bootstrap_id,
            linpreds = linear_form,
            model_form = model_form,
            n_train = length(bootstrap_index),
            n_oob = if (is.null(score)) length(oob_index) else score$n,
            rmse = if (is.null(score)) NA_real_ else score$rmse,
            mae = if (is.null(score)) NA_real_ else score$mae,
            log_rmse = if (is.null(score)) NA_real_ else score$log_rmse,
            log_bias = if (is.null(score)) NA_real_ else score$log_bias,
            status = if (is.null(score)) "error" else "ok",
            error = result$error,
            selected = FALSE
          )
        }
      }
      
      bootstrap_rows <- do.call(rbind, bootstrap_rows)
      valid_forms <- which(bootstrap_rows$status == "ok")
      
      if (length(valid_forms)) {
        # Lowest OOB log-RMSE wins; exact ties favor linear.
        winner <- valid_forms[order(
          bootstrap_rows$log_rmse[valid_forms],
          !bootstrap_rows$linpreds[valid_forms]
        )[1]]
        bootstrap_rows$selected[winner] <- TRUE
      }
      
      candidate_rows[[candidate_row_index]] <- bootstrap_rows
      candidate_row_index <- candidate_row_index + 1L
    }
  }
}

if (!length(candidate_rows)) {
  stop("No bootstrap evaluations were completed.")
}

candidate_metrics <- do.call(rbind, candidate_rows)
row.names(candidate_metrics) <- NULL
utils::write.csv(
  candidate_metrics,
  file.path(output_dir, "BH_bootstrap_candidate_metrics.csv"),
  row.names = FALSE
)

baseline = candidate_metrics[grepl("bottomT_mean",candidate_metrics$variable) &
                  grepl("linear",candidate_metrics$model_form),] |> fgroup_by(substock,variable,model_form) |> 
                  fsummarise(log_rmse = median(log_rmse))

ggplot(candidate_metrics,aes(x = variable,y = log_rmse,color = model_form)) + 
  geom_hline(data = baseline,aes(yintercept = log_rmse),lty = 2,color = "black") + 
  geom_boxplot(fill = NA) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~substock,scale = "free_x")


# -------------------------- #
# Select one covariate per substock
# -------------------------- #

selected_models <- list()
selection_rows <- list()
holdout_metric_rows <- list()
holdout_prediction_rows <- list()
partial_effects_by_substock <- list()

for (substock in expected_substocks) {
  substock_metrics <- candidate_metrics[
    candidate_metrics$substock == substock, ,
    drop = FALSE
  ]
  winners <- substock_metrics[
    substock_metrics$selected & substock_metrics$status == "ok", ,
    drop = FALSE
  ]
  
  if (!nrow(winners)) {
    stop("No successful OOB bootstrap winners for ", substock, ".")
  }
  
  variable_summary <- do.call(rbind, lapply(
    split(winners, winners$variable),
    function(x) {
      data.frame(
        variable = x$variable[1],
        n_wins = nrow(x),
        median_winning_oob_log_rmse = median(x$log_rmse),
        mean_winning_oob_log_rmse = mean(x$log_rmse),
        linear_win_fraction = mean(x$linpreds),
        stringsAsFactors = FALSE
      )
    }
  ))
  
  variable_summary <- variable_summary[
    order(variable_summary$median_winning_oob_log_rmse,
          variable_summary$variable),
    , drop = FALSE
  ]
  selected_variable <- variable_summary$variable[1]
  
  selection_rows[[substock]] <- data.frame(
    substock = substock,
    variable_summary,
    selected = variable_summary$variable == selected_variable,
    row.names = NULL
  )
  
  selected_winners <- winners[
    winners$variable == selected_variable, ,
    drop = FALSE
  ]
  linear_wins <- sum(selected_winners$linpreds)
  mars_wins <- sum(!selected_winners$linpreds)
  
  if (linear_wins > mars_wins) {
    final_linpreds <- TRUE
  } else if (mars_wins > linear_wins) {
    final_linpreds <- FALSE
  } else {
    # Compare OOB log-RMSE spread on bootstraps where both forms succeeded.
    candidates <- substock_metrics[
      substock_metrics$variable == selected_variable &
        substock_metrics$status == "ok" &
        is.finite(substock_metrics$log_rmse), ,
      drop = FALSE
    ]
    
    linear <- candidates[candidates$linpreds, c("bootstrap", "log_rmse")]
    mars <- candidates[!candidates$linpreds, c("bootstrap", "log_rmse")]
    paired <- merge(linear, mars, by = "bootstrap",
                    suffixes = c("_linear", "_mars"))
    
    if (nrow(paired) < 2L) {
      warning("Too few paired bootstraps to compare spread for ",
              substock, "; choosing linear.")
      final_linpreds <- TRUE
    } else {
      linear_spread <- IQR(paired$log_rmse_linear)
      mars_spread <- IQR(paired$log_rmse_mars)
      # Lower IQR wins; equal spread favors linear.
      final_linpreds <- linear_spread <= mars_spread
    }
  }
  
  message("Selected for ", substock, ": ", selected_variable,
          " (linpreds = ", final_linpreds, ")")
  
  # Recompute the split and BH bounds for this substock.
  stock_for_split <- unique(
    dat[dat$substock == substock, c("year", "rec", ssb_column), drop = FALSE]
  )
  valid_stock <- is.finite(stock_for_split$year) &
    is.finite(stock_for_split$rec) &
    is.finite(stock_for_split[[ssb_column]]) &
    stock_for_split$rec > 0 &
    stock_for_split[[ssb_column]] > 0
  stock_for_split <- stock_for_split[valid_stock, , drop = FALSE]
  
  available_years <- sort(unique(stock_for_split$year))
  holdout_years <- tail(available_years, holdout_n_years)
  first_holdout_year <- min(holdout_years)
  
  training_stock <- stock_for_split[
    stock_for_split$year < first_holdout_year, , drop = FALSE
  ]
  bh_bounds <- get_bh_bounds(
    substock = substock,
    rec = training_stock$rec,
    ssb = training_stock[[ssb_column]]
  )
  
  
  # Fit final EMSRR using all usable training observations for this covariate.
  selected_data <- stock_data <- dat[
    dat$substock == substock & dat$variable == selected_variable,
    c("year", "rec", ssb_column, "value"),
    drop = FALSE
  ]
  selected_valid <- complete.cases(selected_data) &
    is.finite(selected_data$year) &
    is.finite(selected_data$rec) &
    is.finite(selected_data[[ssb_column]]) &
    is.finite(selected_data$value) &
    selected_data$rec > 0 &
    selected_data[[ssb_column]] > 0
  selected_data <- selected_data[selected_valid, , drop = FALSE]
  selected_training <- selected_data[
    selected_data$year < first_holdout_year, , drop = FALSE
  ]
  selected_holdout <- selected_data[
    selected_data$year %in% holdout_years, , drop = FALSE
  ]
  
  if (anyDuplicated(selected_training$year) ||
      anyDuplicated(selected_holdout$year)) {
    stop("Duplicate years found for selected covariate ", selected_variable,
         " in ", substock, ".")
  }
  if (nrow(selected_holdout) < 1L) {
    stop("No usable holdout observations for ", substock, " / ",
         selected_variable, ".")
  }
  
  final_fit <- tryCatch({
    bh_model <- fit_bh(
      rec = selected_training$rec,
      ssb = selected_training[[ssb_column]],
      bounds = bh_bounds
    )
    baseline <- predict_bh(bh_model, selected_training[[ssb_column]])
    log_residuals <- log(selected_training$rec) - log(baseline)
    mars_model <- fit_mars_residual(
      log_residuals = log_residuals,
      covariate = selected_training$value,
      linpreds = final_linpreds
    )
    
    list(
      substock = substock,
      variable = selected_variable,
      linpreds = final_linpreds,
      model_form = if (final_linpreds) "linear" else "MARS",
      bh_model = bh_model,
      mars_model = mars_model,
      n_training = nrow(selected_training)
    )
  }, error = function(e) {
    stop("Final EMSRR fit failed for ", substock, " / ",
         selected_variable, ": ", conditionMessage(e))
  })
  
  # Provide the predictor names and adapter expected by calc.partial.effects().
  final_fit$ssb_name <- ssb_column
  final_fit$covariate_name <- selected_variable
  final_fit$vars <- c(ssb_column, selected_variable)
  
  # Calculate and plot partial effects on training data only; the holdout remains
  # reserved for external performance evaluation.
  partial_x <- data.frame(
    c(selected_training[[ssb_column]],selected_holdout[[ssb_column]]),
    c(selected_training$value,selected_holdout$value),
    check.names = FALSE
  )
  names(partial_x) <- c(ssb_column, selected_variable)
  
  partial_effects <- calc.partial.effects(
    model.obj = final_fit,
    # use all data (train + holdout) for partial dependence
    y = c(selected_training$rec,selected_holdout$rec),
    x = partial_x,
    predict.func = predict_emsrr_model
  )
  
  partial_plot <- partial_effects$plot +
    ggplot2::ggtitle(paste(substock, "-", selected_variable))
  
  ggplot2::ggsave(
    filename = file.path(
      output_dir,
      paste0("BH_partial_effects_", substock, ".png")
    ),
    plot = partial_plot,
    width = 8,
    height = 5,
    dpi = 300
  )
  
  partial_effects_by_substock[[substock]] <- partial_effects
  
  # Fit the SSB-only benchmark using all usable stock training years.
  stock_rows <- dat[dat$substock == substock, , drop = FALSE]
  stock_by_year <- unique(stock_rows[, c("year", "rec", ssb_column)])
  stock_valid <- complete.cases(stock_by_year) &
    is.finite(stock_by_year$year) &
    is.finite(stock_by_year$rec) &
    is.finite(stock_by_year[[ssb_column]]) &
    stock_by_year$rec > 0 &
    stock_by_year[[ssb_column]] > 0
  stock_by_year <- stock_by_year[stock_valid, , drop = FALSE]
  
  baseline_training <- stock_by_year[
    stock_by_year$year < first_holdout_year, , drop = FALSE
  ]
  baseline_holdout <- stock_by_year[
    stock_by_year$year %in% holdout_years, , drop = FALSE
  ]
  baseline_fit <- fit_bh(
    rec = baseline_training$rec,
    ssb = baseline_training[[ssb_column]],
    bounds = bh_bounds
  )
  
  emsrr_prediction <- predict_emsrr(
    bh_model = final_fit$bh_model,
    mars_model = final_fit$mars_model,
    ssb = selected_holdout[[ssb_column]],
    covariate = selected_holdout$value
  )
  naive_prediction <- predict_bh(
    baseline_fit,
    selected_holdout[[ssb_column]]
  )
  
  # plot predictions
  # Predictions on the same years used for the selected EMSRR.
  plot_input <- rbind(selected_training, selected_holdout)
  plot_input <- plot_input[order(plot_input$year), , drop = FALSE]
  
  plot_input$period <- ifelse(
    plot_input$year < first_holdout_year, "History", "Holdout"
  )
  plot_input$period <- factor(plot_input$period, levels = c("History", "Holdout"))
  
  # Break lines at missing years and at the history/holdout boundary.
  plot_input$segment <- cumsum(c(
    TRUE,
    diff(plot_input$year) != 1 |
      plot_input$period[-1] != plot_input$period[-nrow(plot_input)]
  ))
  
  emsrr_all <- predict_emsrr(
    bh_model = final_fit$bh_model,
    mars_model = final_fit$mars_model,
    ssb = plot_input[[ssb_column]],
    covariate = plot_input$value
  )
  bh_all <- predict_bh(baseline_fit, plot_input[[ssb_column]])
  
  plot_data <- rbind(
    data.frame(plot_input[, c("year", "period", "segment")],
               series = "Observed", recruitment = plot_input$rec),
    data.frame(plot_input[, c("year", "period", "segment")],
               series = "EMSRR", recruitment = emsrr_all),
    data.frame(plot_input[, c("year", "period", "segment")],
               series = "SSB-only BH", recruitment = bh_all)
  )
  plot_data$series <- factor(
    plot_data$series, levels = c("Observed", "EMSRR", "SSB-only BH")
  )
  
  recruitment_plot <- ggplot(plot_data, aes(x = year, y = recruitment)) +
    annotate(
      "rect",
      xmin = first_holdout_year - 0.5, xmax = max(holdout_years) + 0.5,
      ymin = -Inf, ymax = Inf, fill = "grey90"
    ) +
    geom_line(
      aes(colour = series, linetype = period,
          group = interaction(series, segment)),
      linewidth = 0.8
    ) +
    geom_point(
      data = plot_data[plot_data$series == "Observed", , drop = FALSE],
      colour = "black", size = 1.5
    ) +
    geom_vline(xintercept = first_holdout_year - 0.5, colour = "grey50") +
    scale_colour_manual(values = c(
      "Observed" = "black", "EMSRR" = "#0072B2", "SSB-only BH" = "#D55E00"
    )) +
    scale_linetype_manual(values = c("History" = "solid", "Holdout" = "dashed")) +
    labs(
      title = paste("Recruitment:", substock),
      subtitle = paste("Selected environmental covariate:", selected_variable),
      x = "Year", y = "Recruitment",
      colour = NULL, linetype = NULL
    ) +
    theme_bw()
  
  print(recruitment_plot)
  ggsave(
    filename = file.path(output_dir, paste0("BH_Prediction_rec_ts_", substock, ".png")),
    plot = recruitment_plot, width = 9, height = 5, dpi = 300
  )
  
  observed <- selected_holdout$rec
  emsrr_score <- score_predictions(observed, emsrr_prediction)
  naive_score <- score_predictions(observed, naive_prediction)
  
  holdout_metrics <- rbind(
    cbind(model = "EMSRR", emsrr_score),
    cbind(model = "BH_SSB_only", naive_score)
  )
  holdout_metrics$substock <- substock
  holdout_metrics$selected_variable <- selected_variable
  holdout_metrics$linpreds <- final_linpreds
  holdout_metrics$delta_log_rmse_vs_bh <-
    holdout_metrics$log_rmse -
    holdout_metrics$log_rmse[holdout_metrics$model == "BH_SSB_only"][1]
  
  holdout_metric_rows[[substock]] <- holdout_metrics
  holdout_prediction_rows[[substock]] <- data.frame(
    substock = substock,
    year = selected_holdout$year,
    observed = observed,
    emsrr_predicted = emsrr_prediction,
    bh_ssb_predicted = naive_prediction,
    selected_variable = selected_variable,
    linpreds = final_linpreds
  )
  selected_models[[substock]] <- final_fit
}

# -------------------------- #
# Save results
# -------------------------- #

selection_summary <- do.call(rbind, selection_rows)
holdout_metrics <- do.call(rbind, holdout_metric_rows)
holdout_predictions <- do.call(rbind, holdout_prediction_rows)

utils::write.csv(
  selection_summary,
  file.path(output_dir, "BH_variable_selection_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  holdout_metrics,
  file.path(output_dir, "BH_holdout_metrics.csv"),
  row.names = FALSE
)
utils::write.csv(
  holdout_predictions,
  file.path(output_dir, "BH_holdout_predictions.csv"),
  row.names = FALSE
)
saveRDS(selected_models, file.path(output_dir, "BH_selected_emsrr_models.rds"))

saveRDS(
  partial_effects_by_substock,
  file.path(output_dir, "BH_partial_effects.rds")
)

print(selection_summary)
print(holdout_metrics)