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
# Run from the repository root after data_02_join_SST_with_substock.R.
# Input:  ./data/Cod_Env_data.rda (object env_stk)
# Output: ./results/model_01_fit_EMSRR/
# ----------------------------------------------------------------------------- #

rm(list = ls())
library(earth)

# -------------------------- #
# Configuration
# -------------------------- #

input_file <- "./data/Cod_Env_data.rda"
output_dir <- "./results/model_01_fit_EMSRR"
ssb_column <- "ssb_lag1"
expected_substocks <- c("Southern", "Northwestern", "Viking")
holdout_n_years <- 10L
n_boot <- 1000L
seed <- 2026L
earth_pmethod <- "exhaustive"

if (!file.exists(input_file)) {
  stop("Input data not found: ", input_file,
       ". Run data_02_join_SST_with_substock.R first.")
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

# Fit R = a * S / (S + b) by minimizing squared log-scale residuals.
# For any fixed b, the minimizing log(a) has a closed-form solution, leaving
# only a one-dimensional optimization over log(b).
fit_bh <- function(rec, ssb) {
  keep <- is.finite(rec) & is.finite(ssb) & rec > 0 & ssb > 0
  rec <- rec[keep]
  ssb <- ssb[keep]
  
  if (length(rec) < 3L) {
    stop("At least three positive recruitment/SSB observations are required.")
  }
  
  log_b_bounds <- log(c(min(ssb) * 1e-6, max(ssb) * 1e6))
  if (any(!is.finite(log_b_bounds)) ||
      log_b_bounds[1] >= log_b_bounds[2]) {
    stop("Could not construct finite Beverton-Holt parameter bounds.")
  }
  
  objective <- function(log_b) {
    b <- exp(log_b)
    log_ratio <- log(ssb) - log(ssb + b)
    log_a <- mean(log(rec) - log_ratio)
    sum((log(rec) - log_a - log_ratio)^2)
  }
  
  opt <- stats::optimize(objective, interval = log_b_bounds)
  b <- exp(opt$minimum)
  log_ratio <- log(ssb) - log(ssb + b)
  a <- exp(mean(log(rec) - log_ratio))
  
  list(a = a, b = b)
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
    linpreds = linpreds,
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
            ssb = bootstrap_data[[ssb_column]]
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
      
      # Lower OOB log-RMSE wins. If exactly tied, the first row (linpreds=FALSE)
      # wins deterministically.
      if (length(valid_forms)) {
        winner <- valid_forms[
          which.min(bootstrap_rows$log_rmse[valid_forms])
        ]
        bootstrap_rows$selected[winner] <- TRUE
      }
      
      candidate_rows[[candidate_row_index]] <- bootstrap_rows
      candidate_row_index <- candidate_row_index + 1L
    }
  }http://127.0.0.1:33392/graphics/plot_zoom_png?width=2170&height=1092
}

if (!length(candidate_rows)) {
  stop("No bootstrap evaluations were completed.")
}

candidate_metrics <- do.call(rbind, candidate_rows)
row.names(candidate_metrics) <- NULL
utils::write.csv(
  candidate_metrics,
  file.path(output_dir, "bootstrap_candidate_metrics.csv"),
  row.names = FALSE
)

baseline = candidate_metrics[grepl("SST_mean",candidate_metrics$variable) &
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
    # Break a frequency tie by the lower median OOB error for the form.
    selected_candidates <- substock_metrics[
      substock_metrics$variable == selected_variable &
        substock_metrics$status == "ok", ,
      drop = FALSE
    ]
    form_medians <- tapply(
      selected_candidates$log_rmse,
      selected_candidates$linpreds,
      median
    )
    final_linpreds <- as.logical(names(which.min(form_medians)))
  }
  
  message("Selected for ", substock, ": ", selected_variable,
          " (linpreds = ", final_linpreds, ")")
  
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
      ssb = selected_training[[ssb_column]]
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
    ssb = baseline_training[[ssb_column]]
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
  file.path(output_dir, "variable_selection_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  holdout_metrics,
  file.path(output_dir, "holdout_metrics.csv"),
  row.names = FALSE
)
utils::write.csv(
  holdout_predictions,
  file.path(output_dir, "holdout_predictions.csv"),
  row.names = FALSE
)
saveRDS(selected_models, file.path(output_dir, "selected_emsrr_models.rds"))

print(selection_summary)
print(holdout_metrics)