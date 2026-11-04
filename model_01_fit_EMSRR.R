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

rm(list = ls())

# load packages
library(earth)

# -------------------------- #
# load the fitting functions
# -------------------------- #

source("./functions/MO.fitness.func.BevertonHolt_MARS.R")
source("./functions/MO.fitness.func.envHockeyStick_MARS.R")
source("./functions/calc.partial.effects.R")

# helper function to score predictions
score_predictions <- function(observed, predicted) {
  data.frame(
    n = length(observed),
    rmse = sqrt(mean((observed - predicted)^2)),
    mae = mean(abs(observed - predicted)),
    log_rmse = sqrt(mean((log(observed) - log(predicted))^2)),
    log_bias = mean(log(predicted) - log(observed))
  )
}

input_file <- "./data/Cod_Env_data.rda"
output_dir <- "./results/model_01_fit_EMSRR"
ssb_column <- "ssb_lag1"
expected_substocks <- c("Southern", "Northwestern", "Viking")
holdout_n_years <- 10L

if (!file.exists(input_file)) {
  stop("Input data not found: ", input_file,
       ". Run data_02_join_SST_with_substock.R first.")
}

loaded_data <- new.env(parent = emptyenv())
loaded_objects <- load(input_file, envir = loaded_data)
if (!"env_stk" %in% loaded_objects) {
  stop("Expected object 'env_stk' was not found in ", input_file, ".")
}
if (!is.data.frame(loaded_data$env_stk)) {
  stop("Object 'env_stk' must be a data frame.")
}

dat <- loaded_data$env_stk 

if (anyNA(dat$year) || anyNA(dat$substock) || any(!nzchar(dat$substock))) {
  stop("env_stk contains missing/invalid years or substock labels.")
}

# use only complete cases
dat = dat[complete.cases(dat),]

# Define the holdout from calendar years before removing incomplete rows, so a
# missing SST/recruitment value cannot silently move the evaluation period.
available_years <- sort(unique(dat$year))
if (length(available_years) <= holdout_n_years) {
  stop("Need more than ", holdout_n_years, " distinct years to fit and evaluate.")
}
holdout_years <- tail(available_years, holdout_n_years)
first_holdout_year <- min(holdout_years)


model_objects <- setNames(vector("list", length(expected_substocks)),
                          expected_substocks)
prediction_rows <- list()
metric_rows <- list()
parameter_rows <- list()

# ----------------------------------------------------------------------- #
# fit a model per substock
# use bootstrap out-of-sample predictions to chose the best model for each
# substock
# ------------------------------------------------------------------------ #

for (substock in expected_substocks) {
  stock_data <- dat[dat$substock == substock, , drop = FALSE]
  training_data <- stock_data[stock_data$year < first_holdout_year, , drop = FALSE]
  holdout_data <- stock_data[stock_data$year %in% holdout_years, , drop = FALSE]

  # fit a model for each SST time series seperately
  covars_unique = unique(training_data$variable)
  length(covars_unique)
  n_boot = 50
  k = 1; metrics_test_list = list()
  for(var in covars_unique){
    
   tmp = training_data[training_data$variable == var,c("rec","ssb_lag1","value")]
  names(tmp) = c("rec","ssb_lag1",var)
  # generate bootstrap samples
  for(n in 1:n_boot){
    cat(var,"|",n,"\n")
    indx_boot = sample(1:nrow(tmp),size = nrow(tmp),replace = T)
    train_in = tmp[indx_boot,]
    test_oos = tmp[-unique(indx_boot),]  
    
    # fit model
    bh_fit <- fit.BevertonHolt.MARS(rec = train_in$rec,
                                    ssb = train_in$ssb_lag1,
                                    env.Vars = train_in[,var,drop = F],linpreds = T,
                                    Nrestarts = 1000,thresh = 0.1,pmethod = "exhaustive")
    # calculate out-of-sample metrics
    test_pred = data.frame(obs = test_oos$rec,
                           pred = predict.BevertonHolt.MARS(bh_fit,newdata = test_oos))
    # derive metrics
    metrics_test = score_predictions(observed = test_pred$obs,predicted = test_pred$pred)
    metrics_test$variable = var
    metrics_test$boot_sample = n
    # store
    metrics_test_list[[k]] = metrics_test
    k = k+1
  }
  }
  metrics_test_var = rowbind(metrics_test_list)
  metrics_test_var$substock = substock

  ggplot(metrics_test_var,aes(x = variable,y = log_rmse)) + 
    geom_boxplot()
  ggplot(metrics_test_var,aes(x = variable,y = log_bias)) + 
    geom_boxplot()
  
  
  
  
}


# for testing later on
tmp = training_data[training_data$variable == "SST_Southern_EOT8","value",drop = F]
names(tmp) = "covar"
bh_fit <- fit.BevertonHolt.MARS(rec = training_data$rec[training_data$variable == "SST_Southern_EOT8"],
                                ssb = training_data$ssb_lag1[training_data$variable == "SST_Southern_EOT8"],
                                env.Vars = tmp,linpreds = T,
                                Nrestarts = 1000,thresh = 0.1,pmethod = "exhaustive")
bh_fit$vars = c("ssb_lag1","covar")
tmp = dat[dat$variable == "SST_Southern_EOT8",c("ssb_lag1","value")]
names(tmp)[2] = "covar"
calc.partial.effects(bh_fit,y = dat$rec[dat$variable == "SST_Southern_EOT8"],
                     x = tmp,predict.func = predict.BevertonHolt.MARS)



