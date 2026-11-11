# ==================================================================== #
# Prepare cod substock models, stock data and environmental inputs
#
# Refit the selected models on common-footprint history, then prepare
# historical and projected data. History and projections remain separate.
# Both environmental series use a one-year lag.
# ==================================================================== #

rm(list = ls())

library(FLCore)
library(earth)
library(minqa)

# ------------------------- #
# paths and model selection
# ------------------------- #

# loop over the two different SR-types
types = c("BH","SEG","SEG_trunc_ICES")

for(sr.type in types){
#sr.type = "BH" # change to "SEG" for the segmented models

if(!sr.type %in% c("BH","SEG","SEG_trunc_ICES")){
  stop("sr.type must be BH, SEG or SEG_trunc_ICES.")
}

model.file = file.path(
  "./data/EMSRRs",paste0("Cod.EMSRR_",sr.type,".RData")
)
selection.file = ifelse(sr.type %in% c("BH","SEG"),file.path(
  "./results/model_01_fit_EMSRR",
  paste0(sr.type,"_selected_emsrr_models.rds")
),file.path(
  "./results/model_01_fit_EMSRR_SEG_ICES",
  paste0("SEG","_selected_emsrr_models.rds")
)) 
fit.data.file = "./data/Cod_Env_data.rda"
stock.file = "./data/stock/cod46a7d20_FLstockObjects_assessment_2026.Rdata"
env.file = "./data/Env/processed/bottomT_projections/bottomT_historical_projection_summaries.rds"

refit.seed = 2026L
refit.maxfun = 500L # as in model_01

# ---------------- #
# load the inputs
# ---------------- #

model.env = new.env(parent = emptyenv())
load(model.file,envir = model.env)

if(!exists("Cod.EMSRR",envir = model.env,inherits = FALSE)){
  stop("No Cod.EMSRR object in ",model.file)
}
models = model.env$Cod.EMSRR

# The selected-model RDS retains the chosen linear/MARS form.
selected.models = readRDS(selection.file)

fit.env = new.env(parent = emptyenv())
load(fit.data.file,envir = fit.env)

if(!exists("env_stk",envir = fit.env,inherits = FALSE)){
  stop("No env_stk object in ",fit.data.file)
}
fit.data = fit.env$env_stk

stock.env = new.env(parent = emptyenv())
stock.names = load(stock.file,envir = stock.env)

if(length(stock.names) != 1L){
  stop("Expected one stock-list object in ",stock.file)
}
cod = get(stock.names,envir = stock.env)

env = readRDS(env.file)
substocks = c("Southern","Northwestern","Viking")

if(!is.list(models) || !setequal(names(models),substocks)){
  stop("Expected one fitted model for each of the three substocks.")
}
if(!is.list(selected.models) ||
   !setequal(names(selected.models),substocks)){
  stop("Expected one saved selection for each of the three substocks.")
}
if(!is.data.frame(fit.data) ||
   !all(c("year","substock","rec","ssb_lag1","variable","value") %in%
        names(fit.data))){
  stop("env_stk is missing the model fitting columns.")
}
if(!is.list(env) || !is.data.frame(env$historical_common) ||
   !is.list(env$projections) ||
   !all(c("rcp45","rcp85") %in% names(env$projections))){
  stop("Environmental file must contain historical_common and both RCP scenarios.")
}

# -------------------------- #
# extract the stock histories
# -------------------------- #

# Follow data_02's naming and recruitment-age lag, but match SSB by
# calendar year rather than assuming rows are consecutive.
get.stock = function(stk){
  
  substk = sub(".*\\b(\\w+)\\s+sub-stock\\b.*","\\1",
               stk@name,perl = TRUE)
  substk = if(grepl("Northwest",substk)) "Northwestern" else
    if(grepl("South",substk)) "Southern" else substk
  
  r = rec(stk)
  s = ssb(stk)
  
  age = suppressWarnings(as.numeric(dimnames(r)$age))
  r.year = suppressWarnings(as.integer(dimnames(r)$year))
  s.year = suppressWarnings(as.integer(dimnames(s)$year))
  
  if(length(age) != 1L || is.na(age) || age < 0 ||
     age != floor(age) ||
     length(r.year) != length(as.numeric(r)) ||
     length(s.year) != length(as.numeric(s)) ||
     anyNA(r.year) || anyNA(s.year) ||
     anyDuplicated(r.year) || anyDuplicated(s.year)){
    stop(substk,": expected one recruitment age and one stock observation per year.")
  }
  
  out = data.frame(year = r.year,rec = as.numeric(r))
  out[[paste0("ssb_lag",age)]] = as.numeric(s)[match(r.year-age,s.year)]
  out = out[order(out$year),,drop = FALSE]
  
  list(substock = substk,data = out)
}

stock.list = lapply(cod,get.stock)
stock.ids = vapply(stock.list,`[[`,character(1),"substock")

if(anyDuplicated(stock.ids) || !setequal(stock.ids,substocks)){
  stop("Stock objects must contain exactly Southern, Northwestern and Viking.")
}

stock.list = setNames(lapply(stock.list,`[[`,"data"),stock.ids)

# ---------------- #
# helper functions
# ---------------- #

check.years = function(df,label){
  if(!is.data.frame(df) || !"year" %in% names(df) ||
     !is.numeric(df$year) || anyNA(df$year) ||
     any(!is.finite(df$year)) ||
     any(df$year != floor(df$year)) || anyDuplicated(df$year)){
    stop(label,": expected unique numeric calendar years.")
  }
}

lag.series = function(year,values,lag = 1L){
  values[match(year-lag,year)]
}

# Refit the same L1 log-residual objective and multistart BOBYQA as
# model_01. Reuse its saved full-data bounds, including custom limits.
refit.baseline = function(rec,ssb,bounds,n.starts,type,maxfun){
  
  if(!is.list(bounds) || !all(c("a","b") %in% names(bounds)) ||
     any(vapply(bounds[c("a","b")],function(b)
       length(b) != 2L || any(!is.finite(b)) ||
       any(b <= 0) || b[1] >= b[2],logical(1))) ||
     length(n.starts) != 1L || !is.finite(n.starts) ||
     n.starts < 2L || n.starts != floor(n.starts)){
    stop("Invalid saved baseline bounds or number of starts.")
  }
  
  lower = log(c(a = bounds$a[1],b = bounds$b[1]))
  upper = log(c(a = bounds$a[2],b = bounds$b[2]))
  
  if(type == "BH"){
    
    a0 = max(quantile(rec,0.75))
    b0 = a0/max(quantile(rec/ssb,0.90))
    
    objective = function(par){
      log.pred = par[1]+log(ssb)-log(ssb+exp(par[2]))
      sum(abs(log(rec)-log.pred))
    }
    
  } else {
    
    b0 = median(ssb)
    a0 = exp(median(log(rec)-log(pmin(ssb,b0))))
    
    objective = function(par){
      log.pred = par[1]+pmin(log(ssb),par[2])
      sum(abs(log(rec)-log.pred))
    }
  }
  
  starts = cbind(
    a = abs(rnorm(n.starts-1L,mean = a0,sd = 1.8*a0)),
    b = pmax(bounds$b[1],
             abs(rnorm(n.starts-1L,mean = b0,sd = 1.8*b0)))
  )
  starts = rbind(c(a = a0,b = b0),starts)
  
  starts[,"a"] = pmin(bounds$a[2],pmax(bounds$a[1],starts[,"a"]))
  starts[,"b"] = pmin(bounds$b[2],pmax(bounds$b[1],starts[,"b"]))
  starts = log(starts)
  
  fits = lapply(seq_len(nrow(starts)),function(i){
    tryCatch(
      minqa::bobyqa(starts[i,],objective,lower = lower,upper = upper,
                    control = list(maxfun = maxfun)),
      error = function(e) e
    )
  })
  
  valid = which(vapply(fits,function(x)
    !inherits(x,"error") && is.finite(x$fval) &&
      all(is.finite(x$par)),logical(1)))
  
  if(!length(valid)){
    failed = fits[vapply(fits,inherits,logical(1),what = "error")]
    details = unique(vapply(failed,conditionMessage,character(1)))
    stop("All baseline fits failed or were non-finite. ",
         paste(details,collapse = "; "))
  }
  
  converged = valid[vapply(fits[valid],function(x)
    isTRUE(x$ierr == 0),logical(1))]
  
  if(!length(converged)){
    warning("No baseline fit reported normal convergence; using the best finite fit.")
    converged = valid
  }
  
  best.start = converged[which.min(
    vapply(fits[converged],`[[`,numeric(1),"fval")
  )]
  best = fits[[best.start]]
  
  list(
    params = setNames(exp(best$par),c("a","b")),
    optimization = list(
      objective = best$fval,
      convergence = best$ierr,
      bounds = bounds,
      n_starts = n.starts,
      best_start = best.start
    )
  )
}

# Return scenario -> variable -> wide table, retaining realization IDs.
prepare.proj = function(runs,raw.var,model.var,rcp){
  
  if(!is.list(runs) || !length(runs) || is.null(names(runs)) ||
     any(!nzchar(names(runs))) || anyDuplicated(names(runs)) ||
     any(!grepl("^[0-9]+$",names(runs)))){
    stop(rcp,": realizations must have unique numeric names.")
  }
  
  runs = runs[order(as.numeric(names(runs)))]
  years = NULL
  out = NULL
  
  for(it in names(runs)){
    x = runs[[it]]
    check.years(x,paste(rcp,it))
    
    if(!raw.var %in% names(x) || !is.numeric(x[[raw.var]]) ||
       any(!is.finite(x[[raw.var]]))){
      stop(rcp," realization ",it,": invalid or missing ",raw.var)
    }
    
    x = x[order(x$year),,drop = FALSE]
    
    if(is.null(years)){
      years = x$year
      if(length(years) < 2L || any(diff(years) != 1L)){
        stop(rcp,": projection years must be consecutive.")
      }
      out = data.frame(year = years[-1L])
    } else if(!identical(as.numeric(x$year),as.numeric(years))){
      stop(rcp," realization ",it,": years differ from other realizations.")
    }
    
    # The first projected year has no preceding projection year.
    out[[paste0(model.var,"_it",it)]] =
      lag.series(x$year,x[[raw.var]])[-1L]
  }
  
  out
}

# ----------------------------------------------------------- #
# refit and prepare each substock model with its own inputs
# ----------------------------------------------------------- #

hist = env$historical_common
check.years(hist,"Common-footprint history")
hist = hist[order(hist$year),,drop = FALSE]

Cod.EMSRR.obj = setNames(vector("list",length(substocks)),substocks)
models.common = setNames(vector("list",length(substocks)),substocks)

set.seed(refit.seed)

for(substk in substocks){
  
  model = models[[substk]]
  selected = selected.models[[substk]]
  
  if(length(model$env.vars) != 1L ||
     length(model$SR.model$vars) != 2L ||
     !identical(model$env.vars,model$SR.model$vars[2])){
    stop(substk,": expected a model with one selected environmental variable.")
  }
  
  if(!"ssb_lag1" %in% names(stock.list[[substk]]) ||
     !identical(model$SR.model$vars[1],"ssb_lag1")){
    stop(substk,": model and stock data must use the same SSB age lag.")
  }
  
  var = model$env.vars
  sr.name = if(sr.type == "BH") "BevertonHolt.model" else "HockeyStick.model"
  selected.sr = if(sr.type == "BH") selected$bh_model else selected$seg_model
  
  if(!identical(selected$variable,var) || is.null(selected.sr) ||
     length(selected$linpreds) != 1L || is.na(selected$linpreds) ||
     !is.logical(selected$linpreds) ||
     is.null(model$SR.model[[sr.name]]) ||
     !is.function(model$predict.func)){
    stop(substk,": saved selection and ",sr.type," model are inconsistent.")
  }
  
  # Translate the fitted mean name back to the data_00.d output name.
  raw.var = if(var == paste0("bottomT_mean_",substk)){
    paste0("bottomT_",substk)
  } else {
    var
  }
  
  if(!raw.var %in% names(hist) || !is.numeric(hist[[raw.var]])){
    stop(substk,": historical variable not found: ",raw.var)
  }
  if(length(hist$year) < 2L || any(diff(hist$year) != 1L) ||
     any(!is.finite(hist[[raw.var]]))){
    stop(substk,": historical series must be complete and annual.")
  }
  
  # Historical environmental value for recruitment year t is from t-1.
  hist.var = data.frame(year = hist$year[-1L])
  hist.var[[var]] = lag.series(hist$year,hist[[raw.var]])[-1L]
  
  # Reconstruct the original full-data fit first. This checks that the
  # saved model, selected variable and fitting data belong together.
  d = fit.data[fit.data$substock == substk & fit.data$variable == var,
               c("year","rec","ssb_lag1","value"),drop = FALSE]
  
  d = d[complete.cases(d) & is.finite(d$year) & is.finite(d$rec) &
          is.finite(d$ssb_lag1) & is.finite(d$value) &
          d$rec > 0 & d$ssb_lag1 > 0,,drop = FALSE]
  d = d[order(d$year),,drop = FALSE]
  
  if(nrow(d) < 5L || anyDuplicated(d$year) ||
     !identical(as.numeric(d$year),as.numeric(model$SR.model$years))){
    stop(substk,": original model years do not match the saved fitting data.")
  }
  
  newdata = data.frame(d$ssb_lag1,d$value,check.names = FALSE)
  names(newdata) = model$SR.model$vars
  old.pred = model$predict.func(model$SR.model,newdata)
  
  if(!isTRUE(all.equal(as.numeric(d$rec),
                       as.numeric(model$SR.model$observed),
                       tolerance = 1e-8,check.attributes = FALSE)) ||
     !isTRUE(all.equal(as.numeric(old.pred),
                       as.numeric(model$SR.model$predicted),
                       tolerance = 1e-8,check.attributes = FALSE))){
    stop(substk,": saved original fit does not match the fitting data.")
  }
  
  # Replace only the environmental values; keep the original fit years.
  common.values = hist.var[[var]][match(d$year,hist.var$year)]
  
  if(any(!is.finite(common.values))){
    stop(substk,": common-footprint values are missing for the original fit years.")
  }
  
  cat(substk,"|",var,"| maximum footprint difference:",
      max(abs(d$value-common.values)),"\n")
  
  d$value = common.values
  
  # Refit the baseline and selected MARS form on all original fit years.
  opt = model$SR.model$optimization
  
  if (sr.type == "BH") {
    # Reconstruct the stock data used to set model_01's full-data bounds.
    stock <- unique(fit.data[fit.data$substock == substk,
                             c("year", "rec", "ssb_lag1"), drop = FALSE])
    
    stock <- stock[complete.cases(stock) &
                     is.finite(stock$year) &
                     is.finite(stock$rec) & stock$rec > 0 &
                     is.finite(stock$ssb_lag1) & stock$ssb_lag1 > 0,
                   , drop = FALSE]
    
    if (!nrow(stock)) stop(substk, ": no valid stock data.")
    
    Smax <- max(stock$ssb_lag1)
    expected <- list(
      a = c(0.01 * min(stock$rec), 100 * max(stock$rec)),
      b = c(Smax * (1 / 0.85 - 1), Smax * (1 / 0.75 - 1))
    )
    
    if (!is.list(opt$bounds) ||
        !isTRUE(all.equal(opt$bounds$a, expected$a,
                          tolerance = 1e-8, check.attributes = FALSE)) ||
        !isTRUE(all.equal(opt$bounds$b, expected$b,
                          tolerance = 1e-8, check.attributes = FALSE))) {
      stop(substk, ": saved BH bounds differ from the 75–85% constraints. ",
           "Rerun model_01 and check that model_02 loads the new fit.")
    }
    
    cat(substk, "| BH bounds verified; Smax =", Smax, "\n")
  }
  
  
  if(sr.type %in% c("BH","SEG")){
  refit = refit.baseline(
    d$rec,d$ssb_lag1,opt$bounds,opt$n_starts,
    sr.type,refit.maxfun
  )
  } else if(sr.type == "SEG_trunc_ICES"){
    refit = list(
      params = model$SR.model[[sr.name]]$params,
      optimization = opt
    )
  }
  
  model$SR.model[[sr.name]]$params = refit$params
  baseline = model$SR.model[[sr.name]]$model(refit$params,d$ssb_lag1)
  
  model$SR.model$MARS.model = earth::earth(
    log_residuals ~ x,
    data = data.frame(
      log_residuals = log(d$rec)-log(baseline),
      x = d$value
    ),
    linpreds = selected$linpreds,
    fast.k = 20,
    nprune = 20,
    pmethod = "exhaustive"
  )
  
  model$SR.model$years = d$year
  model$SR.model$observed = d$rec
  model$SR.model$optimization = refit$optimization
  
  newdata = data.frame(d$ssb_lag1,d$value,check.names = FALSE)
  names(newdata) = model$SR.model$vars
  
  model$SR.model$predicted = model$predict.func(model$SR.model,newdata)
  
  if(any(!is.finite(model$SR.model$predicted)) ||
     any(model$SR.model$predicted <= 0)){
    stop(substk,": refit produced invalid recruitment predictions.")
  }
  
  model$sd.residuals =
    sd(log(d$rec)-log(model$SR.model$predicted))
  
  model$footprint = "common"
  model$env.lag = 1L
  model$selection.footprint = "original historical"
  model$holdout.metrics.footprint = "original historical"
  models.common[[substk]] = model
  
  # Prepare both RCP scenarios without joining them to history.
  proj = lapply(c("rcp45","rcp85"),function(rcp){
    setNames(list(prepare.proj(env$projections[[rcp]],
                               raw.var,var,rcp)),var)
  })
  names(proj) = c("rcp45","rcp85")
  
  Cod.EMSRR.obj[[substk]] = list(
    stock.data = stock.list[[substk]],
    model = model,
    env.variables = list(
      historical = hist.var,
      climate.projections = proj
    )
  )
  
  # Confirm that the packaged inputs reproduce the full-data refit.
  check = merge(stock.list[[substk]],hist.var,by = "year")
  check = check[match(model$SR.model$years,check$year),,drop = FALSE]
  
  if(anyNA(check$year) || any(!is.finite(check$rec)) ||
     any(!is.finite(check$ssb_lag1)) ||
     any(!is.finite(check[[var]]))){
    stop(substk,": model refit years cannot be reconstructed from inputs.")
  }
  
  if(!isTRUE(all.equal(check$rec,
                       as.numeric(model$SR.model$observed),
                       tolerance = 1e-8,check.attributes = FALSE))){
    stop(substk,": recruitment differs from the refit input.")
  }
  
  predicted = model$predict.func(model$SR.model,check)
  
  if(!isTRUE(all.equal(as.numeric(predicted),
                       as.numeric(model$SR.model$predicted),
                       tolerance = 1e-8,check.attributes = FALSE))){
    stop(substk,": packaged input predictions differ from the saved model fit.")
  }
  
  cat(substk,"|",var,"|",nrow(check),"refit years verified\n")
}

# ---------------------- #
# save the prepared data
# ---------------------- #

dir.create("./model",showWarnings = FALSE,recursive = TRUE)

# Save the common-footprint models separately as well as in the
# combined object. The original model_01 files remain unchanged.
Cod.EMSRR = models.common

save(
  Cod.EMSRR,
  file = file.path("./data/EMSRRs",
                   paste0("Cod.EMSRR_",sr.type,"_common.RData"))
)

save(
  Cod.EMSRR.obj,
  file = file.path("./model",
                   paste0("Cod.EMSRR_",sr.type,"_model_and_data.RData"))
)
} # end of file