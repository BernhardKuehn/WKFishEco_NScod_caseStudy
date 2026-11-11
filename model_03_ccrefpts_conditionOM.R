# ------------------------------------------------------------- #
#################################################################
# >>>>>>>>>> IMPORTANT NOTE FOR THIS TO RUN PROPERLY: <<<<<<<<< 
# make sure 'FLasher' and the 'mse' package are NOT loaded!!!!
# as those will provide conflicts with the same type of classes
# defined in 'mseEMSRR' and 'FLasherEMSRR'
#################################################################
# ------------------------------------------------------------- #

# Creating operating models incorporating explicit environmental effects
# on stock recruitment:
# - Current, RCP4.5 & RCP8.5 as climate projections with environmental uncertainty
# - autocorrelated recruitment deviances as stochastic noise
# - bootstrapped SRR model versions to represent parameter uncertainty (not included for the moment)

# remove all existing objects ----------------
rm(list = ls())

library(FLCore)
library(FLasherEMSRR)
library(mseEMSRR)
library(stringr)

# 1. required functions ------------------------------------------------------
source("./functions/utilities.R")

# function to calculate autocorrelation of residuals 
calc.rho = function(emsrr.obj,boot.mods){
  
  input.dat = merge(emsrr.obj$stock.data,emsrr.obj$env.variables$historical)
  input.dat = input.dat[complete.cases(input.dat),]
  resids = lapply(boot.mods, function(x){
    pred = emsrr.obj$model$predict.func(x,newdata = input.dat)
    obs = input.dat$rec
    return(log(pred)-log(obs))}) # use log-space
  # get AR1-rho coefficient
  rho = sapply(resids,function(x) acf(x,plot=F)$acf[2])
  return(rho)
}


# 2. load look-up tables ------------------------------------------------------

load(file = "./data/eqsim_lut.Rdata")
load(file = "./data/BRPs.Rdata")

# FINAL year
fy <- 2099

# NUMBER of iterations (use only 100 for the moment)
it <- 100 #500

# -------------------------------------------------- #
# specify deviance type
deviance.type = "lognormal_ar1" # either 'lognormal' or 'lognormal_ar1'
# specify if bootstrapped models should be used
model.type = "single" # either 'single' or 'boot' (not for the moment)

# -------------------------------------------------- #
# define periods with a constant environment
periods = c("None",2030,2040,2050)
# define climate scenarios
climScen = c("noCC","rcp45","rcp85")

# combine
cc = expand.grid(climScen = climScen,
                 period = periods)
# exclude a few combinations
cc = cc[!(cc$climScen %in% c("rcp45","rcp85") & cc$period == "None" |
          cc$climScen == "noCC" & cc$period != "None"),]
print(cc)

# load stock file ----------------

# Load and identify the three assessment stocks
stock.file = "./data/stock/cod46a7d20_FLstockObjects_assessment_2026.Rdata"

load(stock.file,verbose = T)
class(stks)
names(stks) = sapply(1:length(stks),function(x) str_extract(stks[[x]]@name,"Southern|South|Northwestern|Northwest|Viking"))
names(stks) = dplyr::recode(names(stks),"South" = "Southern","Northwest" = "Northwestern")

# check for name matching
if(!all(names(stks) %in% rownames(BRPs))){
  stop("Stk names need to match the rownames in BRPs!")
}

# 3. loop through substocks ---------
srr.type = c("BH","SEG","SEG_trunc_ICES")

for(srt in srr.type){ # loop through srr types
  
  # load SRR obj with data
  emsrr.obj.nm = load(file.path("./model/",
                                paste0(paste("Cod.EMSRR",srt,"model_and_data",sep = "_"),".RData")))
  emsrr.obj_all = get(emsrr.obj.nm)
  
for(substock in rownames(BRPs)){

  # get substock from stock object
  cat("Building OM for:",substock,"...","\n")

  # load stock obj
  stock = stks[[substock]]
  
  DIMS <- dims(stock)
  RANGE <- range(stock)
  
  # INTERMEDIATE year
  iy <- RANGE["maxyear"] + 1
  
  # DATA year
  dy <- iy - 1
  
  # define ref-pts
  BRPs[substock,]  
  refpts <- FLPar(BRPs[substock,], params = c("Flim", "Fpa", "Fmsy", "lFmsy", "uFmsy", "Fmgt", "Blim", "Bpa", "Btrigger"), units = "NA")
  
  set.seed(987)
  
  # project stock
  stockProj = propagate(stock,it)
  stockProj = stf(stockProj,fy-dy)
  
  ## - Stock-recruitment relationship(s) (SRR)
  # ------------------- #
  # prepare EMSRR setup
  
  # select per substock
  emsrr.obj = emsrr.obj_all[[substock]]
  
  # get model & predict function
  mod = emsrr.obj$model$SR.model
  predict.mod = emsrr.obj$model$predict.func
  
  # get shorter stock-name abbreviation
  # Unique names also allow BH and SEG files to be loaded together
  stk.nm = paste(substock,srt,sep = "_")
  
  # define the predict function to match the name of the stock
  assign(paste(stk.nm,"rec.model",sep = "."),emsrr.obj$model$SR.model)
  assign(paste(stk.nm,"predict.func",sep = "."),emsrr.obj$model$predict.func)
  env.vars = emsrr.obj$model$env.vars
  
  # create a function to make predictions
  # (make it as generic as possible)
  
  if(model.type == "single"){
  
  # 1. create empty function
  sr_predict = function(ssb,params){}
  # 2. define (overwrite) the body of the function  
  e <- substitute(
    {
    stk <- stk.nm
    newdat <- data.frame(ssb = as.numeric(ssb),
                         t(params))
    # look in the global environment outside the function (dangerous, careful here!!!)
    predict.func <- get(paste(stk,"predict.func",sep = "."),envir = sys.parent(0))
    rec.model <- get(paste(stk,"rec.model",sep = "."),envir = sys.parent(0))
    
    colnames(newdat) <- rec.model$vars
    # predict recruitment
    pred <- predict.func(rec.model, newdata = newdat)
    pred <- FLQuant(as.vector(pred))
    return(pred)
  }, list(stk.nm = stk.nm)) # explicitly write the name of the stock in the function body
  
  # 4. put into body of the function
  body(sr_predict,envir = environment(sr_predict)) <- e #as.call(c(as.name("{"), e))
  
  } else if(model.type == "boot"){
    
    # add the names of the environmental covars to the model object
    EMSRRs_boot[[substock]]$srr_fit_boot = lapply(EMSRRs_boot[[substock]]$srr_fit_boot, function(x){
      vars = emsrr.obj$model$SR.model$vars
      x$vars = vars
      return(x)
    })
    
    # define a list of models as "rec.model" from the bootstrapped obj
    assign(paste(stk.nm,"rec.model",sep = "."),EMSRRs_boot[[substock]]$srr_fit_boot[1:it])
    
    
    # 1. create empty function
    sr_predict = function(ssb,params){}
    # 2. define (overwrite) the body of the function  
    e <- substitute(
      {
        stk <- stk.nm
        newdat <- data.frame(ssb = as.numeric(ssb),
                             t(params[-length(params)]))
        # last parameter specifies the iteration
        it = params[length(params)]
        # look in the global environment outside the function (dangerous, careful here!!!)
        predict.func <- get(paste(stk,"predict.func",sep = "."),envir = sys.parent(0))
        # specifies a list
        rec.model <- get(paste(stk,"rec.model",sep = "."),envir = sys.parent(0))
        # select the model of the specific iteration
        rec.model <- rec.model[[it]]
        colnames(newdat) <- rec.model$vars
        # predict recruitment
        pred <- predict.func(rec.model, newdata = newdat)
        pred <- FLQuant(as.vector(pred))
        return(pred)
      }, list(stk.nm = stk.nm)) # explicitly write the name of the stock in the function body
    
    # 4. put into body of the function
    body(sr_predict,envir = environment(sr_predict)) <- e #as.call(c(as.name("{"), e))
    
    
    
  } else {
    stop("model.type must be either 'single' or 'boot'!")
  }
  
  sr_predict
  
  # define a specific name for the stock
  assign(paste(stk.nm,"sr_predict",sep = "."),
         sr_predict)
  
  # prepare projection with EMSRR
  srMod = as.FLSR(stockProj)
  it = dim(stockProj)[6]
  
  # define model as a formula
  model(srMod) = as.formula(paste("rec",paste0(paste(stk.nm,"sr_predict",sep="."),"(ssb,params)"),sep = "~"))
  # fill params (covar in this case)
  first.yr = range(srMod)[["minyear"]]
  last.yr = range(srMod)[["maxyear"]]
  
  # ------------------------------------------------- #
  # create covars based on period and climScen
  
  # replace the EMSRRs by a "no change climate signal" (aka. RCP4.5 without trend)
  # remove the median signal from it.
  
    tmp = emsrr.obj
    tmp$env.variables$climate.projections$noCC <- tmp$env.variables$climate.projections$rcp45
    env_nms = names(tmp$env.variables$climate.projections$rcp45)
    
    tmp$env.variables$climate.projections$noCC = lapply(env_nms, function(env) {
      x = tmp$env.variables$climate.projections$rcp45[[env]]
      med = apply(x[, -1], 1, median)
      yrs = x$year
      lag = tmp$model$env.lag
      if(length(lag) != 1L || is.na(lag) || lag != 1L){
        stop("Expected the one-year environmental lag saved by model_02.")
      }
      # set the offset of the noCC run to the HISTORICAL mean over 2006 - 2020
      x_hist = tmp$env.variables$historical[,c("year",env)]
      # offset for the noCC run defined as the historical overlapping period
      offset = mean(x_hist[x_hist$year %in% ((2006:2020)+lag),env])
      df = data.frame(year = x[, 1], x[, -1] - med + offset) 
      # correction for the intercept
      return(df)
    })
    names(tmp$env.variables$climate.projections$noCC) = env_nms
    # calc.mean over respective periods (2030,2040,2050)
    lst_means = list()
    for(r in c("rcp45","rcp85")){
      lst_means[[r]] = list()
      for(p in an(periods[!(periods %in% "None")])){
        tmp_rcp = tmp$env.variables$climate.projections[[r]]
        lst_means[[r]][[ac(p)]] = lapply(tmp_rcp,function(x){
          iter.period_mean = apply(x[x$year >= p & x$year < p+10,-1],2,mean)
          tot.period_mean = mean(iter.period_mean)
          # remove trend
          med = apply(x[,-1],1,median)
          #shift to level of the period average
          df = data.frame(year = x[,1],x[,-1] - med + tot.period_mean)
          return(df)})
        # put into original environmental data list as "new" scenarios
        tmp$env.variables$climate.projections[[paste(r,p,sep = "_")]] <- lst_means[[r]][[ac(p)]]
      }
    }
    # assign
    assign(x = "emsrr.obj",value = tmp)
  
  # --------------------------------------------------------------- #
  # deviances (lognormal distributed iid)
  # !************************************ !
  # might need correction later, e.g. to allow autocorrelation... 
  
  set.seed(42)
  # get cv of the fits
  if(model.type == "boot"){
  stop("Bootstrapped EMSRRs not defined yet...")
  } else if(model.type == "single"){
    cv = rep(emsrr.obj$model$sd.residuals,it)
  } else {
    stop("Specify 'model.type' to be either 'single' or 'boot'!")
  }
  
  
  if(deviance.type == "lognormal"){
    for(i in 1:it){
      if(i == 1){
        srdevs = FLQuant(rlnorm(n = length(dy:fy),meanlog = 0,
                                sdlog = cv[i]), #emsrr.obj$model$sd.residuals),
                         dimnames = list(year = dy:fy),iter = it)
      } else {
        iter(srdevs,i) <- FLQuant(rlnorm(n = length(dy:fy),meanlog = 0,
                                         sdlog = cv[i]), #emsrr.obj$model$sd.residuals),
                                  dimnames = list(year = dy:fy))
      }
    }
    
  } else if(deviance.type == "lognormal_ar1"){
    
    # ---------------------------------------- #
    # create autocorrelated resids
    # calculate AR1 coefficient rho!
    # - by predicting with the bootstrapped model on the original data and calc. deviances
    # - fit AR1 model
    if(model.type == "boot"){
    rho = calc.rho(emsrr.obj,
                   boot.mods = EMSRRs_boot[[substock]]$srr_fit_boot)
    } else if(model.type == "single"){
      rho = rep(calc.rho(emsrr.obj,
                     boot.mods = list(mod)),it)
    } else {
      stop("Specify 'model.type' to be either 'single' or 'boot'!")
    }
    # generate autocorrelated errors
    srdevs <- rlnormar1(sdlog=cv[1:it], 
                        rho = rho[1:it], 
                        years = seq(dy, fy), 
                        bias.correct = F,
                        n = it)
    # plot(srdevs,srdevs2)
  } else {
    stop("Please specify 'deviance.type' either as 'lognormal' or 'lognormal_ar1'!")
  }
  
  dim(srdevs)
  
  # ----------------------------------------- #
  # create different oms based on the selected RCP & period
  
  for(ii in 1:nrow(cc)){
  
  rcp = ac(cc$climScen[ii])
  period = ac(cc$period[ii])
  rcp_period = ifelse(rcp == "noCC",rcp,paste(rcp,period,sep = "_"))
  
  cat("Creating om for: ",">>",ii,"|",rcp_period,"<<"," climate and period combination...","\n")
    
  # fill covar
  # mean of all realisations
  
  env.proj = lapply(emsrr.obj$env.variables$climate.projections[[rcp_period]],
                    function(x){
                      years = x$year
                      tx = t(x[,-which(names(x) %in% "year")])
                      df = split(tx,f = 1:nrow(tx))
                      df = lapply(df,function(d){names(d) <-years; return(d)})
                      # extent to range up to first yr
                      yrs_extent = first.yr:max(years)
                      df_ext = lapply(df,function(d){
                        tmp = vector("numeric",length = length(yrs_extent))
                        names(tmp) = yrs_extent
                        tmp[yrs_extent %in% ac(years)] <- d
                        return(tmp)})
                      # write out
                      return(df_ext)})
  
  
  covar.flq = FLQuants(lapply(env.proj,
                              function(x){
                                # create multi-iter object
                                itl = length(x)
                                tmp.flq = FLQuant(x[[1]],dimnames = list(year = names(x[[1]])),iter = itl)
                                for(i in 1:itl)
                                  iter(tmp.flq,i) <- FLQuant(x[[i]],dimnames = list(year = names(x[[i]])))
                                return(tmp.flq)
                              }))
  
  # convert to array
  env = emsrr.obj$model$env.vars
  
  covar.array = list()
  for(i in 1:length(env)){
    tmp = array(as.numeric(covar.flq[[env[i]]][,ac(first.yr:last.yr),]),
                dim = c(1,length(first.yr:last.yr),it),
                dimnames = list(params = env[i],
                                year = first.yr:last.yr,
                                iter = 1:it)) 
    covar.array[[i]] = tmp
  }
  covar.array = DescTools::Abind(covar.array,along = 1)
  dim(covar.array)
  names(dimnames(covar.array)) = c("params","year","iter")
  
  if(model.type == "boot"){
    # add additional iteration as parameter
    it.array = array(matrix(1:it,
                            nrow = dim(covar.array)[2],
                            ncol = dim(covar.array)[3],
                            byrow = T),dim = c(1,dim(covar.array)[-1]),
                     dimnames = list("iteration",dimnames(covar.array)[[2]],dimnames(covar.array)[[3]]))
  covar.array = DescTools::Abind(covar.array,it.array,along = 1)  
  names(dimnames(covar.array)) = list("params","year","iter")
  }
  
  params(srMod) = FLPar(covar.array)
  srMod
  
  # -------------------------------------------------------------------------- #
  
  # BUILD FLom
  om <- mseEMSRR::FLom(stock=propagate(stock,it),
                  sr = srMod,
                  refpts=refpts,
                  deviances = srdevs)
  #om@sr@params # looks good        

  # SETUP om future: average of last 3 years **
  om <- fwdWindow(om, end=fy)
  
  # RESAMPLE forecast year slots with Eqsim definitions
  # catchability / selectivity
  uiter <- seq(it)
  if(!eqsim_lut[[substock]]$sel_const){ # sample years if not const mean value used
    for(i in uiter){
      # prepare the selectivity & landings fraction based on sel.yrs
      samp.yrs <- sample(x = eqsim_lut[[substock]]$sel_years, size = length(iy:fy), replace = T)
      om@stock@harvest[, ac(iy:fy),,,,i] <- om@stock@harvest[,ac(samp.yrs),,,,i]
      
      om@stock@landings.n[, ac(iy:fy),,,,i] <- om@stock@landings.n[,ac(samp.yrs),,,,i] / om@stock@catch.n[,ac(samp.yrs),,,,i]
      om@stock@discards.n[, ac(iy:fy),,,,i] <- 1 - om@stock@landings.n[, ac(iy:fy),,,,i]
    }
  }else{ # const mean value used
    for(i in uiter){
      mean.yrs <- eqsim_lut[[substock]]$sel_years
      om@stock@harvest[, ac(iy:fy),,,,i] <- apply(om@stock@harvest[,ac(mean.yrs),,,,i], MARGIN = c(1,3,4,5,6), mean, na.rm = TRUE)
      
      om@stock@landings.n[, ac(iy:fy),,,,i] <- apply(om@stock@landings.n[,ac(mean.yrs),,,,i] / om@stock@catch.n[,ac(mean.yrs),,,,i],MARGIN = c(1,3,4,5,6), mean, na.rm = TRUE)
      om@stock@discards.n[, ac(iy:fy),,,,i] <- 1 - om@stock@landings.n[, ac(iy:fy),,,,i]
      
    }
  }
  
  # ind weight (stock, landings, discards)
  # and biological slots
  if(!eqsim_lut[[substock]]$bio_const){ # sample years if not const mean value used
    for(i in uiter){
      samp.yrs <- sample(x = eqsim_lut[[substock]]$bio_years, size = length(iy:fy), replace = T)
      
      om@stock@mat[, ac(iy:fy),,,,i] <- om@stock@mat[,ac(samp.yrs),,,,i]
      om@stock@m[, ac(iy:fy),,,,i] <- om@stock@m[,ac(samp.yrs),,,,i]
      om@stock@stock.wt[, ac(iy:fy),,,,i] <- om@stock@stock.wt[,ac(samp.yrs),,,,i]
      
      om@stock@landings.wt[, ac(iy:fy),,,,i] <- om@stock@landings.wt[,ac(samp.yrs),,,,i]
      om@stock@discards.wt[, ac(iy:fy),,,,i] <- om@stock@discards.wt[,ac(samp.yrs),,,,i]
      om@stock@catch.wt[, ac(iy:fy),,,,i] <- om@stock@catch.wt[,ac(samp.yrs),,,,i]
      
    }
  }else{ # const mean value used
    for(i in uiter){
      mean.yrs <- eqsim_lut[[substock]]$bio_years
      
      om@stock@mat[, ac(iy:fy),,,,i] <- apply(om@stock@mat[,ac(mean.yrs),,,,i], MARGIN = c(1,3,4,5,6), mean, na.rm = TRUE)
      om@stock@m[, ac(iy:fy),,,,i] <- apply(om@stock@m[,ac(mean.yrs),,,,i], MARGIN = c(1,3,4,5,6), mean, na.rm = TRUE)
      
      om@stock@stock.wt[, ac(iy:fy),,,,i] <- apply(om@stock@stock.wt[,ac(mean.yrs),,,,i], MARGIN = c(1,3,4,5,6), mean, na.rm = TRUE)
      om@stock@landings.wt[, ac(iy:fy),,,,i] <- apply(om@stock@landings.wt[,ac(mean.yrs),,,,i], MARGIN = c(1,3,4,5,6), mean, na.rm = TRUE)
      om@stock@discards.wt[, ac(iy:fy),,,,i] <- apply(om@stock@discards.wt[,ac(mean.yrs),,,,i], MARGIN = c(1,3,4,5,6), mean, na.rm = TRUE)
      om@stock@catch.wt[, ac(iy:fy),,,,i] <- apply(om@stock@catch.wt[,ac(mean.yrs),,,,i], MARGIN = c(1,3,4,5,6), mean, na.rm = TRUE)
      
    }
  }
  
  
  # SET stochastic rec dy
  rec(stock(om))[, ac(dy)] <- rec(om)[1, ac(dy)] * srdevs[, ac(dy)]
  
  # PROJECT forward for iy assumption
  ctrl <- FLasherEMSRR::fwdControl( 
    data.frame(
      year = c(iy),
      value = c(1),
      quant = c("f"),
      relYear = c(dy)                               
    )
  )
  # project 1yr with forward (calls FLasherEMSRR)
  om <- FLasherEMSRR::fwd(om, control = ctrl)
  
  # build specific oms
  assign(paste("om",rcp,period,sep = "_"),om)
  
  }
  
  # F and SSB deviances
  shortcut_devs2 <- function(om, Fcv = 0.212, Fphi = 0.423, SSBcv = 0, SSBphi = 0){
    devs <- FLQuants(
      F = ar1rlnorm(rho = Fphi, years = dimnames(om)$year, iters = dims(om)$iter, 
                    meanlog = 0, sdlog = Fcv, bias.correct = F), 
      SSB = ar1rlnorm(rho = 0, years = dimnames(om)$year, iters = dims(om)$iter, 
                      meanlog = 0, sdlog = SSBcv, bias.correct = F))
    return(devs)
  }
  sdevs <- shortcut_devs2(om, Fcv=eqsim_lut[[substock]]$Fcv, 
                          Fphi=eqsim_lut[[substock]]$Fphi, SSBcv=eqsim_lut[[substock]]$SSBcv)
  
  
  # - CONSTRUCT iem, implementation error module w/10% noise **
  # should be based on Fcv?? ****
  iem <- FLiem(method=noise.iem,
               args = list(noise=rlnorm(it, rec(om)[,,,,,1] %=% 0, 0.1)))
  
  
  # - SAVE
  
  save(list = c(
        # export different oms, based on the specified RCP & period
       grep(paste(paste0("om_",paste(cc$climScen,cc$period,sep = "_")),collapse = "|"),ls(),value = T),
       "iem", "sdevs", "iy",
       # export functions related to the EMSRR
       grep(paste(substock,srt,sep = "_"),grep("sr_predict|rec.model|predict.func",ls(),value = T),value = T)
       ), 
       file = file.path("model", paste0("Cod_",substock,"_EMSRR_",srt, "_ccrefpts_input.rda")), compress="xz")
  
  } # end of substock loop
} # end of sr-type loop

