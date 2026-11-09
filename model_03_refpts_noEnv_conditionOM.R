# remove all existing objects ----------------
rm(list = ls())

library(FLCore)
library(FLasher)
library(mse)
library(FLSRTMB)
library(collapse)
library(stringr)

# 1. required functions ------------------------------------------------------
source("functions/utilities.R")

# segmented regression with fixed breakpoint at Blim (Bloss for type 6?)
segreg3 <- function(ab, ssb) log(ifelse(ssb >= BlimBloss, ab$a * BlimBloss, ab$a * ssb))
segreg <- function(params,ssb)  ifelse(ssb <= params$b, params$a * ssb, params$a * params$b)
bh <- function(params,ssb) params$a * ssb/(params$b + ssb)

# function to calculate a bootstrapped rho-estime for recruitment deviations
rho_boot = function(srr.pars,ssb,rec,srr.model.func){
  require(pautocorr)
  # vectorized moving block bootstrap from the pautocorr package, but returning the raw bootstrapped correlations
  mbb = function (acf.est, pacf.est, ahat, data, l, a1, a2, boot, lagmax) 
  {
    acf.y <- matrix(NA, nrow = boot, ncol = lagmax)
    #pacf.y <- matrix(NA, nrow = boot, ncol = lagmax)
    n.t <- length(data)
    k <- floor(n.t/l)
    num.b <- n.t - l + 1
    seq.b <- seq(1, num.b, by = 1)
    odr <- seq.b %x% t(rep(1, l))
    add <- rep(1, num.b) %x% t(0:(l - 1))
    odr <- odr + add
    for (i in 1:boot) {
      for (j in 1:50) {
        temp <- as.vector(sample(seq.b, k + 1, replace = TRUE))
        odri <- t(odr[temp, ])
        pair.id.M <- pair.VMB(odri, lagmax, n.t)
        tmp.acf <- corr.VMB(pair.id.M, data)
        tmp.pacf <- Durbin_Levinson_rho(tmp.acf)
        if (sum(abs(tmp.pacf > 1) + abs(tmp.pacf < (-1)), 
                na.rm = TRUE) == 0) {
          break
        }
      }
      acf.y[i, ] <- tmp.acf
      #pacf.y[i, ] <- tmp.pacf
    }
    return(acf.y)
  }
  
  # srr.pars ->  data.frame of the form
  # iter, par1, par2, sigmaR
  it = max(srr.pars$iter)
  
rho = list()
for(ii in 1:it){
  pars = srr.pars[srr.pars$iter == ii,!(names(srr.pars) %in% c("iter"))]
  mod = srr.model.func[[pars$mod]]

  rec_pred = mod(ssb = ssb,
                 params = pars)
  rec_ts = rec[-1]
  rec_ts_pred = rec_pred[-length(rec_pred)]
  # calc residuals
  resids = log(rec_pred[-length(rec_pred)]) - log(as.numeric(rec(stkSRyrs))[-1])
  # moving block bootstrap
  x = as.matrix(resids)
  nT = nrow(x)
  blocksize = round(nT^(1/3))
  lagmax = 1
  alpha = 0.05
  ahat <- JN.VMB(x, lagmax+1, blocksize)
  acf.est = acf(x,plot = F,lag.max = lagmax)$acf[-1]
  # block bootstrap
  rho[[ii]] = mbb(acf.est = acf.est,pacf.est = NULL,ahat = ahat,data = x,l = blocksize,a1 = alpha/2,
                  a2 = 1 - alpha/2,boot = 1E5/it, # adjust the number of inner bootstrap samples for rho depending on the number of outer iterations
                  lagmax = lagmax)[,1]
}
rho = do.call(c,rho)

# store a sample from the bootstrapped density (from https://stackoverflow.com/questions/17087312/generate-a-random-number-from-a-density-object-or-more-broadly-from-a-set-of-nu)
dens <- density(na.omit(rho))
rho_sample <- sample(x = dens$x, it, prob = dens$y, replace=TRUE) + rnorm(it, 0, dens$bw)
#srpars.boot$rho = rho_sample
# return rho_sample
return(rho_sample)
}

# 2. load look-up tables ------------------------------------------------------

load(file = "data/eqsim_lut.Rdata")
load(file = "data/BRPs.Rdata")

# FINAL year
fy <- 2099

# NUMBER of iterations (for the moment only 100 iterations due to time constraints)
it <- 100 #500

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


# 3. loop through stocks ---------

for(substock in rownames(BRPs)){

  # load stock object
  cat("Building OM for:",substock,"...","\n")
  
  # load stock obj
  stock = stks[[substock]]
  DIMS <- dims(stock)
  RANGE <- range(stock)
  
  # INTERMEDIATE year
  iy <- RANGE["maxyear"] + 1
  
  # DATA year
  dy <- iy - 1

  
  BRPs[substock,]  
  refpts <- FLPar(BRPs[substock,], params = c("Flim", "Fpa", "Fmsy", "lFmsy", "uFmsy", "Fmgt", "Blim", "Bpa", "Btrigger"), units = "NA")
  
  set.seed(987)
  
  # - Stock-recruitment relationship(s) **
  recrLag <- RANGE["min"] # rec_lag between ssb and rec
  srr_years_lag <- (min(eqsim_lut[[substock]]$srr_years) - recrLag) : max(eqsim_lut[[substock]]$srr_years)
  stkSRyrs <- stock[,unique(ac(pmax(min(an(dimnames(stock)$year)),srr_years_lag)))]
  SRmods <- eqsim_lut[[substock]]$models
  SRmods <- tolower(SRmods)
  SRmodsAdj <- replace(SRmods, SRmods == "segregfix", "segreg")
  
  # FIT models
  if("segregfix" %in% SRmods){
    brkpt = ifelse(substock %in% c("Northwestern","Viking"),as.numeric(ssb(stkSRyrs)[,ac(2006),]),BRPs[substock, "Blim"])
    fits <- srrTMB(as.FLSRs(stkSRyrs, models=SRmodsAdj), 
      spr0 = mean(spr0y(stkSRyrs)), Blim = brkpt)
  }else{
    fits <- srrTMB(as.FLSRs(stkSRyrs, models=SRmodsAdj), 
      spr0=mean(spr0y(stkSRyrs)))
  }
  
  
  # PLOT
  plotsrs(fits)
  # plotsrts(fits, relative = TRUE)

  # BOOTSTRAP and SELECT model by largest logLik **
  if("segregfix" %in% SRmods){
    brkpt = ifelse(substock %in% c("Northwestern","Viking"),as.numeric(ssb(stkSRyrs)[,ac(2006),]),BRPs[substock, "Blim"])
    srpars <- bootstrapSR(stkSRyrs, iters=it,
      models=SRmodsAdj, method="best", Blim = brkpt,spr0 = mean(spr0y(stkSRyrs)))
  }else{
    srpars <- bootstrapSR(stkSRyrs, iters=it,
      models=SRmodsAdj, method="best",spr0 = mean(spr0y(stkSRyrs)))
  }
  
  dir.create("./data/SRRs",showWarnings = F)
  save(fits, srpars, file=file.path("data","SRRs", paste0("Cod_",substock, "_SRR_bootstrap.rda")), compress="xz")
  
  
  # - CONSTRUCT OM
  
  # GENERATE future deviances: lognormal autocorrelated **
  if(eqsim_lut[[substock]]$rhologRec){
    # since the bootstrap is not working correctly, do it myself
    #rhologRec <- srpars$rho
    
    srpars_df = as.data.frame(srpars)
    srpars_df = pivot(srpars_df,ids = "iter",names = "params",how = "wider")
    srpars_df$iter = an(ac(srpars_df$iter))
    srpars_df = srpars_df[,c("iter","a","b","m")]
    srpars_df$mod <- ifelse(srpars_df$m == 1,"bh",
                            ifelse(srpars_df$m == 2,"ricker",
                              ifelse(srpars_df$m == 3,"segreg")))
    # print
    cat("# -------------------------------------------------------------------- #","\n")
    cat("Running autocorrelation routine with",">>>",unique(srpars_df$mod),"<<<","selected as SRR.","\n")
    cat("# --------------------------------------------------------------------- #","\n")
    
    srpars_df$m <- NULL
    ssb_rec_df = data.frame(rec = dplyr::lag(as.numeric(rec(stkSRyrs)),RANGE[["min"]]),
                            ssb = as.numeric(ssb(stkSRyrs)))
    
    srpars$rho = rhologRec = rho_boot(srr.pars = srpars_df,ssb = ssb_rec_df$ssb,rec = ssb_rec_df$rec,
                                      srr.model.func = list(segreg = segreg,bh = bh))
    # plot the rho distribution over the bootstrap samples
    hist(rhologRec,main = paste("Autocorr.",substock))
    
  }else{
    rhologRec <- srpars$rho
    rhologRec[] <- 0
  }
  
  srdevs <- rlnormar1(sdlog=srpars$sigmaR, 
    rho = rhologRec, 
    years = seq(dy, fy), 
    bias.correct = F)
  
  #plot(srdevs)

  
  # BUILD FLom
  om <- FLom(stock=propagate(stock, it), refpts=refpts, model='mixedsrr',
    params=srpars, deviances=srdevs)
  
  # SETUP om future: average of last 3 years **
  om <- fwdWindow(om, end=fy)
  
  # RESAMPLE forecast year slots with Eqsim definitions
  # catchability / selectivity
  uiter <- seq(it)
  if(!eqsim_lut[[substock]]$sel_const){ # sample years if not const mean value used
    for(i in uiter){
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
  ctrl <- fwdControl( 
    data.frame(
      year = c(iy),
      value = c(1),
      quant = c("f"),
      relYear = c(dy)                               
    )
  )
  om <- fwd(om, control = ctrl)

  
  
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
  rec_zeros_dimNms = dimnames(rec(om))
  rec_zeros_dimNms$iter = "1"
  rec_zeros = FLQuant(0,dim = c(dim(rec(om))[-6],1),
                      dimnames = rec_zeros_dimNms)
  
  iem <- FLiem(method=noise.iem,
    args = list(noise=rlnorm(it,rec_zeros, 0.1)))
  
  
  # - SAVE
  
  save(om, iem, sdevs, iy, file = file.path("model", paste0("Cod_",substock, "_refpts_input.rda")), compress="xz")

}


