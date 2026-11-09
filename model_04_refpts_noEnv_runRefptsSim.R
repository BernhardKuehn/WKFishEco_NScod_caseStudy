# remove all existing objects ----------------
rm(list = ls())

library(FLCore)
library(FLasher)
library(mse)
library(parallel)

# 0. create output directory --------------------------

dir.create(file.path("model", "refpts"), recursive = TRUE, showWarnings = FALSE)

# 1. required functions ------------------------------------------------------
source("./functions/utilities.R")

# segmented regression with fixed breakpoint at Blim (Bloss for type 6?)
segreg3 <- function(ab, ssb) log(ifelse(ssb >= BlimBloss, ab$a * BlimBloss, ab$a * ssb))

# parallelized functions to run MSE with/without HCR

# MSE with HRC
parFun_noHCR <- function(x){
  library(mse)
  library(FLasher)
  
  Ftarg <- info$Ftarg[x]
  
  # super fast, but I think EQsim still considers SSBcv, Fcv, and Fphi when applying no HCR
  # runFtarg <- fwd(om, control = fwdControl(year = seq(iy, fy), quant = "fbar", value = Ftarg))
  
  # thus, run MSE, but adapt HCR to Btrigger = 0
  mseargs <- list(iy=iy, fy=fy, data_lag=1, management_lag=1, frq=1)
  
  # SETUP standard ICES advice rule
  arule <- mpCtrl(list(
    
    # (est)imation method: shortcut.sa + SSB deviances
    est = mseCtrl(method=shortcut.sa,
                  args=list(SSBdevs=sdevs$SSB)),
    
    # hcr: hockeystick (fbar ~ ssb | lim, trigger, target, min)
    hcr = mseCtrl(method=hockeystick.hcr,
                  args=list(lim=0, trigger=0, target=Ftarg,
                            min=0, metric="ssb", output="fbar")),
    
    # (i)mplementation (sys)tem: tac.is (C ~ F) + F deviances
    isys = mseCtrl(method=tac.is,
                   args=list(recyrs= 10, fmin=0, Fdevs=sdevs$F)) # 10-yr geom mean
  ))
  
  # plot HCR
  # plot_hockeystick.hcr(arule$hcr, labels=c(lim="Blim", trigger="MSYBtrigger",
  #                                          min="", target="Ftarget")) +
  #   xlab(expression(hat(SSB))) + ylab(expression(bar(F)))
  
  # - RUN applying ICES advice rule
  system.time(
    advice <- mp(om, iem=iem, ctrl=arule, args=mseargs,parallel = FALSE)
  )
  
  # extract info
  res <- list(
    rec = rec(advice@om@stock),
    fbar = fbar(advice@om@stock),
    ssb = ssb(advice@om@stock),
    landings = landings(advice@om@stock),
    catch = catch(advice@om@stock)
  )
  
  return(res)
  
} # end of parFun_noHCR

# ------------------------------------------- #

# MSE with HRC
parFun_HCR <- function(x,Btrigger){
  library(mse)
  library(FLasher)
  
  Ftarg <- info$Ftarg[x]
  
  mseargs <- list(iy=iy, fy=fy, data_lag=1, management_lag=1, frq=1)
  
  # SETUP standard ICES advice rule
  arule <- mpCtrl(list(
    
    # (est)imation method: shortcut.sa + SSB deviances
    est = mseCtrl(method=shortcut.sa,
                  args=list(SSBdevs=sdevs$SSB)),
    
    # hcr: hockeystick (fbar ~ ssb | lim, trigger, target, min)
    # use lim? ***
    hcr = mseCtrl(method=hockeystick.hcr,
                  args=list(lim=0, trigger=Btrigger, target=Ftarg,
                            min=0, metric="ssb", output="fbar")),
    
    # (i)mplementation (sys)tem: tac.is (C ~ F) + F deviances
    isys = mseCtrl(method=tac.is,
                   args=list(recyrs= 10, fmin=0, Fdevs=sdevs$F)) # 10-yr geom mean
  ))
  
  # plot HCR
  # plot_hockeystick.hcr(arule$hcr, labels=c(lim="Blim", trigger="MSYBtrigger",
  #                                          min="", target="Ftarget")) +
  #   xlab(expression(hat(SSB))) + ylab(expression(bar(F)))
  
  # - RUN applying ICES advice rule
  system.time(
    advice <- mp(om, iem=iem, ctrl=arule, args=mseargs,parallel = FALSE)
  )
  
  # extract info
  res <- list(
    rec = rec(advice@om@stock),
    fbar = fbar(advice@om@stock),
    ssb = ssb(advice@om@stock),
    landings = landings(advice@om@stock),
    catch = catch(advice@om@stock)
  )
  
  return(res)
  
} # end of parFun with HCR

# 2. load look-up tables ------------------------------------------------------

load(file = "data/eqsim_lut.Rdata")
load(file = "data/BRPs.Rdata")

BRPs <- BRPs[rownames(BRPs) %in% c("COD-NS", "HER", "POK","PLE-NS","WHG-NS"),]

# NUMBER of evaluation years
ey <- 20

# NUMBER of cores
ncores <- 20 

newBRPs_bench <- BRPs*NA

for(stk_name in rownames(BRPs)){
  # load stock object
  
  load(file.path("model", paste0("Cod_",substock, "_refpts_input.rda")), verbose = T)
  
  DIMS <- dims(om@stock)
  RANGE <- range(om@stock)
  
  # DATA year
  dy <- iy - 1
  
  # FINAL year
  fy <- RANGE["maxyear"]
  
  # NUMBER of iterations
  it <- DIMS$iter
  

  # no hcr ------------
  
  # conduct forecasts across Ftarget levels without HCR and later with HCR
  info <- data.frame(Ftarg = seq(0, 1, by = 0.02),
                     stk_name = stk_name)
  
  info$scen <- seq(nrow(info))
  head(info)
  dim(info)
  # -------- #
  # simulate
  # -------- #
  
  # no hcr ------------
  
  # Set-up cluster and run
  clusterType <- ifelse(Sys.info()["sysname"] == "Windows", "PSOCK", "FORK")
  # ARGS <- ls() # to pass anything from the main environment (called below by clusterExport)
  # only pass relevant ARGs to avoid memory cluttering
  ARGS =  c("om","iy","fy","info","sdevs","iem")
  cl <- parallel::makeCluster(ncores, type=clusterType) # set-up cluster
  nn <- split(info$scen, info$scen) # indices for each iteration (used as "x" in the parFun)
  parallel::clusterExport(cl, varlist = ARGS, envir = environment()) # export info to the cluster
  t1 = Sys.time()
  RES <- parLapply(cl, nn, parFun_noHCR) # runs the function through each index in nn (index passed to "x" in parFun)
  stopCluster(cl) # stops cluster
  t2 = Sys.time()
  time.taken1 = t2-t1
  
  
  df <- lapply(as.list(seq(info$Ftarg)), FUN = function(x){
    tmp1 <- as.data.frame(RES[[x]]$ssb)
    tmp1$slot <- "ssb"
    tmp2 <- as.data.frame(RES[[x]]$landings)
    tmp2$slot <- "landings"
    tmp3 <- as.data.frame(RES[[x]]$catch)
    tmp3$slot <- "catch"
    tmp4 <- as.data.frame(RES[[x]]$rec)
    tmp4$slot <- "rec"
    tmp5 <- as.data.frame(RES[[x]]$fbar)
    tmp5$slot <- "fbar"
    tmp <- rbind(tmp1, tmp2, tmp3, tmp4, tmp5)
    tmp$Ftarg <- info$Ftarg[x]
    return(tmp)
  })
  df <- do.call("rbind", df)
  
  # Fmsy
  aggL <- aggregate(data ~ iter + Ftarg, data = df, subset = slot == "landings" & year %in% ((fy-ey):fy), FUN = median)
  #plot(data ~ Ftarg, data = aggL, pch = ".", cex = 5, col = adjustcolor(1,0.2))
  splL <- smooth.spline(x = aggL$Ftarg, y = aggL$data)
  newdat <- data.frame(Ftarg = seq(0,max(info$Ftarg), len = 1001))
  newdat$landings <- predict(splL, x = newdat$Ftarg)$y
  #lines(landings ~ Ftarg, newdat)
  (MSY <- newdat$landings[which.max(newdat$landings)])
  (Fmsy <- newdat$Ftarg[which.max(newdat$landings)])
  FmsyCI <- range(newdat$Ftarg[which(newdat$landings > MSY*0.95)]) # range of Ftarget +/- Fmsy
  hitL <- which(newdat$Ftarg==FmsyCI[1])
  hitU <- which(newdat$Ftarg==FmsyCI[2])
  (lFmsy <- newdat$Ftarg[hitL])
  (uFmsy <- newdat$Ftarg[hitU])
  
  
  # Bmsy, Btrigger
  aggBmed <- aggregate(data ~ iter + Ftarg, data = df, subset = slot == "ssb" & year %in% ((fy-ey):fy), FUN = median)
  aggBlow <- aggregate(data ~ iter + Ftarg, data = df, subset = slot == "ssb" & year %in% ((fy-ey):fy), FUN = quantile, prob = 0.05)
  #plot(data ~ Ftarg, data = aggBmed, pch = ".", cex = 5, col = adjustcolor(1,0.2))
  splBmed <- smooth.spline(x = aggBmed$Ftarg, y = aggBmed$data)
  splBlow <- smooth.spline(x = aggBlow$Ftarg, y = aggBlow$data)
  newdat$Bmed <- predict(splBmed, x = newdat$Ftarg)$y
  newdat$Blow <- predict(splBlow, x = newdat$Ftarg)$y
  #lines(Bmed ~ Ftarg, newdat)
  #lines(Blow ~ Ftarg, newdat, lty = 2)
  (Bmsy <- newdat$Bmed[which.max(newdat$landings)])

  iters = unique(df$iter)
  ssb_at_Fmsy = sapply(iters,function(x){
    tmp = aggBmed[aggBmed$iter == x,]
    #interpol.F = smooth.spline(x= tmp$Ftarg,y = tmp$data) 
    #ssb_at_Fmsy = predict(interpol.F,x = Fmsy)$y
    interpol.F = approxfun(x= tmp$Ftarg,y = tmp$data)
    ssb_at_Fmsy = interpol.F(v = Fmsy)
    return(ssb_at_Fmsy)
  })
  
  # calculate Btrigger as the 5% quantile
  Btrig = quantile(ssb_at_Fmsy,0.05)
  
  # with hcr  ------------
  
  # Set-up cluster and run
  clusterType <- ifelse(Sys.info()["sysname"] == "Windows", "PSOCK", "FORK")
  # ARGS <- ls() # to pass anything from the main environment (called below by clusterExport)
  ARGS =  c("om","iy","fy","info","sdevs","iem","Btrig")
  cl <- parallel::makeCluster(ncores, type=clusterType) # set-up cluster
  nn <- split(info$scen, info$scen) # indices for each iteration (used as "x" in the parFun)
  parallel::clusterExport(cl, varlist = ARGS, envir = environment()) # export info to the cluster
  t1 = Sys.time()
  RES2 <- parLapply(cl, nn, parFun_HCR,Btrigger = Btrig) # runs the function through each index in nn (index passed to "x" in parFun)
  stopCluster(cl) # stops cluster
  t2 = Sys.time()
  time.taken2 = t2-t1
    
  df2 <- lapply(as.list(seq(info$Ftarg)), FUN = function(x){
    tmp1 <- as.data.frame(RES2[[x]]$ssb)
    tmp1$slot <- "ssb"
    tmp2 <- as.data.frame(RES2[[x]]$landings)
    tmp2$slot <- "landings"
    tmp3 <- as.data.frame(RES2[[x]]$catch)
    tmp3$slot <- "catch"
    tmp4 <- as.data.frame(RES2[[x]]$rec)
    tmp4$slot <- "rec"
    tmp5 <- as.data.frame(RES2[[x]]$fbar)
    tmp5$slot <- "fbar"
    tmp <- rbind(tmp1, tmp2, tmp3, tmp4, tmp5)
    tmp$Ftarg <- info$Ftarg[x]
    return(tmp)
  })
  df2 <- do.call("rbind", df2)
  
  
  # Bmsy, Btrigger
  aggBlow <- aggregate(data ~ iter + Ftarg, data = df2, subset = slot == "ssb" & year %in% ((fy-ey):fy), FUN = quantile, prob = 0.05)
  splBlow <- smooth.spline(x = aggBlow$Ftarg, y = aggBlow$data)
  newdat$Blow <- predict(splBlow, x = newdat$Ftarg)$y
  
  # Fpa
  newdat$safe <- newdat$Blow > Bmsy*0.3 # compare to a newly estimated Blim
  (Fpa <- newdat$Ftarg[which(!newdat$safe)[1]])
  
  newBRPs_bench[stk_name, "Fpa"] <- Fpa
  newBRPs_bench[stk_name, "Fmsy"] <- Fmsy
  newBRPs_bench[stk_name, "FmsyL"] <- lFmsy
  newBRPs_bench[stk_name, "FmsyH"] <- uFmsy
  newBRPs_bench[stk_name, "Btrigger"] <- Btrig
  newBRPs_bench[stk_name, "Blim"] <- Bmsy*0.3
  
  # collate results of all oms:
  df$climScen = df2$climScen = "benchmark"
  df$period = df2$period = "None"
  df$stock = df2$stock = stk_name
  
  # track time
  tt = paste(round(time.taken1+time.taken2,digits = 3),units(time.taken1+time.taken2))
  names(tt) = paste(stk_name,"om","benchmark",sep = "_")
  
  # collate results of all oms:
  df$climScen = df2$climScen = "benchmark"
  df$period = df2$period = "None"
  df$stock = df2$stock = stk_name
  
  # track time
  tt = paste(round(time.taken1+time.taken2,digits = 3),units(time.taken1+time.taken2))
  names(tt) = paste(stk_name,"om","benchmark",sep = "_")
  
  # save results
  save(df, df2,tt, file = file.path("model","refpts", paste0(stk_name,"_SRR","_benchmark", "_model.rda")), compress="xz")

}

save(newBRPs_bench, file = file.path("model", paste0("newBRPs","_SRR","_benchmark", "_model.rda")), compress="xz")

