# ------------------------------------------------------------- #
#################################################################
# >>>>>>>>>> IMPORTANT NOTE FOR THIS TO RUN PROPERLY: <<<<<<<<< 
# make sure 'FLasher' and the 'mse' package are NOT loaded!!!!
# as those will provide conflicts with the same type of classes
# defined in 'mseEMSRR' and 'FLasherEMSRR'
#################################################################
# ------------------------------------------------------------- #

# remove all existing objects ----------------
rm(list = ls())

library(FLCore)
library(FLasherEMSRR)
library(mseEMSRR)
library(parallel)

# 1. required functions ------------------------------------------------------
source("utilities.R")

# parallelized functions to run MSE with/without HCR

# MSE with HRC
parFun_noHCR <- function(x){
  library(mseEMSRR)
  library(FLasherEMSRR)
  
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
                   args=list(recyrs=10, fmin=0, Fdevs=sdevs$F)) # use last 10 years for geometric mean rec. assumption
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
  library(mseEMSRR)
  library(FLasherEMSRR)
  
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
                   args=list(recyrs=10, fmin=0, Fdevs=sdevs$F)) # 10-yr geometric mean for short term forecast
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

# NUMBER of evaluation years (last years of the time period)
ey <- 20

# NUMBER of cores
ncores <- 20

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
# -------------------------------------------------- #

# run for BH and Segreg
srr.type = c("BH","SEG")

newBRPs_list <- list() # to store the new ref-pts
ii = 1
for(stk_name in rownames(BRPs)){ # substock loop
  for(srt in srr.type){
    # load stock object
  
  # stk_name <- "PLE-NS"
  om.obj.nms = load(file.path("model", paste0("Cod_",substock,"_EMSRR_",srt, "_ccrefpts_input.rda")), verbose = T)
  
  # BRPs with consideration of climate scenarios
  tmp_BRPs = cbind(BRPs[stk_name,]*NA,cc)
  tmp_BRPs$srr.type = srt # specify srr type in the refpts table
  # run for different oms
  # df_oms = df2_oms =list()
  # tt = vector()
  for(n in 1:nrow(cc)){
    
    # get the respective om
    rcp = cc$climScen[n]
    period = cc$period[n]
  
    cat("Running:",stk_name,"|",srt,"|",paste("om",rcp,period,sep = "_"),"...","\n")
    om = get(paste("om",rcp,period,sep = "_"))
    
    DIMS <- dims(om@stock)
    RANGE <- range(om@stock)
    
    # DATA year
    dy <- iy - 1
    
    # FINAL year
    fy <- RANGE["maxyear"]
    
    # NUMBER of iterations
    it <- DIMS$iter
    
    # conduct forecasts across Ftarget levels without HCR and later with HCR
    info <- data.frame(Ftarg = seq(0, 1, by = 0.02),
                       stk_name = stk_name)
    ## for testing use a courser grid
    #  info <- data.frame(Ftarg = seq(0, 1, by = 0.5),
    #                     stk_name = stk_name)
    info$scen <- seq(nrow(info))
    head(info)
    nrow(info)
    
    # -------- #
    # simulate
    # -------- #
    
  # no hcr ------------
  
  # Set-up cluster and run
  clusterType <- ifelse(Sys.info()["sysname"] == "Windows", "PSOCK", "FORK")
  # ARGS <- ls() # to pass anything from the main environment (called below by clusterExport)
  # only pass relevant ARGs to avoid memory cluttering
  ARGS =  c("om","iy","fy","info","sdevs","iem",
            grep(paste(stk_name,srt,sep = "_"),grep("sr_predict|rec.model|predict.func",ls(),value = T),value = T))
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
  df <- rbindlist(df)
  
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
  #(Btrig <- newdat$Blow[which.max(newdat$landings)])
  #abline(h = Btrig,col = "red",lty = 2)
  #abline(v = Fmsy,col = "red",lty = 2)
  
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
  ARGS =  c("om","iy","fy","info","sdevs","iem","Btrig",
            grep(sub("-.*","",stk_name),grep("sr_predict|rec.model|predict.func",ls(),value = T),value = T))
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
  df2 <- rbindlist(df2)
  
  # Bmsy, Btrigger
  aggBlow <- aggregate(data ~ iter + Ftarg, data = df2, subset = slot == "ssb" & year %in% ((fy-ey):fy), FUN = quantile, prob = 0.05)
  splBlow <- smooth.spline(x = aggBlow$Ftarg, y = aggBlow$data)
  newdat$Blow <- predict(splBlow, x = newdat$Ftarg)$y
  
  # Fpa
  newdat$safe <- newdat$Blow > Bmsy*0.3 # compare to a newly estimated Blim
  (Fpa <- newdat$Ftarg[which(!newdat$safe)[1]])
  
  tmp_BRPs[tmp_BRPs$climScen == rcp &  tmp_BRPs$period == period, "Fpa"] <- Fpa
  tmp_BRPs[tmp_BRPs$climScen == rcp &  tmp_BRPs$period == period, "Fmsy"] <- Fmsy
  tmp_BRPs[tmp_BRPs$climScen == rcp &  tmp_BRPs$period == period, "FmsyL"] <- lFmsy
  tmp_BRPs[tmp_BRPs$climScen == rcp &  tmp_BRPs$period == period, "FmsyH"] <- uFmsy
  tmp_BRPs[tmp_BRPs$climScen == rcp &  tmp_BRPs$period == period, "Btrigger"] <- Btrig
  tmp_BRPs[tmp_BRPs$climScen == rcp &  tmp_BRPs$period == period, "Blim"] <- Bmsy*0.3

  # collate results of all oms:
  df$climScen = rcp
  df$period = period
  df$stock = stk_name
  df2$climScen = rcp
  df2$period = period
  df2$stock = stk_name
  
  # track time
  tt = paste(round(time.taken1+time.taken2,digits = 3),units(time.taken1+time.taken2))
  names(tt) = paste(stk_name,"om",rcp,period,sep = "_")
  
  # save individual oms, otherwise files get too big
  save(df, df2,tt, file = file.path("model","refpts", paste(stk_name,srt,rcp,period,"EMSRR","model.rda",sep = "_")), compress="xz")
  
  # remove to free memory
  rm(RES,RES2,df,df2)
  
  }
  # # bind together
  # df = rbindlist(df_oms)
  # df2 = rbindlist(df2_oms)
  
  # store ref-pts changes
  tmp_BRPs$stock = stk_name
  newBRPs_list[[ii]] = tmp_BRPs
  ii = ii + 1 
  # free memory and remove om-objects
  rm(list = om.obj.nms)
  }
}
newBRPs = rbindlist(newBRPs_list)
# write out ref-pt changes
save(newBRPs,newBRPs_list, file = file.path("model", paste0("newBRPs","_EMSRR", "_model.rda")), compress="xz")
