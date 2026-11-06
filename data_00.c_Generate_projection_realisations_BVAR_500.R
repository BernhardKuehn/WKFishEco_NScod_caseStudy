# =============================================== #
# Generate new artificial time series based
# on autocorrelation and interdependencies between/within 
# spatial fields
# using the BVAR-procedure 
# last change: 2026-11-06
# author: Bernhard Kuehn
# prepare 500 spatial fields
# =============================================== #

rm(list = ls())

library(raster)
library(marmalaid)
library(BVAR)
library(terra)
library(corrplot)
library(tsDyn)

# create and define directory to store output
store.dir = "./data/Env/POLCOMS-ERSEM_biasRM_BVARsims/"
dir.create(store.dir,showWarnings = F)

# -------------------------- #
# load bias corrected fields
# -------------------------- #

# Bottom temp
load("./data/Env/POLCOMS-ERSEM_biasRM/POLCOM_ERSEM_BottomTemp_raster_meanBiasRM.RData")

# --------------------- #
# load helper functions
# --------------------- #

# function for fitting and simulating a spatial BVAR
source("./functions/spatial.BVAR.R")

apply.rasterFunc.on.nestedList = function(x,fn){
  lapply(x,function(xx){
    if(is.list(xx)){
      apply.rasterFunc.on.nestedList(xx,fn = fn)
    } else if(class(xx) %in% c("RasterLayer","RasterBrick","RasterStack")){
      fn(xx) # apply function on raster
    } else if(!is.list(xx) && !(class(xx) %in% c("RasterLayer","RasterBrick","RasterStack"))) {
      xx
    }
  })
  
}


# -------- #
# analysis
# -------- #

clim.proj = c("rcp45","rcp85")
# climate projection loop
for(ii in 1:length(clim.proj)){
  
  rcp = clim.proj[ii]
  
  # -------------#
  ### fit BVAR ###
  
  # put raster in list
  RCP_rst.list = list(bottomT = POLCOM_ERSEM_BottomTemp_raster_meanBiasRM[[rcp]])
                      
  time_rcp = sub("X","",names(POLCOM_ERSEM_BottomTemp_raster_meanBiasRM[[rcp]]))
  
  # generate anomalies (otherwise the seasonal signal will overrule the subsequent PCA analysis)
  RCP_anom_rst.list = lapply(RCP_rst.list,function(x){
    if(class(x) == "list"){
      # simplify
      x = lapply(x,function(xx){
        xx.mat = data.frame(coordinates(xx),as.matrix(xx))
        xx = rasterFromXYZ(xx.mat[complete.cases(xx.mat[,-c(1:2)]),])
        return(xx)})
      # calc anomalies
      out = lapply(x,fieldAnomaly.raster,
                   time = as.Date(time_rcp,format = "%Y.%m.%d"),
                   level = "month")
    } else {
      # simplify raster (when it contains NAs at the edges)
      x.mat = data.frame(coordinates(x),as.matrix(x))
      x = rasterFromXYZ(x.mat[complete.cases(x.mat[,-c(1:2)]),])
      # calc anomalies
      out = fieldAnomaly.raster(x,
                                time = as.Date(time_rcp,format = "%Y.%m.%d"),
                                level = "month")
    }
    return(out)})
  
  
  
  
  # fit spatial BVAR
  RNGkind("Mersenne-Twister")
  set.seed(12345) # set seed to get consistent BVAR (which can be reproduced)
  BVAR_rcp = fit.spatial.BVAR(raster = RCP_anom_rst.list,
                              time_raster = as.Date(time_rcp,format = "%Y.%m.%d"),
                              PC.thresh = 0.9,min.PCs = 5, lags = 1,
                              n_draw = 20000, n_burn = 10000, n_thin = 1,plot = T)
  # check nr. of PCs that got selected with this threshold
  BVAR_rcp$nr_pcs_selected
  # check size of the BVAR
  sum(as.integer(BVAR_rcp$nr_pcs_selected[1,]))*BVAR_rcp$BVAR_model$meta$lags
  
  
  # check convergence of the BVAR
  run_mcmc <- as.mcmc(BVAR_rcp$BVAR_model, vars = "lambda")
  geweke.diag(run_mcmc) # ok (between critical values of 1.96)
  geweke.plot(run_mcmc) # looks ok
  BVAR_rcp$diagnostics
  
  # ----------------- #
  # write out and save
  assign(x = paste("BVAR",rcp,"bottomT",sep = "_"),
         value = BVAR_rcp)
  save(list = paste("BVAR",rcp,"bottomT",sep = "_"),
       file = file.path("./data/",paste0(paste("BVAR",rcp,"bottomT",sep = "_"),".RData")))
  
  # ----------------------------- #
  ### generate new realisations ###
  # 500 new ones
  # problem generating all these 500 in one step
  # just try to do it in batches of 25 (works better and you dont run so fast into out of memory errors)
  
  # also write out in Batch sizes....
  nits = 500
  
  run.in.batch = TRUE
  if(run.in.batch == T){
    batch.size = 25
    n.batch = nits/batch.size
    sims = list()
    for(nb in 1:n.batch){
      cat("Calculating batch",paste0(nb,"/",n.batch),"...","\n")
      
      # when in parallel mode, set it to a RNGkind that is compatible with parallel computation
      RNGkind("L'Ecuyer-CMRG")
      set.seed(42*nb)
      BVAR_RCP_sims = sim.spatial.BVAR(BVAR_rcp,start.val.type = "mean",N.sim = batch.size,L = length(time_rcp),return.raster = T,
                                       parallel = T,ncores = 4)
      # get the sims and store them
      sims[[1]] = BVAR_RCP_sims$sims
      # count the current iteration/replication
      it = (nb-1)*batch.size+(1:batch.size)
      # free disc-space
      rm(BVAR_RCP_sims)
      
      # restructure:
      # concetinate and bring to old format (to have all simulations of 1 variable together)
      sims_restructured = list()
      for(s in seq(sims)){
        tmp = sims[[s]]
        names.vars = names(tmp)
        if(s == 1){
          sims_restructured = tmp 
        } else{
          for(v in names.vars){
            sims_restructured[[v]] = c(sims_restructured[[v]],tmp[[v]]) 
          }
        }
      }
      sims = sims_restructured
      
      # check if I got unwanted dependencies/correlations in there due to seed dependencies
      for(v in seq(sims)){
        tmp = sims[[v]]
        df = list()
        for(iter in seq(tmp)){
          if(class(tmp[[iter]]) == "list"){
            df[[iter]] = do.call(cbind,lapply(tmp[[iter]],function(x) terra::global(x,mean,na.rm = T))) 
          } else {
            df[[iter]] = terra::global(tmp[[iter]],mean,na.rm = T)  
          }
        }
        df = do.call(cbind,df)
        cors = cor(df)
        diag(cors) = 0
        colnames(cors) = rownames(cors) = 1:nrow(cors)
        par(mfrow = c(1,2))
        corrplot(cors,method = "square",order = "original",title = names(sims)[v])
        boxplot(as.vector(cors))
      }
      # --------------------------------- #
      # output simulation files and store
      # --------------------------------- #
      
      # BVAR calculations have been done on seasonal corrected data, however season needs to be added back now
      # again....
      
      dir.create(store.dir)
      for(i in seq(sims)){
        tmp = sims[[i]]
        var = names(sims)[i]
        # add back the seasonal component
        # - get the anomaly
        anom = RCP_anom_rst.list[[var]]
        # - get means
        if(class(anom) != "list"){
          seasonal.means = attributes(anom)$grid.means
        } else{
          seasonal.means = lapply(anom,function(x) attributes(x)$grid.means)
        }
        
        
        dir.create(paste0(store.dir,var))
        for(k in 1:length(tmp)){
          cat(var,"|",k,"...","\n")
          if(class(tmp[[k]]) == "list"){
            for(n in seq(tmp[[k]])){
              
              # add seasonal signal (add slot grid.means to attributes)
              attributes(tmp[[k]][[n]])$grid.means <- seasonal.means[[n]]
              tmp.seasonal <- fieldAnomaly.to.original(tmp[[k]][[n]],
                                                       time = as.Date(names(tmp[[k]][[n]])),
                                                       level = "month")
              fl.nm = paste("POLCOM_ERSEM",rcp,
                            ifelse(n ==1,paste(var,"u",sep = "."),
                                   paste(var,"v",sep = ".")),it[k],sep = "_")
              cat("Writing out:",">>",fl.nm,"<<","\n")
              terra::writeRaster(tmp.seasonal,
                                 filename = file.path(paste0(store.dir,var),
                                                      paste0(fl.nm,".tiff")),
                                 overwrite = T)
            }
          } else{
            # add seasonal signal
            attributes(tmp[[k]])$grid.means <- seasonal.means
            tmp.seasonal <- fieldAnomaly.to.original(tmp[[k]],
                                                     time = as.Date(names(tmp[[k]])),
                                                     level = "month")
            fl.nm = paste("POLCOM_ERSEM",rcp,var,it[k],sep = "_")
            cat("Writing out:",">>",fl.nm,"<<","\n")
            terra::writeRaster(tmp.seasonal,
                               filename = file.path(paste0(store.dir,var),
                                                    paste0(fl.nm,".tiff")),
                               overwrite = T)
          }
        } # end of inner write-out loop
      } # end of outer write-out loop
    } # end of batch loop
  } else{
    cat(paste("Generating",nits,"number of realisations is hardly possible in sequential mode!!! \n
        Please rethink your choice!","\n"))
  }
} # end of rcp loop

# ---------------------------------------------------------------------------------------------------------------- #