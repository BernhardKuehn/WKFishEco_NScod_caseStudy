# ================================================================================ #
# Prepare projections of env. signals for their use in the reference pt estimation
# author: Bernhard Kuehn
# last change: 2026-11-06
# =========================================================== # 

# ------------------------------------------------------------------ #
# steps: 
# 1. load proj. fields
# 2. load spatial-mean/clustered pattern or regions for spatial avg. 
# 3. process proj. the same way that the historical fields have been processed
# 4. project fields onto historical fields
# 5. extract the relevant variables needed for the different stocks
# 6. write out and store
# ------------------------------------------------------------------ #

rm(list = ls())

library(raster)
library(terra)
library(marmalaid)
library(stringr)
library(DescTools)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(Rcpp2doParallel)

# ---------------- #
# helper functions
# ---------------- # 

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

# # fast matrix multiplication in C++
# cppFunction(depends = c("RcppArmadillo","RcppEigen"),
#             'SEXP matMult_rcpp(Eigen::MatrixXd A, Eigen::MatrixXd B){
#             Eigen::MatrixXd C = A * B;
#             return Rcpp::wrap(C);
#             }')


# ---------------------- #
# load historical fields
# ---------------------- #

# load original fields
load("./data/Env/processed/Spatial_bottomT_yr_by_cod_substock.rda")

# load stored clusters
load("./data/Env/processed/SpatialCluster_bottomT_yr_by_cod_substock.rda")


env.hist.dat = list(bottomT = bottomT_yr_by_substock_spatial)

# ------------------------------------------------ #
# read in BVAR-generated files and perform overlay
# ------------------------------------------------ #

dir.proj = "./data/Env/POLCOMS-ERSEM_biasRM_BVARsims/"

# load a file from each subdirectory
sub.dir.proj = list.files(dir.proj)

# in order for this to work the folder names should match the historical file names
if(!all(names(env.hist.dat) %in% sub.dir.proj)){
  stop("Please make sure that the filenames match the list names, otherwise correct matching cannot be quaranted!")
} else {
  # continue from here and perform overlay
  
  # --------------------- #
  # generate overlay mask
  # --------------------- #
  extent_min_list = list()
  overlay.masks = list()
  for(i in seq(sub.dir.proj)){
    tmp.dir = file.path(dir.proj,sub.dir.proj[i])
    tmp = rast(file.path(tmp.dir,list.files(tmp.dir)[1]))
    
    # match
    var = sub.dir.proj[i]
    
    # eof?
    if(class(env.hist.dat[[var]]) == "list"){
      hist = sum(!is.na(brick(env.hist.dat[[var]][[1]])))
    } else{
      hist = sum(!is.na(brick(env.hist.dat[[var]])))
    }
    proj = brick(sum(!is.na(tmp)))
    
    # cut to same extent
    xy_min.max = rbind(as.vector(extent(hist)),as.vector(extent(proj)))
    extent_min = extent(c(max(xy_min.max[,1]),
                          min(xy_min.max[,2]),
                          max(xy_min.max[,3]),
                          min(xy_min.max[,4])))
    # store the extent_min for latter use
    extent_min_list[[i]] = extent_min
    names(extent_min_list)[i] = var
    # crop
    hist = crop(hist,extent_min)
    proj = crop(proj,extent_min)
    
    # create overlay 
    over.hist.proj = overlay(x = hist,
                             y = proj,fun = sum)
    mask.hist.proj = over.hist.proj !=cellStats(over.hist.proj,max)
    
    # store mask
    overlay.masks[[i]] = mask.hist.proj
    names(overlay.masks)[i] = var
  }
  par(mfrow = c(2,2),mar = c(1,1,2,1))
  lapply(overlay.masks,image)
  
  # -------------------------------- #
  # process env. fields through mask
  # -------------------------------- #
  
  #env.hist.dat_season_mask = env.hist.dat_season
  env.hist.dat_mask = env.hist.dat
  for(i in seq(env.hist.dat)){
    hist = env.hist.dat[[i]]
    var = names(env.hist.dat)[i]
    print(var)
    # get corresponding overlay mask
    mask.o = overlay.masks[[var]]
    if(class(hist) == "list"){
      tmp.out = apply.rasterFunc.on.nestedList(hist,
                                               fn = function(x) mask(crop(x,extent(mask.o)),mask = mask.o,
                                                                     maskvalue = 1,inverse = F))
    } else {
      tmp.out = mask(crop(hist,extent(mask.o)),mask = mask.o,
                     maskvalue = 1,inverse = F)
    }
    
    # store 
    env.hist.dat_mask[[i]] = tmp.out
  }
  
  # -------------- #
  # calc anomalies
  # -------------- #
  
  env.hist.dat_anom = lapply(env.hist.dat_mask,function(x){
    rst.anom = list()
    if(class(x) == "list"){
      time.rst = as.Date(sub("X","",names(x[[1]])),format = "%Y.%m.%d")
      rst.anom = lapply(x, function(xx)
        fieldAnomaly.raster(rst1 = xx,
                            time = time.rst,level = "month"))
    } else {
      time.rst = as.Date(sub("X","",names(x)),format = "%Y.%m.%d")
      rst.anom = fieldAnomaly.raster(rst1 = x,
                                     time = time.rst,level = "month")
      
    }
    return(rst.anom)
  })
  
  
  # ----------------------- #
  # calculate seasonal avg.
  # ----------------------- #
  
  # apply on list
  env.hist.dat_season_mask = lapply(env.hist.dat_anom,
                               fn_seasonal,ss = seasons)
  
  
  # ------------------------- #
  # calculate historical EOFs
  # ------------------------- #
  
  env.hist_eofs = list()
  for(i in seq(env.hist.dat_season_mask)){
    tmp = env.hist.dat_season_mask[[i]]
    var = names(env.hist.dat_season_mask)[i]
    if(var == "currents"){
      # since currents enters the spatial.PCA function as list-obj,
      # the 'apply.raster.on.nested.list' function does not work here
      out = list()
      for(j in seq(tmp)){
        out[[j]] = spatial.PCA(x = tmp[[j]],center = T,scale. = F,
                    spatial.extent = extent(tmp[[j]][[1]]),plot = F)
        names(out)[j] = names(tmp)[j]
      }
    } else{
     out =  apply.rasterFunc.on.nestedList(tmp,
                                     fn = function(x) spatial.PCA(x = x,center = T,scale. = F,
                                                                  spatial.extent = extent(x),plot = F))
      
    }
    env.hist_eofs[[i]] = out
    names(env.hist_eofs)[i] = var
  }
  
  # ------------------------------------------------ #
  # added on the 2023-11-30:
  # check flip of PCs of the 'new' cropped version
  # with the 'old' uncropped one
  # ------------------------------------------------ #
  vars = lut.eof.pattern$new
  flip.vars = list()
  for(vv in seq(vars)){
  var = vars[vv]
  old.var = lut.eof.pattern$old[lut.eof.pattern$new == var]
  old.eofs = EOF.pattern$EOFs[[old.var]]
  new.eofs = env.hist_eofs[[var]]
  seasons.eofs = names(old.eofs)
  flip.list = list()
  for(z in 1:length(seasons.eofs)){
    seas = seasons.eofs[z]
    # check the first 9 PCs
    pcs.old = old.eofs[[seas]]$sp.PCA$x[,1:10]
    pcs.new = new.eofs[[seas]]$sp.PCA$x[,1:10]
    par(mfrow = c(3,4),mar = c(2,2,4,2),oma = c(0,0,2,0))
    
    cors = vector()
    for(pc in 1:10){
      # calc correlations
      check1 = pcs.old[,pc]
      check2 = pcs.new[,pc]
      if(length(check1) != length(check2)){
      # simple cut the longer one
      min.l = pmin(length(check1),length(check2))
      check1 = check1[1:min.l]
      check2 = check2[1:min.l]
      }
      cors[pc] = cor(check1,check2)
      
      plot(check1,type = "l")
      lines(check2,col = "red")
      title(paste0("PC",pc),adj = 0)
    }
    mtext(text = paste(var,seas,sep = "-"),side = 3,line = 0,outer = T,font = 2)
    flip = cors<0.3
    names(flip) = paste0("PC",1:10)
    flip.list[[z]] = flip
  }  
  names(flip.list) = seasons.eofs
  flip.vars[[vv]] = do.call(cbind,flip.list)
  }
  names(flip.vars) = vars 
  # stored matrices showing which PCs are flipped 
  # (need to be unflipped)
  
# ============================================== #
# 2. Project MPIOM data onto historical EOFs
# ============================================== #
  
# sub steps: 
# - use anomalies of proj
# - calc. seasonal avg
# - mask with overlay mask
# - project onto seasonal eofs, store only the first 10 pattern
# nr. PCs 2 (proj) x 4 (vars) x 5 (seasons) x 10 (eofs) x 100 (realisations) = 4E5 time series

# !!! this part is a bit slow...think about running in parallel !!!
  
clim.proj = c("rcp45","rcp85")
# loop
EOF.pattern_proj = vector("list",length(sub.dir.proj))
names(EOF.pattern_proj) = sub.dir.proj


# parallelize outer loop
library(snow)
library(doSNOW)
library(foreach)

cl = makeCluster(4,type = "SOCK")
registerDoSNOW(cl)


EOF.pattern_proj = foreach(i = seq(sub.dir.proj),.combine = c,.inorder = T,
                           .packages = c("raster","terra","marmalaid","stringr","DescTools")) %dopar% {
#for(i in seq(sub.dir.proj)){
  var = sub.dir.proj[i]
  EOF.pattern_proj = vector("list",1)
  names(EOF.pattern_proj) = sub.dir.proj[i]
  # get realisations
  proj.files = list.files(file.path(dir.proj,sub.dir.proj[i]))
  EOF.pattern_proj[[var]] = vector("list",length(clim.proj))
  names(EOF.pattern_proj[[var]]) = clim.proj
  for(rcp in clim.proj){
    # get each climate projection seperately
    proj.rcp.sub = proj.files[grep(rcp,proj.files)]
    # order based on numbering
    proj.rcp.sub = proj.rcp.sub[order(as.numeric(gsub("_","",str_extract(proj.rcp.sub,"_[0-9]+"))))]
    # remove auxiliary files (aux.xml)
    proj.rcp.sub = proj.rcp.sub[!grepl("xml",proj.rcp.sub)]
    # read in and calculate seasonal avg. 
    if(var == "currents"){
      n.iters = length(proj.rcp.sub)/2
      # split proj.rcp.sub
      u.runs = grep(".u",proj.rcp.sub,value = T)
      v.runs = grep(".v",proj.rcp.sub,value = T)
      #uv.runs = list(u = u.runs,v = v.runs)
      } else{
      n.iters = length(proj.rcp.sub)
      }
    store.list = vector("list",length = length(seasons))
    
    for(it in 1:n.iters){
      # do all the calculations in terra and only swith for the seasonal 
      # calculations to raster, cos otherwise it is too slow.... 
      if(var == "currents"){
        
        tmp.rst = lapply(list(u = u.runs[it],v = v.runs[it]),
                        function(x) rast(file.path(dir.proj,var,x)))
        time.rst = as.Date(sub("X","",names(tmp.rst[[1]])),format = "%Y.%m.%d")
        # crop
        tmp.rst = lapply(tmp.rst,function(x) crop(x,extent_min_list[[var]]))
        # mask with overlay
        tmp.rst = lapply(tmp.rst,function(x) mask(x,mask = rast(overlay.masks[[var]]),
                       maskvalue = 1))
        # calc anomalies
        for(ll in 1:length(tmp.rst)){
          tmp.rst[[ll]] = fieldAnomaly.raster(rst1 = tmp.rst[[ll]],
                                        rst2.anomalies = env.hist.dat_anom[[var]][[ll]],
                                        time = time.rst,level = "month")
        }
        # calculate seasonal means
        #tmp.rst = lapply(tmp.rst,brick)
        tmp.rst_seasonal = fn_seasonal(tmp.rst,ss = seasons)
      } else {
        tmp.rst = rast(file.path(dir.proj,var,proj.rcp.sub[it]))
        time.rst = as.Date(sub("X","",names(tmp.rst)),format = "%Y.%m.%d")
        # crop
        tmp.rst = crop(tmp.rst,extent_min_list[[var]])
        # mask with overlay
        tmp.rst = mask(tmp.rst,mask = rast(overlay.masks[[var]]),
                     maskvalue = 1)
        # !!! calculate anomalies!!!!
        # using the means of the historical period (added additional functionality in the marmalaid pkg)
        tmp.rst = fieldAnomaly.raster(rst1 = tmp.rst,
                          rst2.anomalies = env.hist.dat_anom[[var]],
                          time = time.rst,level = "month")
      
      # calculate seasonal means
      #tmp.rst = brick(tmp.rst)
      tmp.rst_seasonal = fn_seasonal(tmp.rst,ss = seasons)
      }
      # project onto historical eofs
      for(ss in 1:length(seasons)){
        cat(var,"|",rcp,"|",names(seasons)[ss],"|",it,"\n")
        #time.hist = env.hist_eofs[[var]][[names(seasons)[ss]]]$time
        tmp.eof = env.hist_eofs[[var]][[names(seasons)[ss]]]$sp.PCA
        if(var == "currents"){
        input = lapply(tmp.rst_seasonal[[names(seasons)[ss]]],
                 function(x){
                   mat = as.matrix(x)
                   mat = mat[complete.cases(mat),]
                   return(mat)
                 })
        input = t(do.call(rbind,input))
          
        } else {
          input = as.matrix(tmp.rst_seasonal[[names(seasons)[ss]]])
          input = input[complete.cases(input),]
          input = t(input)
        }
        # get PCs and select only the first 10
        tmp.pcs = predict(tmp.eof,input)[,1:10]
        # ------------------- #
        # added on 2023-11-30: 
        # correct flipped PCs
        # ------------------- #
        
        # flip PCs if necessary: 
        flip.multiplier = ifelse(flip.vars[[var]][,names(seasons)[ss]] == FALSE,1,-1)
        tmp.pcs = t(t(tmp.pcs) * flip.multiplier)
        
        #time.proj = as.numeric(sub("X","",rownames(tmp.pcs)))
        # # compare
        # par(mfrow = c(3,3),mar = c(2,2,3,2))
        # for(kk in 1:9){
        #   plot(time.hist,tmp.eof$x[,kk],type= "l",xlim = c(min(time.eof),2100))
        #   lines(time.proj,tmp.pcs[,kk],col = "red")
        # }
        
        # store in array of form [time,PC,realisation]
        if(it == 1){
          store.list[[ss]] = tmp.pcs
        } else {
          store.list[[ss]] = Abind(store.list[[ss]],tmp.pcs,along = 3)
        }
      }
      names(store.list) = names(seasons)
    }
    EOF.pattern_proj[[var]][[rcp]] = store.list
  } 
  # pass element to the outer loop of the parallelisation routine
  EOF.pattern_proj
}
stopCluster(cl)
# check the different realisations 

par(mfrow = c(3,3))
for(pc in 1:9){
plot(env.hist_eofs$bottomT$DJF$time,
     env.hist_eofs$bottomT$DJF$sp.PCA$x[,pc],lwd = 2,
     type ="l",xlim = c(1950,2100),ylim = c(-100,100))
# add old PC  
lines(EOF.pattern$EOFs$bottomT$DJF$time,
      EOF.pattern$EOFs$bottomT$DJF$sp.PCA$x[,pc],lwd = 2,col ="limegreen",lty = 2)
  
time.proj = as.numeric(sub("X","",rownames(EOF.pattern_proj$bottomT$rcp45[[1]])))
for(i in 1:100)
lines(time.proj,EOF.pattern_proj$bottomT$rcp45[[1]][,pc,i],col = "red")
}


par(mfrow = c(3,3),oma = c(0,0,1,0),mar = c(1,1,1,1))
for(pc in 1:9){
  plot(env.hist_eofs$bottomT$DJF$time,
       env.hist_eofs$bottomT$DJF$sp.PCA$x[,pc],lwd = 2,
       type ="l",xlim = c(1950,2100),ylim = 2*range(env.hist_eofs$bottomT$DJF$sp.PCA$x[,pc]))
  # add original EOF-pattern
  lines(EOF.pattern$EOFs$bottomT$DJF$time,EOF.pattern$EOFs$bottomT$DJF$sp.PCA$x[,pc],col = "deeppink")
  
  time.proj = as.numeric(sub("X","",rownames(EOF.pattern_proj$bottomT$rcp85[[1]])))
  for(i in 1:100){
    lines(time.proj,EOF.pattern_proj$bottomT$rcp45[[1]][,pc,i],col = "olivedrab")
    lines(time.proj,EOF.pattern_proj$bottomT$rcp85[[1]][,pc,i],col = "red")
  }
  # add mean trajectory
  lines(time.proj,apply(EOF.pattern_proj$bottomT$rcp45$DJF[,pc,],1,mean),col = "turquoise",lwd = 2)
  lines(time.proj,apply(EOF.pattern_proj$bottomT$rcp85$DJF[,pc,],1,mean),col = "tan1",lwd = 2)
  lines(env.hist_eofs$bottomT$DJF$time,
        env.hist_eofs$bottomT$DJF$sp.PCA$x[,pc],lwd = 2)
}
mtext("bottomT DJF",side = 3,line = -1,font = 2,outer = T)

par(mfrow = c(3,3),oma = c(0,0,1,0))
for(pc in 1:9){
  plot(env.hist_eofs$currents$JJA$time,
       env.hist_eofs$currents$JJA$sp.PCA$x[,pc],lwd = 2,
       type ="l",xlim = c(1950,2100),ylim = c(-5,5))
  # add original EOF-pattern
  lines(EOF.pattern$EOFs$Currents$JJA$time,EOF.pattern$EOFs$Currents$JJA$sp.PCA$x[,pc],col = "deeppink")
  
  time.proj = as.numeric(sub("X","",rownames(EOF.pattern_proj$currents$rcp45[[3]])))
  for(i in 1:100){
    lines(time.proj,EOF.pattern_proj$currents$rcp45[[3]][,pc,i],col = "olivedrab")
    lines(time.proj,EOF.pattern_proj$currents$rcp85[[3]][,pc,i],col = "red")
  }
  lines(time.proj,apply(EOF.pattern_proj$currents$rcp45$JJA[,pc,],1,mean),col = "turquoise",lwd = 2)
  lines(time.proj,apply(EOF.pattern_proj$currents$rcp85$JJA[,pc,],1,mean),col = "tan1",lwd = 2)
  lines(env.hist_eofs$currents$JJA$time,
        env.hist_eofs$currents$JJA$sp.PCA$x[,pc],lwd = 2)
}
mtext("Currents JJA",side = 3,line = -1,font = 2,outer = T)


par(mfrow = c(3,3),oma = c(0,0,1,0))
for(pc in 1:9){
  plot(env.hist_eofs$salt$DJF$time,
       env.hist_eofs$salt$DJF$sp.PCA$x[,pc],lwd = 2,
       type ="l",xlim = c(1950,2100),ylim = c(5*range(env.hist_eofs$salt$DJF$sp.PCA$x[,pc])))
  # add original EOF-pattern
  lines(EOF.pattern$EOFs$Salt$DJF$time,EOF.pattern$EOFs$Salt$DJF$sp.PCA$x[,pc],col = "deeppink")
  
  time.proj = as.numeric(sub("X","",rownames(EOF.pattern_proj$currents$rcp45[[3]])))
  for(i in 1:100){
    lines(time.proj,EOF.pattern_proj$salt$rcp45$DJF[,pc,i],col = "olivedrab")
    lines(time.proj,EOF.pattern_proj$salt$rcp85$DJF[,pc,i],col = "red")
  }
  lines(time.proj,apply(EOF.pattern_proj$salt$rcp45$DJF[,pc,],1,mean),col = "turquoise",lwd = 2)
  lines(time.proj,apply(EOF.pattern_proj$salt$rcp85$DJF[,pc,],1,mean),col = "tan1",lwd = 2)
  
  lines(env.hist_eofs$salt$DJF$time,
        env.hist_eofs$salt$DJF$sp.PCA$x[,pc],lwd = 2)
}
mtext("Salinity DJF",side = 3,line = -1,font = 2,outer = T)


par(mfrow = c(3,3),oma = c(0,0,1,0))
for(pc in 1:9){
  plot(env.hist_eofs$salt$YR$time,
       env.hist_eofs$salt$YR$sp.PCA$x[,pc],lwd = 2,
       type ="l",xlim = c(1950,2100),ylim = c(5*range(env.hist_eofs$salt$YR$sp.PCA$x[,pc])))
  # add original EOF-pattern
  lines(EOF.pattern$EOFs$Salt$YR$time,EOF.pattern$EOFs$Salt$YR$sp.PCA$x[,pc],col = "deeppink")
  
  time.proj = as.numeric(sub("X","",rownames(EOF.pattern_proj$currents$rcp45[[3]])))
  for(i in 1:100){
    lines(time.proj,EOF.pattern_proj$salt$rcp45$YR[,pc,i],col = "olivedrab")
    lines(time.proj,EOF.pattern_proj$salt$rcp85$YR[,pc,i],col = "red")
  }
  lines(time.proj,apply(EOF.pattern_proj$salt$rcp45$YR[,pc,],1,mean),col = "turquoise",lwd = 2)
  lines(time.proj,apply(EOF.pattern_proj$salt$rcp85$YR[,pc,],1,mean),col = "tan1",lwd = 2)
  
  lines(env.hist_eofs$salt$YR$time,
        env.hist_eofs$salt$YR$sp.PCA$x[,pc],lwd = 2)
}
mtext("Salinity YR",side = 3,line = -1,font = 2,outer = T)


par(mfrow = c(3,3),oma = c(0,0,1,0))
for(pc in 1:9){
  plot(env.hist_eofs$sst$YR$time,
       env.hist_eofs$sst$YR$sp.PCA$x[,pc],lwd = 2,
       type ="l",xlim = c(1950,2100),ylim = c(5*range(env.hist_eofs$sst$YR$sp.PCA$x[,pc])))
  time.proj = as.numeric(sub("X","",rownames(EOF.pattern_proj$currents$rcp45[[3]])))
  for(i in 1:100){
    lines(time.proj,EOF.pattern_proj$sst$rcp45$YR[,pc,i],col = "olivedrab")
    lines(time.proj,EOF.pattern_proj$sst$rcp85$YR[,pc,i],col = "red")
  }
  lines(time.proj,apply(EOF.pattern_proj$sst$rcp45$YR[,pc,],1,mean),col = "turquoise",lwd = 2)
  lines(time.proj,apply(EOF.pattern_proj$sst$rcp85$YR[,pc,],1,mean),col = "tan1",lwd = 2)
  
  lines(env.hist_eofs$sst$YR$time,
       env.hist_eofs$sst$YR$sp.PCA$x[,pc],lwd = 2)
  
}
mtext("SST YR",side = 3,line = -1,font = 2,outer = T)


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# implement better plotting with ggplot!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# decide on an output format!!!!
# just write out as it is (for the moment)
save(EOF.pattern_proj,file = "./data/EOF.pattern_RCPproj_POLCOM.ERSEM_bottomT_sst_salt_currents_NS.RData")
  
}


