# ========================================= #
# Script to look at and aggregate AHOI data
# historical data 1948 - 2024
# monthly
# updated AHOI runs (v2.3.05) from 2026
# last change: 02-11-26
# author: Bernhard Kuehn
# ========================================= #

rm(list = ls())

library(raster)
library(RNetCDF)
library(marmalaid)
library(stringr)
library(terra)
library(DescTools)
library(foreach)
library(parallel)
library(doParallel)

data_dir = "./data/Env/AHOI/v23.05/"

nc.files = list.files(data_dir,pattern = "*.nc")

# check one nc-file
nc = open.nc(file.path(data_dir,nc.files[1]))
print.nc(nc)

# get lat long
lat = as.numeric(unique(var.get.nc(nc,"Latitude")))
lon = unique(as.numeric(var.get.nc(nc,"Longitude")))
res = c(1/unique(round(1/diff(lon),digits = 1)),1/unique(round(1/diff(lat),digits = 1)))
extent.obj = extent(c(range(lon),range(lat)))
extent.obj = extent.obj + c(-res[1]/2,res[1]/2,-res[2]/2,res[2]/2)

# check out depth
depth = var.get.nc(nc,"Depth")
length(depth)
# 87 depth levels

# extract temp, bottom temp, salinity
var.names = c("Salinity","Temperature","Temperature_Bottom","Salinity_Bottom")
var.names_long = c("SurfaceSalt","SST","BottomTemp","BottomSalt")
lookup.tb = data.frame(var.name = var.names,
                       var.names_long = var.names_long,
                       depth.to.extract = c(10,10,NA,NA))

# function to rotate matrix
rotate.mat = function(mat, direction = "left",times = NULL) {
  if(!is.null(times)){
    for(i in 1:times){
      mat = t(apply(mat,2,rev))
    }
    return(mat)
  }
  if(direction == "left"){
    return(t(apply(mat,2,rev)))  
  } else if(direction == "right"){
    for(i in 1:3){
      mat = rotate.mat(mat,"left")
    }
    return(mat)
  }
}

get.bottom = function(ar,depth.dim){
  # how many depth levels
  dlvl = dim(ar)[depth.dim]
  diml = length(dim(ar))
  time.dim = dim(ar)[-depth.dim][3]
  td = which(dim(ar) == time.dim)
  ll = list()
  for(n in 1:time.dim){
    if(td == 4){
      # get smallest index which is not NA
      index.mat = apply(ar[,,,n],1:2,function(x) sum(!is.na(x)))
    } else if(td == 3){
      index.mat = apply(ar[,,n,],1:2,function(x) sum(!is.na(x)))
    } else{
      stop("Something wrong with the time and depth dimension!")
    }
    index.mat[index.mat == 0] <- NA
    mat = matrix(NA,nrow(ar),ncol(ar))
    for(i in 1:ncol(index.mat)){
      for(j in 1:nrow(index.mat)){
        if(td == 4){
          mat[j,i] = ar[j,i,index.mat[j,i],n]
        } else if(td == 3){
          mat[j,i] = ar[j,i,n,index.mat[j,i]]
        }
      }
    }
    ll[[n]] = mat
    svMisc::progress(n)
  }
  bottom.ar = Abind(ll,along = 3)
  return(bottom.ar)
}

# set values where landmasses are to NA
library(maptools)
library(maps)
library(rgeos)
map.world = map(database = "world",fill =T)
IDs = sapply(strsplit(map.world$names, ":"), function(x) x[1])
land.mask = map2SpatialPolygons(map.world, IDs=IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))

CP <- as(extent(-20, 20, 40, 70), "SpatialPolygons")
proj4string(CP) = CRS(proj4string(land.mask))
land.mask.small = gIntersection(land.mask,CP,byid = T)

# load data
for(n in seq(var.names)){
  
  var = lookup.tb$var.name[n]
  var_long = lookup.tb$var.names_long[n]
  depth.to.avg = lookup.tb$depth.to.extract[n]
  
  print(var)
  # get the corresponding nc-file
  var.files = nc.files
  
  # parallelize
  # register parallel workers
  cl = makeCluster(10)
  registerDoParallel(cl)
  
  var_list =  foreach(i = seq(nc.files),
                      .combine = list,.multicombine = T,.inorder = T,.packages = c("RNetCDF","DescTools","stringr")) %dopar% {
                        
    cat("Opening...","\n")
    nc = open.nc(file.path(data_dir,var.files[i]))
    some.var = var.get.nc(nc,variable = ifelse(var == "s","so",
                                               ifelse(var == "t","thetao",var)))
    dim(some.var)
    cat("Calc Avg. over depth...","\n")
    
    if(!is.na(depth.to.avg)){
      # calc averages over respective depth
      some.var_avg = apply(some.var[,,which(depth<=depth.to.avg)],c(1,2),mean,na.rm = T)
    } else if(is.na(depth.to.avg) & grepl("Bottom",var_long) & length(dim(some.var)) >= 3){
      some.var_avg = get.bottom(some.var,depth.dim = 3)
    } else {
      some.var_avg = some.var
    }
    dim(some.var_avg)
    # get time
    t_extr = str_extract(var.files[i],"[0-9]+")
    time.var = as.Date(paste(substr(t_extr,1,4),substr(t_extr,5,6),"15",sep = "-"))
    
    close.nc(nc)
    # store
    output = list(data = some.var_avg,
                  time = time.var)
    output
    
    
    
  }
  stopCluster(cl)
  # get time ordering right
  # extract data
  var_list = flatten.list(var_list)
  data_list = var_list[which(names(var_list) %in% "data")]
  time.nc_all = as.Date(do.call(c,var_list[which(names(var_list) %in% "time")]),origin = "1970-01-01")
  
  # check if the times align
  order.time= order(time.nc_all)
  time.nc_all = time.nc_all[order.time]
  plot(time.nc_all)
  lines(c(1,length(time.nc_all)),
        c(min(time.nc_all),max(time.nc_all)),col = "red",lwd = 1.5)
  
  
  # bring together
  some.var_array = Abind(data_list[order.time],along = 3)
  dim(some.var_array)
  
  
  # rotate to right direction
  some.var_array = Abind(lapply(asplit(some.var_array,3),function(x) rotate.mat(x,"left",times = 3)),along = 3)
  
  # some.var_rst = brick(some.var_array,
  #                      xmn = extent.obj[1],
  #                      xmx = extent.obj[2],
  #                      ymn = extent.obj[3], 
  #                      ymx = extent.obj[4])
  # image(some.var_rst)
  # maps::map(add = T)
  
  # create new "terra"-raster
  some.var_rst = rast(some.var_array,extent = extent.obj,crs = "EPSG:4326")
  time(some.var_rst) <- as.POSIXlt(time.nc_all)
  
  image(some.var_rst)
  maps::map(add = T)
  names(some.var_rst) = time.nc_all
  
  # store as tiff
  dir.create("./data/Env/processed",showWarnings = F)
  terra::writeRaster(some.var_rst,filename = file.path("./data/Env/processed",paste0(paste("AHOI",
                                                                                var_long,
                                                                                "raster",
                                                                                paste(format(range(time.nc_all),"%Y"),collapse = "_"),sep = "_"),".tif")),
                     overwrite = T)
  # # # don`t cut North Sea
  # # some.var_rst = crop(some.var_rst,extent(c(-5,12,49,62)))
  # image(some.var_rst)
  # maps::map(add = T)
  # names(some.var_rst) = time.var_all
  # 
  # # cut with land-sea mask
  # 
  # some.var_rst = mask(x = some.var_rst,mask = land.mask.small,inverse = T)
  # 
  # some.var_rst = readAll(some.var_rst)
  # assign(paste("AHOI",var_long,"raster",sep = "_"),some.var_rst)
  # 
  # # store
  # save(list = paste("AHOI",var_long,"raster",sep = "_"),
  #      file = file.path("./output",
  #                       paste0(paste("AHOI",var_long,"raster",paste(format(range(time.var_all),"%Y"),collapse = "_"),sep = "_"),".RData")))
  
}
# ------------------------------------------------------------- #