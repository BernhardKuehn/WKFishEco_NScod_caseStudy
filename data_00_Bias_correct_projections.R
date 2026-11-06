# ============================== #
# Bias correction of Temperature
# ============================== #

rm(list = ls())

library(raster)
library(marmalaid)
library(sinkr)
library(reshape2)

# ---------------- #
# helper functions
# ---------------- #

#helper function
min.max.scaler = function(x){
  scale(x,min(x),max(x)-min(x))
}
# remove linear trend
trend.rm = function(x){
  lx = 1:length(x)
  mod = lm(x~ lx)
  coef = mod$coefficients
  #x - predict(lm(x~lx))
  x - lx*coef[[2]]
}
trend.add = function(x.rm,x){
  # get trend from x
  lx = 1:length(x)
  mod = lm(x~ lx)
  #x.rm + predict(mod,newdata = data.frame(1:length(x.rm)))
  coef = mod$coefficients
  x.rm + lx*coef[[2]]
}
# --------------------------------------------------------------- #

# load data
AHOI_BottomTemp_raster = brick("./data/Env/processed/AHOI_BottomTemp_raster_1945_2024.tif")
load("./data/Env/POLCOMS-ERSEM/POLCOM_ERSEM_GreaterNorthSea_rcp45_BottomTemp_raster_2006_2099.RData")
load("./data/Env/POLCOMS-ERSEM/POLCOM_ERSEM_GreaterNorthSea_rcp85_BottomTemp_raster_2006_2099.RData")

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

# combine these to one polygon 
NA.mask = union(land.mask.small)
crs(NA.mask) = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" 

# mask with polygons
POLCOM_ERSEM_rcp45_BottomTemp_raster = mask(POLCOM_ERSEM_GreaterNorthSea_rcp45_BottomTemp_raster,mask = NA.mask,inverse = T,updatevalue = NA)
POLCOM_ERSEM_rcp85_BottomTemp_raster = mask(POLCOM_ERSEM_GreaterNorthSea_rcp85_BottomTemp_raster,mask = NA.mask,inverse = T,updatevalue = NA)
AHOI_BottomTemp_raster = mask(AHOI_BottomTemp_raster,mask = NA.mask,inverse = T,updatevalue = NA)

par(mfrow = c(1,3),mar = c(2,2,2,2))
image(AHOI_BottomTemp_raster)
maps::map(add = T)
image(POLCOM_ERSEM_rcp45_BottomTemp_raster)
maps::map(add = T)
image(POLCOM_ERSEM_rcp45_BottomTemp_raster)
maps::map(add = T)

# resample to same resolution (0.2x0.2)
POLCOM_ERSEM_rcp45_BottomTemp_raster = resample(POLCOM_ERSEM_rcp45_BottomTemp_raster,AHOI_BottomTemp_raster)
POLCOM_ERSEM_rcp85_BottomTemp_raster = resample(POLCOM_ERSEM_rcp85_BottomTemp_raster,AHOI_BottomTemp_raster)


# ==================================== #
# set equal number of landmasses to NA
# ==================================== #

# check the landmasses in the ensemble runs and the AHOI-data
BottomTemp.polcom.rasters = list(rcp45 = POLCOM_ERSEM_rcp45_BottomTemp_raster,
                          rcp85 = POLCOM_ERSEM_rcp85_BottomTemp_raster)

par(mfrow = c(1,2))
for(i in 1:length(BottomTemp.polcom.rasters)){
  test1 = BottomTemp.polcom.rasters[[i]]
  test1[!is.na(test1)] = 1
  test2 = AHOI_BottomTemp_raster
  test2[!is.na(test2)] = 2
  
  rgb.col = col2rgb("red",alpha = T)
  col.transparent = rgb(rgb.col[[1]],rgb.col[[2]],rgb.col[[3]],alpha = 0.5*255,maxColorValue = 255)
  image(test2[[1]],col = "olivedrab")
  image(test1[[1]],col = col.transparent,add = T)
  maps::map(add = T)
  legend("bottomright",fill = c("olivedrab",col.transparent),legend = c("AHOI","POLCOM"),bty = "n")
  title(names(BottomTemp.polcom.rasters)[i],adj = 0)
}
# not exact overlap!!!
# take only grid-cells inot account where all datasets have data!!!

# calculate overlay
overlay.AHOI.POLCOM = overlay(stack(AHOI_BottomTemp_raster[[1]],brick(lapply(BottomTemp.polcom.rasters,function(x) x[[1]]))),
                              fun = sum)
overlay.AHOI.POLCOM[!is.na(overlay.AHOI.POLCOM)] = 1
par(mfrow = c(1,1))
image(overlay.AHOI.POLCOM)
maps::map(add = T)

BottomTemp_polcom_crop = lapply(BottomTemp.polcom.rasters,function(x) mask(x,mask = overlay.AHOI.POLCOM))
AHOI.BottomTemp.crop = mask(AHOI_BottomTemp_raster,mask = overlay.AHOI.POLCOM)


# ================== #
# kill seasonal bias
# ================== #

# perform Bias correction for overlap period
time.ahoi = as.Date(sub("X","",names(AHOI.BottomTemp.crop)),"%Y.%m.%d")
time.polcom = as.Date(sub("X","",names(BottomTemp_polcom_crop$rcp45)),"%Y.%m.%d")
# overlap
time.ahoi_overlap = time.ahoi[format(time.ahoi,"%Y%m") %in% format(time.polcom,"%Y%m")]
time.polcom_overlap = time.polcom[format(time.polcom,"%Y%m") %in% format(time.ahoi,"%Y%m")]


# mean of AHOI 2006:2024 (recent overlap)
indx.AHOI = which(time.ahoi %in% time.ahoi_overlap)
BottomTemp.AHOI_overlap = AHOI.BottomTemp.crop[[indx.AHOI]]

BottomTemp_polcom.bias.rm = BottomTemp_polcom_crop 

# for RCP4.5 and RCP8.5
for(j in seq(BottomTemp_polcom_crop)){
  
  Bias.month = all.month.indx = list()
  for(i in 1:12){
    cat("Month:",month.abb[[i]],"\n")
    indx.month.rcp4.5_8.5 = which(as.numeric(format(time.polcom_overlap,"%m")) == i)
    indx.month.ahoi = which(as.numeric(format(time.ahoi_overlap,"%m")) == i)
    # store all indexes for later
    all.month.indx[[i]] = indx.month.rcp4.5_8.5
    
    # calculate mean
    Bias.month[[i]] = mean(BottomTemp_polcom_crop[[j]][[indx.month.rcp4.5_8.5]] - BottomTemp.AHOI_overlap[[indx.month.ahoi]])
    
    # remove monthly Bias for projection run
    index.rcp4.5_8.5.month.proj = which(as.numeric(format(time.polcom,"%m")) == i) 
    for(jj in index.rcp4.5_8.5.month.proj){
      BottomTemp_polcom.bias.rm[[j]][[jj]] = BottomTemp_polcom_crop[[j]][[jj]] - Bias.month[[i]]
    }
  }
}


# --------------------- #
# check bias correction
# --------------------- #

par(mfrow = c(2,1),mar = c(2,2,3,2))
# original data
plot(time.ahoi,cellStats(AHOI.BottomTemp.crop,mean),
     type = "l",xlim = c(as.Date("1960-01-01"),as.Date("2100-01-01")),col = "black",xlab = "time",ylab = "BottomTemp")
# projections
for(j in 1:2){
  lines(time.polcom,
        cellStats(BottomTemp_polcom_crop[[j]],mean),col = c("olivedrab","red")[j])
}

# bias corrected time series
plot(time.ahoi,cellStats(AHOI.BottomTemp.crop,mean),
     type = "l",xlim = c(as.Date("1960-01-01"),as.Date("2100-01-01")),col = "black",xlab = "time",ylab = "BottomTemp")
# projections
for(j in 1:2){
  lines(time.polcom,
        cellStats(BottomTemp_polcom.bias.rm[[j]],mean),col = c("olivedrab","red")[j])
}
# looks ok!                     

# ------------------------- #
# calculate seasonal means
# ------------------------- 
# ================================ #
# calculate mean over month for seasons
# DJF - Q1
# MAM - Q2
# JJA - Q3
# SON - Q4
# ============================================================== #

library(marmalaid)

month.codes = list(DJF = c(12,1,2),
                   MAM = 3:5,
                   JJA = 6:8,
                   SON = 9:11,
                   YR = 1:12)

BottomTemp_polcom.bias.rm.month = BottomTemp_polcom.month = list()
AHOI.month = list()
for(j in seq(BottomTemp_polcom.bias.rm)){
  
  tmp.proj = BottomTemp_polcom.bias.rm[[j]]
  BottomTemp_polcom.bias.rm.month[[j]] = BottomTemp_polcom.month[[j]] = list()
  
  for(i in 1:length(month.codes)) {
    month = month.codes[[i]]
    # for BottomTemp
    BottomTemp_polcom.bias.rm.month[[j]][[i]] = calc.mean.over.Month(raster = tmp.proj,time = time.polcom,month = month,shiftYear = F)
    BottomTemp_polcom.month[[j]][[i]] = calc.mean.over.Month(raster = BottomTemp_polcom_crop[[j]],time = time.polcom,month = month,shiftYear = F)
    if(j == 1){
      AHOI.month[[i]] = calc.mean.over.Month(raster = AHOI.BottomTemp.crop,time = time.ahoi,month = month,shiftYear = F)
    }
  }
  # rename
  names(BottomTemp_polcom.bias.rm.month[[j]]) = names(BottomTemp_polcom.month[[j]]) = names(month.codes) 
}
names(AHOI.month) = names(month.codes)
names(BottomTemp_polcom.bias.rm.month) = names(BottomTemp_polcom.month) = names(BottomTemp_polcom.bias.rm)

# ------------------- #
# plot bias corrected
# ------------------- #

par(mfrow = c(2,3),oma = c(0,0,2,0))
for(i in names(month.codes)){
  plot(unique(as.numeric(format(time.ahoi,"%Y"))),
       cellStats(AHOI.month[[i]],mean),type = "l",
       col ="black",xlab = "time",ylab = "BottomTemp [?C]",las = 1,xlim = c(1948,2099),ylim = range(cellStats(AHOI.month[[i]],mean))+c(-1,2))
  for(nn in 1:2){  
    #lines(unique(as.numeric(format(BottomTemp_rcp4.5_8.5.ensembles.bias.rm$historical[[1]]$time,"%Y"))),cellStats(tmp.hist[[i]],mean),type = "l",col ="black")
    lines(unique(as.numeric(format(time.polcom,"%Y"))),cellStats(BottomTemp_polcom.bias.rm.month[[nn]][[i]],mean),type = "l",col = c("olivedrab","red")[nn])
    title(i, adj = 0)
  }
  mtext("Bias removed",3,line = 0,outer = T,cex = 1.2,font = 2)
}
# ------------------ #
# not bias corrected
# ------------------ #

par(mfrow = c(2,3))
for(i in names(month.codes)){
  plot(unique(as.numeric(format(time.ahoi,"%Y"))),
       cellStats(AHOI.month[[i]],mean),type = "l",
       col ="black",xlab = "time",ylab = "BottomTemp [?C]",las = 1,xlim = c(1948,2099),ylim = range(cellStats(AHOI.month[[i]],mean))+c(-1,2))
  for(nn in 1:2){  
    #lines(unique(as.numeric(format(BottomTemp_rcp4.5_8.5.ensembles.bias.rm$historical[[1]]$time,"%Y"))),cellStats(tmp.hist[[i]],mean),type = "l",col ="black")
    lines(unique(as.numeric(format(time.polcom,"%Y"))),cellStats(BottomTemp_polcom.month[[nn]][[i]],mean),type = "l",col = c("olivedrab","red")[nn])
    title(i, adj = 0)
  }
  mtext("Raw data",3,line = 0,outer = T,cex = 1.2,font = 2)
}

# --------- #
# write out #
# --------- #

# write times as names in rasterfiles
names(BottomTemp_polcom.bias.rm$rcp45) = names(BottomTemp_polcom.bias.rm$rcp85) = time.polcom

POLCOM_ERSEM_BottomTemp_raster_meanBiasRM = BottomTemp_polcom.bias.rm

POLCOM_ERSEM_BottomTemp_raster_meanBiasRM = lapply(POLCOM_ERSEM_BottomTemp_raster_meanBiasRM,readAll)

dir.create("./data/Env/POLCOMS-ERSEM_biasRM/",showWarnings = F)
save(POLCOM_ERSEM_BottomTemp_raster_meanBiasRM,file = "./data/Env/POLCOMS-ERSEM_biasRM/POLCOM_ERSEM_BottomTemp_raster_meanBiasRM.RData")


