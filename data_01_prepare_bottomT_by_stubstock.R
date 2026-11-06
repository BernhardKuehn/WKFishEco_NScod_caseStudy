# -------------------------------------------------------- #
# Script to: 
# - calculate yearly avg. of bottomT over the 3 substock areas
# -------------------------------------------------------- #

rm(list = ls())

library(raster)
library(terra)
library(sf)
library(ggplot2)
library(readxl)
library(collapse)
library(marmalaid)

# set seed
set.seed(42)

# --------------------------------- #
# get the substock area definitions
# --------------------------------- # 

# load ICES rectangles
sf::read_sf("./data/ICES_stat_rect/StatRec_map_Areas_Full_20170124.shp") -> ICES_rects
#ICES_rects = as.data.frame(ICES_rects)
ICES_rects$area = gsub("\\.","",ICES_rects$Area_27)

# get 6a, 7d, 4 and 3
unique(ICES_rects$area)
ICES_rects_NS = ICES_rects[grep("^(4a|4b|4c|7d|6a|3a20)$",ICES_rects$area),]

land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
bb <- sf::st_bbox(ICES_rects_NS)

# show on a curvature grid
rects <- sf::st_transform(ICES_rects_NS, 4326)
land  <- sf::st_transform(land, 4326)
land  <- sf::st_crop(land, sf::st_bbox(rects))

rotated_crs <- sf::st_crs(
  "+proj=omerc +lat_0=55 +lonc=-5 +alpha=0 +gamma=5 +k=1 +datum=WGS84 +units=m +no_defs"
)

# plot
ggplot() +
  geom_sf(data = rects, aes(fill = area),
          color = "white", linewidth = 0.2) +
  geom_sf(data = land, fill = "grey85", color = "grey60",
          linewidth = 0.15) +
  coord_sf(crs = rotated_crs) +
  geom_sf_text(data = rects,aes(x = stat_x,y = stat_y,label = ICESNAME),size = 2,color = "white") + 
  theme_minimal() +
  labs(fill = "ICES area")

# list those that are specific to the substocks

stk_areas = read_excel("./data/stock/Area_NSCOD3_table_6a_inshore_offshore.xlsx")
head(stk_areas)

rects_substock = merge(rects,stk_areas[,c("Statistical Rectangle","areaPopulation")],by.x = c("ICESNAME"),by.y = "Statistical Rectangle")

ggplot() +
  geom_sf(data = rects_substock, aes(fill = areaPopulation),
          color = "white", linewidth = 0.2) +
  geom_sf(data = land, fill = "grey85", color = "grey60",
          linewidth = 0.15) +
  coord_sf(crs = rotated_crs) +
  geom_sf_text(data = rects,aes(x = stat_x,y = stat_y,label = ICESNAME),size = 2,color = "white") + 
  theme_minimal() +
  labs(fill = "ICES area")

# merge Northwestern parts
rects_substock$substock = ifelse(grepl("Northwestern",rects_substock$areaPopulation),"Northwestern",rects_substock$areaPopulation)
# -------------------------------------------------------------------------------------------------- #

# ----------------------------------------- # 
# load bottomT data and cut with substock areas
# ----------------------------------------- #

# load
bottomT = rast("./data/Env/processed/AHOI_BottomTemp_raster_1945_2024.tif")

# extract time
time_bottomT = as.Date(time(bottomT))

# convert subareas to vector in terra package
areas_v <- terra::project(
  terra::vect(rects_substock[,"substock",drop = F]),
  terra::crs(bottomT)
)

# unique substocks
substock_regions <- unique(areas_v$substock)

# crop by substock
bottomT_by_substock <- setNames(lapply(substock_regions, function(r) {
  region_polygon <- areas_v[areas_v$substock == r, ,drop =F]
  terra::crop(bottomT, region_polygon, mask = TRUE)
}), substock_regions)


par(mfrow = c(1,3))
image(bottomT_by_substock$Southern,main = "Southern")
maps::map(add =T)
image(bottomT_by_substock$Northwestern,main = "Northwestern")
maps::map(add =T)
image(bottomT_by_substock$Viking,main = "Viking")
maps::map(add =T)

# -------------------------------------- #
# calculate yearly avg. over the regions
# -------------------------------------- #

# avg. over the year
bottomT_yr_by_substock = lapply(bottomT_by_substock,function(x){
    calc.mean.over.Month(x,time = time_bottomT,month = 1:12,shiftYear = F)
})

# avg. over the subareas
bottomT_yr_by_substock_ts = lapply(1:length(bottomT_yr_by_substock),function(i){
  x = bottomT_yr_by_substock[[i]]
  time = as.numeric(names(x))
  tmp = brick(x)
  ts = cellStats(tmp,mean)
  df = data.frame(year = time,ts = ts)
  names(df)[2] = paste("bottomT",names(bottomT_yr_by_substock)[i],sep = "_")
  return(df)
})

# join
yr = data.frame(year = 1945:2024)
bottomT_substk = mult.merge("yr","bottomT_yr_by_substock_ts",by = "year")

# check difference
par(mfrow = c(2,1))
plot(bottomT_substk$year,bottomT_substk$bottomT_Southern,type = "l",col = "aquamarine4",ylim = c(8,14))
lines(bottomT_substk$year,bottomT_substk$bottomT_Northwestern,col = "red")
lines(bottomT_substk$year,bottomT_substk$bottomT_Viking,col = "blue")

plot(bottomT_substk$year,scale(bottomT_substk$bottomT_Southern),type = "l",col = "aquamarine4")
lines(bottomT_substk$year,scale(bottomT_substk$bottomT_Northwestern),col = "red")
lines(bottomT_substk$year,scale(bottomT_substk$bottomT_Viking),col = "blue")

cor(bottomT_substk)
# differences are marginal, correlations all above ~ 0.8

# # ---------------------------------------- #
# # alternative: 
# # calculate main pattern per area via EOTs
# # ---------------------------------------- #
# 
# library(fast.EOT)
# 
# EOTs_bottomT_yr_by_substock = lapply(bottomT_yr_by_substock,function(x)
#                                  fast.EOT::fast.eot(x = brick(x),n = 10,
#                                    standardised = T,
#                                    n_cores = 8,verbose = T))
# 
# # check the different EOTs per subarea
# lapply(EOTs_bottomT_yr_by_substock,function(x) plot_eot(x,n = 10))
# 
# # extract temporal pattern
# 
# eots_ts_out = lapply(names(EOTs_bottomT_yr_by_substock),
#        function(x){
#          eot_out = EOTs_bottomT_yr_by_substock[[x]]
#          eot_ts = as.data.frame(sapply(eot_out, function(x) x@eot))
#          names(eot_ts) = paste("bottomT",x,names(eot_ts),sep = "_")
#          # add yrs 
#          eot_ts = data.frame(year = as.numeric(names(bottomT_yr_by_substock[[x]])),
#                              eot_ts)
#          return(eot_ts)
#        })
# 
# # par(mfrow = c(1,1))
# # plot(scale(eot_ts[,1]),type = "l")
# # for(i in 2:10)
# # lines(scale(eot_ts[,i]))
# 
# # join
# eots_bottomT_substk = mult.merge("yr","eots_ts_out",by = "year")

# --------------------------------------------- #
# perform a clustering of similar ocean pattern
# --------------------------------------------- #

source("functions/spatial_clustering.R")
set.seed(42)
bottomT_clusters <- lapply(
  bottomT_yr_by_substock,
  spatial_kmeans,
  k_values = 2:10
)

# Selected k for each substock
sapply(bottomT_clusters, `[[`, "selected_k")

# Example selected map and area-weighted bottomT covariates
par(mfrow = c(1,3))
for(i in 1:length(bottomT_clusters)){
plot(bottomT_clusters[[i]]$selected$cluster_map)
maps::map(add = TRUE)
title(names(bottomT_clusters)[i],adj =0)
}

cluster_ts_out = lapply(1:length(bottomT_clusters),function(i){
                        x = bottomT_clusters[[i]]
                        substk = names(bottomT_clusters)[i]
                        tmp = as.data.frame(x$selected$timeseries)
                        tmp$time = as.numeric(gsub("X","",tmp$time))
                        names(tmp)[names(tmp) == "time"] <- "year"
                        #tmp$substock = substk
                        names(tmp)[grep("Cluster",names(tmp))] = paste("bottomT",substk,
                                                                       names(tmp)[grep("Cluster",names(tmp))],sep = "_") 
                        return(tmp)
                        })

# join
cluster_bottomT_substk = mult.merge("yr","cluster_ts_out",by = "year")

cluster_bottomT_substk_spatial = lapply(bottomT_clusters, function(x) {
  map = x$selected$cluster_map
  if (!raster::inMemory(map)) map = raster::readAll(map)
  map
})

# ------------ #
# write out ts
# ------------ #

save(object = bottomT_substk,file = file.path("./data/Env/processed","bottomT_yr_by_cod_substock.rda"))
#save(object = eots_bottomT_substk,file = file.path("./data/Env/processed","EOT_bottomT_yr_by_cod_substock.rda"))
save(object = cluster_bottomT_substk,file = file.path("./data/Env/processed","Cluster_bottomT_yr_by_cod_substock.rda"))

# -------------------------- # 
# write out the spatial data
# -------------------------- #

save(object = cluster_bottomT_substk_spatial,
     file = file.path("./data/Env/processed","SpatialCluster_bottomT_yr_by_cod_substock.rda"))

# SpatRaster external pointers do not survive save()/load(); store the
# annual fields as in-memory RasterBricks so the spatial product is portable.
bottomT_yr_by_substock_spatial = lapply(bottomT_yr_by_substock, function(x) {
  x = raster::brick(x)
  if (!raster::inMemory(x)) x = raster::readAll(x)
  x
})
save(object = bottomT_yr_by_substock_spatial,
     file = file.path("./data/Env/processed","Spatial_bottomT_yr_by_cod_substock.rda"))
