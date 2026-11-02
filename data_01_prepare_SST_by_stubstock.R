# -------------------------------------------------------- #
# Script to: 
# - calculate yearly avg. of SST over the 3 substock areas
# -------------------------------------------------------- #

library(raster)
library(terra)
library(sf)
library(ggplot2)
library(readxl)
library(collapse)
library(marmalaid)


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
# load SST data and cut with substock areas
# ----------------------------------------- #

# load
sst = rast("./data/Env/processed/AHOI_SST_raster_1945_2024.tif")

# extract time
time_sst = as.Date(time(sst))

# convert subareas to vector in terra package
areas_v <- terra::project(
  terra::vect(rects_substock[,"substock",drop = F]),
  terra::crs(sst)
)

# unique substocks
substock_regions <- unique(areas_v$substock)

# crop by substock
sst_by_substock <- setNames(lapply(substock_regions, function(r) {
  region_polygon <- areas_v[areas_v$substock == r, ,drop =F]
  terra::crop(sst, region_polygon, mask = TRUE)
}), substock_regions)


par(mfrow = c(1,3))
image(sst_by_substock$Southern,main = "Southern")
maps::map(add =T)
image(sst_by_substock$Northwestern,main = "Northwestern")
maps::map(add =T)
image(sst_by_substock$Viking,main = "Viking")
maps::map(add =T)

# -------------------------------------- #
# calculate yearly avg. over the regions
# -------------------------------------- #

# avg. over the year
sst_yr_by_substock = lapply(sst_by_substock,function(x){
    calc.mean.over.Month(x,time = time_sst,month = 1:12,shiftYear = F)
})

# avg. over the subareas
sst_yr_by_substock_ts = lapply(1:length(sst_yr_by_substock),function(i){
  x = sst_yr_by_substock[[i]]
  time = as.numeric(names(x))
  tmp = brick(x)
  ts = cellStats(tmp,mean)
  df = data.frame(year = time,ts = ts)
  names(df)[2] = paste("SST",names(sst_yr_by_substock)[i],sep = "_")
  return(df)
})

# join
yr = data.frame(year = 1945:2024)
sst_substk = mult.merge("yr","sst_yr_by_substock_ts",by = "year")

# check difference
par(mfrow = c(2,1))
plot(sst_substk$year,sst_substk$SST_Southern,type = "l",col = "aquamarine4",ylim = c(8,14))
lines(sst_substk$year,sst_substk$SST_Northwestern,col = "red")
lines(sst_substk$year,sst_substk$SST_Viking,col = "blue")

plot(sst_substk$year,scale(sst_substk$SST_Southern),type = "l",col = "aquamarine4")
lines(sst_substk$year,scale(sst_substk$SST_Northwestern),col = "red")
lines(sst_substk$year,scale(sst_substk$SST_Viking),col = "blue")

cor(sst_substk)
# differences are marginal, correlations all above 0.77

# --------- #
# write out
# --------- #

save(object = sst_substk,file = file.path("./data/Env/processed","SST_yr_by_cod_substock.rda"))










