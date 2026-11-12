# -------------------------------- #
# plot the clusters of temperature
# -------------------------------- #

rm(list = ls())

library(terra)
library(raster)
library(wesanderson)

load(file.path("./data/Env/processed","SpatialCluster_bottomT_yr_by_cod_substock.rda"))

pal = c(wesanderson::wes_palette("GrandBudapest1",n = 4),wesanderson::wes_palette("Darjeeling1",n = 5))

par(mfrow = c(1,3))
plot(cluster_bottomT_substk_spatial$Northwestern,main = "Northwest",col = pal)
maps::map(add = T)
plot(cluster_bottomT_substk_spatial$Southern,main = "Southern",col = pal)
maps::map(add = T)
plot(cluster_bottomT_substk_spatial$Viking,main = "Viking",col = pal)
maps::map(add = T)

