# ------------------------------------------------------------------------- #
# Script to bring env. data together with the stock-recruitment time series
# ------------------------------------------------------------------------- #

rm(list = ls())

library(readxl)
library(FLCore)
library(ggplot2)
library(stringr)


# --------- #
# load data
# --------- #

# load SST
load("./data/Env/processed/SST_yr_by_cod_substock.rda")
# # load SST EOTs
# load("./data/Env/processed/EOT_SST_yr_by_cod_substock.rda")
# load clustering results
load("./data/Env/processed/Cluster_SST_yr_by_cod_substock.rda")

# load stock data
tmp.nm = load("./data/stock/cod46a7d20_FLstockObjects_assessment_2026.Rdata")
cod = get(tmp.nm)

# ------------------- #
# extract rec and ssb
# ------------------- #

get_stk_data = function(stk_obj){
  substk =  sub(".*\\b(\\w+)\\s+sub-stock\\b.*", "\\1", stk_obj@name, perl = TRUE)
  r = rec(stk_obj)
  ssb = ssb(stk_obj)
  yr = seq(stk_obj@range["minyear"],stk_obj@range["maxyear"])
  lag = as.numeric(dimnames(r)$age)
  df = data.frame(year = yr,rec = as.numeric(r),
                  ssb_lag = dplyr::lag(as.numeric(ssb),n = lag))
  names(df)[grepl("ssb",names(df))] = paste0("ssb_","lag",lag)
  df$substock = substk
  return(df)
}

cod_substks_lst = lapply(cod,get_stk_data)
# join
cod_substks = do.call(rbind,cod_substks_lst)

# ggplot(cod_substks,aes(x = ssb_lag1,y = rec,color = substock)) + 
#   geom_point() + 
#   geom_smooth() + 
#   facet_wrap(~substock) + 
#   theme_bw()

# -------------------------------- #
# join with the environmental data
# -------------------------------- #

# rename SST mean
names(sst_substk)[grep("SST",names(sst_substk))] <- 
  gsub("SST","SST_mean",names(sst_substk)[grep("SST",names(sst_substk))])

# # join mean and EOTs
# eot_mean_sst_substk = merge(eots_sst_substk,sst_substk,by = "year")

# join mean and clustering
clustering_mean_sst_substk = merge(cluster_sst_substk,sst_substk,by = "year")


# # lag SST
# eot_mean_sst_substk_lag = eot_mean_sst_substk
# eot_mean_sst_substk_lag[,-1] = apply(eot_mean_sst_substk[,-1],2,function(x) dplyr::lag(x,1))
# # bring to long format
# eot_mean_sst_substk_long = pivot(eot_mean_sst_substk_lag,ids = "year")
# eot_mean_sst_substk_long$substock = str_extract(pattern = "Northwestern|Southern|Viking",
#                                                 string = eot_mean_sst_substk_long$variable)
# 
# unique(cod_substks$substock)
# unique(eot_mean_sst_substk_long$substock)
# cod_substks$substock = ifelse(grepl("Northwest",cod_substks$substock),"Northwestern",
#                               ifelse(grepl("South",cod_substks$substock),"Southern",cod_substks$substock))
# env_stk = merge(cod_substks,eot_mean_sst_substk_long[,c("year","variable","value","substock")],by = c("year","substock"))
# head(env_stk)

# lag SST
clustering_mean_sst_substk_lag = clustering_mean_sst_substk
clustering_mean_sst_substk_lag[,-1] = apply(clustering_mean_sst_substk[,-1],2,function(x) dplyr::lag(x,1))
# bring to long format
clustering_mean_sst_substk_long = pivot(clustering_mean_sst_substk_lag,ids = "year")
clustering_mean_sst_substk_long$substock = str_extract(pattern = "Northwestern|Southern|Viking",
                                                string = clustering_mean_sst_substk_long$variable)

unique(cod_substks$substock)
unique(clustering_mean_sst_substk_long$substock)
cod_substks$substock = ifelse(grepl("Northwest",cod_substks$substock),"Northwestern",
                              ifelse(grepl("South",cod_substks$substock),"Southern",cod_substks$substock))
env_stk = merge(cod_substks,clustering_mean_sst_substk_long[,c("year","variable","value","substock")],by = c("year","substock"))
head(env_stk)



# check relationship between recruitment and SST_lag1
ggplot(env_stk,aes(x = value,y = rec,color = substock,group = substock)) + 
  geom_point() + 
  geom_smooth() + 
  facet_wrap(~variable,scale = "free_x")


# --------- #
# write out 
# --------- #

save(env_stk,file = file.path("./data","Cod_Env_data.rda"))


