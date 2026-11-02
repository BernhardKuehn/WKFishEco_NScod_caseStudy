# ------------------------------------------------------------------------- #
# Script to bring env. data together with the stock-recruitment time series
# ------------------------------------------------------------------------- #

rm(list = ls())

library(readxl)
library(FLCore)
library(ggplot2)

# --------- #
# load data
# --------- #

# load SST
load("./data/Env/processed/SST_yr_by_cod_substock.rda")

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

# lag SST
sst_substk_lag = sst_substk
sst_substk_lag[,-1] = apply(sst_substk[,-1],2,function(x) dplyr::lag(x,1))
# bring to long format
sst_substk_long = pivot(sst_substk_lag,ids = "year",names = list("variable","SST_lag1"))
sst_substk_long$substock = gsub(".*_","",sst_substk_long$variable)

unique(cod_substks$substock)
unique(sst_substk_long$substock)
cod_substks$substock = ifelse(grepl("Northwest",cod_substks$substock),"Northwestern",
                              ifelse(grepl("South",cod_substks$substock),"Southern",cod_substks$substock))
env_stk = merge(cod_substks,sst_substk_long[,c("year","SST_lag1","substock")],by = c("year","substock"))
head(env_stk)


# check relationship between recruitment and SST_lag1
ggplot(env_stk,aes(x = SST_lag1,y = rec,color = substock,group = substock)) + 
  geom_point() + 
  geom_smooth()


# --------- #
# write out 
# --------- #

save(env_stk,file = file.path("./data","Cod_Env_data.rda"))


