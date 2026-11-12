# ------------------------ #
# Prepare output for F_eco
# ------------------------ # 

rm(list = ls())

library(collapse)

# load refpts
load("./data/BRPs.Rdata")
load("./model/newBRPs_EMSRR_model.rda")
load("./model/newBRPs_SRR_benchmark_model.rda")
load("./model/newBRPs_EMSRR_truncICES_model.rda")

BRPs$stock = rownames(BRPs)
BRPs$srr.type = "SEG"
BRPs$analysis_type = "SRR_eqsim"

newBRPs$analysis_type = "EMSRR_FLR"
newBRPs_bench$analysis_type = "SRR_FLR"
newBRPs_bench$srr.type = "SEG"
newBRPs_bench$stock = rownames(newBRPs_bench)

newBRPs_truncICES$analysis_type = "EMSRR_FLR_ICES.based.SRR"

BRPs_all = rowbind(BRPs,newBRPs,newBRPs_truncICES,newBRPs_bench,fill = T)

# select Fmsy, FmsyL, FmsyU
out = BRPs_all[,c("stock","Fmsy","FmsyL","FmsyH","srr.type","climScen","period","analysis_type")]

# no BH
out = out[out$srr.type != "BH",]
out$period = dplyr::recode(out$period,"None" = "Current")

out$Note = ifelse(out$analysis_type == "SRR_eqsim","Those are the official ICES refpts from the 2026 assessment.",
              ifelse(out$analysis_type == "SRR_FLR","Those are the refpts calculated in FLR, following the benchmark/2026 assessment procedure as close as possible.",
                ifelse(out$analysis_type == "EMSRR_FLR","Those are the refpts calculated in FLR with an environmental influence on recruitment.",
                  ifelse(out$analysis_type == "EMSRR_FLR_ICES.based.SRR","Those are the refpts calculated in FLR with an environmental influence on rec, \nbut with the R-SSB relationship fixed at the benchmark relationship.",NA))))


ggplot(out[out$climScen %in% c("None",NA,"noCC"),],aes(x = analysis_type,y = Fmsy)) +
  geom_point() + 
  facet_wrap(~stock)


# write out
write.table(out,"./results/Cod_substock_Fmsy_ccrefpts.csv",col.names = colnames(out),
            row.names = F,sep = ",")

