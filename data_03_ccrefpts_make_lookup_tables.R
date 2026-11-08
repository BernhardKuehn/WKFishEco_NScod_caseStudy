# Northern shelf cod lookup tables -----------------------------------------
#
# EqSim settings: ICES 2026 reference-point scripts
# Reference points: ICES Advice 2026, cod.27.46a7d20, Table 4
#
# Substocks: Southern, Northwestern and Viking


# 1. EqSim assumptions -----------------------------------------------------

# srr_years describes the stock input window. Because recruitment is
# lagged relative to SSB, the report describes the recruitment period
# for these fits as 1998-2025.
#
# For Northwestern and Viking, ICES uses a fixed-breakpoint Segreg model
# for FMSY/Fpa, but a separate full-period Segreg fit for Blim.
# The fixed-breakpoint fit is described as type 5/6 in the ICES scripts.
#
# rhologRec and SSBcv were not explicitly supplied in the ICES scripts;
# the values below are defaults in the inspected EqSim implementation.
# SSBcv is not the assessment uncertainty used to calculate Bpa.

eqsim_lut = list(
  
  list(stock = "Southern",
       srr_years = 1997:2025, srr_type = 2,
       bio_years = 2016:2025, sel_years = 2023:2025,
       bio_const = FALSE, sel_const = FALSE,
       models = "Segreg", rhologRec = TRUE,
       Fcv = 0.212, Fphi = 0.423, SSBcv = 0,
       blim_srr_years = 1997:2025, blim_srr_type = 2,
       blim_models = "Segreg"),
  
  list(stock = "Northwestern",
       srr_years = 1997:2025, srr_type = NA_integer_,
       bio_years = 2016:2025, sel_years = 2023:2025,
       bio_const = FALSE, sel_const = FALSE,
       models = "SegregFix", rhologRec = TRUE,
       Fcv = 0.212, Fphi = 0.423, SSBcv = 0,
       blim_srr_years = 1983:2025, blim_srr_type = 2,
       blim_models = "Segreg", bloss_year = 2006,
       notes = "FMSY fit: type 5/6; breakpoint fixed at 2006 SSB (Bloss)."),
  
  list(stock = "Viking",
       srr_years = 1997:2025, srr_type = NA_integer_,
       bio_years = 2016:2025, sel_years = 2023:2025,
       bio_const = FALSE, sel_const = FALSE,
       models = "SegregFix", rhologRec = TRUE,
       Fcv = 0.212, Fphi = 0.423, SSBcv = 0,
       blim_srr_years = 1983:2025, blim_srr_type = 2,
       blim_models = "Segreg", bloss_year = 2006,
       notes = "FMSY fit: type 5/6; breakpoint fixed at 2006 SSB (Bloss).")
)

names(eqsim_lut) = vapply(eqsim_lut, function(x) x$stock, character(1))

if(!dir.exists("data")) dir.create("data", recursive = TRUE)
save(eqsim_lut, file = "data/eqsim_lut.Rdata")


# 2. Biological reference points ------------------------------------------
#
# Values from Table 4 of ICES Advice 2026, cod.27.46a7d20.
# Biomass reference points are in tonnes; F is for ages 2-4.
# MSY Btrigger is equal to Bpa for all three substocks.
# Flim and Fmgt are not given in the supplied table.

BRPs = data.frame(
  Flim = c(NA_real_, NA_real_, NA_real_),
  Fpa = c(0.395, 0.442, 0.409),
  Fmsy = c(0.243, 0.189, 0.198),
  FmsyL = c(0.154, 0.108, 0.121),
  FmsyH = c(0.381, 0.325, 0.345),
  Fmgt = c(NA_real_, NA_real_, NA_real_),
  Blim = c(14433, 23422, 14027),
  Bpa = c(21549, 32547, 19491),
  Btrigger = c(21549, 32547, 19491),
  row.names = c("Southern", "Northwestern", "Viking")
)

save(BRPs, file = "data/BRPs.Rdata")