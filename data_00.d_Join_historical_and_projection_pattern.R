# ================================================================================ #
# Prepare historical and projected bottom-temperature time series
#
# 1. Recalculate the historical substock means and cluster time series
# 2. Compare these with the time series saved by data_01
# 3. Find the cells shared by history and all RCP realizations
# 4. Calculate annual means for the projected fields using the historical clusters
# 5. Compare history and projections, including their overlapping years
#
# History and projections remain separate; no lag or splicing is applied here.
# Run from the repository root after data_01, data_00.b and data_00.c.
# ================================================================================ #

rm(list = ls())

library(raster)
library(terra)

# ---------------- #
# helper functions
# ---------------- #

# Read the monthly dates stored in the projection raster or its layer names.
get.proj.time = function(x,filename){
  time.rst = terra::time(x)
  if(length(time.rst) != nlyr(x) || anyNA(time.rst)){
    time.names = gsub("\\.","-",sub("^X","",names(x)))
    if(any(!grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$",time.names))){
      stop(filename,": no usable monthly dates in the raster.")
    }
    time.rst = as.Date(time.names,format = "%Y-%m-%d")
  } else {
    time.rst = as.Date(time.rst)
  }
  if(anyNA(time.rst)) stop(filename,": invalid monthly dates.")
  year = as.integer(format(time.rst,"%Y"))
  month = as.integer(format(time.rst,"%m"))
  month.index = year*12L + month
  if(anyDuplicated(month.index) || any(diff(month.index) != 1L) ||
     any(table(year) != 12L)){
    stop(filename,": expected consecutive, complete years of monthly data.")
  }
  return(time.rst)
}

# Apply the same spatial summaries to annual historical and projected values.
# The substock mean is unweighted (as in data_01); cluster means use the
# raster::area weights used by spatial_kmeans on the original cluster map.
calc.substock.ts = function(fields,years,region.info,common = TRUE){
  out = data.frame(year = years)
  for(substk in names(region.info)){
    info = region.info[[substk]]
    cells = if(common) info$common else info$historical
    vals = fields[[substk]]
    out[[paste0("bottomT_",substk)]] =
      colMeans(vals[cells,,drop = FALSE],na.rm = !common)
    for(cluster in info$cluster.ids){
      index = which(cells & !is.na(info$cluster) & info$cluster == cluster)
      if(!length(index)) stop(substk,": no cells remain in Cluster",cluster)
      weights = info$area[index]/sum(info$area[index])
      out[[paste0("bottomT_",substk,"_Cluster",cluster)]] =
        drop(crossprod(weights,vals[index,,drop = FALSE]))
    }
  }
  return(out)
}

# ---------------------- #
# load historical fields
# ---------------------- #

dir.hist = "./data/Env/processed"
dir.proj = "./data/Env/POLCOMS-ERSEM_biasRM_BVARsims/bottomT"
dir.out = "./data/Env/processed/bottomT_projections"

load(file.path(dir.hist,"Spatial_bottomT_yr_by_cod_substock.rda"))
load(file.path(dir.hist,"SpatialCluster_bottomT_yr_by_cod_substock.rda"))
load(file.path(dir.hist,"bottomT_yr_by_cod_substock.rda"))
load(file.path(dir.hist,"Cluster_bottomT_yr_by_cod_substock.rda"))

if(!is.list(bottomT_yr_by_substock_spatial) ||
   !length(bottomT_yr_by_substock_spatial) ||
   is.null(names(bottomT_yr_by_substock_spatial)) ||
   any(!nzchar(names(bottomT_yr_by_substock_spatial))) ||
   anyDuplicated(names(bottomT_yr_by_substock_spatial)) ||
   !is.list(cluster_bottomT_substk_spatial) ||
   is.null(names(cluster_bottomT_substk_spatial)) ||
   anyDuplicated(names(cluster_bottomT_substk_spatial)) ||
   !setequal(names(bottomT_yr_by_substock_spatial),names(cluster_bottomT_substk_spatial))){
  stop("Historical fields and cluster maps must have the same unique substock names.")
}
if(!is.data.frame(bottomT_substk) || !is.data.frame(cluster_bottomT_substk) ||
   !"year" %in% names(bottomT_substk) || !"year" %in% names(cluster_bottomT_substk) ||
   anyNA(bottomT_substk$year) || anyNA(cluster_bottomT_substk$year) ||
   anyDuplicated(bottomT_substk$year) || anyDuplicated(cluster_bottomT_substk$year)){
  stop("Saved historical time series must have unique, non-missing years.")
}

# Store the historical values, cluster IDs and cell areas in the same order.
# A cluster map includes only cells with a complete historical time series.
region.info = hist.fields = list()
hist.extents = list()
years.hist = NULL
for(substk in names(bottomT_yr_by_substock_spatial)){
  hist = tryCatch({
    x = terra::rast(bottomT_yr_by_substock_spatial[[substk]])
    nlyr(x)
    x
  },error = function(e){
    stop(substk,": cannot restore historical raster (",conditionMessage(e),
         "). Regenerate it with the updated data_01 script.",call. = FALSE)
  })
  clusters = terra::rast(cluster_bottomT_substk_spatial[[substk]])
  if(nlyr(clusters) != 1L ||
     !isTRUE(terra::compareGeom(hist[[1]],clusters,crs = TRUE,ext = TRUE,
                                 rowcol = TRUE,res = TRUE,stopOnError = FALSE))){
    stop(substk,": historical field and cluster map do not share a grid.")
  }
  if(!isTRUE(terra::is.lonlat(hist))) stop(substk,": expected a longitude/latitude grid.")
  year.names = sub("^X","",names(hist))
  if(any(!grepl("^[0-9]{4}$",year.names))){
    stop(substk,": historical raster layers must be named by year.")
  }
  years = as.integer(year.names)
  if(anyDuplicated(years) || is.unsorted(years,strictly = TRUE)){
    stop(substk,": historical years must be unique and increasing.")
  }
  if(is.null(years.hist)) years.hist = years
  if(!identical(years,years.hist)) stop(substk,": historical years differ among substocks.")

  if(length(region.info)){
    reference = terra::rast(region.info[[1]]$template)
    same = terra::compareGeom(hist[[1]],reference,crs = TRUE,ext = FALSE,
                              rowcol = FALSE,res = TRUE,stopOnError = FALSE)
    offset = (terra::origin(hist)-terra::origin(reference))/terra::res(reference)
    if(!isTRUE(same) || any(abs(offset-round(offset)) > 1e-6)){
      stop(substk,": historical substocks do not share the same grid.")
    }
  }
  values.hist = terra::values(hist,mat = TRUE)
  cluster.ids = terra::values(clusters,mat = FALSE)
  ids = sort(unique(cluster.ids[!is.na(cluster.ids)]))
  if(!length(ids) || any(!is.finite(ids) | ids < 1 | ids != floor(ids))){
    stop(substk,": cluster IDs must be positive integers.")
  }
  cell.area = raster::getValues(raster::area(raster::raster(hist[[1]])))
  complete = rowSums(is.finite(values.hist)) == ncol(values.hist) &
    is.finite(cell.area) & cell.area > 0
  if(any(!is.na(cluster.ids) & !complete)){
    stop(substk,": cluster map includes cells without complete historical data.")
  }
  region.info[[substk]] = list(historical = rowSums(is.finite(values.hist)) > 0,
                               common = complete,cluster = cluster.ids,
                               cluster.ids = ids,area = cell.area,
                               xy = terra::xyFromCell(hist,seq_len(ncell(hist))),
                               template = raster::raster(hist[[1]]))
  hist.fields[[substk]] = values.hist
  hist.extents[[substk]] = as.vector(terra::ext(hist))
}

# ---------------------------------------------------- #
# recalculate and compare the historical time series
# ---------------------------------------------------- #

hist.saved = merge(bottomT_substk,cluster_bottomT_substk,by = "year",all = TRUE)
hist.full = calc.substock.ts(hist.fields,years.hist,region.info,common = FALSE)
variables = setdiff(names(hist.full),"year")
if(!all(variables %in% names(hist.saved))){
  stop("Saved historical time series are missing: ",
       paste(setdiff(variables,names(hist.saved)),collapse = ", "))
}

# A full-footprint mismatch points to different inputs or processing in data_01.
# The common-footprint comparison is added after finding the projection coverage.
hist.comparison = lapply(variables,function(var){
  years = sort(unique(c(hist.saved$year,hist.full$year)))
  saved = hist.saved[[var]][match(years,hist.saved$year)]
  full = hist.full[[var]][match(years,hist.full$year)]
  data.frame(year = years,variable = var,saved = saved,recomputed_full = full,
             full_minus_saved = full-saved)
})
hist.comparison = do.call(rbind,hist.comparison)
if(any(!is.finite(hist.comparison$saved) | !is.finite(hist.comparison$recomputed_full))){
  warning("Saved and recalculated historical values do not all match by year.")
}
if(any(abs(hist.comparison$full_minus_saved) > 1e-6,na.rm = TRUE)){
  warning("Full-footprint historical values differ from the saved time series.")
}

# ------------------------------------------------------- #
# find the footprint shared by history and both RCP runs
# ------------------------------------------------------- #

# Select and order the files written by data_00.c.
clim.proj = c("rcp45","rcp85")
if(!dir.exists(dir.proj)) stop("Missing projection directory: ",dir.proj)
all.files = list.files(dir.proj,pattern = "\\.tiff?$",full.names = TRUE,ignore.case = TRUE)
file.pattern = "^POLCOM_ERSEM_(rcp45|rcp85)_bottomT_([0-9]+)\\.tiff?$"
if(any(!grepl(file.pattern,basename(all.files),ignore.case = TRUE))){
  stop("Unexpected projection TIFF filename(s): ",
       paste(basename(all.files[!grepl(file.pattern,basename(all.files),
                                       ignore.case = TRUE)]),collapse = ", "))
}
proj.files = setNames(vector("list",length(clim.proj)),clim.proj)
for(rcp in clim.proj){
  files = all.files[grepl(paste0("^POLCOM_ERSEM_",rcp,"_bottomT_[0-9]+\\.tiff?$"),
                          basename(all.files),ignore.case = TRUE)]
  if(!length(files)) stop("No projection files found for ",rcp)
  ids = as.integer(sub(file.pattern,"\\2",basename(files),ignore.case = TRUE))
  if(anyNA(ids) || anyDuplicated(ids)) stop("Invalid or duplicate realization numbers for ",rcp)
  proj.files[[rcp]] = setNames(files[order(ids)],as.character(sort(ids)))
}

# Scan every realization. Retain only cells with complete historical data
# and finite values in every projection month, in both scenarios.
xy.ext = do.call(rbind,hist.extents)
hist.bounds = terra::ext(min(xy.ext[,1]),max(xy.ext[,2]),
                         min(xy.ext[,3]),max(xy.ext[,4]))
reference = terra::rast(region.info[[1]]$template)
months.rcp = list()
for(rcp in clim.proj){
  for(it in names(proj.files[[rcp]])){
    filename = proj.files[[rcp]][[it]]
    cat("Checking footprint:",basename(filename),"\n")
    proj = terra::rast(filename)
    same = terra::compareGeom(proj[[1]],reference,crs = TRUE,ext = FALSE,
                              rowcol = FALSE,res = TRUE,stopOnError = FALSE)
    offset = (terra::origin(proj)-terra::origin(reference))/terra::res(reference)
    if(!isTRUE(same) || any(abs(offset-round(offset)) > 1e-6)){
      stop(filename,": CRS, resolution or grid origin differs from history.")
    }
    time.proj = get.proj.time(proj,filename)
    months = format(time.proj,"%Y-%m")
    if(is.null(months.rcp[[rcp]])){
      months.rcp[[rcp]] = months
    } else if(!identical(months.rcp[[rcp]],months)){
      stop(filename,": months differ among ",rcp," realizations.")
    }
    proj = terra::crop(proj,hist.bounds)
    valid = terra::app(is.finite(proj),"sum") == nlyr(proj)
    valid.values = terra::values(valid,mat = FALSE)
    for(substk in names(region.info)){
      index = terra::cellFromXY(valid,region.info[[substk]]$xy)
      present = valid.values[index]
      region.info[[substk]]$common = region.info[[substk]]$common &
        !is.na(present) & present == 1
    }
  }
}

# Report how much of each substock and cluster is retained.
coverage = list()
common.masks = list()
for(substk in names(region.info)){
  info = region.info[[substk]]
  groups = c(list(mean = info$historical),
             setNames(lapply(info$cluster.ids,function(id)
               !is.na(info$cluster) & info$cluster == id),
               paste0("Cluster",info$cluster.ids)))
  for(group in names(groups)){
    n.hist = sum(groups[[group]])
    n.common = sum(groups[[group]] & info$common)
    if(!n.common) stop(substk," ",group,": no cells on the common footprint.")
    coverage[[length(coverage)+1L]] =
      data.frame(substock = substk,group = group,historical_cells = n.hist,
                 common_cells = n.common,retained_fraction = n.common/n.hist)
  }
  common.masks[[substk]] =
    raster::setValues(info$template,ifelse(info$common,1,NA_real_))
}
coverage = do.call(rbind,coverage)
print(coverage)
if(any(coverage$retained_fraction < 1)){
  warning("The common footprint excludes historical cells; inspect the coverage and comparison.")
}

# Recalculate history on that same footprint and compare with the saved series.
hist.common = calc.substock.ts(hist.fields,years.hist,region.info)
hist.comparison$recomputed_common = mapply(function(year,var){
  hist.common[[var]][match(year,hist.common$year)]
},hist.comparison$year,hist.comparison$variable)
hist.comparison$common_minus_saved =
  hist.comparison$recomputed_common-hist.comparison$saved

# ------------------------------------------------------------ #
# annual projected means using the historical areas and clusters
# ------------------------------------------------------------ #

bottomT.proj = setNames(vector("list",length(clim.proj)),clim.proj)
for(rcp in clim.proj){
  bottomT.proj[[rcp]] = setNames(vector("list",length(proj.files[[rcp]])),
                                 names(proj.files[[rcp]]))
  for(it in names(proj.files[[rcp]])){
    filename = proj.files[[rcp]][[it]]
    cat("Aggregating:",basename(filename),"\n")
    proj = terra::crop(terra::rast(filename),hist.bounds)
    time.proj = get.proj.time(proj,filename)
    year.proj = as.integer(format(time.proj,"%Y"))

    # Average the twelve monthly values per cell, as in data_01.
    proj.yr = terra::tapp(proj,year.proj,"mean")
    proj.values = terra::values(proj.yr,mat = TRUE)
    proj.fields = lapply(region.info,function(info){
      index = terra::cellFromXY(proj.yr,info$xy)
      proj.values[index,,drop = FALSE]
    })
    for(substk in names(region.info)){
      if(any(!is.finite(proj.fields[[substk]][region.info[[substk]]$common,,drop = FALSE]))){
        stop(filename,": missing data on the common footprint for ",substk)
      }
    }
    bottomT.proj[[rcp]][[it]] =
      calc.substock.ts(proj.fields,sort(unique(year.proj)),region.info)
  }
}

# ----------------------------------------------------- #
# compare the historical and projected temporal patterns
# ----------------------------------------------------- #

# Compare the median and 10-90% range of the realizations against history.
# Bias correction does not imply that individual stochastic runs match history.
overlap = list()
for(rcp in clim.proj){
  runs = bottomT.proj[[rcp]]
  overlap.years = intersect(hist.common$year,runs[[1]]$year)
  if(!length(overlap.years)){
    warning("Historical and ",rcp," time series have no overlapping years.")
    next
  }
  for(var in variables){
    index = match(overlap.years,runs[[1]]$year)
    values = matrix(vapply(runs,function(x) x[[var]][index],
                           numeric(length(overlap.years))),nrow = length(overlap.years))
    qs = t(apply(values,1,quantile,probs = c(0.1,0.5,0.9),names = FALSE))
    observed = hist.common[[var]][match(overlap.years,hist.common$year)]
    overlap[[length(overlap)+1L]] =
      data.frame(scenario = rcp,variable = var,year = overlap.years,
                 historical = observed,projection_p10 = qs[,1],
                 projection_median = qs[,2],projection_p90 = qs[,3],
                 median_minus_historical = qs[,2]-observed)
  }
}
if(length(overlap)){
  overlap = do.call(rbind,overlap)
} else {
  overlap = data.frame(scenario = character(),variable = character(),year = integer(),
                       historical = numeric(),projection_p10 = numeric(),
                       projection_median = numeric(),projection_p90 = numeric(),
                       median_minus_historical = numeric())
}

# Plot saved vs. recalculated history and the projections for every variable.
if(!dir.exists(dir.out) && !dir.create(dir.out,recursive = TRUE)){
  stop("Could not create output directory: ",dir.out)
}
pdf(file.path(dir.out,"historical_projection_comparison.pdf"),width = 13,height = 5)
par(mfrow = c(1,2),mar = c(4,4,3,1))
cols = c(rcp45 = "olivedrab",rcp85 = "red")
for(var in variables){
  check = hist.comparison[hist.comparison$variable == var,]
  ylim.hist = range(c(check$saved,check$recomputed_full,check$recomputed_common),finite = TRUE)
  plot(check$year,check$saved,type = "l",col = "grey40",lty = 2,
       ylim = ylim.hist,xlab = "Year",ylab = "Annual bottom temperature",main = var)
  lines(check$year,check$recomputed_full,col = "blue")
  lines(check$year,check$recomputed_common,lwd = 2)
  legend("topleft",c("Saved","Recalculated, full","Recalculated, common"),
         col = c("grey40","blue","black"),lty = c(2,1,1),bty = "n",cex = 0.8)

  ranges = lapply(bottomT.proj,function(runs){
    years = runs[[1]]$year
    values = matrix(vapply(runs,function(x) x[[var]],numeric(length(years))),
                    nrow = length(years))
    qs = t(apply(values,1,quantile,probs = c(0.1,0.5,0.9),names = FALSE))
    data.frame(year = years,low = qs[,1],median = qs[,2],high = qs[,3])
  })
  xlim = range(c(hist.common$year,unlist(lapply(ranges,`[[`,"year"))))
  ylim = range(c(hist.common[[var]],unlist(lapply(ranges,function(x) c(x$low,x$high)))),
               finite = TRUE)
  plot(hist.common$year,hist.common[[var]],type = "n",xlim = xlim,ylim = ylim,
       xlab = "Year",ylab = "Annual bottom temperature",main = var)
  for(rcp in clim.proj){
    x = ranges[[rcp]]
    polygon(c(x$year,rev(x$year)),c(x$low,rev(x$high)),border = NA,
            col = adjustcolor(cols[rcp],alpha.f = 0.16))
    lines(x$year,x$median,col = cols[rcp],lwd = 2)
  }
  lines(hist.common$year,hist.common[[var]],lwd = 2)
  legend("topleft",c("History, common","RCP4.5 median (10-90%)",
                     "RCP8.5 median (10-90%)"),
         col = c("black",cols),lty = 1,lwd = 2,bty = "n",cex = 0.8)
}
dev.off()

# ------------------- #
# write out the data
# ------------------- #

# Each projected realization is an unlagged wide table with the same variable
# names as the saved historical series. History and projections are not joined.
bottomT_projection_results = list(historical_full = hist.full,
                                  historical_common = hist.common,
                                  historical_comparison = hist.comparison,
                                  coverage = coverage,common_masks = common.masks,
                                  projections = bottomT.proj,overlap = overlap,
                                  projection_files = proj.files)
saveRDS(bottomT_projection_results,file.path(dir.out,"bottomT_historical_projection_summaries.rds"))
write.csv(hist.comparison,file.path(dir.out,"historical_comparison.csv"),row.names = FALSE)
write.csv(coverage,file.path(dir.out,"common_footprint_coverage.csv"),row.names = FALSE)
write.csv(overlap,file.path(dir.out,"historical_projection_overlap.csv"),row.names = FALSE)
