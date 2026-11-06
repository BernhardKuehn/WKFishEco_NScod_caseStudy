# ==================================== #
# wrapper function to fit a BVAR
# Bayesian Vector autoregressive model
# to spatial data
# - used to capture dependencies in
#   the climate projections of various
#   oceanographic variables
# - generate new "data" that follow a 
#   certain autocorrelation and covariance
#   structure found in the data
# steps:
# - perform EOF analysis per variable
#   to have a lower dimensional space
# - remove nonlinear, monotone trend
# - fit BVAR
# - simulate new trajectories with VAR
# - generate new spatial fields from PCs
# ==================================== #

# 1. function to fit BVAR model on raster files
fit.spatial.BVAR = function(raster,
                            time_raster,
                            PC.thresh = 0.8,min.PCs = 5,
                            lags = 1, n_draw = 10000, n_burn = 5000, n_thin = 1,
                            lambda = bv_lambda(mode = 0.2, sd = 0.4, min = 0.0001, max = 5),
                            alpha = bv_alpha(mode = 2, sd = 0.25, min = 1, max = 3),
                            metropolis.h = bv_metropolis(scale_hess = 0.005, adjust_acc = TRUE,
                                                         acc_lower = 0.25, acc_upper = 0.35, acc_change = 0.02),
                            plot = T
){
  library(raster)
  library(marmalaid)
  library(BVAR)
  library(coda)
  
  # helper function to remove nonlinear trend (via scam)
  
  nonlinear.Trend = function(time,ts,...){
    library(mgcv)
    library(scam)
    data = data.frame(time,ts)
    # check linear trend
    lm.trend = lm(ts ~ time,data = data)$coefficients[[2]]
    if(lm.trend >= 0){
      # monotone increasing p-spline
      formula = as.formula(paste("ts","s(time,bs = 'mpi')",sep = " ~ "))
    } else{
      # monotone decreasing p-spline
      formula = as.formula(paste("ts","s(time,bs = 'mpd')",sep = " ~ "))
    }
    #formula = as.formula(paste("ts","s(time,bs = 'ps',m = c(3,2))",sep = " ~ "))
    #mod = gam(formula,data = data[complete.cases(data),],family = gaussian(link = "identity"))
    mod = scam(formula,data = data[complete.cases(data),],family = gaussian(link = "identity"),...)
    
    return(predict(mod,newdata = data.frame(time = time)))
  }
  
  
  # 1. decompose via EOF-analysis
  EOF_list = lapply(raster,function(x){
    if(class(x) %in% c("Raster","RasterBrick","RasterStack")){
      extent.x = extent(x)
    } else {
      extent.x = extent(x[[1]])
    }
    eof = spatial.PCA(x,center = T,scale. = F,
                      spatial.extent = extent.x)
    return(eof)})
  # 2. select limited nr. of PCs
  cum.var = sapply(EOF_list,
                   function(x) max(c(min.PCs,which(cumsum(x$sp.PCA$sdev^2/sum(x$sp.PCA$sdev^2)) >= PC.thresh)[1])))
  df.expl.var = sapply(EOF_list,
                    function(x) cumsum(x$sp.PCA$sdev^2/sum(x$sp.PCA$sdev^2)))
  expl.var = sapply(1:length(cum.var),function(x) df.expl.var[cum.var[x],x])
  
  PCs_select = do.call(cbind,lapply(1:length(EOF_list),
                                    function(x) EOF_list[[x]]$sp.PCA$x[,1:cum.var[[x]]]))
  # 3. remove non-linear trend
  PCs_trend = apply(PCs_select,2,function(x) nonlinear.Trend(time = 1:length(x),ts = x,optimizer = "efs"))
  PCs_select_trendRM = PCs_select - PCs_trend
  
  # additional scaling step to make the generation of new realisations easier and more stable
  PCs_select_trendRM_scaled = scale(PCs_select_trendRM)
  
  
  
  # 4. fit BVAR model
  # specify prior
  mn <- bv_minnesota(lambda = lambda,
                     alpha = alpha,var = 1e07)
  priors <- bv_priors(hyper = "auto", mn = mn)
  # metropolis-hastings
  metropolis.h <- metropolis.h
  
  #fit BVAR
  mar = bvar(data = PCs_select_trendRM_scaled,
             lags = lags,
             n_draw = n_draw, n_burn = n_burn, n_thin = n_thin,
             priors = priors, mh = metropolis.h,verbose = T)
  # plot 
  if(plot == T){
    plot(mar)
  }
  # diagnostics
  BIC.mar = BIC(mar)
  AIC.mar = AIC(mar)
  # check convergence
  run_mcmc <- as.mcmc(mar)
  # test for convergence Geweke`s convergence diagnostic
  convergence = apply(run_mcmc,2,geweke.diag)
  df.convergence = data.frame(parameter = names(convergence),
                              z.score = sapply(convergence,function(x) x$z))
  df.convergence$p.value = 2*pnorm(-abs(df.convergence$z.score),
                                   mean = 0, sd = 1, lower.tail = TRUE)
  # add sign
  df.convergence$sign = ifelse(df.convergence$p.value <= 0.05,"*",F)
  df.convergence$converged = ifelse(df.convergence$p.value <= 0.05,F,T)
  # different ordering
  df.convergence = df.convergence[,c("parameter","converged","z.score","p.value","sign")]
  
  # 5. store everything to be used in further functions
  bvar_summary = summary(mar)
  coefficients = bvar_summary$coef
  vcov = bvar_summary$vcov
  
  # 6. read out as list
  out = list(BVAR_model = mar,
             BVAR_coefficients = coefficients,
             BVAR_vcov = vcov,
             diagnostics = list(convergence.mcmc = df.convergence,
                                fit = c(AIC = AIC.mar,
                                        BIC = BIC.mar)),
             EOFs = EOF_list,
             nr_pcs_selected = data.frame(rbind(no.PCs = as.character(cum.var),
                                                cum.expl.var = round(100*expl.var,digits = 2))),
             y = PCs_select_trendRM_scaled,
             scale = list(center = attributes(PCs_select_trendRM_scaled)[["scaled:center"]],
                          scale = attributes(PCs_select_trendRM_scaled)[["scaled:scale"]]),
             Trend = PCs_trend)
  return(out)
}

# # 1.1 function to calculate further diagnostics that 
# # compare the quality of simulated VAR processes
# diagnostics.spatial.BVAR = function(BVAR.fit,N.sim){
#   
#   # helper func
#   acf_compare = function(obs.mat,pred.mat,metric){
#     obs.acf = apply(obs.mat,2,function(x) acf(x,plot = F)$acf[-1])
#     pred.acf = apply(pred.mat,2,function(x) acf(x,plot = F)$acf[-1])
#     m = vector()
#     for(i in 1:ncol(obs.acf)){
#       m[i] = metric(data.frame(obs = obs.acf[,i],
#                                pred = pred.acf[,i]))
#     }
#     return(m)
#   }
#   ccf_compare = function(obs.mat,pred.mat,metric){
#     combs = combn(1:ncol(obs.mat),2)
#     
#     obs.ccf = apply(combs,2,function(i) ccf(obs.mat[,i[1]],
#                                             obs.mat[,i[2]],plot = F)$acf)
#     pred.ccf = apply(combs,2,function(i) ccf(pred.mat[,i[1]],
#                                              pred.mat[,i[2]],plot = F)$acf)
#     
#     m = vector()
#     for(i in 1:ncol(obs.ccf)){
#       m[i] = metric(data.frame(obs = obs.ccf[,i],
#                                pred = pred.ccf[,i]))
#     }
#     return(m)
#   }
#   
#   # 1. get coefficients and var-cov-matrix for a BVAR model
#   coefs = BVAR.fit$BVAR_coefficients
#   vcov = BVAR.fit$BVAR_vcov
#   lags = BVAR.fit$BVAR_model$meta$lags
#   y = BVAR.fit$y
#   trend = BVAR.fit$Trend
#   
#   if(is.null(L)){
#     L = nrow(y)
#   }
#   # 2. simulate/produce new realisations
#   acf.match = vector()
#   for(i in 1:N.sim){
#     # sample starting values
#     index.start = sample(1:nrow(y),1,replace = T)
#     start.val = t(y[index.start,])
#     # simulate VAR
#     tmp = VAR.sim(B = t(coefs),n = L,lag = lags,starting = start.val)
#     # calculate diagnostics in terms of acf, ccf
#     m1 = acf_compare(obs.mat = y,pred.mat = tmp,metric = function(x) cor(x[,1],x[,2]))
#     
#     m1_random = acf_compare(obs.mat = y,pred.mat = apply(y,2,function(x) rnorm(1:nrow(y),
#                                                                                mean = mean(x),sd(x))),metric = function(x) cor(x[,1],x[,2]))
#     par(mfrow = c(1,1))
#     plot(m1,type = "l",ylim = range(c(m1,m1_random)))
#     
#     for(i in 1:100){
#       m1_random = acf_compare(obs.mat = y,pred.mat = apply(y,2,function(x) rnorm(1:nrow(y),
#                                                                                  mean = mean(x),sd(x))),metric = function(x) cor(x[,1],x[,2]))
#       
#     points(m1_random,col = "red",type = "l")
#     }
#     abline(v = cumsum(BVAR.fit$nr_pcs_selected),lty =2)
#     sum(m1)
#     sum(m1_random)
#     
#     
#     m2 = ccf_compare(obs.mat = y,pred.mat = tmp,metric = RMSE)
#     par(mfrow = c(1,2))
#     plot(m1)
#     plot(m2)
#     median(m)
#     acf.match[i] = sum(m) # collect
#     
#     ## add Trend back on 
#     #tmp = tmp + trend
#   }
#   
#   
#   # 3. write out/compress
#   
#   
#   
#   
# }

# 2. function to generate new realisations from this model
# !!!! Needs parallelization in order to be effective !!!!!
sim.spatial.BVAR = function(BVAR.fit,N.sim,L = NULL,start.val.type = "mean",return.raster = F,parallel = F,ncores = NULL){
  # function arguments
  # - BVAR.fit = object coming from the "fit.spatial.BVAR" function
  # - N.sim = number of simulations
  # L = time series length
  # start.val.type should be either mean or first10
  # return.raster = if a raster file should be returned, otherwise only the matrices
  
  library(Rcpp)
  library(Rcpp2doParallel)
  library(tsDyn)
  library(raster)
  library(terra)
  if(parallel == T){
    library(parallel)
    library(doSNOW)
    library(foreach)
  }
  
  # helper func
  split.by.col = function(df,f){
    list.spl = split.data.frame(t(df),f)
    list.spl = lapply(list.spl,t)
    return(list.spl)
  }
  
  # 1. get coefficients and var-cov-matrix for a BVAR model
  coefs = BVAR.fit$BVAR_coefficients
  vcov = BVAR.fit$BVAR_vcov
  lags = BVAR.fit$BVAR_model$meta$lags
  scale.params = BVAR.fit$scale # scale parameters to bring back to original scale
  y = BVAR.fit$y
  split.vec = rep(1:ncol(BVAR.fit$nr_pcs_selected),as.integer(BVAR.fit$nr_pcs_selected[1,]))
  trend = BVAR.fit$Trend
  if(is.null(L)){
    L = nrow(y)
  }
  # 2. simulate/produce new realisations
  vars_sim = list()
  for(i in 1:N.sim){
    # sample starting values (from the first 10% of data)
    if(start.val.type == "first10%"){
      index.start = sample(1:ceiling(0.1*nrow(y)),1,replace = T)
      # start index representing the one with the largest lag
      start.val = matrix(y[index.start:(index.start+(lags-1)),],
                         nrow = lags,ncol = ncol(y))
    } else if(start.val.type == "mean"){
      start.val = matrix(apply(y[1:ceiling(0.1*nrow(y)),],2,mean),
                         nrow = lags,ncol = ncol(y))
      
    } else if(start.val.type == "first"){
      start.val = matrix(y[1:(lags-1),],
                         nrow = lags,ncol = ncol(y))
      
      } else{
      stop("start.val.type needs to be either `mean`, `first` or `first10%`!")
    }
        # simulate VAR
    tmp = VAR.sim(B = t(coefs),n = L,lag = lags,starting = start.val,
                  varcov = vcov)
    # rescale
    tmp = t(t(tmp)*scale.params[["scale"]] + scale.params[["center"]])
    # add Trend back on 
    vars_sim[[i]] = tmp + trend
  }
  
  if(parallel == T){
    # --- parallel mode --- #
    if(is.null(ncores)){
      ncores = detectCores(logical = F)-1
    }
    
    # set up progress bar
    ntasks <- length(vars_sim)
    pb <- txtProgressBar(max = ntasks, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    
    # set up parallel mode
    cl = makeCluster(ncores,type = "SOCK")
    registerDoSNOW(cl)
    #registerDoParallel(cl)
    
    opts <- list(progress=progress)
    
    # split and add further PCs
    par.out = foreach(n = 1:length(vars_sim),.combine = "list",
                      .multicombine = T,.packages = c("raster","Rcpp","Rcpp2doParallel","terra"),.options.snow = opts) %dopar% {
                        spl = split.by.col(vars_sim[[n]],f = split.vec)
                        
                        # use C++ for fast matrix multiplication
                        # (problem: needs to be defined on each worker, rather than globally for all...
                        # or put into a R-package....)
                        cppFunction(depends = c("RcppArmadillo","RcppEigen"),
                                    'SEXP matMult_rcpp(Eigen::MatrixXd A, Eigen::MatrixXd B){
                                     Eigen::MatrixXd C = A * B;
                                     return Rcpp::wrap(C);
                                     }')
                        
                        
                        #vars_sim_Covar = list()
                        vars_sim_Covar = vector("list",ncol(BVAR.fit$nr_pcs_selected))
                        for(i in 1:ncol(BVAR.fit$nr_pcs_selected)){
                          # add other PCs
                          PCs.sim = BVAR.fit$EOFs[[i]]$sp.PCA$x
                          PCs.sim[,1:ncol(spl[[i]])] = spl[[i]]
                          # build spatial field
                          Mat.sp.sim = matMult_rcpp(PCs.sim,t(BVAR.fit$EOFs[[i]]$sp.PCA$rotation))
                          Mat.sp.sim = scale(Mat.sp.sim,center = -BVAR.fit$EOFs[[i]]$sp.PCA$center,scale = F)
                          
                          if(return.raster == F){
                            vars_sim_Covar[[i]] = Mat.sp.sim
                            cat(names(BVAR.fit$nr_pcs_selected)[i],"|",paste("iteration:",n),"\n")
                          } else {
                            # get coordinates and create raster
                            ll = length(BVAR.fit$EOFs[[i]]$raster)
                            if(class(BVAR.fit$EOFs[[i]]$raster) == "list" & ll > 1){
                              rst.sim = vector("list",ll)
                              spl.index = rep(1:ll,each = ncol(Mat.sp.sim)/ll)
                              Mat.sp.sim_split = split.by.col(Mat.sp.sim,spl.index)
                              for(ii in 1:ll){
                                # use terra as it is faster
                                tmp.mat = terra::as.data.frame(BVAR.fit$EOFs[[i]]$raster[[ii]])
                                coords = coordinates(BVAR.fit$EOFs[[i]]$raster[[ii]])
                                coords = coords[complete.cases(tmp.mat),]
                                #rst.sim[[ii]] = rasterFromXYZ(cbind(coords,t(Mat.sp.sim_split[[ii]])))
                                rst.sim[[ii]] = terra::rast(cbind(coords,t(Mat.sp.sim_split[[ii]])),
                                                            type = "xyz",crs = "EPSG:4326")
                                names(rst.sim[[ii]]) = BVAR.fit$EOFs[[i]]$time
                                # load raster into workspace (pack)
                                 rst.sim[[ii]] = terra::wrap(rst.sim[[ii]])
                              }
                            } else {
                              tmp.mat = terra::as.data.frame(BVAR.fit$EOFs[[i]]$raster)
                              coords = coordinates(BVAR.fit$EOFs[[i]]$raster)
                              coords = coords[complete.cases(tmp.mat),]
                              #rst.sim = rasterFromXYZ(cbind(coords,t(Mat.sp.sim)))
                              rst.sim = terra::rast(cbind(coords,t(Mat.sp.sim)),
                                                    type = "xyz",crs = "EPSG:4326")
                              names(rst.sim) = BVAR.fit$EOFs[[i]]$time
                            # load raster into workspace (pack)
                             rst.sim = terra::wrap(rst.sim)
                             }
                            cat(names(BVAR.fit$nr_pcs_selected)[i],"|",paste("iteration:",n),"\n")
                            # store in nested list
                            vars_sim_Covar[[i]] = rst.sim
                          }
                        }
                        # return from parallel calculations
                        vars_sim_Covar
                      }
    close(pb)
    stopCluster(cl)
    # bring to same structure as in the serial mode
    vars.sim_rst = vector("list",ncol(BVAR.fit$nr_pcs_selected))
    for(n in 1:length(par.out)){
      tmp = par.out[[n]]
      for(i in 1:length(tmp)){
        vars.sim_rst[[i]][[n]] = tmp[[i]]
      }
    }
    names(vars.sim_rst) = names(BVAR.fit$nr_pcs_selected)
    # and unpack the wrapped rasters coming from the nodes
    vars.sim_rst = rapply(vars.sim_rst,f = rast,how = "list")
    
  } else {
    # use C++ for fast matrix multiplication
    cppFunction(depends = c("RcppArmadillo","RcppEigen"),
                'SEXP matMult_rcpp(Eigen::MatrixXd A, Eigen::MatrixXd B){
    Eigen::MatrixXd C = A * B;
    return Rcpp::wrap(C);
    }')
    # --- serial mode --- #  
    # split and add further PCs
    vars.sim_rst = vector("list",ncol(BVAR.fit$nr_pcs_selected))
    for(n in 1:length(vars_sim)){
      spl = split.by.col(vars_sim[[n]],f = split.vec)
      for(i in 1:ncol(BVAR.fit$nr_pcs_selected)){
        # add other PCs
        PCs.sim = BVAR.fit$EOFs[[i]]$sp.PCA$x
        PCs.sim[,1:ncol(spl[[i]])] = spl[[i]]
        # build spatial field
        Mat.sp.sim = matMult_rcpp(PCs.sim,t(BVAR.fit$EOFs[[i]]$sp.PCA$rotation))
        Mat.sp.sim = scale(Mat.sp.sim,center = -BVAR.fit$EOFs[[i]]$sp.PCA$center,scale = F)
        
        if(return.raster == F){
          vars.sim_rst[[i]][[n]] = Mat.sp.sim
          cat(names(BVAR.fit$nr_pcs_selected)[i],"|",paste("iteration:",n),"\n")
        } else {
          # get coordinates and create raster
          ll = length(BVAR.fit$EOFs[[i]]$raster)
          if(class(BVAR.fit$EOFs[[i]]$raster) == "list" & ll > 1){
            rst.sim = vector("list",ll)
            spl.index = rep(1:ll,each = ncol(Mat.sp.sim)/ll)
            Mat.sp.sim_split = split.by.col(Mat.sp.sim,spl.index)
            for(ii in 1:ll){
              # use terra as it is faster than raster
              tmp.mat = terra::as.data.frame(BVAR.fit$EOFs[[i]]$raster[[ii]])
              coords = coordinates(BVAR.fit$EOFs[[i]]$raster[[ii]])
              coords = coords[complete.cases(tmp.mat),]
              #rst.sim[[ii]] = rasterFromXYZ(cbind(coords,t(Mat.sp.sim_split[[ii]])))
              rst.sim[[ii]] = terra::rast(cbind(coords,t(Mat.sp.sim_split[[ii]])), 
                                          type = "xyz",crs = "EPSG:4326")
              names(rst.sim[[ii]]) = BVAR.fit$EOFs[[i]]$time
            }
          } else {
            tmp.mat = terra::as.data.frame(BVAR.fit$EOFs[[i]]$raster)
            coords = coordinates(BVAR.fit$EOFs[[i]]$raster)
            coords = coords[complete.cases(tmp.mat),]
            #rst.sim = rasterFromXYZ(cbind(coords,t(Mat.sp.sim)))
            rst.sim = terra::rast(cbind(coords,t(Mat.sp.sim)),
                                  type = "xyz",crs = "EPSG:4326")
            names(rst.sim) = BVAR.fit$EOFs[[i]]$time
          }
          cat(names(BVAR.fit$nr_pcs_selected)[i],"|",paste("iteration:",n),"\n")
          # store in nested list
          vars.sim_rst[[i]][[n]] = rst.sim
        }
      }
    }
    names(vars.sim_rst) = names(BVAR.fit$nr_pcs_selected)
  }
  # 3. write out/compress
  # return
  out.list = list(type = ifelse(return.raster == TRUE,"raster","matrix"),
                  sims = vars.sim_rst)
  return(out.list)
}




