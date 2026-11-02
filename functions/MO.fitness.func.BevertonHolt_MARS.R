# ==================================================== #
# Multiobjective fitness function for BevertonHolt
# combined with a MARS model
# ==================================================== #

MO.fitness.func.BevertonHolt_MARS = function(feat,x,y,
                                       seed,
                                       metric.func,
                                       pmethod = "none",
                                       interaction.degree = 1,
                                       fast.k = 10,
                                       reps = 20,
                                       folds = 5,
                                       holdout = list(include = FALSE,
                                                      size = NA), # option to include a holdout set directly 
                                       # in the search and track its evolution
                                       Nrestarts = 10,
                                       force.SSB = T,
                                       write.logs = F,
                                       sink.folder = "./tmp",
                                       sink.file.name = "train.test.out.csv",
                                       only.obs.vs.pred = FALSE,
                                       return.scaled.metric = TRUE,...){
  
  library(foreach)
  library(optimx)
  library(minqa)
  
  if(exists(".Random.seed")){
    old.seed = .Random.seed
  } else{
    set.seed(NULL)
    old.seed = .Random.seed
  }
  on.exit(assign(".Random.seed",value = old.seed,envir = .GlobalEnv))
  
  if(write.logs == T){
    
    if(!dir.exists(sink.folder)){
     dir.create(sink.folder)
    }
  # create log-file
    if(!file.exists(file.path(sink.folder,sink.file.name))){
      file.create(file.path(sink.folder,sink.file.name))
      df = data.frame(vars = NA,
                    train = NA,
                    test = NA,
                    holdout = NA)
      df = df[-1,]
      write.table(df,file = file.path(sink.folder,sink.file.name),
                append = T, sep = ";",row.names = F,col.names = T)
    }
  }
    
  if(force.SSB == TRUE){
    indx.ssb = grep("ssb",colnames(x))
    feat[indx.ssb] = 1
  }
  
  # ----------------- #
  # fit BevertonHolt model
  # ----------------- # 
  
  bevertonHolt = function(par,x){
    par[["a"]]*x/(x + par[["b"]])
    #par[["a"]]/(1 + par[["b"]]/x) # reparameterised version
  }
  
  opt.BH = function(par,y,x){
    out = bevertonHolt(par,x)
    return(sqrt(sum(log(out) - log(y))^2))
  }
  
  
  if(force.SSB == TRUE){
    # select subset, but without ssb
    x.subset = x[,which(feat == 1)[-indx.ssb],drop = FALSE]
  }else{
    x.subset = x[,which(feat == 1),drop = FALSE]
  }
  
  if(holdout$include == TRUE & !is.na(holdout$size)){
    index.holdout = (1+length(y)-holdout$size):length(y)
    
    # holdout
    x.holdout = x.subset[index.holdout,,drop = F]
    y.holdout = y[index.holdout]
    data.holdout = data.frame(rec = y.holdout,
                              ssb = x[index.holdout,grep("ssb",names(x))])
    
    # train + test
    y.new = y[-index.holdout]
    x.subset.new = x.subset[-index.holdout,,drop = F]
    data = data.frame(rec = y.new,
                      ssb = x[-index.holdout,grep("ssb",names(x))])
  } else if(holdout$include == TRUE & is.na(holdout$size)){
    stop("Please specify size of the holdout set! \n Second argument of 'holdout' needs to be numeric!")
  } else if(holdout$include == FALSE){
    data = data.frame(rec = y,
                      ssb = x[,grep("ssb",names(x))])
    y.new = y
    x.subset.new = x.subset
  }
  
  set.seed(seed)
  cc.folds = caret::createMultiFolds(y = y.new,k = folds,times = reps)  
  
  # calculate fraction of features
  frac.feat = sum(feat)/length(feat)
  # control for empty set
  if(frac.feat == 0){
    return(c(metric = Inf,params =  frac.feat))
  }
  
  # perform cv
  data.list = list()
  for(k in 1:interaction.degree){ # choose degree of interaction
    i = 1:length(cc.folds)
    out = foreach::foreach(i = i,.packages = c("caret","earth","pracma","optimx","minqa"),.multicombine = TRUE) %do% {
      
      # ------------------------- #
      # 1. fit BevertonHolt model
      # ------------------------- #
      
      # define starting parameters (follow the formulation done in FLR)
      a <- max(quantile(data$rec[cc.folds[[i]]], 0.75, na.rm = TRUE))
      b <- a/max(quantile(data$rec[cc.folds[[i]]]/data$ssb[cc.folds[[i]]], 0.9, na.rm = TRUE))
      
      # sample for restarts
      a = rnorm(Nrestarts,a,sd = 0.1*a)
      b = rnorm(Nrestarts,b,sd = 0.1*b)
      
      # # define starting parameters
      # a = rnorm(Nrestarts,mean(data$rec[cc.folds[[i]]]),sd = 0.1*mean(data$rec[cc.folds[[i]]]))
      # b = rnorm(Nrestarts,min(data$ssb[cc.folds[[i]]]),sd = 0.1*min(data$ssb[cc.folds[[i]]]))
      
      par_mat = cbind(a = a,b = b)
      
      # multistarting
      BH <- multistart(parmat =  par_mat,
                       fn = opt.BH,
                       y = data$rec[cc.folds[[i]]],
                       x = data$ssb[cc.folds[[i]]],
                       method = "bobyqa",
                       lower = c(a = 0, b = 0.01*min(data$ssb[cc.folds[[i]]])), # constraint to avoid too steep solutions
                       upper = c(a = Inf,b = Inf), control=list(trace=0))
      # get only solutions that converged
      pars = BH[BH$convergence == 0,c("a","b")]
      # use the median of the parameter evaluation
      pars = apply(pars,2,median,na.rm = T)
      
      if(pracma::isempty(x.subset)){
        # predicted (in original scale)
        pred.BevertonHolt.test = bevertonHolt(par = pars,x = data$ssb[-cc.folds[[i]]]) 
        pred.BevertonHolt.train = bevertonHolt(par = pars,x = data$ssb[cc.folds[[i]]])
        
        obs.test = y.new[-cc.folds[[i]]]
        obs.train = y.new[cc.folds[[i]]]
        
        out.test = data.frame(pred = as.numeric(pred.BevertonHolt.test),obs = obs.test,fold.index = i)
        out.train = data.frame(pred = as.numeric(pred.BevertonHolt.train),obs = obs.train,fold.index = i)
        
        # combine
        out = list(train = out.train,test = out.test)
      } else {
        
        # ------------------------------ #
        # 2. get residuals to fit MARS
        # ------------------------------ #
        
        # get residuals
        resids.rec = log(data$rec[cc.folds[[i]]])-log(bevertonHolt(par = pars,x = data$ssb[cc.folds[[i]]]))
        
        
        # create train & test
        y.train = resids.rec
        x.train = x.subset.new[cc.folds[[i]],,drop = FALSE]
        x.test = x.subset.new[-cc.folds[[i]],,drop = FALSE]
        
        xy.train = data.frame(y.train = y.train,x.train)
        
        # MARS model with earth pkg
        MARS.out = earth::earth(y.train~.,data = xy.train,fast.k = fast.k,pmethod = pmethod,degree = k,...)
        
        # if a folds yields an intercept-only model, return the original SRR model
        if(length(MARS.out$coefficients) == 1){
          pred.BevertonHolt.test = bevertonHolt(par = pars,x = data$ssb[-cc.folds[[i]]]) 
          pred.BevertonHolt.train = bevertonHolt(par = pars,x = data$ssb[cc.folds[[i]]]) 
        } else{
          # predicted (in original scale)
          pred.BevertonHolt.test = exp(log(bevertonHolt(par = pars,x = data$ssb[-cc.folds[[i]]])) + 
                                         as.numeric(stats::predict(MARS.out,x.test)))
          pred.BevertonHolt.train = exp(log(bevertonHolt(par = pars,x = data$ssb[cc.folds[[i]]])) + 
                                         as.numeric(stats::predict(MARS.out,x.train)))
        }
        obs.test = y.new[-cc.folds[[i]]]
        obs.train = y.new[cc.folds[[i]]]
        
        out.test = data.frame(pred = as.numeric(pred.BevertonHolt.test),obs = obs.test,fold.index = i)
        out.train = data.frame(pred = as.numeric(pred.BevertonHolt.train),obs = obs.train,fold.index = i)
        
        out = list(train = out.train,test = out.test)
      }
      out
    }
    
    if(holdout$include == TRUE){
      # calculate a model for the whole dataset & evaluate on the holdout set
      # fit BH
      # ------------------------- #
      # 1. fit BevertonHolt model
      # ------------------------- #
      
      # define starting parameters (follow the formulation done in FLR)
      a <- max(quantile(data$rec, 0.75, na.rm = TRUE))
      b <- a/max(quantile(data$rec/data$ssb, 0.9, na.rm = TRUE))
      
      # sample for restarts
      a = rnorm(Nrestarts,a,sd = 0.1*a)
      b = rnorm(Nrestarts,b,sd = 0.1*b)
      
      # # define starting parameters
      # a = rnorm(Nrestarts,mean(data$rec),sd = 0.1*mean(data$rec))
      # b = rnorm(Nrestarts,min(data$ssb),sd = 0.1*min(data$ssb))
      
      par_mat = cbind(a = a,b = b)
      
      # multistarting
      BH <- multistart(parmat =  par_mat,
                       fn = opt.BH,
                       y = data$rec,
                       x = data$ssb,
                       method = "bobyqa",
                       lower = c(a = 0, b = 0.01*min(data$ssb)), # constraint to avoid too steep solutions
                       upper = c(a = Inf,b = Inf), control=list(trace=0))
      # get only solutions that converged
      pars = BH[BH$convergence == 0,c("a","b")]
      # use the median of the parameter evaluation
      pars = apply(pars,2,median,na.rm = T)
      
      if(pracma::isempty(x.subset.new)){
        # predicted (in original scale)
        pred.BevertonHolt.holdout = bevertonHolt(par = pars,x = data.holdout$ssb) 
      } else{
      # get residuals
      resids.rec = log(data$rec)-log(bevertonHolt(par = pars,x = data$ssb))
      
      xy = data.frame(y = resids.rec,x.subset.new)
      
      # MARS model with earth pkg
      MARS.out = earth::earth(y~.,data = xy,fast.k = fast.k,pmethod = pmethod,degree = k,...)
      
      # evaluate
      # if it yields an intercept-only model, return the original SRR model
      if(length(MARS.out$coefficients) == 1){
        pred.BevertonHolt.holdout = bevertonHolt(par = pars,x = data.holdout$ssb) 
      } else{
        # predicted (in original scale)
        pred.BevertonHolt.holdout = exp(log(bevertonHolt(par = pars,x = data.holdout$ssb)) + 
                                          as.numeric(stats::predict(MARS.out,x.holdout)))
        }
      }
      out.holdout = data.frame(pred = as.numeric(pred.BevertonHolt.holdout),
                               obs = y.holdout)
    } else{
      out.holdout = data.frame(pred = NA,obs = NA)
    }
    
    # combine all obs and pred to one
    data.test = data.table::rbindlist(lapply(out,function(x) x$test))
    data.train = data.table::rbindlist(lapply(out,function(x) x$train))
    
    data.list[[k]] = list(test = data.test,
                          train = data.train,
                          holdout = out.holdout)
  }
  # option to only return the predicted vs. observed CV results
  if(only.obs.vs.pred == TRUE | is.null(metric.func)){
    return(data.list)
    }
  # calculate metric over all folds and repetitions (to minimize the metric)
  metric.test = sapply(data.list,function(x) metric.func(x$test))
  metric.train = sapply(data.list,function(x) metric.func(x$train))
  metric.holdout = sapply(data.list,function(x) metric.func(x$holdout))
  
  # choose the lowest
  test.metric = min(metric.test)
  index.test.metric = which(test.metric == min(metric.test))
  train.metric = metric.train[index.test.metric]
  holdout.metric = metric.holdout[index.test.metric]
  
  # write out to disk
  if(write.logs == TRUE){
  index.feat = paste(which(feat == 1),collapse = "~")
  df = data.frame(vars = index.feat,
                  train= train.metric,
                  test = test.metric,
                  holdout = holdout.metric)
  write.table(df,file = file.path(sink.folder,sink.file.name),
              append = T, sep = ";",row.names = F,col.names = F)
  }
  
  # mean
  if(return.scaled.metric == TRUE){
    # return a metric scaled between 0 and 1
    naive.metric = metric.func(data = data.frame(pred = mean(y,na.rm = TRUE),obs = y))
    rel.metric = test.metric/naive.metric
    # control for predictions, which are worse than the mean
    rel.metric = ifelse(rel.metric>1,1,rel.metric)
  } else{
    rel.metric = test.metric
  }
  return.out = c(metric = rel.metric,params =  frac.feat)
  return(return.out)
}

# ------------------------------------------------------------ #
# env. Hockeystick + MARS for env. Vars

fit.BevertonHolt.MARS = function(rec,ssb,env.Vars,Nrestarts = 10,...){
  require(optimx)
  require(minqa)
  # df that holds rec & ssb
  data = data.frame(rec = rec,
                    ssb = ssb)
  
  # ------------------------ #
  # 1. fit BevertonHolt model
  # ------------------------ #
  
  bevertonHolt = function(par,x){
    par[["a"]]*x/(x + par[["b"]])
    #par[["a"]]/(1 + par[["b"]]/x) # reparameterised version
  }
  
  opt.BH = function(par,y,x){
    out = bevertonHolt(par,x)
    return(sqrt(sum(log(out) - log(y))^2))
  }

  
  # fit BevertonHolt
  # define starting parameters (follow the formulation done in FLR)
  a <- max(quantile(data$rec, 0.75, na.rm = TRUE))
  b <- a/max(quantile(data$rec/data$ssb, 0.9, na.rm = TRUE))
  
  # sample for restarts
  a = rnorm(Nrestarts,a,sd = 0.1*a)
  b = rnorm(Nrestarts,b,sd = 0.1*b)
  
  # # define starting parameters
  # a = rnorm(Nrestarts,mean(data$rec),sd = 0.1*mean(data$rec))
  # b = rnorm(Nrestarts,min(data$ssb),sd = 0.1*min(data$ssb))
  
  par_mat = cbind(a = a,b = b)
  
  # multistarting
  BH <- multistart(parmat =  par_mat,
                   fn = opt.BH,
                   y = data$rec,
                   x = data$ssb,
                   method = "bobyqa",
                   lower = c(a = 0, b = 0.01*min(data$ssb)), # constraint to avoid too steep solutions
                   upper = c(a = Inf,b = Inf), control=list(trace=0))
  # get only solutions that converged
  pars = BH[BH$convergence == 0,c("a","b")]
  # use the median of the parameter evaluation
  pars = apply(pars,2,median,na.rm = T)
  
  # -------------------------------------- #
  # 2. get residuals to fit env. Model to
  # -------------------------------------- #
  
  # get residuals
  resids.rec = log(data$rec)-log(bevertonHolt(par = pars,x = data$ssb))
  
  if(pracma::isempty(env.Vars) | all(is.na(env.Vars)) | is.null(env.Vars)){
    predicted = bevertonHolt(par = pars,x = data$ssb)
    residual.model = NA
    
  } else{
    # ----------------------- #
    # 3. fit MARS to residuals
    # ----------------------- #
    data.resids = data.frame(resids = resids.rec,env.Vars)
    residual.model = earth::earth(resids~.,data = data.resids,...)
    
    # check if residual model is an intercept only model
    if(length(residual.model$coefficients) == 1){
      # if it is only an intercept only model return just the fitted BevertonHolt
      predicted = bevertonHolt(par = pars,x = data$ssb)
    } else{
      predicted = exp(log(bevertonHolt(par = pars,x = data$ssb)) + as.numeric(predict(residual.model)))
      
    }
    
    
  }
  # return stuff in a list
  return.obj = list(predicted = predicted,
                    optimx.output = BH,
                    BevertonHolt.model = list(model = bevertonHolt,
                                              params = pars),
                    MARS.model = residual.model
  )
  
  return(return.obj)
}

# predict function for "fitted model obj"
predict.BevertonHolt.MARS = function(model.obj,newdata){
  library(earth)
  names(newdata)[grep("ssb",names(newdata))] <- "ssb" # rename ssb to match with BevertonHolt model
  BevertonHolt.mod = model.obj$BevertonHolt.model$model
  BH.pars = model.obj$BevertonHolt.model$params
  # check if MARS model is empty or an intercept-only model
  if(all(is.na(model.obj$MARS.model)) || length(model.obj$MARS.model$coefficients) == 1){
    predict.MARSOut = 0
  } else{
    predict.MARSOut = as.numeric(predict(model.obj$MARS.model,newdata))
  }
  predict.all = exp(log(BevertonHolt.mod(par = BH.pars,x = newdata$ssb))+predict.MARSOut)
  
  return(predict.all)
}




