# ==================================================== #
# Multiobjective fitness function for env. Hockeystick
# combined with a MARS model
# ==================================================== #

MO.fitness.func.envHockeyStick_MARS = function(feat,x,y,
                                               seed,
                                               metric.func,
                                               pmethod = "none",
                                               interaction.degree = 1,
                                               fast.k = 10,
                                               reps = 20,
                                               folds = 5,
                                               Nrestars = 10,
                                               holdout = list(include = FALSE,
                                                              size = NA), # option to include a holdout set directly 
                                               # in the search and track its evolution
                                               force.SSB = T,
                                               write.logs = FALSE,
                                               sink.folder = "./tmp",
                                               sink.file.name = "train.test.out.csv",
                                               only.obs.vs.pred = FALSE,
                                               return.scaled.metric = TRUE,...){
  
  library(foreach)
  library(optimx)
  library(minqa)
  library(SimDesign)
  
  if(exists(".Random.seed")){
    old.seed = .Random.seed
  } else{
    set.seed(NULL)
    old.seed = .Random.seed
  }
  on.exit(assign(".Random.seed",value = old.seed,envir = .GlobalEnv))
  
  if(write.logs == TRUE){
    if(!dir.exists(sink.folder)) {
      dir.create(sink.folder)
    }
    if (!file.exists(file.path(sink.folder, sink.file.name))) {
      file.create(file.path(sink.folder, sink.file.name))
      df = data.frame(vars = NA,
                      train = NA,
                      test = NA,
                      holdout = NA)
      df = df[-1, ]
      write.table(df,
                  file = file.path(sink.folder, sink.file.name),
                  append = T,sep = ";",row.names = F,col.names = T)
    }
  }
  
  if(force.SSB == TRUE){
    indx.ssb = grep("ssb",colnames(x))
    feat[indx.ssb] = 1
  }
  
  # ---------------------------- #
  # fit global Hockeystick model
  # ---------------------------- # 
  # Hockey stick model
  HS.recruitment = function(par,ssb){
    return(ifelse(c(ssb) <= par[["b"]], par[["a"]] * c(ssb), par[["a"]] * par[["b"]]))
  }
  opt.HS = function(par,y,x){
    out = HS.recruitment(par,x)
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
    out = foreach::foreach(i = i,.packages = c("caret","earth","pracma","optimx","minqa","SimDesign"),.multicombine = TRUE) %do% {
      
      # ------------------------ #
      # 1. fit Hockeystick model
      # ------------------------ #
      
      # MLE <- optim(fn= HS.likelihood,
      #              data = data[cc.folds[[i]],],par=start.params,method="Nelder-Mead")
  
      a = rnorm(Nrestarts,mean(data$rec[cc.folds[[i]]])/(2.5*min(data$ssb[cc.folds[[i]]])),sd = 0.1*mean(data$rec[cc.folds[[i]]])/(2.5*min(data$ssb[cc.folds[[i]]])))
      b = seq(min(data$ssb[cc.folds[[i]]]),0.9*max(data$ssb[cc.folds[[i]]]),length.out = Nrestarts)
      
      par_mat = expand.grid(a = a,b = b)
      #par_mat = cbind(a = a,b = b)
      
      par_mat$fitness = apply(par_mat,1,function(x) opt.HS(par = x,y = data$rec[cc.folds[[i]]],x = data$ssb[cc.folds[[i]]]))
      #plot(par_mat$a,par_mat$fitness)
      #plot(par_mat$b,par_mat$fitness)
      # get Nrestarts number of best params
      par_mat = as.matrix(par_mat[order(par_mat$fitness),][1:pmin(Nrestarts,pmin(Nrestarts,10)),c("a","b")])
      # sample around these good starting values
      par_mat[,"a"] = par_mat[,"a"] + rnorm(10,0,0.1*min(par_mat[,"a"]))
      par_mat[,"b"] = par_mat[,"b"] + rnorm(10,0,0.1*min(par_mat[,"b"]))
      par_mat[,"b"] = ifelse(par_mat[,"b"]< min(data$ssb[cc.folds[[i]]]),min(data$ssb[cc.folds[[i]]]),par_mat[,"b"])
      
      HS <- SimDesign::quiet(multistart(parmat =  par_mat,
                       fn = opt.HS,
                       y = data$rec[cc.folds[[i]]],
                       x = data$ssb[cc.folds[[i]]],
                       method = "bobyqa",
                       lower = c(a = 0, b = min(data$ssb[cc.folds[[i]]])), # constraint to avoid too steep solutions
                       upper = c(a = Inf,b = Inf), control=list(trace=0)))
    
      
      pars = data.frame(HS[which(HS$value< quantile(HS$value,0.25)),c("a","b")])
      pars = apply(pars,2,median,na.rm = T)

      if(pracma::isempty(x.subset.new)){
        # predicted (in original scale)
        pred.HS.test = HS.recruitment(par = pars,ssb = data$ssb[-cc.folds[[i]]])
        pred.HS.train = HS.recruitment(par = pars,ssb = data$ssb[cc.folds[[i]]])
        
        obs.test = y.new[-cc.folds[[i]]]
        obs.train = y.new[cc.folds[[i]]]
        out.test = data.frame(pred = as.numeric(pred.HS.test),obs = obs.test,fold.index = i)
        out.train = data.frame(pred = as.numeric(pred.HS.train),obs = obs.train,fold.index = i)
        # combine
        out = list(train = out.train,test = out.test)
      } else {
        
        # ------------------------------ #
        # 2. get residuals to fit MARS
        # ------------------------------ #
        
        # get residuals
        resids.rec = log(data$rec[cc.folds[[i]]])-log(HS.recruitment(pars,ssb = data$ssb[cc.folds[[i]]]))
        
        
        # create train & test
        y.train = resids.rec
        x.train = x.subset.new[cc.folds[[i]],,drop = FALSE]
        x.test = x.subset.new[-cc.folds[[i]],,drop = FALSE]
        
        xy.train = data.frame(y.train = y.train,x.train)
        
        # MARS model with earth pkg
        MARS.out = earth::earth(y.train~.,data = xy.train,fast.k = fast.k,pmethod = pmethod,degree = k,...)
        
        # if a folds yields an intercept-only model, return the original SRR model
        if(length(MARS.out$coefficients) == 1){
          pred.HS.test = HS.recruitment(pars,ssb = data$ssb[-cc.folds[[i]]])
          pred.HS.train = HS.recruitment(pars,ssb = data$ssb[cc.folds[[i]]])
        } else{
        # predicted (in original scale)
        pred.HS.test = exp(log(HS.recruitment(pars,ssb = data$ssb[-cc.folds[[i]]])) +
                             as.numeric(stats::predict(MARS.out,x.test)))
        pred.HS.train = exp(log(HS.recruitment(pars,ssb = data$ssb[cc.folds[[i]]])) +
                             as.numeric(stats::predict(MARS.out,x.train)))
        }
        obs.test = y.new[-cc.folds[[i]]]
        obs.train = y.new[cc.folds[[i]]]
        
        out.test = data.frame(pred = as.numeric(pred.HS.test),obs = obs.test,fold.index = i)
        out.train = data.frame(pred = as.numeric(pred.HS.train),obs = obs.train,fold.index = i)
        
        out = list(train = out.train,test = out.test)
      }
      out
    }
    
    if(holdout$include == TRUE){
      # calculate a model for the whole dataset & evaluate on the holdout set
      # fit HS
      a = rnorm(Nrestarts,mean(data$rec)/(2.5*min(data$ssb)),sd = 0.1*mean(data$rec)/(2.5*min(data$ssb)))
      b = seq(min(data$ssb),0.9*max(data$ssb),length.out = Nrestarts)
      
      par_mat = expand.grid(a = a,b = b)
      #par_mat = cbind(a = a,b = b)
      
      par_mat$fitness = apply(par_mat,1,function(x) opt.HS(par = x,y = data$rec,x = data$ssb))
      # get Nrestarts number of best params
      par_mat = as.matrix(par_mat[order(par_mat$fitness),][1:pmin(Nrestarts,pmin(Nrestarts,10)),c("a","b")])
      # sample around these good starting values
      par_mat[,"a"] = par_mat[,"a"] + rnorm(10,0,0.1*min(par_mat[,"a"]))
      par_mat[,"b"] = par_mat[,"b"] + rnorm(10,0,0.1*min(par_mat[,"b"]))
      par_mat[,"b"] = ifelse(par_mat[,"b"]< min(data$ssb),min(data$ssb),par_mat[,"b"])
      
      HS <- SimDesign::quiet(multistart(parmat =  par_mat,
                                        fn = opt.HS,
                                        y = data$rec,
                                        x = data$ssb,
                                        method = "bobyqa",
                                        lower = c(a = 0, b = min(data$ssb)), # constraint to avoid too steep solutions
                                        upper = c(a = Inf,b = Inf), control=list(trace=0)))
      
      
      pars = data.frame(HS[which(HS$value< quantile(HS$value,0.25)),c("a","b")])
      pars = apply(pars,2,median,na.rm = T)
      # get residuals
      resids.rec = log(data$rec)-log(HS.recruitment(pars,ssb = data$ssb))
      xy = data.frame(y = resids.rec,x.subset.new)
      
      if(pracma::isempty(x.subset.new)){
        # predicted (in original scale)
        pred.HS.holdout = HS.recruitment(par = pars,ssb = data.holdout$ssb)
        out.holdout = data.frame(pred = as.numeric(pred.HS.holdout),
                                 obs = y.holdout)
      } else {
        # MARS model with earth pkg
        MARS.out = earth::earth(y~.,data = xy,fast.k = fast.k,pmethod = pmethod,degree = k,...)
        
        # evaluate
        # if it yields an intercept-only model, return the original SRR model
        if(length(MARS.out$coefficients) == 1){
          pred.HS.holdout = HS.recruitment(par = pars,ssb = data.holdout$ssb)
        } else{
          # predicted (in original scale)
          pred.HS.holdout = exp(log(HS.recruitment(pars,ssb = data.holdout$ssb)) +
                               as.numeric(stats::predict(MARS.out,x.holdout)))
        }
        out.holdout = data.frame(pred = as.numeric(pred.HS.holdout),
                                 obs = y.holdout)
      }
    } else {
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

fit.env.Hockeystick.MARS = function(rec,ssb,Nrestarts,env.Vars,...){
  
  library(optimx)
  library(minqa)
  library(SimDesign)
  
  # df that holds rec & ssb
  data = data.frame(rec = rec,
                    ssb = ssb)
  #params = start.params
  
  # Hockey stick model
  HS.recruitment = function(par,ssb){
    return(ifelse(c(ssb) <= par[["b"]], par[["a"]] * c(ssb), par[["a"]] * par[["b"]]))
  }

  opt.HS = function(par,y,x){
    out = HS.recruitment(par,x)
    return(sqrt(sum((log(out) - log(y))^2)))
  }

  # ------------------------ #
  # 1. fit Hockeystick model
  # ------------------------ #
  
  # define starting parameters
  a = rnorm(Nrestarts,mean(data$rec)/(2.5*min(data$ssb)),sd = 0.1*mean(data$rec)/(2.5*min(data$ssb)))
  b = seq(min(data$ssb),0.9*max(data$ssb),length.out = Nrestarts)
  
  par_mat = expand.grid(a = a,b = b)
  #par_mat = cbind(a = a,b = b)
  
  par_mat$fitness = apply(par_mat,1,function(x) opt.HS(par = x,y = data$rec,x = data$ssb))
  #plot(par_mat$a,par_mat$fitness)
  #plot(par_mat$b,par_mat$fitness)
  # get Nrestarts number of best params
  par_mat = as.matrix(par_mat[order(par_mat$fitness),][1:pmin(Nrestarts,pmin(Nrestarts,10)),c("a","b")])
  # sample around these good starting values
  par_mat[,"a"] = par_mat[,"a"] + rnorm(10,0,0.1*min(par_mat[,"a"]))
  par_mat[,"b"] = par_mat[,"b"] + rnorm(10,0,0.1*min(par_mat[,"b"]))
  par_mat[,"b"] = ifelse(par_mat[,"b"]< min(data$ssb),min(data$ssb),par_mat[,"b"])
  
  HS <- SimDesign::quiet(multistart(parmat =  par_mat,
                                    fn = opt.HS,
                                    y = data$rec,
                                    x = data$ssb,
                                    method = "bobyqa",
                                    lower = c(a = 0, b = min(data$ssb)), # constraint to avoid too steep solutions
                                    upper = c(a = Inf,b = Inf), control=list(trace=0)))
  
  
  pars = data.frame(HS[which(HS$value< quantile(HS$value,0.25)),c("a","b")])
  pars = apply(pars,2,median,na.rm = T)
  
  # plot(data$ssb,data$rec)
  # for(i in 1:nrow(pars)){
  #   lines(seq(0,3E6,length.out = 100),HS.recruitment(par = pars[i,],ssb = seq(0,3E6,length.out = 100)))
  # }
  # 
 
  # lines(seq(0,3E6,length.out = 100),HS.recruitment(par = pars,ssb = seq(0,3E6,length.out = 100)),
  #       col = "red",lwd = 2)
  # ------------------------------ #
  # 2. get residuals to fit glm to
  # ------------------------------ #
  
  # get residuals
  resids.rec = log(data$rec)-log(HS.recruitment(pars,ssb = data$ssb))
  
  if(pracma::isempty(env.Vars) | all(is.na(env.Vars)) | is.null(env.Vars)){
    predicted = HS.recruitment(pars,ssb = data$ssb)
    residual.model = NA
    
  } else{
  # ----------------------- #
  # 3. fit MARS to residuals
  # ----------------------- #
  data.resids = data.frame(resids = resids.rec,env.Vars)
  residual.model = earth::earth(resids~.,data = data.resids,...)
  
  
  # check if residual model is an intercept only model
  if(length(residual.model$coefficients) == 1){
    # if it is only an intercept only model return just the fitted Ricker
    predicted = HS.recruitment(pars,ssb = data$ssb)
      } else {
    predicted = exp(log(HS.recruitment(pars,ssb = data$ssb)) +
                      as.numeric(predict(residual.model)))
  
    }
  }
  # return stuff in a list
  return.obj = list(predicted = predicted,
                    HS.params = pars,
                    optimx.output = HS,
                    HockeyStick.model = list(model = HS.recruitment,
                                              params = pars),
                    MARS.model = residual.model)
  
  return(return.obj)
}

fit.env.Hockeystick.bagMARS = function(rec,ssb,start.params,env.Vars,...){
  
  # df that holds rec & ssb
  data = data.frame(rec = rec,
                    ssb = ssb)
  params = start.params
  
  # Hockey stick model
  HS.recruitment = function(ssb,a,b){
    return(ifelse(c(ssb) <= b, a * c(ssb), a * b))
  }
  
  # likelihood function
  HS.likelihood = function(params,data){
    a = params[["a"]]
    b = params[["b"]]
    ssb.name = names(data)[grep("ssb",names(data))]
    ssb = data[,ssb.name]
    rec = log(data$rec)
    
    loglAR1.flcore = function(obs, hat, rho = 0) {
      # calculates likelihood for AR(1) process
      n <- length(obs)
      rsdl <- (obs[-1] - rho * obs[-n] - hat[-1] + rho * hat[-n])
      s2 <- sum(rsdl^2, na.rm = T)
      s1 <- s2
      if (!all(is.na(rsdl[1]))) 
        s1 <- s1 + (1 - rho^2) * (obs[1] - hat[1])^2
      sigma2 <- sum((obs - hat)^2)
      n <- length(obs[!is.na(obs)])
      sigma2.a <- (1 - rho^2) * sigma2
      res <- (log(1/(2 * pi)) - n * log(sigma2.a) + log(1 - rho^2) - s1/(2 * sigma2.a))/2
      if (!is.finite(res)) res <- -1e+100
      return(res)
    }
    
    hat = log(ifelse(c(ssb) <= b, a * c(ssb), a * b))
    return(-loglAR1.flcore(obs = rec,hat = hat,rho = 0))
    #sum((rec - ifelse(c(ssb) <= b, a * c(ssb), a * b))^2)
  }
  
  # ------------------------ #
  # 1. fit Hockeystick model
  # ------------------------ #
  
  MLE <- optim(fn=HS.likelihood,
               data = data,par=params,method="Nelder-Mead")
  
  # ------------------------------ #
  # 2. get residuals to fit glm to
  # ------------------------------ #
  
  # get residuals
  resids.rec = log(data$rec)-log(HS.recruitment(ssb = data$ssb,
                                                a = MLE$par[["a"]],
                                                b = MLE$par[["b"]]))
  
  if(pracma::isempty(env.Vars) | all(is.na(env.Vars)) | is.null(env.Vars)){
    predicted = HS.recruitment(ssb = data$ssb,
                               a = MLE$par[["a"]],
                               b = MLE$par[["b"]])
    residual.model = NA
    
  } else{
    # ----------------------- #
    # 3. fit MARS to residuals
    # ----------------------- #
    data.resids = data.frame(resids = resids.rec,env.Vars)
    residual.model = caret::bagEarth(resids~.,data = data.resids,...)
    
    predicted = exp(log(HS.recruitment(ssb = data$ssb,
                                       a = MLE$par[["a"]],
                                       b = MLE$par[["b"]])) + as.numeric(predict(residual.model)))
    
  }
  # return stuff in a list
  return.obj = list(predicted = predicted,
                    HS.params = c(a = MLE$par[["a"]],b = MLE$par[["b"]]),
                    MLE.obj = MLE,
                    MARS.model = residual.model
  )
  
  return(return.obj)
}



# predict function for "fitted model obj"
predict.HockeyStick.MARS = function(model.obj,newdata){
  library(earth)
  ssb = newdata[,grep("ssb",names(newdata))]
  a = model.obj$HS.params[["a"]]
  b = model.obj$HS.params[["b"]]
  predict.HS = ifelse(c(ssb) <= b, a * c(ssb), a * b)
  if(all(is.na(model.obj$MARS.model)) || length(model.obj$MARS.model$coefficients) == 1){
    predict.MARSOut = 0
  } else{
    predict.MARSOut = as.numeric(predict(model.obj$MARS.model,newdata))
  }
  predict.all = exp(log(predict.HS)+predict.MARSOut)
  
  return(predict.all)
}




