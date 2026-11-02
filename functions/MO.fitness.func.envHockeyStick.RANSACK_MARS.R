# ==================================================== #
# Multiobjective fitness function for env. Hockeystick
# combined with a MARS model
# ==================================================== #

MO.fitness.func.envHockeyStick.robust_MARS = function(feat,x,y,
                                               seed,
                                               metric.func,
                                               HS.g = 0.1, RANSAC.n = 100,
                                               RANSAC.thresh = NULL,
                                               pmethod = "none",
                                               interaction.degree = 1,
                                               fast.k = 10,
                                               reps = 20,
                                               folds = 5,
                                               holdout = list(include = FALSE,
                                                              size = NA), # option to include a holdout set directly 
                                               # in the search and track its evolution
                                               force.SSB = T,
                                               write.logs = F,
                                               sink.folder = "./tmp",
                                               sink.file.name = "train.test.out.csv",
                                               only.obs.vs.pred = FALSE,
                                               return.scaled.metric = TRUE,...){
  
  library(foreach)
  library(StockRecruitSET)
  
  if(exists(".Random.seed")){
    old.seed = .Random.seed
  } else{
    set.seed(NULL)
    old.seed = .Random.seed
  }
  on.exit(assign(".Random.seed",value = old.seed,envir = .GlobalEnv))
  
  if(write.logs == TRUE){
    if(!dir.exists(sink.folder)){
      dir.create(sink.folder)
    }
  
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
  
  # ---------------------------- #
  # fit global Hockeystick model
  # ---------------------------- # 
  # continous Hockey stick model
  contHS.recruitment = function(pars,S,g){
    pars[[1]]*(S + sqrt(pars[[2]]^2 + (g^2)/4) - sqrt((S - pars[[2]])^2 + (g^2)/4))
  }
  # function to fit a HS through 2 points
  fast.HS.fitting = function(dat){
    if(dim(dat)[1] != 2 & dim(dat)[2] != 2)
      stop("Input must be of the dimension 2x2")
    
    x = dat[,which(names(dat) %in% c("x","ssb","SSB","S"))]
    y = dat[,which(names(dat) %in% c("y","rec","R","Rec"))]
    
    mean.y = mean(y)
    x.min = min(x)
    
    p1 = c(x[x == x.min],y[x == x.min])
    # get slope between zero and lowest point
    slope = as.numeric(p1[2]/p1[1])
    # determine crossing point between line and mean
    x.break = mean.y/slope
    # check if x.min is larger than x.break
    x.break = ifelse(x.break > x.min,x.break,x.min)
    bp = c(x.break,mean.y)
    # translate to model parameters of continous HS fit
    a = 0.5*bp[2]/bp[1]
    b = bp[1]
    return(c(a = a,b = b))
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
  
  if(is.null(RANSAC.thresh)){
    nn.dist = abs(median(data$rec) - data$rec)
    RANSAC.thresh = median(nn.dist)
  }
  
  # perform cv
  data.list = list()
  for(k in 1:interaction.degree){ # choose degree of interaction
    i = 1:length(cc.folds)
    out = foreach::foreach(i = i,.packages = c("caret","earth","pracma","optimx","minqa","SimDesign"),.multicombine = TRUE) %do% {
      #for(i in 1:length(cc.folds)){
      data.train  = data[cc.folds[[i]],]
      fit.accepted  = list(); n = 1
      run = T; counter = 1
      max.time = pmax(500,5*RANSAC.n)
      while(run == T){
       # define a criterion that ensures convergence
       # but not so high runtimes...
       # run maximum 500 or 5*RANSAC.n times or collect RANSAC.n inlier solutions
        counter = counter + 1
        if(n == RANSAC.n | counter > max.time){
          run = F
        }
      #for(ii in 1:RANSAC.n){
        # 1. select minimum sample
        
        p2 = data.train[sample(1:nrow(data.train),size = 2),]
        
        # fast fitting
        pars1 = fast.HS.fitting(p2)
        
        # get residuals to minimal model
        y.head = contHS.recruitment(pars1[1:2],
                                S = data.train$ssb,
                                g = HS.g)
        
        a = abs(data.train$rec - y.head)
        inliers = a<RANSAC.thresh
        # accept
        if(sum(inliers)>= 0.5*nrow(data.train)){
          p.inlier = data.train[inliers,]
          
          # fast SR fit implementation approach would be
          out = fitSRCurve(S = p.inlier$ssb,R = p.inlier$rec,
                           shape = "contHockey",g = HS.g) 
          pars1 = out$env$last.par[1:2]
          
          # store accepted fits
          fit.accepted[[n]] = data.frame(a = pars1[1],
                                         b = pars1[2],
                                         no.inlier = sum(inliers))
          n = n+1
        } else {
          next
        }
      }
      fit.accepted.df = data.frame(do.call(rbind,fit.accepted))
      pars.best = fit.accepted.df[fit.accepted.df$no.inlier == max(fit.accepted.df$no.inlier),]
      if(length(dim(pars.best)) == 2){
        pars.best = apply(pars.best,2,median)[1:2]
      }
      
      # epmty set without any covariates
      if(pracma::isempty(x.subset)){
        
        # predicted (in original scale)
        pred.HS.test = contHS.recruitment(pars = exp(pars.best),
                                          S = data$ssb[-cc.folds[[i]]],g = HS.g)
        pred.HS.train = contHS.recruitment(pars = exp(pars.best),
                                          S = data$ssb[cc.folds[[i]]],g = HS.g)
        
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
        resids.rec = log(data$rec[cc.folds[[i]]])-log(contHS.recruitment(pars = exp(pars.best),
                                                                         S = data$ssb[cc.folds[[i]]],g = HS.g))
        
        
        # create train & test
        y.train = resids.rec
        x.train = x.subset.new[cc.folds[[i]],,drop = FALSE]
        x.test = x.subset.new[-cc.folds[[i]],,drop = FALSE]
        
        xy.train = data.frame(y.train = y.train,x.train)
        
        # MARS model with earth pkg
        MARS.out = earth::earth(y.train~.,data = xy.train,fast.k = fast.k,pmethod = pmethod,degree = k)#,...)
        
        # if a folds yields an intercept-only model, return the original SRR model
        if(length(MARS.out$coefficients) == 1){
          pred.HS.test = contHS.recruitment(pars = exp(pars.best),
                                            S = data$ssb[-cc.folds[[i]]],
                                            g = HS.g)
          pred.HS.train = contHS.recruitment(pars = exp(pars.best),
                                            S = data$ssb[cc.folds[[i]]],
                                            g = HS.g)
        } else{
        # predicted (in original scale)
        pred.HS.test = exp(log(contHS.recruitment(pars = exp(pars.best),
                                                  S = data$ssb[-cc.folds[[i]]],
                                                  g = HS.g)) +
                             as.numeric(stats::predict(MARS.out,x.test)))
        
        pred.HS.train = exp(log(contHS.recruitment(pars = exp(pars.best),
                                                  S = data$ssb[cc.folds[[i]]],
                                                  g = HS.g)) +
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
    
    # combine all obs and pred to one
    data.test = data.table::rbindlist(lapply(out,function(x) x$test))
    data.train = data.table::rbindlist(lapply(out,function(x) x$train))
    
    
    if(holdout$include == TRUE){
      # calculate a model for the whole dataset & evaluate on the holdout set
      # fit HS
      data.all  = data
      fit.accepted  = list(); n = 1
      run = T; counter = 1
      max.time = pmax(500,5*RANSAC.n)
      while(run == T){
        # define a criterion that ensures convergence
        # but not so high runtimes...
        # run maximum 500 or 5*RANSAC.n times or collect RANSAC.n inlier solutions
        counter = counter + 1
        if(n == RANSAC.n | counter > max.time){
          run = F
        }
        #for(ii in 1:RANSAC.n){
        # 1. select minimum sample
        
        p2 = data.all[sample(1:nrow(data.all),size = 2),]
        
        # fast fitting
        pars1 = fast.HS.fitting(p2)
        
        # get residuals to minimal model
        y.head = contHS.recruitment(pars1[1:2],
                                    S = data.all$ssb,
                                    g = HS.g)
        
        a = abs(data.all$rec - y.head)
        inliers = a<RANSAC.thresh
        # accept
        if(sum(inliers)>= 0.5*nrow(data.all)){
          p.inlier = data.all[inliers,]
          
          # fast SR fit implementation approach would be
          out = fitSRCurve(S = p.inlier$ssb,R = p.inlier$rec,
                           shape = "contHockey",g = HS.g) 
          pars1 = out$env$last.par[1:2]
          
          # store accepted fits
          fit.accepted[[n]] = data.frame(a = pars1[1],
                                         b = pars1[2],
                                         no.inlier = sum(inliers))
          n = n+1
        } else {
          next
        }
      }
      fit.accepted.df = data.frame(do.call(rbind,fit.accepted))
      pars.best = fit.accepted.df[fit.accepted.df$no.inlier == max(fit.accepted.df$no.inlier),]
      if(length(dim(pars.best)) == 2){
        pars.best = apply(pars.best,2,median)[1:2]
      }
      
      # epmty set without any covariates
      if(pracma::isempty(x.subset)){
        
      pred.HS.holdout = contHS.recruitment(pars = exp(pars.best),
                                        S = data.holdout$ssb,g = HS.g)
      
      } else{
      # get residuals
      resids.rec = log(data$rec)-log(contHS.recruitment(pars = exp(pars.best),
                                                        S = data$ssb,g = HS.g))
      xy = data.frame(y = resids.rec,x.subset.new)
      
      # MARS model with earth pkg
      MARS.out = earth::earth(y~.,data = xy,fast.k = fast.k,pmethod = pmethod,degree = k)#,...)
      
      # evaluate
      # if it yields an intercept-only model, return the original SRR model
      if(length(MARS.out$coefficients) == 1){
        pred.HS.holdout = contHS.recruitment(pars = exp(pars.best),
                                          S = data.holdout$ssb,
                                          g = HS.g)
      } else{
        # predicted (in original scale)
        pred.HS.holdout = exp(log(contHS.recruitment(pars = exp(pars.best),
                                                         S = data.holdout$ssb,
                                                         g = HS.g)) + as.numeric(stats::predict(MARS.out,x.holdout)))
        }
      }
    out.holdout = data.frame(pred = as.numeric(pred.HS.holdout),
                               obs = y.holdout)
    } else{
      out.holdout = data.frame(pred = NA,obs = NA)
    }
    
    
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

fit.env.HS.robust.MARS = function(rec,ssb,RANSAC.n,env.Vars,RANSAC.thresh = NULL,HS.g = 0.1,...){
  
  library(StockRecruitSET)
  
  # df that holds rec & ssb
  data = data.frame(rec = rec,
                    ssb = ssb)
  
  if(is.null(RANSAC.thresh)){
    nn.dist = abs(median(data$rec) - data$rec)
    RANSAC.thresh = median(nn.dist)
  }
  
  # continous Hockey stick model
  contHS.recruitment = function(pars,S,g){
    pars[[1]]*(S + sqrt(pars[[2]]^2 + (g^2)/4) - sqrt((S - pars[[2]])^2 + (g^2)/4))
  }
  # function to fit a HS through 2 points
  fast.HS.fitting = function(dat){
    if(dim(dat)[1] != 2 & dim(dat)[2] != 2)
      stop("Input must be of the dimension 2x2")
    
    x = dat[,which(names(dat) %in% c("x","ssb","SSB","S"))]
    y = dat[,which(names(dat) %in% c("y","rec","R","Rec"))]
    
    mean.y = mean(y)
    x.min = min(x)
    
    p1 = c(x[x == x.min],y[x == x.min])
    # get slope between zero and lowest point
    slope = as.numeric(p1[2]/p1[1])
    # determine crossing point between line and mean
    x.break = mean.y/slope
    # check if x.min is larger than x.break
    x.break = ifelse(x.break > x.min,x.break,x.min)
    bp = c(x.break,mean.y)
    # translate to model parameters of continous HS fit
    a = 0.5*bp[2]/bp[1]
    b = bp[1]
    return(c(a = a,b = b))
  }

  # ------------------------ #
  # 1. fit Hockeystick model
  # ------------------------ #
  fit.accepted  = list(); n = 1
  for(ii in 1:RANSAC.n){
    # 1. select minimum sample
    
    p2 = data[sample(1:nrow(data),size = 2),]
    
    # fast fitting
    pars1 = fast.HS.fitting(p2)
    
    # get residuals to minimal model
    y.head = contHS.recruitment(pars1[1:2],
                                S = data$ssb,
                                g = HS.g)
    
    a = abs(data$rec - y.head)
    inliers = a<RANSAC.thresh
    # accept
    if(sum(inliers)>= 0.5*nrow(data)){
      p.inlier = data[inliers,]
      
      # fast SR fit implementation approach would be
      out = fitSRCurve(S = p.inlier$ssb,R = p.inlier$rec,
                       shape = "contHockey",g = HS.g) 
      pars1 = out$env$last.par[1:2]
      
      # store accepted fits
      fit.accepted[[n]] = data.frame(a = pars1[1],
                                     b = pars1[2],
                                     no.inlier = sum(inliers))
      n = n+1
    } else {
      next
    }
  }
  fit.accepted.df = data.frame(do.call(rbind,fit.accepted))
  pars.best = fit.accepted.df[fit.accepted.df$no.inlier == max(fit.accepted.df$no.inlier),]
  if(length(dim(pars.best)) == 2){
    pars.best = apply(pars.best,2,median)[1:2]
  }
 pars = pars.best 
 
  # ------------------------------ #
  # 2. get residuals to fit glm to
  # ------------------------------ #
  
  # get residuals
  resids.rec = log(data$rec)-log(contHS.recruitment(pars = exp(pars),S = data$ssb,g = HS.g))
  
  if(pracma::isempty(env.Vars) | all(is.na(env.Vars)) | is.null(env.Vars)){
    predicted = contHS.recruitment(pars = exp(pars),S = data$ssb,g = HS.g)
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
    predicted = contHS.recruitment(pars = exp(pars),S = data$ssb,g = HS.g)
      } else {
    predicted = exp(log(contHS.recruitment(pars = exp(pars),S = data$ssb,g = HS.g)) +
                      as.numeric(predict(residual.model)))
  
    }
  }
  # return stuff in a list
  return.obj = list(predicted = predicted,
                    HS.params = c(pars,g = HS.g),
                    #optimx.output = HS,
                    HockeyStick.model = list(model = contHS.recruitment,
                                              params = c(pars,g = HS.g)),
                    MARS.model = residual.model)
  
  return(return.obj)
}

predict.robust.HS = function(model.obj,newdata){
  library(earth)
  ssb = newdata[,grep("ssb",names(newdata))]
  a = exp(model.obj$HS.params[["a"]])
  b = exp(model.obj$HS.params[["b"]])
  g = model.obj$HS.params[["g"]]
  predict.HS = a*(c(ssb) + sqrt(b^2 + (g^2)/4) - sqrt((c(ssb) - b)^2 + (g^2)/4))
  if(all(is.na(model.obj$MARS.model)) || length(model.obj$MARS.model$coefficients) == 1){
    predict.MARSOut = 0
  } else{
    predict.MARSOut = as.numeric(predict(model.obj$MARS.model,newdata))
  }
  predict.all = exp(log(predict.HS)+predict.MARSOut)
  
  return(predict.all)
}

