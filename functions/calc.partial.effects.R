# ============================= #
# function to calculate partial 
# effects of the fitted models
# ============================= #


calc.partial.effects = function(model.obj,y,x,predict.func) {
  stopifnot(require(pdp))
  stopifnot(require(reshape2))
  stopifnot(require(ggplot2))
  
  
  xy = data.frame(y = y,x)
  # ------------------------ #
  # evaluate partial effects
  # ------------------------ #
  
  vars = model.obj$vars
  pdp_out.all = list()
  median_pdp.all = list()
  n = 1
  for(v in vars){
    if(grepl("ssb",v)){
      pred.grid = data.frame(seq(0,max(xy[,grepl("ssb",names(xy))]),length.out = 100))
      names(pred.grid) = v
      pdp_out = partial(object = model.obj,pred.var = v,train = xy,pred.grid = pred.grid,
                        pred.fun = function(object,newdata) predict.func(object,newdata),plot = F)
    } else {
      # add 20% extrapolation bounds
      pred.grid = data.frame(seq(ifelse(min(xy[,v])<0,1.2*min(xy[,v]),0.8*min(xy[,v])),
                                 ifelse(max(xy[,v])<0,0.8*max(xy[,v]),1.2*max(xy[,v])),length.out = 100))
      names(pred.grid) = v
      pdp_out = partial(object = model.obj,pred.var = v,train = xy,pred.grid = pred.grid,
                        pred.fun = function(object,newdata) predict.func(object,newdata),plot = F)
    }
    by.id = pdp_out[,v] #rep(1:max(pdp_out$yhat.id),nrow(pdp_out)/max(pdp_out$yhat.id))
    median_pdp = aggregate(pdp_out,by = list(by.id),FUN = median)
    median_pdp = median_pdp[,1:3]
    names(median_pdp)[which(names(median_pdp) == "Group.1")] = "yhat.id"
    
    
    pdp_out.all[[n]] = melt(pdp_out,id.vars = c("yhat.id","yhat"))
    median_pdp.all[[n]] = melt(median_pdp,id.vars = c("yhat.id","yhat"))
    n = n+1
  }
  pdp_out = do.call(rbind,pdp_out.all)
  median_pdp = do.call(rbind,median_pdp.all)
  vars_pdp = melt(xy[,c("y",vars),drop = F],id.vars = "y")
  
  pdp.plot = ggplot(pdp_out,aes(x = value,y = yhat)) + 
             geom_line(aes(group = yhat.id)) +
             geom_line(data = median_pdp,aes(x = value,y = yhat),col = "red",linewidth = 2) +
             geom_point(data = vars_pdp,aes(x = value,y = y),pch = 21,fill = "burlywood1",colour = "gray50") +          
             facet_wrap(~variable,scales = "free") + 
             theme_bw()
  
  # return
  return.list = list(pdp = pdp_out,
       median_pdp = median_pdp,
       plot = pdp.plot)
  return(return.list)
}
