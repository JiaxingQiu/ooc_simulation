format_forward <- function(fwd_obj){
  #forward
  model_size <- c()
  x_picked <- c()
  dev <- c()
  aic <- c()
  bic <- c()
  nic <- c()
  nicc <- c()
  cvpred <- c()
  cvdev <- c()
  loopred <- c()
  loodev <- c()
  for(s in seq(1,length(fwd_obj$res_ls))){
    model_size <- c(model_size, s)
    x_picked <- c(x_picked, fwd_obj$res_ls[[s]][['x_pick']])
    dev <- c(dev, fwd_obj$res_ls[[s]][[paste0("s",s)]][[1]]$Deviance)
    aic <- c(aic, fwd_obj$res_ls[[s]][[paste0("s",s)]][[1]]$AIC)
    bic <- c(bic, fwd_obj$res_ls[[s]][[paste0("s",s)]][[1]]$BIC)
    nic <- c(nic, fwd_obj$res_ls[[s]][[paste0("s",s)]][[1]]$NIC)
    nicc <- c(nicc, fwd_obj$res_ls[[s]][[paste0("s",s)]][[1]]$NICc)
    cvpred <- c(cvpred, fwd_obj$res_ls[[s]][[paste0("s",s)]][[1]]$cvpred)
    cvdev <- c(cvdev, fwd_obj$res_ls[[s]][[paste0("s",s)]][[1]]$cvDeviance)
    loopred <- c(loopred, fwd_obj$res_ls[[s]][[paste0("s",s)]][[1]]$loopred)
    loodev <- c(loodev, fwd_obj$res_ls[[s]][[paste0("s",s)]][[1]]$looDeviance)
  }
  df <- data.frame(model_size, x_picked)
  if(length(dev)>0) df$dev <- dev
  if(length(aic)>0) df$aic <- aic
  if(length(bic)>0) df$bic <- bic
  if(length(nic)>0) df$nic <- nic
  if(length(nicc)>0) df$nicc <- nicc
  if(length(cvpred)>0) df$cvpred <- cvpred
  if(length(cvdev)>0) df$cvdev <- cvdev
  if(length(loopred)>0) df$loopred <- loopred
  if(length(loodev)>0) df$loodev <- loodev
  
  return(df)
  
}


format_backward <- function(bwd_obj){
  model_size <- c()
  x_removed <- c()
  dev <- c()
  aic <- c()
  bic <- c()
  nic <- c()
  cvpred <- c()
  cvdev <- c()
  loopred <- c()
  loodev <- c()
  for(s in seq(1,length(bwd_obj$res_ls))){
    model_size <- c(model_size, length(bwd_obj$res_ls)-s+1)
    x_removed <- c(x_removed, bwd_obj$res_ls[[s]][['x_removed']])
    dev <- c(dev, bwd_obj$res_ls[[s]][[paste0("s",length(bwd_obj$res_ls)-s+1)]][[1]]$Deviance)
    aic <- c(aic, bwd_obj$res_ls[[s]][[paste0("s",length(bwd_obj$res_ls)-s+1)]][[1]]$AIC)
    bic <- c(bic, bwd_obj$res_ls[[s]][[paste0("s",length(bwd_obj$res_ls)-s+1)]][[1]]$BIC)
    nic <- c(nic, bwd_obj$res_ls[[s]][[paste0("s",length(bwd_obj$res_ls)-s+1)]][[1]]$NIC)
    cvpred <- c(cvpred, bwd_obj$res_ls[[s]][[paste0("s",length(bwd_obj$res_ls)-s+1)]][[1]]$cvpred)
    cvdev <- c(cvdev, bwd_obj$res_ls[[s]][[paste0("s",length(bwd_obj$res_ls)-s+1)]][[1]]$cvDeviance)
    loopred <- c(loopred, bwd_obj$res_ls[[s]][[paste0("s",length(bwd_obj$res_ls)-s+1)]][[1]]$loopred)
    loodev <- c(loodev, bwd_obj$res_ls[[s]][[paste0("s",length(bwd_obj$res_ls)-s+1)]][[1]]$looDeviance)
  }
  df <- data.frame(model_size, x_removed)
  if(length(dev)>0) df$dev <- dev
  if(length(aic)>0) df$aic <- aic
  if(length(bic)>0) df$bic <- bic
  if(length(nic)>0) df$nic <- nic
  if(length(cvpred)>0) df$cvpred <- cvpred
  if(length(cvdev)>0) df$cvdev <- cvdev
  if(length(loopred)>0) df$loopred <- loopred
  if(length(loodev)>0) df$loodev <- loodev
  
  return(df)
  
}
