modified_stepwise_glm_backward_parallel <- function(df,
                                  y,
                                  x,
                                  c,
                                  maxstep,
                                  eval_ls,
                                  eval_by,
                                  nfold=10, # default 10 fold cv for "cvpred" and "cvDeviance"
                                  family = c("binomial", "gaussian")[1],
                                  free_cores = 2
){
  res_ls <- list()
  tune_ls <- list()
  x_remove <- c() # remove one variable per each step
  x_picked <- x
  if(is.null(maxstep)) maxstep <- length(x)
  for(s in rev(seq(1+length(x)-maxstep,length(x)))){
    
    if(s == length(x)){ 
      rmvd_mat <- combn(x_picked, 0) # full model
    }else{
      rmvd_mat <- combn(x_picked, 1) 
    }
    
    # ---- find optimal variable at this model size (step) ---- 
    numCores <- detectCores() - free_cores  # Leave two cores free for system processes
    registerDoParallel(cores=numCores)
    tune_ls[[length(x)-s+1]] <- foreach(i = c(1:ncol(rmvd_mat)), .packages = c("pROC", "dplyr")) %dopar% {
      # define result list object to return
      tune_score <- NULL
      
      # ---- train ----
      # add "fold" column to original df
      x_sub <- setdiff(x_picked, rmvd_mat[,i])
      f <- "fold"
      fold_idx_df <- data.frame(c=unique(as.character(df[,c])), 
                                fold=cut(seq(1, length(unique(as.character(df[,c])))), 
                                         min(length(unique(as.character(df[,c]))),nfold), 
                                         labels=c(seq(1, nfold))))
      df$c <- df[,c]
      df <- merge(df,fold_idx_df,all.x=TRUE)
      df <- df[,setdiff(colnames(df), c("c")) ]
      df_mdl <- df[complete.cases(df[,c(x_sub,c,f,y)]),c(x_sub,c,f,y)] #dplyr::distinct() # only keep complete rows but must keep duplicates rows for same scale selection 
      # logistic regression
      fml <- paste0(y, "~", paste0(x_sub, collapse = "+"))
      mdl <- glm(fml, family = family, data = df_mdl)
      mdl$c <- df_mdl[,c]
      mdl$x_subset <- x_sub
      # ---- eval ----
      if(eval_by == "NIC") { tune_score <- NIC(mdl, family=family)$nic }
      if(eval_by == "AIC") { tune_score <- NIC(mdl, family=family)$aic }
      if(eval_by == "NICc") { tune_score <- NIC(mdl, family=family)$nicc }
      if(eval_by == "BIC") { tune_score <- BIC(mdl) }
      if(eval_by == "Deviance") { tune_score <- NIC(mdl, family=family)$dev }
      if(eval_by %in% c("cvpred","cvDeviance") ){
        fold_vec <- df_mdl[,f]
      }
      if(eval_by %in% c("loopred","looDeviance") ){
        fold_vec <- df_mdl[,c]
      }
      if(eval_by %in% c("cvpred","cvDeviance","loopred","looDeviance") ){
        # do cv or loo
        y_prob <- c()
        y_true <- c()
        for (f_sub in unique(fold_vec)){
          mdl_sub <- glm(fml, family = family, data = df_mdl[!fold_vec==f_sub,])
          y_prob <- c(y_prob, predict(mdl_sub, type = "response",
                                      newdata = df_mdl[fold_vec==f_sub,]))
          y_true <- c(y_true, df_mdl[fold_vec==f_sub,y])
        }
        if(eval_by %in% c("cvpred","loopred") ){
          if(family=="binomial"){
            tune_score <- as.numeric(pROC::auc(pROC::roc(response = y_true, predictor = y_prob))) # AUC
          }
          if(family=="gaussian"){
            tune_score <- sqrt(sum( (y_true-y_prob)^2 )) # MSE
          }
        }
        if(eval_by %in% c("cvDeviance","looDeviance") ){
          if(family=="binomial"){
            tune_score <- -2*sum(y_true * log(y_prob) + (1 - y_true) * log(1 - y_prob), na.rm = TRUE)
          }
          if(family=="gaussian"){
            s <- sigma(mdl)
            tune_score <- length(y_true)*log(2*pi*s^2) + 1/s^2*sum( (y_true-y_prob)^2 )  
          }
        }
      }
      return(list("tune_score" = tune_score,
                  "x_removed" = rmvd_mat[,i])) # x_pick means x picked to remove
    }
    
    
    # ---- find optimal variable at this model size (step) ----
    tune_vec <- unlist(lapply(tune_ls[[length(x)-s+1]], function(e){return(e$tune_score)}))
    x_removed_vec <- unlist(lapply(tune_ls[[length(x)-s+1]], function(e){return(e$x_removed)}))
    if(is.null(x_removed_vec)) x_removed_vec <- "All"
    if(eval_by %in% c("cvpred","loopred") ) tune_vec <- -tune_vec
    x_removed <- x_removed_vec[which(tune_vec==min(tune_vec))][1]
    at_score <- unlist(lapply(tune_ls[[length(x)-s+1]], function(e){return(e$tune_score)}))[which(tune_vec==min(tune_vec))]
    # ---- calculate other matrices with current setting ----
    x_picked <- setdiff(x_picked, x_removed)
    x_remove <- c(x_remove, x_removed)
    res_ls[[length(x)-s+1]] <- all_subset_glm_parallel(df, y, x_picked, c, 
                                  size=length(x_picked),
                                  eval_ls=eval_ls,
                                  family = family,
                                  free_cores = free_cores)
    res_ls[[length(x)-s+1]][['x_removed']] <- x_removed
    res_ls[[length(x)-s+1]][['at_score']] <- at_score
  }
  
  return(list("res_ls" = res_ls, 
              "tune_ls" = tune_ls))
}
