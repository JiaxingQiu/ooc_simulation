#' All subset model selection for logistic regression in parallel
#' 
#' @param df data.frame to run logistic regression on
#' @param x set of individual predictors to create all possible combination of subsets
#' @param y response variable name
#' @param c cluster variable name, to calculate nic, loocv_devienve, loocv_auc, 10fold_cvpred
#' @param size of the model, list of integer or continuous sequence i.e seq(4,6)
#' @param forward default T using forward step wise, set to F to use backward
#' if not specified, all possible model size will be used which will take longer time
#' @param eval_ls evaluation matrices to use default c("Deviance", "AIC", "NIC", "cvpred", "cvDeviance","loopred", "looDeviance"), can choose only a subset
#' @returns 


library(foreach)
library(doParallel)
modified_stepwise_glm_parallel <- function(df,
                           y,
                           x,
                           c,
                           forward = T, 
                           maxstep = NULL,
                           eval_ls = c("Deviance", "AIC", "BIC", "NIC", "NICc", 
                                     "cvpred", "cvDeviance",
                                     "loopred", "looDeviance")[1:6], # loopred = looauc if binomial, loopred = loomse if gaussian 
                           eval_by = c("Deviance", "AIC", "BIC", "NIC", "NICc", 
                                       "cvpred", "cvDeviance",
                                       "loopred", "looDeviance")[4],
                           nfold=10, # default 10 fold cv for "cvpred" and "cvDeviance"
                           family = c("binomial", "gaussian")[1],
                           free_cores = 2
){
  
  res_ls <- list()
  tune_ls <- list()
  if(is.null(maxstep)) maxstep <- length(x)
  # forward step wise
  if(forward){
    x_left <- x
    x_picked <- c()
    for(st in seq(1, maxstep)){
      comb_mat <- combn(x_left, 1)
      # ---- find optimal variable at this model size (step) ---- 
      numCores <- detectCores() - free_cores  # Leave two cores free for system processes
      registerDoParallel(cores=numCores)
      tune_ls[[st]] <- foreach(i = c(1:ncol(comb_mat)), .packages = c("pROC", "dplyr")) %dopar% {
        # define result list object to return
        tune_score <- NULL
        
        # ---- train ----
        # add "fold" column to original df
        x_sub <- c(x_picked,comb_mat[,i])
        stopifnot(length(x_sub)==st)
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
              st <- sigma(mdl)
              tune_score <- length(y_true)*log(2*pi*st^2) + 1/st^2*sum( (y_true-y_prob)^2 )  
            }
          }
        }
        return(list("tune_score" = tune_score,
                    "x_pick" = comb_mat[1,i]))
      }
      
      
      # ---- find optimal variable at this model size (step) ----
      tune_vec <- unlist(lapply(tune_ls[[st]], function(e){return(e$tune_score)}))
      x_pick_vec <- unlist(lapply(tune_ls[[st]], function(e){return(e$x_pick)}))
      if(eval_by %in% c("cvpred","loopred") ) tune_vec <- -tune_vec
      x_pick <- x_pick_vec[which(tune_vec==min(tune_vec))][1]
      at_score <- unlist(lapply(tune_ls[[st]], function(e){return(e$tune_score)}))[which(tune_vec==min(tune_vec))][1]
      # ---- calculate other matrices with current setting ----
      x_picked <- c(x_picked, x_pick)
      x_left <- setdiff(x_left, x_pick)
      res_ls[[st]] <- all_subset_glm_parallel(df, y, x_picked, c, 
                                    size=length(x_picked),
                                    eval_ls=eval_ls,
                                    family = family,
                                    free_cores = free_cores)
      res_ls[[st]][['x_pick']] <- x_pick
      res_ls[[st]][['at_score']] <- at_score
    }
  }
  
  # backward step wise
  if(!forward){
    back_obj <- modified_stepwise_glm_backward_parallel(df, y, x, c, maxstep, eval_ls, eval_by, nfold, family)
    res_ls <- back_obj$res_ls
    tune_ls <- back_obj$tune_ls
  }
  

  final_obj <- list("info"=list("forward"=forward,
                                "eval_by"=eval_by),
                    "res_ls" = res_ls,
                    "tune_ls" = tune_ls )
  return(final_obj)
  
}
