#' All subset model selection for logistic regression in parallel
#' 
#' @param df data.frame to run logistic regression on
#' @param x set of individual predictors to create all possible combination of subsets
#' @param y response variable name
#' @param c cluster variable name, to calculate nic, loocv_devienve, loocv_auc, 10fold_cvpred
#' @param size of the model, list of integer or continuous sequence i.e seq(4,6)
#' if not specified, all possible model size will be used which will take longer time
#' @param eval_ls evaluation matrices to use default c("Deviance", "AIC", "NIC", "cvpred", "cvDeviance","loopred", "looDeviance"), can choose only a subset
#' @returns 


library(foreach)
library(doParallel)
all_subset_glm_parallel <- function(df,
                           y,
                           x,
                           c,
                           size=NULL,
                           eval_ls=c("Deviance", "AIC", "BIC", "NIC", "NICc",
                                   "cvpred", "cvDeviance",
                                   "loopred", "looDeviance")[1:6],
                           nfold=10, # default 10 fold cv for "cvpred" and "cvDeviance"
                           family = c("binomial", "gaussian")[1],
                           free_cores = 2){
  
  if(is.null(size)) size <- seq(1,length(x))
  if(length(size)>0) size <- size[size<=length(x)]
  if(length(size)==0) return("invalid size")
  res_ls <- list()
  for(s in size){
    # Generate all combinations of size s in a matrix
    comb_mat <- combn(x, s)
    # parallel train logistic regression on each subset
    numCores <- detectCores() - free_cores  # Leave two cores free for system processes
    registerDoParallel(cores=numCores)
    res_ls[[paste0("s",s)]] <- foreach(i = c(1:ncol(comb_mat)), .packages = c("pROC", "dplyr")) %dopar% {
      # define result list object to return
      res <- list(Deviance = NULL,
                  AIC=NULL,
                  BIC=NULL,
                  NIC=NULL,
                  NICc=NULL,
                  cvpred=NULL,
                  cvDeviance=NULL,
                  loopred = NULL,
                  looDeviance = NULL)
      
      # ---- train ----
      # add "fold" column to original df
      x_sub <- comb_mat[,i]
      f <- "fold"
      fold_idx_df <- data.frame(c=unique(as.character(df[,c])), 
                                fold=cut(seq(1, length(unique(as.character(df[,c])))), 
                                         min(length(unique(as.character(df[,c]))),nfold), 
                                         labels=c(seq(1, nfold))))
      df$c <- df[,c]
      df <- merge(df,fold_idx_df,all.x=TRUE)
      df <- df[,setdiff(colnames(df), c("c")) ]
      df_mdl <- df[complete.cases(df[,c(x_sub,c,f,y)]),c(x_sub,c,f,y)] # only keep complete rows (keep duplicates)
      # logistic regression
      fml <- paste0(y, "~", paste0(x_sub, collapse = "+"))
      mdl <- glm(fml, family = family, data = df_mdl)
      mdl$c <- df_mdl[,c]
      mdl$x_subset <- x_sub
      
      
      
      # ---- eval ----
      if( length(intersect(c("NIC","NICc","AIC","BIC","Deviance"),eval_ls))>0){
        ics <- NIC(mdl, family=family)
        res$AIC <- ics$aic
        res$NIC <- ics$nic
        res$NICc <- ics$nicc
        res$Deviance <- ics$dev
        res$BIC <- BIC(mdl)
      }
      eval_extra <- setdiff(eval_ls, c("NIC","NICc","AIC","BIC","Deviance"))
      if(length(eval_extra)>0){
        for(e in eval_extra){
          if(e %in% c("cvpred","cvDeviance") ){
            fold_vec <- df_mdl[,f]
          }
          if(e %in% c("loopred","looDeviance") ){
            fold_vec <- df_mdl[,c]
          }
          # do cv or loo
          y_prob <- c()
          y_true <- c()
          for (f_sub in unique(fold_vec)){
            mdl_sub <- glm(fml, family = family, data = df_mdl[!fold_vec==f_sub,])
            y_prob <- c(y_prob, predict(mdl_sub, type = "response",
                                        newdata = df_mdl[fold_vec==f_sub,]))
            y_true <- c(y_true, df_mdl[fold_vec==f_sub,y])
          }
          if(family=="binomial"){
            cvpred <- as.numeric(pROC::auc(pROC::roc(response = y_true, predictor = y_prob))) # AUC
            cvDeviance <- -2*sum(y_true * log(y_prob) + (1 - y_true) * log(1 - y_prob), na.rm = TRUE)
          }
          if(family=="gaussian"){
            cvpred <- sqrt(sum( (y_true-y_prob)^2 )) # MSE
            s <- sigma(mdl)
            cvDeviance <- length(y_true)*log(2*pi*s^2) + 1/s^2*sum( (y_true-y_prob)^2 )  
          }
          
          if(e %in% c("cvpred","cvDeviance") ){
            res$cvpred <- cvpred
            res$cvDeviance <- cvDeviance
          }
          if(e %in% c("loopred","looDeviance") ){
            res$loopred <- cvpred
            res$looDeviance <- cvDeviance
          }
        }
      }
      
      res$mdl <- mdl
      return(res)
    }
  }
  return(res_ls)
}
