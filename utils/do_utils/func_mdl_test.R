mdl_test <- function(
  y_true,
  y_prob,
  threshold=0.5
){
  
  stopifnot(!is.na(threshold))
  # create score dataframe
  res_df <- data.frame(message="external validation initiated")
  
  # ---- test data is engineered the same way to have valid scores ----
  y_true <- as.numeric(y_true)
  y_prob <- as.numeric(y_prob)
  stopifnot(length(y_prob)==length(y_true))
  
  df_tmp <- data.frame(y_true, y_prob)
  df_tmp <- df_tmp[complete.cases(df_tmp),]# only keep complete finite values
  df_tmp <- df_tmp[which(is.finite(df_tmp$y_true) & is.finite(df_tmp$y_prob)),]
  y_true <- df_tmp$y_true 
  y_prob <- df_tmp$y_prob
  
  y_prob <- round(as.numeric(as.character(y_prob)),6)
  # if(!all(y_prob>=0 & y_prob<=1)){
  #   y_prob <- (y_prob-min(y_prob))/(max(y_prob)-min(y_prob))
  # }
  # stopifnot(all(y_prob>=0 & y_prob<=1)) # not a must though 
  y_pred <- as.numeric(ifelse(y_prob<=threshold,0,1))
  
  logloss <- NA
  tryCatch({
    logloss = MLmetrics::LogLoss(y_prob, y_true)
  },error=function(e){
    print(e)
  })
  AUROC <- NA
  tryCatch({
    AUROC = ifelse(is.na(MLmetrics::AUC(y_prob, y_true)), round(pROC::auc(pROC::roc(y_true, round(y_prob,6) )),6), MLmetrics::AUC(y_prob, y_true))
  },error=function(e){
    print(e)
  })
  AUPRC <- NA
  tryCatch({
    AUPRC = MLmetrics::PRAUC(y_pred, y_true)
  },error=function(e){
    print(e)
  })
  accuracy <- NA
  tryCatch({
    accuracy = MLmetrics::Accuracy(y_pred,y_true)
  },error=function(e){
    print(e)
  })
  f1score <- NA
  tryCatch({
    f1score = F1_Score(y_true = y_true, y_pred = y_pred)
  },error=function(e){
    print(e)
  })
  res_df <- data.frame(
    data = "test data",
    logloss = logloss,
    AUROC = AUROC,
    AUPRC = AUPRC,
    accuracy = accuracy,
    f1score = f1score)
 
  return(list("res_df"=res_df))
}
