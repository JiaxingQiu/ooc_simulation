library(glmmLasso)

lasso_x_select_cluster <- function(
    data, 
    y_col,
    cluster_col,# cluster colname
    x_cols_nonlin_rcs5=c(),
    x_cols_nonlin_rcs4=c(),
    x_cols_nonlin_rcs3=c(),
    x_cols_linear=c(), 
    x_cols_fct=c(),
    x_cols_tag=c(),
    family = c("binomial")[1],
    dict_data=NULL, # dictionary table is optional
    lambda_seq = NULL, # seq(500,0,by=-5) external specified lasso lambda value sequence of lambdas to tune by if given
    tune_by = c("AIC", "BIC", "auc", "logloss", "misclass")[1],
    foldid_col=NULL
){
  
  # ---- Description ----
  # dictionary oriented (do) predictor variables selection
  
  # ---- Arguments ----
  # data: a dataframe object with essential dictionary information as attributes on each column
  # dict: dictionary for corresponding data frame input
  # x_cols: pre-selected predictor column names (x)
  
  # ---- Value ----
  # data_final: a dataframe object with updated redundency attributes on each column
  # dict_final: updated dictionary for correspondingdata_final
  
  print("--- response distribution ---")
  print(table(data[,y_col]))
  
  # ---- pre-processing ----
  # linear > rcs3 > rcs4 > rcs5
  # variable name population space
  if(!is.null(dict_data)) {
    input_cols_choices <- intersect(colnames(data), dict_data$varname[which(dict_data$type=="num")])
  } else {
    input_cols_choices <- colnames(data)
  }
  # organize model formula string
  x_cols_linear <- intersect(x_cols_linear, input_cols_choices)
  x_cols_nonlin_rcs3 <- setdiff(x_cols_nonlin_rcs3, x_cols_linear) # remove linear from rcs3
  x_cols_nonlin_rcs3 <- intersect(x_cols_nonlin_rcs3, input_cols_choices)# make sure rcs 3 all come from numeric columns 
  x_cols_nonlin_rcs4 <- setdiff(setdiff(x_cols_nonlin_rcs4, x_cols_linear),x_cols_nonlin_rcs3)
  x_cols_nonlin_rcs4 <- intersect(x_cols_nonlin_rcs4, input_cols_choices)
  x_cols_nonlin_rcs5 <- setdiff(setdiff(setdiff(x_cols_nonlin_rcs5, x_cols_linear),x_cols_nonlin_rcs3),x_cols_nonlin_rcs4)
  x_cols_nonlin_rcs5 <- intersect(x_cols_nonlin_rcs5, input_cols_choices) # now we successfully split the rcs knots groups
  x_cols <- unique(c(x_cols_linear, x_cols_nonlin_rcs3, x_cols_nonlin_rcs4, x_cols_nonlin_rcs5,x_cols_tag))
  
  data_org <- data
  #data[,x_cols] <- scale(data_org[,x_cols], center = standardize, scale = standardize)
  data[,x_cols_tag] <- data_org[,x_cols_tag] # overwrite dummy columns, they should not be scaled
  
  
  # add transformations to predictor variables and save grouping info in a dataframe object
  x_col_df_all <- data.frame()
  # --- rcs5 ---
  for (col in x_cols_nonlin_rcs5){
    success <- FALSE
    tryCatch({
      # additional columns
      data[,col] <- as.numeric(rcs(data[,col],5)[,1])
      data[,paste0(col,"'")] <- as.numeric(rcs(data[,col],5)[,2])
      data[,paste0(col,"''")] <- as.numeric(rcs(data[,col],5)[,3])
      data[,paste0(col,"'''")] <- as.numeric(rcs(data[,col],5)[,4])
      
      # append this varname and group info as a data frame to overall column info
      x_col_df <- data.frame(
        x_colname = c(col, paste0(col,"'"), paste0(col,"''"),paste0(col,"'''")),
        x_group = paste0(col,"_rcs5")
      )
      x_col_df_all <- bind_rows(x_col_df_all, x_col_df)
      success <- TRUE
    },error=function(e){
      print(paste0("failed to assign rcs 4 for ", col, ", using linear instead"))
    })
    if(!success){
      x_cols_linear <- c(x_cols_linear, col)
      x_cols_nonlin_rcs5 <- setdiff(x_cols_nonlin_rcs5, col)
    }
  }
  
  # --- rcs4 ---
  for (col in x_cols_nonlin_rcs4){
    success <- FALSE
    tryCatch({
      # additional columns
      data[,col] <- as.numeric(rcs(data[,col],4)[,1])
      data[,paste0(col,"'")] <- as.numeric(rcs(data[,col],4)[,2])
      data[,paste0(col,"''")] <- as.numeric(rcs(data[,col],4)[,3])
      
      # append this varname and group info as a data frame to overall column info
      x_col_df <- data.frame(
        x_colname = c(col, paste0(col,"'"), paste0(col,"''")),
        x_group = paste0(col,"_rcs4")
      )
      x_col_df_all <- bind_rows(x_col_df_all, x_col_df)
      success <- TRUE
    },error=function(e){
      print(paste0("failed to assign rcs 4 for ", col, ", using linear instead"))
    })
    if(!success){
      x_cols_linear <- c(x_cols_linear, col)
      x_cols_nonlin_rcs4 <- setdiff(x_cols_nonlin_rcs4, col)
    }
  }
  
  # --- rcs3 ---
  for (col in x_cols_nonlin_rcs3){
    success <- FALSE
    tryCatch({
      # additional columns
      data[,col] <- as.numeric(rcs(data[,col],3)[,1])
      data[,paste0(col,"'")] <- as.numeric(rcs(data[,col],3)[,2])
      
      # append this varname and group info as a data frame to overall column info
      x_col_df <- data.frame(
        x_colname = c(col, paste0(col,"'")),
        x_group = paste0(col,"_rcs3")
      )
      x_col_df_all <- bind_rows(x_col_df_all, x_col_df)
      success <- TRUE
    },error=function(e){
      print(paste0("failed to assign rcs 3 for ", col, ", using linear instead"))
    })
    if(!success){
      x_cols_linear <- c(x_cols_linear, col)
      x_cols_nonlin_rcs3 <- setdiff(x_cols_nonlin_rcs3, col)
    }
  }
  
  # --- linear ---
  for (col in x_cols_linear){
    success <- FALSE
    tryCatch({
      # append this varname and group info as a data frame to overall column info
      x_col_df <- data.frame(
        x_colname = c(col),
        x_group = paste0(col,"_linear")
      )
      x_col_df_all <- bind_rows(x_col_df_all, x_col_df)
      success <- TRUE
    },error=function(e){
      print(paste0("failed to assign linear for ", col, ", variable removed"))
    })
    if(!success) x_cols_linear <- setdiff(x_cols_linear, col)
  }
  
  # --- fct ---
  for (col in x_cols_fct){
    # dummy the variable
    levels <- unique(as.character(data[,col]))
    for(l in levels){
      l_fix <- gsub("[^[:alnum:]]","_",l)
      data[,paste0(col,"___",l_fix)] <- ifelse(data[,col]==l,1,0)
      # append this varname and group info as a data frame to overall column info
      x_col_df <- data.frame(
        x_colname = c(paste0(col,"___",l_fix)),
        x_group = paste0(col,"_level")
      )
      x_col_df_all <- bind_rows(x_col_df_all, x_col_df)
    }
  }
  # --- tag ---
  for (col in x_cols_tag){
    # append this varname and group info as a data frame to overall column info
    x_col_df <- data.frame(
      x_colname = c(col),
      x_group = paste0(col,"_tag01")
    )
    x_col_df_all <- bind_rows(x_col_df_all, x_col_df)
  }
  # store modeling data in the same format
  x <- data.matrix(data[complete.cases(data[,c(x_col_df_all$x_colname,y_col)]),c(x_col_df_all$x_colname)])
  y <- data.matrix(data[complete.cases(data[,c(x_col_df_all$x_colname,y_col)]),y_col])
  colnames(y) <- y_col
  y <- ifelse(y==1, 1, -1)
  
  
  
  # ---- engineer input arguments ----
  # tune_by
  if(!tune_by%in%c("AIC", "BIC", "auc", "logloss", "misclass")) tune_by <- "BIC"
  # lambda_seq
  if(is.null(lambda_seq)) lambda_seq <- exp(seq(log(460),log(1e-4),length.out = 50)) # c(seq(10000,2000,by=-1000),seq(1000,200,by=-100),seq(100,20,by=-5),seq(10,2,by=-0.5),seq(1,0.01,by=-0.1))
  
  
  
  # ---- lasso regression tuning ----
  coeff_lambda_trace <- data.frame()
  score_lambda_trace <- data.frame()
  ## loop over lambda grid
  for(j in 1:length(lambda_seq)){
    tryCatch({
      # tune by AIC / BIC (always run)
      all_df <- data
      all_df$cluster <- as.factor(rownames(all_df))
      all_df <- all_df[complete.cases(all_df[,c(x_col_df_all$x_colname,y_col,"cluster")]),c(x_col_df_all$x_colname,y_col,"cluster")]
      # train overall model on data 
      lasso_all <- NULL
      lasso_all <- glmmLasso::glmmLasso(fix = formula(paste0(y_col," ~ ",paste0(paste0("`",x_col_df_all$x_colname,"`"), collapse = " + "))),
                                        rnd = list(cluster=~1),
                                        data = all_df,
                                        lambda=lambda_seq[j],
                                        family = binomial(link="logit"),
                                        control=list(print.iter=TRUE))
                                        
      Delta.start <- as.matrix(t(as.numeric( lasso_all$Deltamatrix[lasso_all$conv.step,] )))
      Q.start <- as.numeric( lasso_all$Q_long[[lasso_all$conv.step+1]] )
      # get coefficient
      coeff_lambda_df <- as.data.frame(t(as.data.frame(lasso_all$coefficients)))
      rownames(coeff_lambda_df) <- NULL
      # skip current lambda if all coefficients become zero
      stopifnot(!rowSums( coeff_lambda_df[,setdiff(colnames(coeff_lambda_df), "(Intercept)")], na.rm = TRUE)==0)
      coeff_lambda_df$lambda <- lambda_seq[j]
      coeff_lambda_trace <- bind_rows(coeff_lambda_trace, coeff_lambda_df)
      score_lambda_df <- data.frame(lambda = lambda_seq[j], 
                                    AIC = lasso_all$aic, 
                                    BIC = lasso_all$bic)
      # tune by AUC / logloss 
      if (tune_by %in% c("auc", "logloss", "misclass")){
        # ------ manually 10 fold cross validation scores ---------
        valid_yhat_df_all <- data.frame()
        score_final_cv <- data.frame()
        for(i in unique(data[,foldid_col])){
          try({
            # prepare suitable data format to model
            validset <- data[which(data[,foldid_col]==i), ]
            trainset <- data[which(!data[,foldid_col]==i), ]
            
            train_df <- as.data.frame(trainset)
            train_df$cluster <- as.factor(rownames(train_df))
            valid_df <- as.data.frame(validset)
            valid_df$cluster <- as.factor(rownames(valid_df))
            train_df <- train_df[complete.cases(train_df[,c(x_col_df_all$x_colname,y_col,"cluster")]),]
            valid_df <- valid_df[complete.cases(valid_df[,c(x_col_df_all$x_colname,y_col,"cluster")]),]
            
            
            # train current model on train_df
            lasso_fold <- NULL
            lasso_fold <- glmmLasso::glmmLasso(fix = formula(paste0(y_col," ~ ",paste0(paste0("`",x_col_df_all$x_colname,"`"), collapse = " + "))),
                                               rnd = list(cluster=~1),
                                               data = train_df,
                                               lambda=lambda_seq[j],
                                               family = binomial(link="logit"),
                                               control=list(print.iter=TRUE,
                                                            start=Delta.start,
                                                            q_start=Q.start)
                                               )
            # predict on left out validation df 
            y_prob_valid <- predict(lasso_fold, valid_df) # probability scale
            valid_yhat_df <- data.frame( 
              rowid = rownames( validset[complete.cases(validset[,c(x_col_df_all$x_colname,y_col)]),x_col_df_all$x_colname] ),
              y_true = as.numeric(  ifelse(valid_df[,y_col]>0,1,0) ),
              y_prob =  as.numeric( y_prob_valid) ,
              fold = i
            )
            valid_yhat_df_all <- bind_rows(valid_yhat_df_all, valid_yhat_df)
            
          },TRUE)
        }
        
        # cross-validated labeling (yhat)
        cv_final_scores <- mdl_test(y_true = valid_yhat_df_all$y_true,
                                    y_prob = valid_yhat_df_all$y_prob,
                                    threshold = mean(data[,y_col],na.rm=TRUE))
        score_final_cv <- cv_final_scores$res_df
        score_final_cv$data <- "cv_yhat_df"
        score_final_cv$AUROC <- ifelse(score_final_cv$AUROC>0.5, score_final_cv$AUROC, 1-score_final_cv$AUROC)
        score_lambda_df <- data.frame(lambda = lambda_seq[j], 
                                      auc = score_final_cv$AUROC,
                                      logloss = score_final_cv$logloss,
                                      misclass = 1-score_final_cv$accuracy,
                                      AIC = lasso_all$aic, 
                                      BIC = lasso_all$bic)
      }
      # get score traces
      score_lambda_trace <- bind_rows(score_lambda_trace, score_lambda_df)
      print(paste0("keep lamdba -- ", lambda_seq[j]))
    },error=function(e){
      print(paste0("skip lamdba -- ", lambda_seq[j]))
      print(e)
    })
  }
  
  
  
  # ---- final optimal model with optimized lambda ----
  lambda_opt <- score_lambda_trace$lambda[which.min(score_lambda_trace[,tune_by])][1]
  if(tune_by=="auc") lambda_opt <- score_lambda_trace$lambda[which.max(score_lambda_trace[,tune_by])][1]
  all_df <- data
  all_df$cluster <- as.factor(rownames(all_df))
  all_df <- all_df[complete.cases(all_df[,c(x_col_df_all$x_colname,y_col,"cluster")]),]
  lasso_optimal <- glmmLasso::glmmLasso(fix = formula(paste0(y_col," ~ ",paste(paste0("`",x_col_df_all$x_colname,"`"), collapse = " + "))),
                                        rnd = list(cluster=~1),
                                        data = all_df,
                                        lambda=lambda_opt,
                                        family = binomial(link="logit"))
  
  
  lasso_optimal$x <- x
  lasso_optimal$y <- ifelse(y<0,0,1)
  lasso_optimal$group_info <- x_col_df_all
  lasso_optimal$c_info <- cluster_col
  lasso_optimal$lambda.opt <- lambda_opt
  
  # ---- final 10 fold evaluation ----
  valid_yhat_df_all <- data.frame()
  score_final_cv <- data.frame()
  for(i in unique(data[,foldid_col])){
    try({
      # prepare suitable data format to model
      validset <- data[which(data[,foldid_col]==i), ]
      trainset <- data[which(!data[,foldid_col]==i), ]
      train_df <- as.data.frame(trainset)
      train_df$cluster <- as.factor(rownames(train_df))
      valid_df <- as.data.frame(validset)
      valid_df$cluster <- as.factor(rownames(valid_df))
      train_df <- train_df[complete.cases(train_df[,c(x_col_df_all$x_colname,y_col,"cluster")]),]
      valid_df <- valid_df[complete.cases(valid_df[,c(x_col_df_all$x_colname,y_col,"cluster")]),]
      # train current model on train_df
      lasso_fold <- NULL
      lasso_fold <- glmmLasso::glmmLasso(fix = formula(paste0(y_col," ~ ",paste(paste0("`",x_col_df_all$x_colname,"`"), collapse = " + "))),
                                         rnd = list(cluster=~1),
                                         data = train_df,
                                         lambda=lambda_opt,
                                         family = binomial(link="logit"))
     
      y_prob_valid <- predict(lasso_fold, valid_df) # probability scale
      valid_yhat_df <- data.frame( 
        rowid = rownames( validset[complete.cases(validset[,c(x_col_df_all$x_colname,y_col)]),x_col_df_all$x_colname] ),
        #cluster = validset[complete.cases(validset[,c(x_col_df_all$x_colname,y_col)]),cluster_col],
        y_true = as.numeric(  ifelse(valid_df[,y_col]>0,1,0) ),
        y_prob =  as.numeric( y_prob_valid) ,
        fold = i
      )
      valid_yhat_df_all <- bind_rows(valid_yhat_df_all, valid_yhat_df)
    },TRUE)
  }
  # cross-validated labeling (yhat)
  cv_final_scores <- mdl_test(y_true = valid_yhat_df_all$y_true,
                              y_prob = valid_yhat_df_all$y_prob,
                              threshold = mean(data[,y_col],na.rm=TRUE))
  score_final_cv <- cv_final_scores$res_df
  score_final_cv$data <- "cv_yhat_df"
  score_final_cv$AUROC <- ifelse(score_final_cv$AUROC>0.5, score_final_cv$AUROC, 1-score_final_cv$AUROC)
  
  
  # ---- return ----
  x_select_mdls_cluster <- list( lasso_coefs = coeff_lambda_trace,
                                 lasso_trace = score_lambda_trace[,c("lambda",tune_by)],
                                 lasso_trace_bonus = score_lambda_trace,
                                 lasso_optimal = lasso_optimal,
                                 cv_yhat_df = valid_yhat_df_all,
                                 score_final_cv = score_final_cv)
  
  return(x_select_mdls_cluster)
  
}



