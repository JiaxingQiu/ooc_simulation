library(gglasso)

lasso_x_select_group <- function(
  data, 
  y_col,
  x_cols_nonlin_rcs5=c(),
  x_cols_nonlin_rcs4=c(),
  x_cols_nonlin_rcs3=c(),
  x_cols_linear=c(), 
  x_cols_fct=c(),
  x_cols_tag=c(),
  family = c("binomial", "multinomial", "gaussian")[1],
  dict_data=NULL, # dictionary table is optional
  lambda=c("auto","1se","min")[1],
  lambda_value = NULL, # external specified lasso lambda value
  tune_by = c("misclass", 
              "logloss",
              "auc"),
  lambda_seq = NULL,
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
  
  # ---- preprocess ----
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
  
  x <- data.matrix(data[complete.cases(data[,c(x_col_df_all$x_colname,y_col)]),x_col_df_all$x_colname])
  y <- data.matrix(data[complete.cases(data[,c(x_col_df_all$x_colname,y_col)]),y_col])
  colnames(y) <- y_col
  y <- ifelse(y==1, 1, -1)
  #y <- as.factor(y)
  v.group <- as.numeric(as.factor(x_col_df_all$x_group))
  
  
  # ---- fix arguments (tune_by, lamdba_seq, foldid_vector cannot be empty) ----
  # --- tune_by ---
  if (tune_by%in%c("auc","misclass","logloss")){
    if (tune_by=="logloss") tune_by <- "loss"
  }else{
    # default
    tune_by <- "misclass"
  }
  stopifnot(tune_by%in%c("auc", "misclass", "loss", "AIC", "BIC"))
  # lambda_seq
  if(is.null(lambda_seq)){
    tryCatch({
      # Friedman, Hastie & Tibshirani (2010) 'strategy 
      yt <- ifelse(y>0,1,0)
      lambda_max <- max(abs(colSums(x*as.numeric(yt),na.rm=TRUE)), na.rm = TRUE)/dim(x)[1]
      epsilon <- .00001
      K <- 300
      lambda_seq <- round(exp(seq(log(lambda_max), log(lambda_max*epsilon),length.out = K)), digits = 10)
      print("Lambda sequence by Friedman, Hastie & Tibshirani (2010) strategy")
    },error=function(e){
      print("fail to generate lambda sequence manually, using default one in package")
      print(e)
    })
  }
  if(is.null(lambda_seq)){
    # get package default lambda sequence
    lasso_trace <- gglasso::gglasso(x=x, 
                                    y=y, 
                                    group=v.group, 
                                    lambda = NULL,
                                    loss="logit")
    lambda_seq <- lasso_trace$lambda
  }
  stopifnot(!is.null(lambda_seq))
  # --- foldid_vector ---
  foldid_vector <- NULL
  if(!is.null(foldid_col)){
    print("tuning lambda by external cv fold id")
    foldid_vector <- as.numeric(data[complete.cases(data[,c(x_col_df_all$x_colname,y_col)]),foldid_col])
  }else{
    print("tuning lambda by internal cv fold id")
    NR <- nrow(data[complete.cases(data[,c(x_col_df_all$x_colname,y_col)]),])
    foldid_vector <- c(1:NR)%/%ceiling(NR/10)+1
    foldid_vector[which(foldid_vector<1)] <- 1
    foldid_vector[which(foldid_vector>10)] <- 10
  }
  stopifnot(!is.null(foldid_vector))
  # ---- lasso regression tuning ----
  if(tune_by%in%c("misclass", "loss") ){
    # cv by package
    lasso_cv <- gglasso::cv.gglasso(x=x, 
                                    y=y, 
                                    group=v.group, 
                                    loss="logit", # deviance function for regular logistic function
                                    pred.loss = tune_by,#"misclass", #"misclass"
                                    lambda = lambda_seq,
                                    foldid = foldid_vector )
  }else if(tune_by%in%c("auc") ){
    getmin_joy <- function(lambda,cvm,cvsd){
      cvmin=min(cvm,na.rm=TRUE)
      idmin=cvm<=cvmin
      lambda.min=max(lambda[idmin],na.rm=TRUE)
      idmin=match(lambda.min,lambda)
      semin=(cvm+cvsd)[idmin]
      idmin=cvm<=semin
      lambda.1se=max(lambda[idmin],na.rm=TRUE)
      list(lambda.min=lambda.min,lambda.1se=lambda.1se)
    }
    gglasso.object <- gglasso::gglasso(x=x, 
                                       y=y, 
                                       group=v.group, 
                                       lambda = lambda_seq,
                                       loss="logit")
    lambda_seq <- gglasso.object$lambda
    foldid <- foldid_vector
    nfolds <- n_distinct(foldid)
    outlist <- as.list(seq(nfolds))
    ### fit the nfold models and store them
    for (i in seq(nfolds)) {
      which <- foldid == i
      y_sub <- y[!which]
      outlist[[i]] <- gglasso::gglasso(x = x[!which, , drop = FALSE], 
                                       y = y_sub, 
                                       group = v.group, 
                                       lambda = lambda_seq,
                                       loss="logit")
      
    }
    fun <- paste("gglasso_cv", class(gglasso.object)[[2]], sep = "_")
    cvstuff <- do.call(fun, list(outlist, lambda_seq, x, y, foldid, tune_by))
    cvm <- cvstuff$cvm
    cvsd <- cvstuff$cvsd
    cvname <- cvstuff$name
    out <- list(lambda = lambda_seq, cvm = cvm, cvsd = cvsd, cvupper = cvm + cvsd, 
                cvlo = cvm - cvsd, name = cvname, gglasso.fit = gglasso.object)
    lamin <- getmin_joy(lambda_seq, 1-cvm, cvsd)
    obj <- c(out, as.list(lamin))
    class(obj) <- "cv.gglasso"
    lasso_cv <- obj
  }
  
  # collect penalization trace
  lasso_trace <- gglasso::gglasso(x=x, 
                                  y=y, 
                                  group=v.group, 
                                  lambda = lambda_seq,
                                  loss="logit")
  
  #  get optimal lambda
  opt_lambda <- lasso_cv$lambda.1se
  if (lasso_cv$cvm[which(lasso_cv$lambda==lasso_cv$lambda.min)] == min(lasso_cv$cvm) ){
    opt_lambda <- lasso_cv$lambda.min
  }
  if(lambda %in% c("min","1se")){
    if(lambda=="min") {
      opt_lambda <- lasso_cv$lambda.min
    }else{
      opt_lambda <- lasso_cv$lambda.1se
    }
  }
  if(!is.null(lambda_value)){
    opt_lambda <- lambda_value
  }
  # ---- final optimal model with optimized lambda ----
  lasso_optimal <- gglasso::gglasso(x=x, 
                                    y=y, 
                                    group=v.group, 
                                    loss="logit",
                                    lambda = opt_lambda)
 
  lasso_optimal$x <- x
  lasso_optimal$y <- ifelse(y<0,0,1)
  lasso_optimal$group_info <- x_col_df_all
  lasso_optimal$opt_lambda <- opt_lambda
  
  # ------ manually 10 fold cross validation scores ---------
  valid_yhat_df_all <- data.frame()
  valid_yhat_df_permu_all <- data.frame()
  for(i in unique(data[,foldid_col])){
    try({
      validset <- data[which(data[,foldid_col]==i), ]
      trainset <- data[which(!data[,foldid_col]==i), ]
      
      train_x <- data.matrix(trainset[complete.cases(trainset[,c(x_col_df_all$x_colname,y_col)]),x_col_df_all$x_colname])
      train_y <- data.matrix(trainset[complete.cases(trainset[,c(x_col_df_all$x_colname,y_col)]),y_col])
      train_y <- ifelse(train_y==1, 1, -1)
      
      valid_x <- data.matrix(validset[complete.cases(validset[,c(x_col_df_all$x_colname,y_col)]),x_col_df_all$x_colname])
      valid_y <- data.matrix(validset[complete.cases(validset[,c(x_col_df_all$x_colname,y_col)]),y_col])
      valid_y <- ifelse(valid_y==1, 1, -1)
      lasso_fold_optimal <- NULL
      lasso_fold_optimal <- gglasso::gglasso(x=train_x, 
                                             y=train_y, 
                                             group=v.group, 
                                             loss="logit",
                                             lambda = opt_lambda)
     
      y_prob_valid <- exp(predict(lasso_fold_optimal, newx = valid_x, type="link")) / (1 + exp(predict(lasso_fold_optimal, newx = valid_x, type="link")))
      valid_yhat_df <- data.frame( 
        rowid = rownames( validset[complete.cases(validset[,c(x_col_df_all$x_colname,y_col)]),x_col_df_all$x_colname] ),
        y_true = as.numeric( ifelse(valid_y<0,0,1) ),
        y_prob =  as.numeric( y_prob_valid) ,
        fold = i
      )
      valid_yhat_df_all <- bind_rows(valid_yhat_df_all, valid_yhat_df)
      
      # predict on left out validation set with permuted predictors one by one
      for(x_group in unique(x_col_df_all$x_group) ){
        valid_x_permu <- valid_x
        # shuffle
        for(x_col in x_col_df_all$x_colname[which(x_col_df_all$x_group==x_group)] ){
          valid_x_permu[,x_col] <- sample(valid_x_permu[,x_col], size=length(valid_x_permu[,x_col]))
        }
        # make new prediction
        y_prob_valid_permu <- exp(predict(lasso_fold_optimal, newx = valid_x_permu, type="link")) / (1 + exp(predict(lasso_fold_optimal, newx = valid_x_permu, type="link")))
        valid_yhat_df_permu <- data.frame( 
          rowid = rownames( validset[complete.cases(validset[,c(x_col_df_all$x_colname,y_col)]),x_col_df_all$x_colname] ),
          y_true = as.numeric(ifelse(valid_y<0,0,1) ),
          y_prob =  as.numeric(y_prob_valid_permu),
          fold = i,
          permu_x = x_group
        )
        valid_yhat_df_permu_all <- bind_rows(valid_yhat_df_permu_all, valid_yhat_df_permu)
      }
      
    },TRUE)
  }
  
  # cross-validated labeling (yhat)
  cv_final_scores <- mdl_test(y_true = valid_yhat_df_all$y_true,
                              y_prob = valid_yhat_df_all$y_prob,
                              threshold = mean(data[,y_col],na.rm=TRUE))
  score_final_cv <- cv_final_scores$res_df
  score_final_cv$data <- "cv_yhat_df"
  score_final_cv$AUROC <- ifelse(score_final_cv$AUROC>0.5, score_final_cv$AUROC, 1-score_final_cv$AUROC)
  
  # cross-validated scores by permutation importance
  score_final_cv_permu <- data.frame()
  for(permu_x in unique(valid_yhat_df_permu_all$permu_x) ){
    s <- mdl_test(y_true = valid_yhat_df_permu_all$y_true[which(valid_yhat_df_permu_all$permu_x==permu_x)],
                  y_prob = valid_yhat_df_permu_all$y_prob[which(valid_yhat_df_permu_all$permu_x==permu_x)],
                  threshold = mean(data[,y_col],na.rm=TRUE))
    s <- s$res_df
    s$data <- paste0("permutate ",permu_x)
    score_final_cv_permu <- bind_rows(score_final_cv_permu, s)
  }
  
  x_select_mdls_grouped <- list( lasso_cv = lasso_cv,
                                 lasso_trace = lasso_trace,
                                 lasso_optimal = lasso_optimal,
                                 cv_yhat_df = valid_yhat_df_all,
                                 score_final_cv = score_final_cv,
                                 cv_yhat_df_permu = valid_yhat_df_permu_all, # permutation importance cv version
                                 score_final_cv_permu = score_final_cv_permu,
                                 tune_by = tune_by)
  
  
  # ---- return ----
  return(x_select_mdls_grouped)
  
}



