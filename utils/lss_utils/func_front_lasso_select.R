front_lasso_select <- function(
  data,
  dict_data,
  y_label, 
  cluster_label,
  x_labels_linear = c(),
  x_labels_nonlin_rcs5 = c(), 
  x_labels_nonlin_rcs4 = c(),
  x_labels_nonlin_rcs3 = c(),
  x_labels_fct = c(), 
  x_labels_tag = c(), 
  x_labels = unique(c(x_labels_linear,x_labels_nonlin_rcs5,x_labels_nonlin_rcs4,x_labels_nonlin_rcs3,x_labels_fct,x_labels_tag)), 
  # --- engineer ---
  trim_by_label=NULL,
  trim_vec = c(-Inf, Inf), # trim relative time [from, to)
  trim_first = TRUE, # TRUE = trim the data at the beginning before other engineering; FALSE = trim the data as the last step
  time_unit = 1, # the increment scale of relative time
  pctcut_num_labels = c(),
  pctcut_num_vec = c(0.1, 99.9),
  pctcut_num_coerce=TRUE,
  filter_tag_labels=c(),
  winsorizing=FALSE,
  aggregate_per=c("cluster_trim_by_unit", "cluster","row")[2], # should set to be cluster wise or time unit wise
  aggregate_conditioned_on_labels = c(),
  imputation=c("None","Mean", "Median", "Zero")[1],
  imputeby_median = c(),
  imputeby_zero = c(),
  imputeby_mean = c(), 
  impute_per_cluster=FALSE,
  standardize_df = "None", # can be c("None", "Standard(msd)", "Robust(IQR)", "Percentile", a data.frame(varname=c(), center=c(), scale=c())) 
  # --- local ---
  lambda=c("auto","min","1se")[1],
  lambda_value = NULL,
  lasso_by = c("none", "group", "cluster")[2],
  tune_by = c("auc", # area under the ROC curve, two-class logistic regression only
              "misclass", # misclassification error = 1-accuracy
              "logloss", # a.k.a "deviance"
              "AIC", 
              "BIC")[1],
  lambda_seq=NULL,
  trim_ctrl = TRUE,
  test_data=NULL,
  y_map_func=c("fold_risk", "probability", "log_odds")[1],
  y_map_max=3,
  return_performance = TRUE,
  return_effect_plots = FALSE,
  fold_idx_df_ex=NULL # data.frame(cluster_col = c(), fold = c())
){
  
  x_select_mdls <- NULL
  x_select_mdls_grouped <- NULL
  x_select_mdls_cluster <- NULL
  
  # ---- pre-process ----
  x_cols_linear <- dict_data$varname[which(dict_data$label%in%x_labels_linear&dict_data$type=="num")]# linear numeric columns #dict_data$mlrole=="input"&
  x_cols_nonlin_rcs5 <- dict_data$varname[which(dict_data$label%in%x_labels_nonlin_rcs5&dict_data$type=="num")]
  x_cols_nonlin_rcs4 <- dict_data$varname[which(dict_data$label%in%x_labels_nonlin_rcs4&dict_data$type=="num")]
  x_cols_nonlin_rcs3 <- dict_data$varname[which(dict_data$label%in%x_labels_nonlin_rcs3&dict_data$type=="num")]
  x_cols_fct <- dict_data$varname[which(dict_data$label%in%x_labels_fct & dict_data$type=="fct" & dict_data$unit!="tag01")]
  x_cols_tag <- dict_data$varname[which(dict_data$label%in%x_labels_tag & dict_data$type=="fct" & dict_data$unit=="tag01")]
  y_col_tag <- dict_data$varname[which(dict_data$label==y_label & dict_data$type=="fct" & dict_data$unit=="tag01" )]
  y_col_num <- dict_data$varname[which(dict_data$label==y_label & dict_data$type=="num" )]
  y_col <- union(y_col_tag, y_col_num)
  fct_cols <- unique(c(x_cols_fct, x_cols_tag, y_col_tag)) # group columns together by typr for data engineering
  num_cols <- unique(c(x_cols_linear, x_cols_nonlin_rcs3,x_cols_nonlin_rcs4,x_cols_nonlin_rcs5, y_col_num))
  cluster_col <- dict_data$varname[which(dict_data$label==cluster_label)]
  rel_time_col <- dict_data$varname[which(dict_data$label==trim_by_label)]
  trim_by_col <- dict_data$varname[which(dict_data$label==trim_by_label)]
  pctcut_num_cols <- dict_data$varname[which(dict_data$label%in%pctcut_num_labels)]
  filter_tag_cols <- dict_data$varname[which(dict_data$label%in%filter_tag_labels)]
  aggregate_conditioned_on_cols <- intersect(colnames(data), dict_data$varname[which(dict_data$label %in% aggregate_conditioned_on_labels)])
  
  # ---- create standardize_df from input data ----
  # standardize_df can be can be c("None", "Standard(msd)", "Robust(IQR)", "Percentile", a data.frame(varname=c(), center=c(), scale=c())) 
  
  if( paste0(as.character(standardize_df),collapse = "")=="None" ){ # standardize_df is NULL
    print("--- Not standardize ---")
    standardize_df <- NULL
  }
  else if(paste0(as.character(standardize_df),collapse = "")=="Standard(msd)"){
    print("--- Make standardize_df of Standard(msd) from internal engineered dataset ---")
    standardize_df <- NULL
    print("---- 1. engineer internal data without standardization----")
    tryCatch({
      if(trim_ctrl){
        data_in <- engineer(data = data,
                            num_cols = num_cols,
                            fct_cols = fct_cols,
                            cluster_col = cluster_col,
                            trim_by_col = trim_by_col,
                            trim_min = trim_vec[1],
                            trim_max = trim_vec[2],
                            trim_first = trim_first,
                            trim_step_size = time_unit,
                            pctcut_num_cols = pctcut_num_cols,
                            pctcut_num_vec = pctcut_num_vec,
                            pctcut_num_coerce = pctcut_num_coerce,
                            filter_tag_cols = filter_tag_cols,
                            imputation = imputation,
                            imputeby_median = imputeby_median,
                            imputeby_zero = imputeby_zero,
                            imputeby_mean = imputeby_mean, 
                            impute_per_cluster = impute_per_cluster,
                            winsorizing = winsorizing,
                            aggregate_per = aggregate_per,
                            aggregate_conditioned_on_cols = aggregate_conditioned_on_cols,
                            standardize_df = standardize_df) # standardize_df is null here
      }else{
        if (all(unique(as.character(data[,y_col])) %in% c(1,0,NA))){
          data_event <- engineer(data = data[which(data[,y_col]==1),],
                                 num_cols = num_cols,
                                 fct_cols = fct_cols,
                                 cluster_col = cluster_col,
                                 trim_by_col = trim_by_col,
                                 trim_min=trim_vec[1],
                                 trim_max=trim_vec[2],
                                 trim_first = trim_first,
                                 trim_step_size = time_unit,
                                 pctcut_num_cols = pctcut_num_cols,
                                 pctcut_num_vec = pctcut_num_vec,
                                 pctcut_num_coerce = pctcut_num_coerce,
                                 filter_tag_cols = filter_tag_cols,
                                 imputation = imputation,
                                 imputeby_median = imputeby_median,
                                 imputeby_zero = imputeby_zero,
                                 imputeby_mean = imputeby_mean, 
                                 impute_per_cluster = impute_per_cluster,
                                 winsorizing = winsorizing,
                                 aggregate_per = aggregate_per,
                                 aggregate_conditioned_on_cols=aggregate_conditioned_on_cols,
                                 standardize_df = standardize_df)# standardize_df is null here
          data_cntrl <- engineer(data = data[which(data[,y_col]==0),],
                                 num_cols = num_cols,
                                 fct_cols = fct_cols,
                                 cluster_col = cluster_col,
                                 trim_by_col = trim_by_col,
                                 trim_min=-Inf,
                                 trim_max=Inf,
                                 trim_step_size = time_unit,
                                 trim_keepna = TRUE,
                                 pctcut_num_cols = pctcut_num_cols,
                                 pctcut_num_vec = pctcut_num_vec,
                                 pctcut_num_coerce = pctcut_num_coerce,
                                 filter_tag_cols = filter_tag_cols,
                                 imputation = imputation,
                                 imputeby_median = imputeby_median,
                                 imputeby_zero = imputeby_zero,
                                 imputeby_mean = imputeby_mean, 
                                 impute_per_cluster = impute_per_cluster,
                                 winsorizing = winsorizing,
                                 aggregate_per = aggregate_per,
                                 aggregate_conditioned_on_cols=aggregate_conditioned_on_cols,
                                 standardize_df = standardize_df)# standardize_df is null here
          data_in <- bind_rows(data_cntrl, data_event)
        }
      }
      data_in <- assign.dict(data_in, dict_data)
    },error=function(e){
      print("Error!")
      print(e)
    })
    print("---- 2. generate standardization dataframe ----")
    tryCatch({
      standardize_df <- data.frame(varname = num_cols, 
                                   center=apply(data_in[,num_cols],2,mean,na.rm=TRUE),
                                   scale=apply(data_in[,num_cols],2,sd,na.rm=TRUE))
    },error=function(e){
      print("Error!")
      print(e)
    })
  }
  else if(paste0(as.character(standardize_df),collapse = "")=="Robust(IQR)"){
    print("--- Make standardize_df of Robust(IQR) from internal engineered dataset ---")
    standardize_df <- NULL
    print("---- 1. engineer internal data without standardization----")
    tryCatch({
      if(trim_ctrl){
        data_in <- engineer(data = data,
                            num_cols = num_cols,
                            fct_cols = fct_cols,
                            cluster_col = cluster_col,
                            trim_by_col = trim_by_col,
                            trim_min = trim_vec[1],
                            trim_max = trim_vec[2],
                            trim_first = trim_first,
                            trim_step_size = time_unit,
                            pctcut_num_cols = pctcut_num_cols,
                            pctcut_num_vec = pctcut_num_vec,
                            pctcut_num_coerce = pctcut_num_coerce,
                            filter_tag_cols = filter_tag_cols,
                            imputation = imputation,
                            imputeby_median = imputeby_median,
                            imputeby_zero = imputeby_zero,
                            imputeby_mean = imputeby_mean, 
                            impute_per_cluster = impute_per_cluster,
                            winsorizing = winsorizing,
                            aggregate_per = aggregate_per,
                            aggregate_conditioned_on_cols=aggregate_conditioned_on_cols,
                            standardize_df = standardize_df) # standardize_df is null here
      }else{
        if (all(unique(as.character(data[,y_col])) %in% c(1,0,NA))){
          data_event <- engineer(data = data[which(data[,y_col]==1),],
                                 num_cols = num_cols,
                                 fct_cols = fct_cols,
                                 cluster_col = cluster_col,
                                 trim_by_col = trim_by_col,
                                 trim_min=trim_vec[1],
                                 trim_max=trim_vec[2],
                                 trim_first = trim_first,
                                 trim_step_size = time_unit,
                                 pctcut_num_cols = pctcut_num_cols,
                                 pctcut_num_vec = pctcut_num_vec,
                                 pctcut_num_coerce = pctcut_num_coerce,
                                 filter_tag_cols = filter_tag_cols,
                                 imputation = imputation,
                                 imputeby_median = imputeby_median,
                                 imputeby_zero = imputeby_zero,
                                 imputeby_mean = imputeby_mean, 
                                 impute_per_cluster = impute_per_cluster,
                                 winsorizing = winsorizing,
                                 aggregate_per = aggregate_per,
                                 aggregate_conditioned_on_cols=aggregate_conditioned_on_cols,
                                 standardize_df = standardize_df)# standardize_df is null here
          data_cntrl <- engineer(data = data[which(data[,y_col]==0),],
                                 num_cols = num_cols,
                                 fct_cols = fct_cols,
                                 cluster_col = cluster_col,
                                 trim_by_col = trim_by_col,
                                 trim_min=-Inf,
                                 trim_max=Inf,
                                 trim_step_size = time_unit,
                                 trim_keepna = TRUE,
                                 pctcut_num_cols = pctcut_num_cols,
                                 pctcut_num_vec = pctcut_num_vec,
                                 pctcut_num_coerce = pctcut_num_coerce,
                                 filter_tag_cols = filter_tag_cols,
                                 imputation = imputation,
                                 imputeby_median = imputeby_median,
                                 imputeby_zero = imputeby_zero,
                                 imputeby_mean = imputeby_mean, 
                                 impute_per_cluster = impute_per_cluster,
                                 winsorizing = winsorizing,
                                 aggregate_per = aggregate_per,
                                 aggregate_conditioned_on_cols=aggregate_conditioned_on_cols,
                                 standardize_df = standardize_df)# standardize_df is null here
          data_in <- bind_rows(data_cntrl, data_event)
        }
      }
      data_in <- assign.dict(data_in, dict_data)
    },error=function(e){
      print("Error!")
      print(e)
    })
    print("---- 2. generate standardization dataframe ----")
    tryCatch({
      standardize_df <- data.frame(varname = num_cols, 
                                   center=apply(data_in[,num_cols],2,median,na.rm=TRUE),
                                   scale=(apply(data_in[,num_cols],2,quantile,0.75,na.rm=TRUE)-apply(data_in[,num_cols],2,quantile,0.25,na.rm=TRUE)))
    },error=function(e){
      print("Error!")
      print(e)
    })
  }
  else if(paste0(as.character(standardize_df),collapse = "")=="Percentile"){
    print("--- Convert engineered data to percentile ---")
    standardize_df <- NULL
    tryCatch({
      data <- engineer(data = data,
                       num_cols = num_cols,
                       fct_cols = fct_cols,
                       cluster_col = cluster_col,
                       trim_by_col = trim_by_col,
                       trim_min = trim_vec[1],
                       trim_max = trim_vec[2],
                       trim_first = trim_first,
                       trim_step_size = time_unit,
                       pctcut_num_cols = pctcut_num_cols,
                       pctcut_num_vec = pctcut_num_vec,
                       pctcut_num_coerce = pctcut_num_coerce,
                       filter_tag_cols = filter_tag_cols,
                       imputation = imputation,
                       imputeby_median = imputeby_median,
                       imputeby_zero = imputeby_zero,
                       imputeby_mean = imputeby_mean, 
                       impute_per_cluster = impute_per_cluster,
                       winsorizing = winsorizing,
                       aggregate_per = aggregate_per,
                       aggregate_conditioned_on_cols=aggregate_conditioned_on_cols,
                       standardize_df = NULL) # standardize_df is null here
      for(num_col in num_cols){
        tryCatch({
          data[[num_col]] <- est_pctl(data[[num_col]])
        },error=function(e){
          print(paste0("skip percentiling",num_col))
          print(e)
        })
      }
    },error=function(e){
      print("Error!")
      print(e)
    })
  }
  
  # until now standardize_df can be NULL or a dataframe 
  
  
  # ---- engineering ----
  data_in <- NULL
  data_inorg <- NULL
  data_ex <- NULL
  data_exorg <- NULL
  # --- prepare original / no-engineering train and validation dataset (internal)  ---
  print("---- prepare original / no-engineering train and validation dataset (internal)  ----")
  tryCatch({
    data_inorg <- engineer(data = data,
                           num_cols = num_cols,
                           fct_cols = fct_cols,
                           cluster_col = cluster_col,
                           trim_by_col = trim_by_col,
                           standardize_df = standardize_df) # minimum data engineering = using default settings
    data_inorg <- assign.dict(data_inorg, dict_data)
  },error=function(e){
    print("Error!")
    print(e)
  })
  # --- prepare engineered training and validation dataset (internal) ---
  print("---- prepare engineered train and validation dataset (internal) ----")
  tryCatch({
    if(trim_ctrl){
      data_in <- engineer(data = data,
                          num_cols = num_cols,
                          fct_cols = fct_cols,
                          cluster_col = cluster_col,
                          trim_by_col = trim_by_col,
                          trim_min = trim_vec[1],
                          trim_max = trim_vec[2],
                          trim_first = trim_first,
                          trim_step_size = time_unit,
                          pctcut_num_cols = pctcut_num_cols,
                          pctcut_num_vec = pctcut_num_vec,
                          pctcut_num_coerce = pctcut_num_coerce,
                          filter_tag_cols = filter_tag_cols,
                          imputation = imputation,
                          imputeby_median = imputeby_median,
                          imputeby_zero = imputeby_zero,
                          imputeby_mean = imputeby_mean, 
                          impute_per_cluster = impute_per_cluster,
                          winsorizing = winsorizing,
                          aggregate_per = aggregate_per,
                          aggregate_conditioned_on_cols=aggregate_conditioned_on_cols,
                          standardize_df = standardize_df)
    }else{
      if (all(unique(as.character(data[,y_col])) %in% c(1,0,NA))){
        data_event <- engineer(data = data[which(data[,y_col]==1),],
                               num_cols = num_cols,
                               fct_cols = fct_cols,
                               cluster_col = cluster_col,
                               trim_by_col = trim_by_col,
                               trim_min=trim_vec[1],
                               trim_max=trim_vec[2],
                               trim_first = trim_first,
                               trim_step_size = time_unit,
                               pctcut_num_cols = pctcut_num_cols,
                               pctcut_num_vec = pctcut_num_vec,
                               pctcut_num_coerce = pctcut_num_coerce,
                               filter_tag_cols = filter_tag_cols,
                               imputation = imputation,
                               imputeby_median = imputeby_median,
                               imputeby_zero = imputeby_zero,
                               imputeby_mean = imputeby_mean, 
                               impute_per_cluster = impute_per_cluster,
                               winsorizing = winsorizing,
                               aggregate_per = aggregate_per,
                               aggregate_conditioned_on_cols=aggregate_conditioned_on_cols,
                               standardize_df = standardize_df)
        data_cntrl <- engineer(data = data[which(data[,y_col]==0),],
                               num_cols = num_cols,
                               fct_cols = fct_cols,
                               cluster_col = cluster_col,
                               trim_by_col = trim_by_col,
                               trim_min=-Inf,
                               trim_max=Inf,
                               trim_step_size = time_unit,
                               trim_keepna = TRUE,
                               pctcut_num_cols = pctcut_num_cols,
                               pctcut_num_vec = pctcut_num_vec,
                               pctcut_num_coerce = pctcut_num_coerce,
                               filter_tag_cols = filter_tag_cols,
                               imputation = imputation,
                               imputeby_median = imputeby_median,
                               imputeby_zero = imputeby_zero,
                               imputeby_mean = imputeby_mean, 
                               impute_per_cluster = impute_per_cluster,
                               winsorizing = winsorizing,
                               aggregate_per = aggregate_per,
                               aggregate_conditioned_on_cols=aggregate_conditioned_on_cols,
                               standardize_df = standardize_df)
        data_in <- bind_rows(data_cntrl, data_event)
      }
    }
    data_in <- assign.dict(data_in, dict_data)
  },error=function(e){
    print("Error!")
    print(e)
  })
  # --- prepare original / no-engineering test dataset (external)  ---
  print("---- prepare original / no-engineering test dataset (external)  ----")
  tryCatch({
    if (!is.null(test_data)){# if test_data given
      data_exorg <- engineer(data = test_data,
                             num_cols = num_cols,
                             fct_cols = fct_cols,
                             cluster_col = cluster_col,
                             trim_by_col = trim_by_col,
                             standardize_df = standardize_df)
      data_exorg <- assign.dict(data_exorg, dict_data)
    }
  },error=function(e){
    print("Error!")
    print(e)
  })
  
  # --- prepare engineered test dataset (external)  ---
  print("---- prepare engineered test dataset (external)  ----")
  tryCatch({
    if (!is.null(test_data)){
      if(trim_ctrl){
        data_ex <- engineer(data = test_data,
                            num_cols = num_cols,
                            fct_cols = fct_cols,
                            cluster_col = cluster_col,
                            trim_by_col = trim_by_col,
                            trim_min=trim_vec[1],
                            trim_max=trim_vec[2],
                            trim_first = trim_first,
                            trim_step_size = time_unit,
                            pctcut_num_cols = pctcut_num_cols,
                            pctcut_num_vec = pctcut_num_vec,
                            pctcut_num_coerce = pctcut_num_coerce,
                            filter_tag_cols = filter_tag_cols,
                            imputation = imputation,
                            imputeby_median = imputeby_median,
                            imputeby_zero = imputeby_zero,
                            imputeby_mean = imputeby_mean, 
                            impute_per_cluster = impute_per_cluster,
                            winsorizing = winsorizing,
                            aggregate_per = aggregate_per,
                            aggregate_conditioned_on_cols=aggregate_conditioned_on_cols,
                            standardize_df = standardize_df)
      }else{
        if (all(unique(as.character(test_data[,y_col])) %in% c(1,0,NA))){
          data_event <- engineer(data = test_data[which(test_data[,y_col]==1),],
                                 num_cols = num_cols,
                                 fct_cols = fct_cols,
                                 cluster_col = cluster_col,
                                 trim_by_col = trim_by_col,
                                 trim_min=trim_vec[1],
                                 trim_max=trim_vec[2],
                                 trim_first = trim_first,
                                 trim_step_size = time_unit,
                                 pctcut_num_cols = pctcut_num_cols,
                                 pctcut_num_vec = pctcut_num_vec,
                                 pctcut_num_coerce = pctcut_num_coerce,
                                 filter_tag_cols = filter_tag_cols,
                                 imputation = imputation,
                                 imputeby_median = imputeby_median,
                                 imputeby_zero = imputeby_zero,
                                 imputeby_mean = imputeby_mean, 
                                 impute_per_cluster = impute_per_cluster,
                                 winsorizing = winsorizing,
                                 aggregate_per = aggregate_per,
                                 aggregate_conditioned_on_cols=aggregate_conditioned_on_cols,
                                 standardize_df = standardize_df)
          data_cntrl <- engineer(data = test_data[which(test_data[,y_col]==0),],
                                 num_cols = num_cols,
                                 fct_cols = fct_cols,
                                 cluster_col = cluster_col,
                                 trim_by_col = trim_by_col,
                                 trim_min=-Inf,
                                 trim_max=Inf,
                                 trim_step_size = time_unit,
                                 trim_keepna = TRUE,
                                 pctcut_num_cols = pctcut_num_cols,
                                 pctcut_num_vec = pctcut_num_vec,
                                 pctcut_num_coerce = pctcut_num_coerce,
                                 filter_tag_cols = filter_tag_cols,
                                 imputation = imputation,
                                 imputeby_median = imputeby_median,
                                 imputeby_zero = imputeby_zero,
                                 imputeby_mean = imputeby_mean, 
                                 impute_per_cluster = impute_per_cluster,
                                 winsorizing = winsorizing,
                                 aggregate_per = aggregate_per,
                                 aggregate_conditioned_on_cols=aggregate_conditioned_on_cols,
                                 standardize_df = standardize_df)
          data_ex <- bind_rows(data_cntrl, data_event)
        }
      }
      data_ex <- assign.dict(data_ex, dict_data)
    }
  },error=function(e){
    print("Error!")
    print(e)
  })
  
  
  # ---- add cv fold id columns ----
  if(is.null(fold_idx_df_ex)){
    fold_idx_df_ex <- data.frame(cluster_col=unique(data_in[,cluster_col]), 
                            fold=ceiling(c(1:n_distinct(data_in[,cluster_col]))/round(n_distinct(data_in[,cluster_col])/10) ) )
    fold_idx_df_ex$fold <- sample(fold_idx_df_ex$fold, length(fold_idx_df_ex$fold), replace = FALSE)
    fold_idx_df_ex$fold[which(fold_idx_df_ex$fold>10)] <- 10
    stopifnot(n_distinct(fold_idx_df_ex$fold)<=10)
  }else{
    print("using external fold indicator instead!")
  }
  colnames(fold_idx_df_ex)[which(colnames(fold_idx_df_ex)=="cluster_col")] <- cluster_col
  fold_idx_df_ex <- dplyr::distinct(fold_idx_df_ex) # keep unique cluster - fold mapping
  data_in <- merge(data_in[,setdiff(colnames(data_in),"fold")], fold_idx_df_ex, all.x=TRUE)
  
  
  # ---- lasso_x_select ----
  if(length(y_col_tag)>0 ){
    family <- "binomial"
  } else if(length(y_col_num) > 0){
    family <- "gaussian"
  }
  print("---- lasso_x_select ----")
  tryCatch({
    stopifnot(lasso_by=="none")
    x_select_mdls <- lasso_x_select(data = data_in,
                                    y_col = y_col,
                                    x_cols_nonlin_rcs3 = x_cols_nonlin_rcs3,
                                    x_cols_nonlin_rcs4 = x_cols_nonlin_rcs4,
                                    x_cols_nonlin_rcs5 = x_cols_nonlin_rcs5,
                                    x_cols_linear = x_cols_linear, 
                                    x_cols_fct = x_cols_fct,
                                    x_cols_tag = x_cols_tag,
                                    family = family,
                                    dict_data = dict_data,
                                    lambda = lambda,
                                    lambda_value = lambda_value,
                                    tune_by = tune_by,
                                    lambda_seq = lambda_seq,
                                    foldid_col = "fold")
  },error=function(e){
    print("Not run lasso_x_select")
    print(e)
  })
  
  
  # ---- lasso_x_select_group ----
  print("---- lasso_x_select_group ----")
  tryCatch({
    stopifnot(lasso_by=="group")
    x_select_mdls_grouped <- lasso_x_select_group(data = data_in,
                                                  y_col = y_col,
                                                  x_cols_nonlin_rcs3 = x_cols_nonlin_rcs3,
                                                  x_cols_nonlin_rcs4 = x_cols_nonlin_rcs4,
                                                  x_cols_nonlin_rcs5 = x_cols_nonlin_rcs5,
                                                  x_cols_linear = x_cols_linear, 
                                                  x_cols_fct = x_cols_fct,
                                                  x_cols_tag = x_cols_tag,
                                                  family = family,
                                                  dict_data = dict_data,
                                                  lambda = lambda,
                                                  lambda_value = lambda_value,
                                                  tune_by = tune_by,
                                                  lambda_seq = lambda_seq,
                                                  foldid_col = "fold")
    
  },error=function(e){
    print("Not run lasso_x_select_group")
    print(e)
  })

  
  # ---- lasso_x_select_cluster ----
  print("---- lasso_x_select_cluster ----")
  tryCatch({
    stopifnot(lasso_by=="cluster")
    x_select_mdls_cluster <- lasso_x_select_cluster(data = data_in,
                                                  y_col = y_col,
                                                  cluster_col = cluster_col,
                                                  family = family,
                                                  x_cols_nonlin_rcs3 = x_cols_nonlin_rcs3,
                                                  x_cols_nonlin_rcs4 = x_cols_nonlin_rcs4,
                                                  x_cols_nonlin_rcs5 = x_cols_nonlin_rcs5,
                                                  x_cols_linear = x_cols_linear, 
                                                  x_cols_fct = x_cols_fct,
                                                  x_cols_tag = x_cols_tag,
                                                  dict_data = dict_data,
                                                  tune_by = tune_by,
                                                  lambda_seq = lambda_seq,
                                                  foldid_col = "fold" )
  },error=function(e){
    print("Not run lasso_x_select_cluster")
    print(e)
  })
  
  # ---- inference ----
  if(lasso_by=="none") lasso_obj <- list(tune_trace=x_select_mdls$cv_mdls$lasso_cv,
                                         coef_trace=x_select_mdls$trace_mdls$lasso_trace,
                                         opt_model=x_select_mdls$optimal_mdls$lasso_optimal)
  if(lasso_by=="group") lasso_obj <- list(tune_trace=x_select_mdls_grouped$lasso_cv,
                                          coef_trace=x_select_mdls_grouped$lasso_trace,
                                          opt_model=x_select_mdls_grouped$lasso_optimal)
  if(lasso_by=="cluster") lasso_obj <- list(tune_trace=x_select_mdls_cluster$lasso_trace,
                                            coef_trace=x_select_mdls_cluster$lasso_coefs,
                                            opt_model=x_select_mdls_cluster$lasso_optimal)
  
  infer_obj <- lss_infer(lasso_obj, lasso_by)
  
  
  
  
  
  # ---- performances ----
  cv_scores_all_final <- NULL
  if(lasso_by=="none") {
    mdl_obj <- x_select_mdls$optimal_mdls$lasso_optimal
    cv_scores_all_final <- x_select_mdls$score_final_cv_permu
    cv_scores_all_final_none <- x_select_mdls$score_final_cv
    cv_scores_all_final_none$data <- "none"
    cv_scores_all_final <- bind_rows(cv_scores_all_final, cv_scores_all_final_none)
  }
  if(lasso_by=="group") {
    mdl_obj <- x_select_mdls_grouped$lasso_optimal
    cv_scores_all_final <- x_select_mdls_grouped$score_final_cv_permu
    cv_scores_all_final_none <- x_select_mdls_grouped$score_final_cv
    cv_scores_all_final_none$data <- "none"
    cv_scores_all_final <- bind_rows(cv_scores_all_final, cv_scores_all_final_none)
  } 
  if(lasso_by=="cluster") {
    mdl_obj <- x_select_mdls_cluster$lasso_optimal
  }
  if(!is.null(cv_scores_all_final)){
    cv_scores_all_final$data <- gsub("_rcs[0-9]+","",gsub("_linear","",gsub("permutate ","",cv_scores_all_final$data)))
    colnames(cv_scores_all_final)[which(colnames(cv_scores_all_final)=="data")] <- "removed_variable"
    # special case of LASSO, remove zeroed out coefs
    opt_df <- infer_obj$opt_model_df
    selected_vars <- opt_df$varname[which(!opt_df$coef==0)]
    permu_vars <- c()
    for(v in cv_scores_all_final$removed_variable){
      if(any(grepl(v, selected_vars))){
        permu_vars <- c(permu_vars,v)
      }
    }
    cv_scores_all_final <- cv_scores_all_final[which(cv_scores_all_final$removed_variable%in%c("none",permu_vars)),]
  }
  
  
  print("----- lss_perform in -----")
  lss_perform_in <- NULL
  tryCatch({
    stopifnot(return_performance)
    lss_perform_in <- lss_perform(
      mdl_obj = mdl_obj, # a grouped lasso regression model object
      df = data_in, # a dataset to test out performance on
      y_map_func = y_map_func, # response type
      y_map_max = y_map_max, # response upper cutoff
      rel_time_col=rel_time_col,
      return_effect_plots=return_effect_plots,
      lasso_by = lasso_by,
      cv_scores_all_final = cv_scores_all_final
    )
  },error=function(e){
    print("Error!")
    print(e)
  })
  
  print("----- lss_perform inorg -----")
  lss_perform_inorg <- NULL
  tryCatch({
    stopifnot(return_performance)
    lss_perform_inorg <- lss_perform(
      mdl_obj = mdl_obj , # a grouped lasso regression model object
      df = data_inorg, # a dataset to test out performance on
      y_map_func = y_map_func, # response type
      y_map_max = y_map_max, # response upper cutoff
      rel_time_col=rel_time_col,
      return_effect_plots=return_effect_plots,
      lasso_by = lasso_by
    )
  },error=function(e){
    print("Error!")
    print(e)
  })
  
  print("----- lss_perform ex -----")
  lss_perform_ex <- NULL
  tryCatch({
    stopifnot(return_performance)
    lss_perform_ex <- lss_perform(
      mdl_obj = mdl_obj , # a grouped lasso regression model object
      df = data_ex, # a dataset to test out performance on
      y_map_func = y_map_func, # response type
      y_map_max = y_map_max, # response upper cutoff
      rel_time_col=rel_time_col,
      return_effect_plots=return_effect_plots,
      lasso_by = lasso_by
    )
  },error=function(e){
    print("Error!")
    print(e)
  })
  
  print("----- lss_perform exorg -----")
  lss_perform_exorg <- NULL
  tryCatch({
    stopifnot(return_performance)
    lss_perform_exorg <- lss_perform(
      mdl_obj = mdl_obj , # a grouped lasso regression model object
      df = data_exorg, # a dataset to test out performance on
      y_map_func = y_map_func, # response type
      y_map_max = y_map_max, # response upper cutoff
      rel_time_col=rel_time_col,
      return_effect_plots=return_effect_plots,
      lasso_by = lasso_by
    )
  },error=function(e){
    print("Error!")
    print(e)
  })
  # model performance reports on new dataset 
  perform_in_df_hat <- lss_perform_in$df_hat
  perform_inorg_df_hat <- lss_perform_inorg$df_hat
  perform_ex_df_hat <- lss_perform_ex$df_hat
  perform_exorg_df_hat <- lss_perform_exorg$df_hat
  
  perform_in_fitted_eff_plot <- lss_perform_in$fitted_eff_plot
  perform_inorg_fitted_eff_plot <- lss_perform_inorg$fitted_eff_plot
  perform_ex_fitted_eff_plot <- lss_perform_ex$fitted_eff_plot
  perform_exorg_fitted_eff_plot <- lss_perform_exorg$fitted_eff_plot
  
  perform_in_tte_plot <- lss_perform_in$tte_plot
  perform_inorg_tte_plot <- lss_perform_inorg$tte_plot
  perform_ex_tte_plot <- lss_perform_ex$tte_plot
  perform_exorg_tte_plot <- lss_perform_exorg$tte_plot
  
  perform_in_cali_plot <- lss_perform_in$cali_plot
  perform_inorg_cali_plot <- lss_perform_inorg$cali_plot
  perform_ex_cali_plot <- lss_perform_ex$cali_plot
  perform_exorg_cali_plot <- lss_perform_exorg$cali_plot
  
  perform_in_scores_plot <- lss_perform_in$scores_plot
  perform_inorg_scores_plot <- lss_perform_inorg$scores_plot
  perform_ex_scores_plot <- lss_perform_ex$scores_plot
  perform_exorg_scores_plot<- lss_perform_exorg$scores_plot
  
  perform_in_scores_tbl <- lss_perform_in$scores_all_final
  perform_inorg_scores_tbl <- lss_perform_inorg$scores_all_final
  perform_ex_scores_tbl <- lss_perform_ex$scores_all_final
  perform_exorg_scores_tbl <- lss_perform_exorg$scores_all_final
  
  perform_in_tradeoff_plot <- lss_perform_in$tradeoff_plot
  perform_inorg_tradeoff_plot <- lss_perform_inorg$tradeoff_plot
  perform_ex_tradeoff_plot <- lss_perform_ex$tradeoff_plot
  perform_exorg_tradeoff_plot <- lss_perform_exorg$tradeoff_plot
  
  return( list(mdl_obj = list(x_select_mdls = x_select_mdls,
                              x_select_mdls_grouped = x_select_mdls_grouped,
                              x_select_mdls_cluster = x_select_mdls_cluster),
               infer_obj = list(tune_trace_plot = infer_obj$tune_trace_plot,
                                coef_trace_plot = infer_obj$coef_trace_plot,
                                opt_model_df = infer_obj$opt_model_df,
                                lambda_zero_coef = infer_obj$lambda_zero_coef,
                                lasso_coef_trace = infer_obj$lasso_coef ),
               perform_obj = list(perform_in_df_hat = perform_in_df_hat,
                                  perform_inorg_df_hat = perform_inorg_df_hat,
                                  perform_ex_df_hat = perform_ex_df_hat,
                                  perform_exorg_df_hat = perform_exorg_df_hat, 
                                  perform_in_fitted_eff_plot = perform_in_fitted_eff_plot,
                                  perform_inorg_fitted_eff_plot = perform_inorg_fitted_eff_plot,
                                  perform_ex_fitted_eff_plot =perform_ex_fitted_eff_plot,
                                  perform_exorg_fitted_eff_plot = perform_exorg_fitted_eff_plot,
                                  perform_in_tte_plot = perform_in_tte_plot,
                                  perform_inorg_tte_plot = perform_inorg_tte_plot,
                                  perform_ex_tte_plot = perform_ex_tte_plot,
                                  perform_exorg_tte_plot = perform_exorg_tte_plot,
                                  perform_in_cali_plot = perform_in_cali_plot,
                                  perform_inorg_cali_plot = perform_inorg_cali_plot,
                                  perform_ex_cali_plot =perform_ex_cali_plot,
                                  perform_exorg_cali_plot = perform_exorg_cali_plot,
                                  perform_in_scores_plot = perform_in_scores_plot,
                                  perform_inorg_scores_plot = perform_inorg_scores_plot,
                                  perform_ex_scores_plot = perform_ex_scores_plot,
                                  perform_exorg_scores_plot = perform_exorg_scores_plot,
                                  perform_in_scores_tbl = perform_in_scores_tbl,
                                  perform_inorg_scores_tbl = perform_inorg_scores_tbl,
                                  perform_ex_scores_tbl = perform_ex_scores_tbl,
                                  perform_exorg_scores_tbl = perform_exorg_scores_tbl,
                                  perform_in_tradeoff_plot = perform_in_tradeoff_plot,
                                  perform_inorg_tradeoff_plot = perform_inorg_tradeoff_plot,
                                  perform_ex_tradeoff_plot = perform_ex_tradeoff_plot,
                                  perform_exorg_tradeoff_plot = perform_exorg_tradeoff_plot)
  ))
}





# ##################### not run #############################
# # variable to use
# x_vars_linear = c("demog_mage",
#                   "apgar1min"
# )
# x_vars_tag <- c("baby_gender_factor_Male",
#                 "baby_multiple_factor_Yes")
# x_vars_fct <- c("baby_insurance_factor")
# 
# 
# data = data_ml[which(data_ml$ca_days==7),]
# dict_data = dict_ml
# y_label = dict_ml[which(dict_ml$varname=="primary_outcome_factor_Unfavorable"),"label"]
# cluster_label = dict_ml[which(dict_ml$varname=="m_id"),"label"]
# x_labels_linear = dict_ml$label[which(dict_ml$varname%in%x_vars_linear)]
# x_labels_nonlin_rcs5 = c("Maternal age")
# x_labels_nonlin_rcs4 = c("Gestational Age")
# x_labels_nonlin_rcs3 = c("Birth weight")
# x_labels_fct = dict_ml$label[which(dict_ml$varname%in%x_vars_fct)]
# x_labels_tag = dict_ml$label[which(dict_ml$varname%in%x_vars_tag)]
# x_labels = unique(c(x_labels_linear,x_labels_nonlin_rcs5,x_labels_nonlin_rcs4,x_labels_nonlin_rcs3,x_labels_fct,x_labels_tag))
# # --- engineer ---
# trim_by_label = "Post-menstrual Age"
# trim_vec = c(-Inf, Inf) # trim relative time [from, to)
# time_unit = 7 # the increment scale of relative time
# pctcut_num_labels = c()
# pctcut_num_vec = c(0.1, 99.9)
# pctcut_num_coerce=TRUE
# filter_tag_labels=c()
# impute_per_cluster=FALSE
# winsorizing=FALSE
# aggregate_per=c("row", "cluster_trim_by_unit", "cluster")[1]
# imputation=c("None","Mean", "Median", "Zero")[1]
# imputeby_median = c()
# imputeby_zero = c()
# imputeby_mean = c()
# impute_per_cluster=FALSE
# standardize_df = NULL # data.frame(varname=c(), center=c(), scale=c())
# # --- local ---
# lasso_by = c("none", "group", "cluster")[3]
# lambda=c("auto","min","1se")[1]
# lambda_value = NULL
# lambda_seq=NULL
# tune_by = c("AIC", "BIC", "cvAUC", "cvLogLoss")[1]
# trim_ctrl = TRUE
# standardize = TRUE # set to TRUE to standardize data by internal engineered dataset if standardize_df is also null
# test_data= data_ml[c(20000:40000),]
# y_map_func=c("fold_risk", "probability", "log_odds")[1]
# y_map_max=3
# return_performance = TRUE
# return_effect_plots = TRUE
# cluster_col <- "m_id"
# fold_idx_df_ex=data.frame(cluster_col=unique(data_ml[,cluster_col]), fold= ceiling(c(1:n_distinct(data_ml[,cluster_col]))/ round(n_distinct(data_ml[,cluster_col])/10) ) )
# fold_idx_df_ex$fold[which(fold_idx_df_ex$fold>10)] <- 10

