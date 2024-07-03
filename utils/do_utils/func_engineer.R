engineer <- function(
  data, # dataframe object to engineer
  num_cols = c(), # vector of numeric columns
  fct_cols = c(), # vector of factor columns
  cluster_col, # cluster column
  trim_by_col = NULL,
  #--- trim data ---
  trim_min = -Inf, # [from,
  trim_max = Inf, # to)
  trim_step_size = 1,
  trim_keepna = FALSE, # whether or not to keep NA in relative time variable
  trim_first = TRUE, # TRUE = trim the data at the beginning before other engineering; FALSE = trim the data as the last step
  pctcut_num_cols=c(),
  pctcut_num_vec=c(0.1, 99.9),
  pctcut_num_coerce = TRUE,
  filter_tag_cols=c(),
  #--- decorate data ---
  winsorizing=FALSE,
  aggregate_per=c("row", "cluster_trim_by_unit", "cluster")[1],
  aggregate_conditioned_on_cols = c(),
  imputation=c("None", "Mean", "Median", "Zero")[1], # global imputation strategy for the rest of vars not specified in imputeby_...
  imputeby_median = c(),
  imputeby_zero = c(),
  imputeby_mean = c(), 
  impute_per_cluster=FALSE,
  standardize_df = NULL, # can be c(NULL, a dataframe)
  scaler = c("none", "normal", "robust", "rankit", "percent")[1],
  sample_per_cluster = NULL
){
  
  #### clean up inputs ####
  trim_by_col <- intersect(colnames(data), trim_by_col)
  if(length(trim_by_col)>0){
    trim_step_size <- max(1,as.numeric(trim_step_size),na.rm=TRUE)
    trim_min <- trim_min*trim_step_size
    trim_max <- trim_max*trim_step_size
    if (trim_max<=trim_min){
      trim_min <- -Inf
      trim_max <- Inf
    }
  }
  cluster_col <- intersect(colnames(data), cluster_col)
  num_cols <- intersect(colnames(data), num_cols)
  num_cols <- union(num_cols, trim_by_col)
  num_cols <- union(num_cols, pctcut_num_cols)
  fct_cols <- intersect(colnames(data), fct_cols)
  fct_cols <- union(fct_cols, cluster_col)
  fct_cols <- union(fct_cols, filter_tag_cols)
  imputeby_median_cols <- intersect(num_cols,intersect(colnames(data), imputeby_median))
  imputeby_zero_cols <- intersect(num_cols,intersect(colnames(data), imputeby_zero))
  imputeby_mean_cols <- intersect(num_cols,intersect(colnames(data), imputeby_mean))
  aggregate_conditioned_on_cols <- intersect(fct_cols, aggregate_conditioned_on_cols) # aggregation should be done by combination of groups in these columns
  
  # initiate dataframe to return "data_engineered"
  # in the returned dataframe, only trim_by_col, num_cols, fct_cols, cluster_col will be in the columns
  data <- data[,unique(c(trim_by_col, num_cols, fct_cols, cluster_col)) ]
  data_engineered <- data
  
  print("------ Engineering ------")
  #### coerce data type
  if(length(num_cols)>0){
    for(col in num_cols){
      data[, col] <- as.numeric(data[, col]) # make it numeric
      data[which(is.infinite(data[, col])), col] <- NA
    }
  }
  if(length(trim_by_col)>0){
    data[,trim_by_col] <- as.numeric(data[,trim_by_col])
    data[which(is.infinite(data[, trim_by_col])), trim_by_col] <- NA
  }
  if(length(cluster_col)>0){
    data[,cluster_col] <- as.character(data[,cluster_col])
  }
  if(length(fct_cols)>0){
    # (tag) - max
    tag_cols <- c()
    cat_cols <- c()
    for (col in fct_cols){
      if(all(unique(as.character(data[,col]) )%in%c("0","1",NA) )) {
        tag_cols <- c(tag_cols, col)
        data[, col] <- as.numeric(data[, col]) # make it numeric 0 1
      }else {
        cat_cols <- c(cat_cols, col)
        data[, col] <- as.character(data[, col]) # make it character
      }
    }
  }
  
  #### trim data by trim_by_col(relative time indicator) ####
  if (length(trim_by_col)>0){
    if(trim_first){
      print("trimming data at the beginning")
      if (trim_keepna){
        data <- data[which((data[,trim_by_col]>=trim_min & data[,trim_by_col]<trim_max)|is.na(data[,trim_by_col])),]
      }else{
        data <- data[which(data[,trim_by_col]>=trim_min & data[,trim_by_col]<trim_max),]
      }
    }
  }
  #### cutoff numeric variables by percentiles ####
  print(paste0("--- Cutting following variables by percentiles and keep rows within range --- ",paste0(pctcut_num_cols,collapse = ", ")))
  for(pctcut_num_col in pctcut_num_cols){
    tryCatch({
      quantiles <- quantile( data[,pctcut_num_col], c(as.numeric(pctcut_num_vec[1])/100, as.numeric(pctcut_num_vec[2])/100 ), na.rm =TRUE)
      if(pctcut_num_coerce){
        # if coerce extrema
        data[,pctcut_num_col] <- ifelse(data[,pctcut_num_col] < quantiles[1], quantiles[1], data[,pctcut_num_col])
        data[,pctcut_num_col] <- ifelse(data[,pctcut_num_col] > quantiles[2], quantiles[2], data[,pctcut_num_col])
      }else{
        # otherwise remove extrema
        data[,pctcut_num_col] <- ifelse(data[,pctcut_num_col] < quantiles[1], NA, data[,pctcut_num_col])
        data[,pctcut_num_col] <- ifelse(data[,pctcut_num_col] > quantiles[2], NA, data[,pctcut_num_col])
      }
    },error=function(e){
      print(paste0("--- Skip percentile cutoff for num variable ",pctcut_num_col," ---"))
      print(e)
    })
  }
  data <- data[which(complete.cases(data[,pctcut_num_cols])),]# filter data with complete numeric cutoffs
  
  
  #### winsorize numeric columns except trim by col ####
  if(winsorizing){
    if(length(setdiff(num_cols,trim_by_col))>0){
      print("--- Winsorize all the numeric variables ---")
      data[,setdiff(num_cols,trim_by_col)] <- winsorize(data[,setdiff(num_cols,trim_by_col)])
    }
  }
  
  #### subset data by binary filters ####
  if(length(filter_tag_cols)>0){
    for(col in filter_tag_cols){
      tryCatch({
        print(paste0("--- filter data where ",col,"==1 ---"))
        data <- data[which(as.numeric( data[,col] )==1),]
      },error=function(e){
        print(paste0("--- Skip one-hot filter for tag variable", col, " ---"))
        print(e)
      })
    }
  }
  #### aggregation ####
  if(aggregate_per %in% c("cluster_trim_by_unit", "cluster") ){ 
    print(paste0("--- Aggregate data by ", aggregate_per," ---"))
    df_key <- NULL
    df_tag <- NULL
    df_cat <- NULL
    df_num <- NULL
    if (aggregate_per=="cluster_trim_by_unit"){
      # rounding trim by column as multiple value of trim by step size / time unit
      data[,trim_by_col] <- floor(data[,trim_by_col]/trim_step_size)*trim_step_size
      data$key_col <- paste0(data[,cluster_col], "___", floor(data[,trim_by_col]/trim_step_size))
    }else if (aggregate_per=="cluster"){
      data$key_col <- data[,cluster_col]
    }
    # add conditional groups to key col
    if(length(aggregate_conditioned_on_cols)>0){
      data$key_col <- apply(data[, c("key_col",aggregate_conditioned_on_cols)], 1, paste, collapse = "-" )
    }
    
    # key 
    df_key <- data.frame(key_col = unique(data$key_col), stringsAsFactors = FALSE)
    # num - mean
    if(length(num_cols)>0){
      df_num <- data %>% group_by(key_col) %>% 
        summarise_at(vars(all_of(num_cols)), ~mean(.,na.rm = TRUE)) %>% 
        as.data.frame()
      colnames(df_num) <- c("key_col", num_cols)
    }
    # fct 
    if(length(fct_cols)>0){
      df_tag <- data %>% group_by(key_col) %>% 
        summarise_at(vars(all_of(tag_cols)), ~max(.)) %>%
        mutate_at(vars(all_of(tag_cols)), function(x) ifelse(is.infinite(x), NA, x)) %>%
        as.data.frame()
      colnames(df_tag) <- c("key_col", tag_cols)
      df_cat <- data %>% group_by(key_col) %>% 
        summarise_at(vars(all_of(cat_cols)), ~unique(.)[which(!is.na(unique(.)))][1] ) %>% # use the first not na value
        as.data.frame()
      colnames(df_cat) <- c("key_col", cat_cols)
    }
    data <- df_key
    if (!is.null(df_tag)){
      data <- merge(data, df_tag)
    }
    if (!is.null(df_cat)){
      data <- merge(data, df_cat)
    }
    if (!is.null(df_num)){
      data <- merge(data, df_num)
    }
  }
  data_engineered <- data[, intersect(colnames(data),unique(c("key_col",cluster_col, trim_by_col, num_cols, fct_cols)))]
  
  
  #### imputation ####
  for(num_col in num_cols){ data[which(!is.finite(data[,num_col])), num_col] <- NA } # clean inf to NA
  print("--- Imputation ---")
  data <- data_engineered
  tryCatch({
    if(!impute_per_cluster){
      print("---- by cohort ----")
      # special imputation
      data <- data %>% mutate_at(vars(all_of(imputeby_median_cols)), ~tidyr::replace_na(., median(.,na.rm = TRUE))) %>% as.data.frame()
      data <- data %>% mutate_at(vars(all_of(imputeby_mean_cols)), ~tidyr::replace_na(., mean(.,na.rm = TRUE))) %>% as.data.frame()
      data <- data %>% mutate_at(vars(all_of(imputeby_zero_cols)), ~tidyr::replace_na(., 0)) %>% as.data.frame()
      # global imputation
      if(imputation=="Mean") data <- data %>% mutate_at(vars(all_of(num_cols)), ~tidyr::replace_na(., mean(.,na.rm = TRUE))) %>% as.data.frame()
      if(imputation=="Median") data <- data %>% mutate_at(vars(all_of(num_cols)), ~tidyr::replace_na(., median(.,na.rm = TRUE))) %>% as.data.frame()
      if(imputation=="Zero") data <- data %>% mutate_at(vars(all_of(num_cols)), ~tidyr::replace_na(., 0)) %>% as.data.frame()
    }else{
      print("---- by cluster ----")
      
      # cluster-wise imputation
      if(imputation=="Mean") data %>% group_by(vars(cluster_col)) %>% mutate_at(vars(all_of(num_cols)),  ~tidyr::replace_na(., mean(.,na.rm = TRUE)) ) %>% as.data.frame()
      if(imputation=="Median") data %>% group_by(vars(cluster_col)) %>% mutate_at(vars(all_of(num_cols)), ~ifelse(is.na(.), median(., na.rm = TRUE),.) ) %>% as.data.frame()
      if(imputation=="Zero")  data %>% group_by(vars(cluster_col)) %>% mutate_at(vars(all_of(num_cols)),  ~tidyr::replace_na(., 0)) %>% as.data.frame()
      
      # then global imputation
      if(imputation=="Mean") data <- data %>% mutate_at(vars(all_of(num_cols)), ~tidyr::replace_na(., mean(.,na.rm = TRUE))) %>% as.data.frame()
      if(imputation=="Median") data <- data %>% mutate_at(vars(all_of(num_cols)), ~tidyr::replace_na(., median(.,na.rm = TRUE))) %>% as.data.frame()
      if(imputation=="Zero") data <- data %>% mutate_at(vars(all_of(num_cols)), ~tidyr::replace_na(., 0)) %>% as.data.frame()
    }
  },error=function(e){
    print(e)
    print("Skip imputation due to error")
  })
  
  data_engineered <- data
  
  #### standardization by df ####
  if(!is.null(standardize_df)){
    tryCatch({
      stopifnot(all(c("varname", "center", "scale")%in%colnames(standardize_df)))
      for(col in unique(standardize_df$varname) ){
        tryCatch({
          data_engineered[,col] <- as.numeric(scale(data_engineered[,col], 
                                         center = standardize_df$center[which(standardize_df$varname==col)],
                                         scale = standardize_df$scale[which(standardize_df$varname==col)]))
          data_engineered[which(!is.finite(data_engineered[,col] )),col]<-NA
        },error = function(e){
          print(e)
          print(paste0("skip standardize -- ", col))
        })
      }
    },error=function(e){
      print(e)
      print("skip standardization -- error in input standardize dataframe ")
    })
  }
  
  ### standardization ###
  if(is.null(standardize_df)&(!scaler=="none")){
    for(c in num_cols){
      tryCatch({
        if(scaler=="normal"){
          data_engineered[,c] <- normal_scale(data_engineered[,c])
        }
        if(scaler=="robust"){
          data_engineered[,c] <- robust_scale(data_engineered[,c])
        }
        if(scaler=="rankit"){
          data_engineered[,c] <- rankit(data_engineered[,c])
        }
        if(scaler=="percent"){
          data_engineered[,c] <- est_pctl(data_engineered[,c])
        }
      },error = function(e){
        print(e)
        print(paste0(scaler," scaler failed on", c))
      })
    }
    
  }
  
  if(!is.null(sample_per_cluster)){
    tryCatch({
      data_engineered$cluster <- data_engineered[,cluster_col]
      data_engineered_sampled <- data_engineered %>% group_by(cluster) %>% 
        dplyr::slice_sample(n=sample_per_cluster, replace=TRUE) %>%
        as.data.frame()
      data_engineered <- data_engineered_sampled[,setdiff(colnames(data_engineered_sampled),c("cluster"))]
    },error=function(e){
      print(e)
    })
    
  }
  
  #### trim data by trim_by_col(relative time indicator) if trim_first is FALSE ####
  data <- data_engineered
  if (length(trim_by_col)>0){
    if(!trim_first){
      print("trimming data at the end")
      if (trim_keepna){
        data <- data[which((data[,trim_by_col]>=trim_min & data[,trim_by_col]<trim_max)|is.na(data[,trim_by_col])),]
      }else{
        data <- data[which(data[,trim_by_col]>=trim_min & data[,trim_by_col]<trim_max),]
      }
    }
  }
  data_engineered <- data
  
  
  # ---- return final engineered dataset ----
  return(data_engineered)
}








# ############################### not run ######################################
# 
# data = data_ml # dataframe object to engineer
# num_cols = c("baby_weight_lin") # vector of numeric columns
# fct_cols = c("delivery_mode_factor", "posair_ynunk_drvd___tag_factor_Yes") # vector of factor columns
# cluster_col = "subjectnbr" # cluster column
# trim_by_col = "pma_days"
# #--- trim data ---
# trim_min = 35 # [from,
# trim_max = 37 # to)
# trim_step_size = 7
# trim_keepna = FALSE # whether or not to keep NA in relative time variable
# pctcut_num_cols=c()
# pctcut_num_vec=c(0.1, 99.9)
# pctcut_num_coerce = TRUE
# filter_tag_cols=c()
# #--- decorate data ---
# imputation=c("None","Mean", "Median", "Zero")[1]
# impute_per_cluster=FALSE
# winsorizing=FALSE
# aggregate_per=c("row", "cluster_trim_by_unit", "cluster")[2]


