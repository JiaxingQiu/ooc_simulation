do_init <- function(df, 
                    dict_df, 
                    y_col, 
                    x_cols, 
                    cluster_col, 
                    rcs5_low="70%",
                    rcs4_low="50%",
                    linear_cols=NULL){
  # ---- Description ----
  # dictionary oriented (do) initiate degree of freedom to variables by adjusted spearman rho
  
  # ---- Arguments ----
  # df: a dataframe object with essential dictionary information as attributes on each column
  # dict: dictionary for corresponding data frame input
  # x_cols: pre-selected predictor column names (x)
  # r2: ordinary or adjusted R^2 cutoff for redundancy
  
  # ---- Value ----
  # df_final: a dataframe object with updated redundency attributes on each column
  # dict_final: updated dictionary for correspondingdf_final
  
  # initiate rho object
  rho <- NULL
  
  x_cols_num <- intersect(x_cols,rownames(dict_df[which(dict_df$type=="num"),]))
  x_cols_fct <- intersect(x_cols,rownames(dict_df[which(dict_df$type=="fct"),]))
  if (length(x_cols_num)>0){
    
    # --- calculate quadratic spearman rank for numeric predictors ---
    fml <- formula(paste(y_col," ~ ", paste(x_cols_num,collapse = "+")))
    rho <- spearman2(fml, data=df, p=2)
    rho_df <- data.frame(adj_rho=rho[,"Adjusted rho2"])
    rho_df$cols <- rownames(rho)
    rho_df <- rho_df[order(rho_df$adj_rho, decreasing = TRUE),]
    scaler <- data.frame(q=quantile(rho_df$adj_rho, probs = seq(0, 1, 0.01) , na.rm = TRUE))
    
    for (col in rho_df$cols){
      attr(df[,col], "spearman_rho") <- rho_df$adj_rho[which(rho_df$cols==col)]
    }
    df <- assign.dict(df,dict_df)
    
    # --- assign dof to each predictor ---
    for (col in rho_df[which(rho_df$adj_rho>=scaler[rcs5_low,"q"]),'cols']){
      attr(df[,col],"mdl_term_init") <- "rcs5"
    }
    for (col in rho_df[which(rho_df$adj_rho>=scaler[rcs4_low,"q"] & rho_df$adj_rho<scaler[rcs5_low,"q"]),'cols']){
      attr(df[,col],"mdl_term_init") <- "rcs4"
    }
    for (col in rho_df[which(rho_df$adj_rho<scaler[rcs4_low,"q"]),'cols']){
      attr(df[,col],"mdl_term_init") <- "rcs3"
    }
  }
  if (length(x_cols_fct)>0){
    for (col in x_cols_fct){
      if(n_distinct(as.character( df[,col] ) )>=2) attr(df[,col],"mdl_term_init") <- "factor"
    }
  }
  attr(df[,cluster_col],"mdl_term_init") <- "cluster"
  attr(df[,y_col],"mdl_term_init") <- "y"
  
  
  dict_mdl <- get.dict(df)
  selected_cols <- rownames(dict_mdl[which(dict_mdl$mdl_term_init!=""),])
  df <- df[,selected_cols]
  
  
  # ---- change customized inputs and redundency filtered linear vars ----
  linear_only_cols <- union(linear_cols, dict_df$varname[which(dict_df$redun=="linear")])
  if (length(linear_only_cols)>0){
    new_lin_cols <- intersect(linear_only_cols, colnames(df))
    for (col in new_lin_cols){
      attr(df[,col],"mdl_term_init") <- "linear"
    }
  }
  
  # ---- add missingness fraction as a attribute to columns ----
  for (var in colnames(df)){
    attr(df[,var],"na_frac") <- sum(is.na(df[,var]))/length(df[,var])
  }
  
  df_final <- df
  dict_final <- get.dict(df_final)
  df_final <- assign.dict(df_final, dict_final)
  
  return(list("df_final" = df_final,
              "dict_final" = dict_final,
              "spearman2_rho" = rho))
  
}

