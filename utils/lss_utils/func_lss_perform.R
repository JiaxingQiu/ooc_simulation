# validation / testing data are always standardized or engineered in the same way as the modeling data

lss_perform <- function(
  mdl_obj, # a grouped lasso regression model object
  df, # a dataset to test out performance on
  y_map_func = "fold_risk", # response type
  y_map_max = 7, # response upper cutoff
  rel_time_col=NULL,
  return_effect_plots = TRUE,
  lasso_by = c("none", "group", "cluster")[2],
  cv_scores_all_final=NULL
){
  
  library(stringr)
  library(tidytext)
  library(ggplot2)
  # initiate the return objects
  fitted_eff_plot <- NULL
  scores_plot <- NULL
  scores_all_final <- NULL
  tte_plot <- NULL
  df_hat <- NULL

  # define logit reversing function
  logit2prob <- function(logit){
    odds <- exp(logit)
    prob <- odds / (1 + odds)
    return(prob)
  }
  
  # ------------- y-hat generation ------------
  # global mean of outcome
  base_mean_mdl <- mean(mdl_obj$y, na.rm=TRUE)
  # define the right mapping function for y-axis
  if(y_map_func == "probability"){
    ymap <- function(y_logodds, base_mean=NULL){
      y_prob <- ifelse(logit2prob(y_logodds) > y_map_max, y_map_max,logit2prob(y_logodds))
      return(y_prob)
    }
  }else if(y_map_func == "fold_risk"){
    ymap <- function(y_logodds, base_mean=base_mean_mdl){
      print(paste0("--- baseline responce mean for fold increase /relative risk --- ", base_mean))
      y_fold_risk <- ifelse(logit2prob(y_logodds)/base_mean > y_map_max, y_map_max, logit2prob(y_logodds)/base_mean)
      return(y_fold_risk)
    }
  }else{
    ymap <- function(y_logodds, base_mean=NULL){
      y_logodds <- ifelse(y_logodds > y_map_max, y_map_max, y_logodds)
      return(y_logodds)
    }
  }
  
  # find column name for response variable and predictor variables
  y_col <- colnames(mdl_obj$y)
  x_cols <- colnames(mdl_obj$x)
  Xobj <- reconstructX(lasso_optimal_mdl_obj = mdl_obj,df=df)
  x_cols_raw <- Xobj$x_cols_raw
  new_x <- Xobj$new_x
  
  if(lasso_by%in%c("none", "group")){
    # prediction on dataframe
    df$y_pred  <- ymap( predict(mdl_obj, newx=new_x , type = "link") ) # NA are kept
    df$y_logodds  <- predict(mdl_obj, newx=new_x , type = "link")
    df$y_true <- as.factor(df[,y_col])
    df_hat <- df
  }
  if(lasso_by=="cluster"){
    new_df <- Xobj$new_df
    new_df$cluster <- as.factor(df[,mdl_obj$c_info])
    # prediction on dataframe
    df$y_logodds <- NA
    df[complete.cases(new_x),"y_logodds"] <- log(predict(mdl_obj, newdata=new_df)/(1-predict(mdl_obj, newdata=new_df)))
    df$y_pred  <- ymap( df$y_logodds )
    df$y_true <- as.factor(df[,y_col])
    df_hat <- df
  }
  
  
  # ------------- tradeoff curves ------------
  df_hat$y_true01 <- ifelse(as.numeric(as.character(df[,y_col]))==1, 1, 0)
  curves_obj <- mdl_tradeoffs(y_true=df_hat$y_true01, y_hat=df_hat$y_pred)
  tradeoff_plot <- ggpubr::ggarrange(plotlist=curves_obj, nrow=1)
  
  # ----------- calibration plot ---------------
  df_hat$y_cali_groups <- cut(est_pctl( df_hat$y_pred ), 10)
  df_hat$y_true01 <- ifelse(as.numeric(as.character(df[,y_col]))==1, 1, 0)
  base_mean_obs <- mean( df_hat$y_true01, na.rm=TRUE )
  df_cali <- df_hat %>%
    group_by(y_cali_groups) %>%
    summarise(y_cali_observed = ymap(log(mean(y_true01, na.rm=TRUE)/(1-mean(y_true01, na.rm=TRUE))), base_mean=base_mean_obs),
              y_cali_predicted = mean(y_pred, na.rm=TRUE)) %>%
    as.data.frame()
  cali_plot <- ggplot(df_cali, aes(x=y_cali_predicted, y=y_cali_observed)) +
    geom_point() +
    geom_line(color="grey") +
    geom_abline(intercept=0, slope = 1, size=0.3, linetype="dotted") +
    xlab(paste0("Predicted ", y_map_func))+
    ylab(paste0("Observed ", y_map_func))


  # --------- fitted marginal effects -------
  # reformat plot dataframes by group
  tryCatch({
    stopifnot(return_effect_plots)
    fit_eff_plot_list <- list()
    i=0
    df_plot_num <- data.frame()
    df_plot_fct <- data.frame()
    for (x_col in x_cols_raw ){
      df_plot_pred <- df[,c(x_col, "y_pred", "y_true")]
      colnames(df_plot_pred) <- c("predictor_value", "y_pred", "y_true")
      df_plot_pred$predictor_name <- x_col
      if ( n_distinct(df_plot_pred$predictor_value) < 5 | !is.numeric(df_plot_pred$predictor_value)){
        df_plot_pred$predictor_value <- as.character(df_plot_pred$predictor_value)
        df_plot_fct <- bind_rows(df_plot_fct, df_plot_pred)
      }else{
        df_plot_num <- bind_rows(df_plot_num, df_plot_pred)
      }
    }
    if(nrow(df_plot_num)>0){
      # tuning smoothing method and formula based on data size
      smooth_method <- NULL
      smooth_formula <- NULL
      if(nrow(df_plot_num)<5e+03) {
        smooth_method <- "loess"
      }
      if(nrow(df_plot_num)>=5e+03) {
        smooth_method <- "glm"
        smooth_formula <- "y ~ poly(x, 7)"
      }
      i = i+1
      # fitted version of marginal effect plots
      fit_eff_plot_list[[i]] <- ggplot(df_plot_num, aes(x=predictor_value, y=y_pred, color=y_true)) +
        geom_smooth(method = smooth_method, formula = smooth_formula, span=0.3)+
        facet_wrap(~predictor_name, ncol=3, scales = "free_x") +
        ylab(y_map_func) +
        xlab(NULL) +
        scale_color_discrete(name = y_col)+
        scale_color_manual(values=c( "blue", "orange"))+
        theme(legend.position="top") 
    }
    if(nrow(df_plot_fct)>0){
      i = i+1
      df_plot_fct$predictor_value <- paste0(df_plot_fct$predictor_name,"=",df_plot_fct$predictor_value)
      # fitted version of marginal effect plots
      fit_eff_plot_list[[i]] <- ggplot(df_plot_fct, aes(x=predictor_value, y=y_pred, group=y_true, color=y_true)) +
        geom_boxplot()+
        facet_wrap(~predictor_value, ncol=3, scales = "free_x") + # note !!! use predictor_value (level) here for factors
        ylab(y_map_func) +
        xlab(NULL) +
        scale_color_discrete(name = y_col)+
        scale_color_manual(values=c( "blue", "orange"))+
        theme(legend.position="top",
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank()) 
    }

    fit_eff_plot_list <- fit_eff_plot_list[!sapply(fit_eff_plot_list,is.null)]
    fitted_eff_plot <- ggpubr::ggarrange(plotlist = fit_eff_plot_list, ncol=1, common.legend = TRUE, legend = "top")

  },error=function(e){print(e)})

  # --------- relative time plot ----------
  tryCatch({
    stopifnot(return_effect_plots)
    stopifnot(length(rel_time_col)>0)
    stopifnot(rel_time_col %in% colnames(df))
    
    # zoom in relative time column by the range of event group
    rel_time_min <- min(df[which(df[,y_col]==1),rel_time_col],na.rm=TRUE)
    rel_time_max <- max(df[which(df[,y_col]==1),rel_time_col],na.rm=TRUE)
    df[which(df[,rel_time_col]<rel_time_min | df[,rel_time_col]>=rel_time_max),rel_time_col] <- NA
    df_plot_pred <- data.frame(relative_time = df[,rel_time_col],
                               y_pred = as.numeric(df$y_pred),
                               y_true = as.factor(as.character( df$y_true)))
    df_agg <- df_plot_pred %>% group_by(relative_time, y_true) %>% summarise( y_pred_mean = mean(y_pred, na.rm=TRUE),
                                                                               y_pred_median = median(y_pred, na.rm=TRUE) ) %>% as.data.frame()
    df_agg_plot_all <- data.frame()
    df_agg_plot <- df_agg[,c("relative_time", "y_true", "y_pred_mean")]
    colnames(df_agg_plot) <- c("relative_time", "y_true", "y_pred_stat")
    df_agg_plot$y_pred_stat_type <- "mean"
    df_agg_plot_all <- bind_rows(df_agg_plot_all, df_agg_plot)
    df_agg_plot <- df_agg[,c("relative_time", "y_true", "y_pred_median")]
    colnames(df_agg_plot) <- c("relative_time", "y_true", "y_pred_stat")
    df_agg_plot$y_pred_stat_type <- "median"
    df_agg_plot_all <- bind_rows(df_agg_plot_all, df_agg_plot)
    
    # tuning smoothing method and formula based on data size
    smooth_method <- NULL
    smooth_formula <- NULL
    if(nrow(df_plot_pred)<5e+03) {
      smooth_method <- "loess"
    }
    if(nrow(df_plot_pred)>=5e+03) {
      smooth_method <- "glm"
      smooth_formula <- "y ~ poly(x, 6)"
    }
    tte_plot <- ggplot(data=df_agg_plot_all, aes(x=relative_time, color=y_true)) +
      #geom_smooth(data = df_plot_pred, aes(y=y_pred), method = smooth_method, formula = smooth_formula, span=0.3)+
      geom_smooth(aes(y=y_pred_stat, linetype=y_pred_stat_type), method = "loess", span=0.3)+
      ylab(y_map_func) +
      xlab(rel_time_col) +
      scale_color_discrete(name = y_col)+
      scale_color_manual(values=c( "blue", "orange"))+
      theme(legend.position="top") 
    
  },error=function(e){
    print("Skip relative time (tte) plot")
    print(e)
  })
  
  # ----------- feature permutation importance -------------
  tryCatch({
    if(!is.null(cv_scores_all_final)){
      scores_all_final <- cv_scores_all_final
    }else{
      y_true_org <- as.numeric(as.character( df[,y_col] ) )
      y_prob_org <- as.numeric(logit2prob(predict(mdl_obj, newx = new_x, type="link")))
      y_true <- y_true_org[which(!is.na(y_true_org)&!is.na(y_prob_org))]
      y_prob <- y_prob_org[which(!is.na(y_true_org)&!is.na(y_prob_org))]
      test_obj_raw <- mdl_test(y_true = y_true,
                               y_prob = y_prob,
                               threshold=base_mean_mdl)
      scores_raw <- test_obj_raw$res_df[1,c("logloss", "AUROC", "AUPRC", "accuracy", "f1score")]
      scores_all <- data.frame()
      for (x_col in setdiff(x_cols_raw,rel_time_col)){
        df_shuffled <- df
        # shuffle the predictor
        df_shuffled[,x_col] <- sample(df[,x_col], size=length(df[,x_col]))
        Xobjtest <- reconstructX(mdl_obj, df_shuffled)
        new_x_test = Xobjtest$new_x
        y_true_org <- as.numeric(as.character( df_shuffled[,y_col] ) )
        y_prob_org <- as.numeric(logit2prob(predict(mdl_obj, newx = new_x_test, type="link")))
        y_true <- y_true_org[which(!is.na(y_true_org)&!is.na(y_prob_org))]
        y_prob <- y_prob_org[which(!is.na(y_true_org)&!is.na(y_prob_org))]
        
        test_obj <-  mdl_test(y_true = y_true,
                              y_prob = y_prob,
                              threshold=base_mean_mdl)
        scores <- test_obj$res_df[1,c("logloss", "AUROC", "AUPRC", "accuracy", "f1score")]
        scores$varname <- x_col
        scores_all <- bind_rows(scores_all, scores)
      }
      scores_all <- scores_all[,union("varname", colnames(scores_all))]
      scores_raw$varname <- "none"
      scores_all_final <- bind_rows(scores_all, scores_raw)
      colnames(scores_all_final)[which(colnames(scores_all_final) == "varname")] <- "removed_variable"
    }
    
    # make permutation importance plot
    df_plot_scores_all <- data.frame()
    for (score_col in setdiff(colnames(scores_all_final), "removed_variable") ){
      df_plot_scores <- scores_all_final[,c("removed_variable", score_col)]
      colnames(df_plot_scores) <- c("removed_variable", "score_value")
      df_plot_scores$score_by <- score_col
      df_plot_scores_all <- bind_rows(df_plot_scores_all, df_plot_scores)
    }
    
    scores_plot <- ggplot(data=df_plot_scores_all, aes(x=score_value, y=tidytext::reorder_within(removed_variable, score_value, score_by, sep="_________") )) +
      geom_point()+
      geom_vline(data=df_plot_scores_all[which(df_plot_scores_all$removed_variable=="none"),], aes(xintercept=score_value))+
      tidytext::scale_y_reordered(sep="_________") +
      facet_wrap(~score_by, scales = "free", ncol = 2)+
      theme(axis.title.x=element_blank(),
            axis.title.y=element_blank())

  },error=function(e){print(e)})

  print("test")
  
  
  return(list(
    df_hat = df_hat,
    fitted_eff_plot = fitted_eff_plot,
    scores_plot = scores_plot,
    scores_all_final = scores_all_final,
    cali_plot = cali_plot,
    tte_plot=tte_plot,
    tradeoff_plot = tradeoff_plot
  ))

}




reconstructX <- function(lasso_optimal_mdl_obj, df){
  
  df_org <- df
  #x_cols_org <- names(attr(lasso_optimal_mdl_obj$x,"scaled:center"))
  #df[,x_cols_org] <- scale(df[,x_cols_org], center = attr(lasso_optimal_mdl_obj$x,"scaled:center"), scale = attr(lasso_optimal_mdl_obj$x,"scaled:scale") )
  # overwrite tag columns
  x_cols_tag <- lasso_optimal_mdl_obj$group_info$x_colname[which(endsWith( lasso_optimal_mdl_obj$group_info$x_group, "_tag01") )]
  df[,x_cols_tag] <- df_org[,x_cols_tag]
  
  # find column name for response variable and predictor variables
  y_col <- colnames(lasso_optimal_mdl_obj$y)
  x_cols <- colnames(lasso_optimal_mdl_obj$x)
  x_cols_raw <- c()
  group_info_df <- lasso_optimal_mdl_obj$group_info
  
  # reconstruct rcs variables
  for (xg in unique(group_info_df$x_group) ){
    gname <- unlist(strsplit(xg, "_"))[length(unlist(strsplit(xg, "_")))]
    col <- gsub(paste0("_",gname), "",xg)
    x_cols_raw <- union(x_cols_raw, col)
    if(gname == "rcs3"){
      df[,col] <- as.numeric(rcs(df[,col],3)[,1])
      df[,paste0(col,"'")] <- as.numeric(rcs(df[,col],3)[,2])
    }
    if(gname == "rcs4"){
      df[,col] <- as.numeric(rcs(df[,col],4)[,1])
      df[,paste0(col,"'")] <- as.numeric(rcs(df[,col],4)[,2])
      df[,paste0(col,"''")] <- as.numeric(rcs(df[,col],4)[,3])
    }
    if(gname == "rcs5"){
      df[,col] <- as.numeric(rcs(df[,col],5)[,1])
      df[,paste0(col,"'")] <- as.numeric(rcs(df[,col],5)[,2])
      df[,paste0(col,"''")] <- as.numeric(rcs(df[,col],5)[,3])
      df[,paste0(col,"'''")] <- as.numeric(rcs(df[,col],5)[,4])
    }
    if(gname == "level"){
      levels <- unique(as.character(df[,col]))
      for(l in levels){
        l_fix <- gsub("[^[:alnum:]]","_",l)
        df[,paste0(col,"___",l_fix)] <- ifelse(df[,col]==l,1,0)
      }
    }
  }
  stopifnot(all(x_cols %in% colnames(df)))
  #stopifnot(all(x_cols_org %in% x_cols_raw))
  new_x <- data.matrix(df[,x_cols])
  return(list(new_x = new_x,
              new_df = df,
              x_cols_raw = x_cols_raw))
}



# ######################## not run ############################
# mdl_obj = x_select_mdls_grouped$lasso_optimal # a grouped lasso regression model object
# df = data_in # a dataset to test out performance on
# y_map_func = "fold_risk" # response type
# y_map_max = 7 # response upper cutoff
# rel_time_col=NULL
# x_cols=c() # list of predictors to see the performance




