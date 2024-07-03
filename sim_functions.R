
# generate AR(1) series
generate_ar1 <- function(n, phi, sigma=1) {
  # phi = how strong is the autocorrelation
  e <- rnorm(n, mean = 0, sd = sigma) # Generate white noise
  Y <- rep(0, n) # Placeholder for the AR(1) process
  Y[1] <- e[1] # Initial value
  for (t in 2:n) { Y[t] <- phi * Y[t - 1] + e[t]}
  return(Y)
}


random_removal <- function(c, cluster_na_rate) {
  remaining_indices <- integer(0)
  for (i in unique(c)) {
    cluster_indices <- which(c == i)
    num_remove <- ceiling(length(cluster_indices) * cluster_na_rate[i])
    remove_indices <- sample(cluster_indices, num_remove)
    keep_indices <- setdiff(cluster_indices, remove_indices)
    remaining_indices <- c(remaining_indices, keep_indices)
  }
  return(sort(remaining_indices))
}

# data generating process
generate_data <- function(n_cluster,
                          n_obs_per_cluster,
                          n_ttl_betas, 
                          fix_rdm_ratio,
                          sigma_fix,
                          sigma_rdm_fix_ratio,
                          ar1_phi,
                          na_rate,
                          family=c("binomial","gaussian")[1]) {
  library(dplyr)
  library(MASS)
  # note: intercepts will be allowed in random effects, all random effects will be in fixed effects as well
  
  sigma_x <- 1
  sigma_beta <- sigma_fix
  sigma_b <- sigma_fix * sigma_rdm_fix_ratio
  
  # cluster id
  cid_df <- data.frame(cluster = rep(1:n_cluster, each = n_obs_per_cluster))
  
  # generate fixed effect data
  n_fix_betas <- ceiling(n_ttl_betas*fix_rdm_ratio)
  fix_x_ls <- list()
  for(i in 1:n_fix_betas){
    fix_x_ls[[i]] <- rnorm(n_cluster * n_obs_per_cluster, 0, sigma_x)
  }
  fix_x_df <- as.data.frame(fix_x_ls)
  colnames(fix_x_df) <- paste0("fix",c(1:ncol(fix_x_df)))
  fix_x_df$fix_itc <- 1
  
  # generate random effect data by AR(1)
  rdm_x_df <- NULL
  n_rdm_betas <- n_ttl_betas - n_fix_betas
  if(n_rdm_betas>0){
    rdm_phi_vec_ls <- lapply(as.list(rep(n_rdm_betas, n_cluster)), function(n){runif(n, ar1_phi, ar1_phi+0.2)})  
    rdm_x_df_ls <- lapply(rdm_phi_vec_ls, function(phis){
      phis <- as.list(phis)
      x_ls <- lapply(phis, function(x) generate_ar1(n=n_obs_per_cluster, phi=x)) # AR1(phi)
      # x_ls <- lapply(phis, function(x) generate_ar1(n=n_obs_per_cluster, phi=ar1_phi)) # white noise
      x_df <- data.frame(x_ls)
      colnames(x_df) <- paste0("rdm",c(1:ncol(x_df)))
      return(x_df)   } )
    rdm_x_df <- bind_rows(rdm_x_df_ls)
    rdm_x_df$rdm_itc <- 1
  }
  
  
  # generate fixed effects betas  
  fix_betas_mat <- matrix(rnorm(ncol(fix_x_df), 0, sigma_beta), nrow = 1, ncol=ncol(fix_x_df)) # from multi-gaussian
  
  # generate clustered random effects betas
  if(!is.null(rdm_x_df)){
    # ---- within-clustaer random effects ----
    # # from multi-uniform
    # rdm_betas_mat_center <- runif(ncol(rdm_x_df), -5, 5)
    # rdm_betas_mat <- c()
    # for(i in 1:n_cluster){
    #   rdm_betas_mat_error_i <- rnorm(ncol(rdm_x_df), 0, 0.1)
    #   rdm_betas_mat_i <- t(replicate(n_obs_per_cluster, rdm_betas_mat_center + rdm_betas_mat_error_i))
    #   rdm_betas_mat_i <- matrix(rdm_betas_mat_i, nrow=n_obs_per_cluster, ncol=ncol(rdm_x_df))
    #   rdm_betas_mat <- rbind(rdm_betas_mat, rdm_betas_mat_i)
    # }
    # # from multi gaussian
    rdm_betas_mat <- c()
    rdm_betas_mat_i <- mvrnorm(n_cluster, rep(0,ncol(rdm_x_df)), diag(rep(sigma_b, ncol(rdm_x_df))))
    rdm_betas_mat_i <- matrix(rdm_betas_mat_i, nrow=n_cluster, ncol=ncol(rdm_x_df))
    rdm_betas_mat <- rdm_betas_mat_i[rep(1:nrow(rdm_betas_mat_i), each = n_obs_per_cluster),]
    rdm_betas_mat <- matrix(rdm_betas_mat, nrow=n_cluster*n_obs_per_cluster, ncol=ncol(rdm_x_df))
    
    # ---- fixed effects from random effect variables as well (except intercept) -----
    fix_rdm_x_df <- as.matrix(rdm_x_df[,setdiff(colnames(rdm_x_df),"rdm_itc")], nrow=nrow(rdm_x_df), ncol=length(setdiff(colnames(rdm_x_df),"rdm_itc")))
    colnames(fix_rdm_x_df) <- paste0("fix_",setdiff(colnames(rdm_x_df),"rdm_itc"))
    fix_rdm_betas_mat <- matrix(rnorm(ncol(fix_rdm_x_df), 0, sigma_beta), nrow = 1, ncol=ncol(fix_rdm_x_df)) # from multi-gaussian
    
  }
  
  # return betas and data
  fix_betas_df <- as.data.frame(fix_betas_mat[rep(1:nrow(fix_betas_mat), each = n_cluster * n_obs_per_cluster),])
  colnames(fix_betas_df) <- paste0("beta_",colnames(fix_x_df))
  betas <- fix_betas_df
  data <- fix_x_df
  
  if(!is.null(rdm_x_df)){
    rdm_betas_df <- as.data.frame(rdm_betas_mat)
    colnames(rdm_betas_df) <- paste0("beta_",colnames(rdm_x_df))
    betas <- bind_cols(betas, rdm_betas_df)
    data <- bind_cols(data, rdm_x_df)
    
    fix_rdm_betas_df <-  as.data.frame(fix_rdm_betas_mat[rep(1:nrow(fix_rdm_betas_mat), each = n_cluster * n_obs_per_cluster),])
    colnames(fix_rdm_betas_df) <- paste0("beta_",colnames(fix_rdm_x_df))
    betas <- bind_cols(betas, fix_rdm_betas_df)
    data <- bind_cols(data, fix_rdm_x_df)
  } 
  stopifnot(dim(data)[1] == dim(betas)[1])
  stopifnot(dim(data)[2] == dim(betas)[2])
  
  # returns complete
  c <- as.numeric(cid_df$cluster)
  data <- as.matrix(data)
  betas <- as.matrix(betas)
  l <- as.numeric(rowSums(data * betas))
  if(family=="binomial"){
    p <- 1/(1+exp(-l))
    y <- rbinom(length(p), 1, p)
  }else{
    p <- NULL
    y <- l + rnorm(nrow(data), 0, 2)
  }
  
  # cluster-wise missing at random
  if(na_rate>0 & na_rate<1){
    cluster_na_rate <- rnorm(n_cluster, na_rate, 0.1)
    cluster_na_rate[cluster_na_rate>0.8] <- 0.8
    cluster_na_rate[cluster_na_rate<0] <- 0
    remaining_indices <- random_removal(c, cluster_na_rate)
    data <- data[remaining_indices, ]
    betas <- betas[remaining_indices, ]
    p <- p[remaining_indices]
    y <- y[remaining_indices]
    c <- c[remaining_indices] 
  }else if(na_rate == 1){
    cluster_na_rate <- runif(n_cluster, 0.1, 0.9)
    remaining_indices <- random_removal(c, cluster_na_rate)
    data <- data[remaining_indices, ]
    betas <- betas[remaining_indices, ]
    p <- p[remaining_indices]
    y <- y[remaining_indices]
    c <- c[remaining_indices] 
  }
  
  
  # # y label group by c
  # cluster_means <- tapply(y, c, mean)# Calculate the condition for each cluster
  # cluster_conditions <- ifelse(cluster_means > 0.5, 1, 0)
  # yy <- cluster_conditions[c]# Map the conditions back to each element
  # yy <- as.numeric(yy)
  
  return(list("y"=y,#yy,
              "p"=p,
              "c"=c,
              "data" = data,
              "betas" = betas))
}

fit_eval_glmer <- function(y, c, data, family=c("binomial","gaussian")[1], skiploo = F){
  library(Matrix)
  library(lme4)
  
  df_mdl <- as.data.frame(data)
  df_mdl$y <- y
  df_mdl$c <- as.factor(c)
  
  fix_vars <- grep("^fix\\d+$", colnames(data), value = TRUE)
  rdm_vars <- grep("^rdm\\d+$", colnames(data), value = TRUE)
  # paste0(paste0("(0+",rdm_vars,"|c)"),collapse = "+"), "+ (1|c)"
  #paste0("(",paste0(rdm_vars,collapse = "+"),"+1|c)"
  # fml <- formula(paste0("y~", paste0(c(fix_vars, rdm_vars),collapse = "+"), "+", paste0(paste0("(0+",rdm_vars,"|c)"),collapse = "+"), "+ (1|c)"))
  fml <- formula(paste0("y~", paste0(c(fix_vars, rdm_vars),collapse = "+"), "+", paste0("(",paste0(rdm_vars,collapse = "+"),"+1|c)")))
  
  if(family == "binomial"){
    mdl <- glmer(fml,  data = df_mdl, family = binomial)
  }else if(family == "gaussian"){
    mdl <- lmer(fml, data = df_mdl, REML=F)
  }
  # calculate sigma only based on fixed effects
  yfit <- predict(mdl, newdata = df_mdl, re.form=NA, type = "response")
  RSS <- sum((y-yfit)^2)
  n <- nrow(data)
  # Get the number of fixed effect parameters (p)
  p <- length(fixef(mdl))
  sigma_manual <- sqrt(RSS / (n - p))
  s <- sigma_manual# s <- sigma(mdl)
  
  
  loopred <- NA
  looDeviance <- NA
  if(!skiploo){
    y_pred <- c()
    y_true <- c()
    for (f_sub in unique(df_mdl$c)){
      print(f_sub)
      if(family == "binomial"){
        mdl_sub <- glmer(fml,  data = df_mdl[!df_mdl$c==f_sub,], family = binomial)
      }else if(family == "gaussian"){
        mdl_sub <- lmer(fml, data = df_mdl[!df_mdl$c==f_sub,])
      }
      y_pred_sub <- predict(mdl_sub, newdata = df_mdl[df_mdl$c==f_sub,], re.form=NA, type = "response", allow.new.levels=T)
      y_true_sub <- df_mdl[df_mdl$c==f_sub,"y"]
      y_pred <- c(y_pred, y_pred_sub)
      y_true <- c(y_true, y_true_sub)
    }
    if(family=="binomial"){
      loopred <- as.numeric(pROC::auc(pROC::roc(response = y_true, predictor = y_pred)))
      looDeviance <- -2*sum(y_true * log(y_pred) + (1 - y_true) * log(1 - y_pred), na.rm = TRUE)
    }
    if(family=="gaussian"){
      loopred <- sqrt(sum( (y_true-y_pred)^2 )) # MSE
      looDeviance <- length(y_true)*log(2*pi*s^2) + 1/s^2*sum( (y_true-y_pred)^2 )  
    }
  }
  
  return(list("mdl" = mdl,
              "loopred"= loopred,
              "looDeviance" = looDeviance))
}



fit_eval_glm <- function(y, c, data, family=c("binomial","gaussian")[1]){
  
  
  df_mdl <- as.data.frame(data)
  df_mdl$y <- y
  df_mdl$c <- as.factor(c)
  
  fix_vars <- grep("^fix\\d+$", colnames(data), value = TRUE)
  rdm_vars <- grep("^rdm\\d+$", colnames(data), value = TRUE)
  fml <- formula(paste0("y~", paste0(c(fix_vars, rdm_vars),collapse = "+")))
  mdl <- glm(fml, data = df_mdl, family = family)
  mdl$c <- c
  
  # calculate sigma only based on fixed effects
  yfit <- predict(mdl, newdata = df_mdl, type = "response")
  RSS <- sum((y-yfit)^2)
  n <- nrow(data)
  # Get the number of fixed effect parameters (p)
  p <- length(coef(mdl))
  s <- sqrt(RSS / (n - p))
  
  
  
  # loo auc and loo deviance
  y_pred <- c()
  y_true <- c()
  for (f_sub in unique(df_mdl$c)){
    mdl_sub <- glm(fml, family = family, data = df_mdl[!df_mdl$c==f_sub,])
    y_pred_sub <- predict(mdl_sub, type = "response", newdata = df_mdl[df_mdl$c==f_sub,])
    y_true_sub <- df_mdl[df_mdl$c==f_sub,"y"]
    y_pred <- c(y_pred, y_pred_sub)
    y_true <- c(y_true, y_true_sub)
  }
  if(family=="binomial"){
    loopred <- as.numeric(pROC::auc(pROC::roc(response = y_true, predictor = y_pred)))
    looDeviance <- -2*sum(y_true * log(y_pred) + (1 - y_true) * log(1 - y_pred), na.rm = TRUE)
  }
  if(family=="gaussian"){
    loopred <- sqrt(sum( (y_true-y_pred)^2 )) # MSE
    # s <- sigma(mdl)
    looDeviance <- length(y_true)*log(2*pi*s^2) + 1/s^2*sum( (y_true-y_pred)^2 )  
  }
  
  
  
  return(list("mdl"=mdl,
              "loopred" = loopred,
              "looDeviance" = looDeviance))
}





fit_eval_gee <- function(y, c, data, family=c("binomial","gaussian")[1], skiploo = F) {
  library(geepack)
  library(pROC)
  
  
  df_mdl <- as.data.frame(data)
  df_mdl$y <- y
  df_mdl$c <- as.factor(c)
  
  fix_vars <- grep("^fix\\d+$", colnames(data), value = TRUE)
  rdm_vars <- grep("^rdm\\d+$", colnames(data), value = TRUE)
  fml <- formula(paste0("y ~ ", paste0(c(fix_vars, rdm_vars), collapse = " + ")))
  
  if (family == "binomial") {
    mdl <- geeglm(fml, data = df_mdl, family = binomial, id = c, corstr = "ar1")
  } else if (family == "gaussian") {
    mdl <- geeglm(fml, data = df_mdl, family = gaussian, id = c, corstr = "ar1")
  }
  
  # calculate sigma only based on fixed effects
  yfit <- predict(mdl, newdata = df_mdl, type = "response")
  RSS <- sum((y-yfit)^2)
  n <- nrow(data)
  # Get the number of fixed effect parameters (p)
  p <- length(coef(mdl))
  s <- sqrt(RSS / (n - p))
  
  
  loopred <- NA
  looDeviance <- NA
  
  if (!skiploo) {
    y_pred <- c()
    y_true <- c()
    for (f_sub in unique(df_mdl$c)) {
      print(f_sub)
      if (family == "binomial") {
        mdl_sub <- geeglm(fml, data = df_mdl[!df_mdl$c == f_sub, ], family = binomial, id = c, corstr = "ar1")
      } else if (family == "gaussian") {
        mdl_sub <- geeglm(fml, data = df_mdl[!df_mdl$c == f_sub, ], family = gaussian, id = c, corstr = "ar1")
      }
      y_pred_sub <- predict(mdl_sub, newdata = df_mdl[df_mdl$c == f_sub, ], type = "response")
      y_true_sub <- df_mdl[df_mdl$c == f_sub, "y"]
      y_pred <- c(y_pred, y_pred_sub)
      y_true <- c(y_true, y_true_sub)
    }
    if (family == "binomial") {
      loopred <- as.numeric(pROC::auc(pROC::roc(response = y_true, predictor = y_pred)))
      looDeviance <- -2 * sum(y_true * log(y_pred) + (1 - y_true) * log(1 - y_pred), na.rm = TRUE)
    }
    if (family == "gaussian") {
      loopred <- sqrt(sum((y_true - y_pred)^2)) # MSE
      # s <- sqrt(summary(mdl)$dispersion$Estimate)
      looDeviance <- length(y_true) * log(2 * pi * s^2) + 1/s^2 * sum((y_true - y_pred)^2)
    }
  }
  
  return(list("mdl" = mdl,
              "loopred" = loopred,
              "looDeviance" = looDeviance))
}


calculate_bias <- function(res, m0, m1){
  fix_effect_betas <- c(grep("^beta_fix_itc$", colnames(res$betas), value = TRUE),
                        grep("^beta_fix\\d+$", colnames(res$betas), value = TRUE),
                        grep("^beta_fix_rdm\\d+$", colnames(res$betas), value = TRUE) )
  fix_effect_betas <- unique(res$betas[,fix_effect_betas])
  bias0 <- sum((coef(summary(m0$mdl))[,"Estimate"] - fix_effect_betas)^2)/length(fix_effect_betas)
  bias1 <- sum((coef(m1$mdl) - fix_effect_betas)^2)/length(fix_effect_betas)
  
  return(list(bias0 = bias0,
              bias1 = bias1))
}

calculate_se_accuracy <- function(res, m0, m1){
  fix_effect_betas <- c(grep("^beta_fix_itc$", colnames(res$betas), value = TRUE),
                        grep("^beta_fix\\d+$", colnames(res$betas), value = TRUE),
                        grep("^beta_fix_rdm\\d+$", colnames(res$betas), value = TRUE) )
  fix_effect_betas <- unique(res$betas[,fix_effect_betas])
  se_ratio0 <- mean(abs(coef(summary(m0$mdl))[,"Estimate"] - fix_effect_betas)/ sqrt(diag(vcov(m0$mdl)))) #coef(summary(m0))[,"Std. Error"])
  se_ratio1 <- mean(abs(coef(m1$mdl) - fix_effect_betas)/sqrt(diag(m1$vcov)))
  
  return(list(se_ratio0 = se_ratio0,
              se_ratio1 = se_ratio1))
  
}


generate_overfit <- function(res, nd = 5, type = c("poly","rcs", "random")[1]){
  res$data_org <- res$data
  # for each of the randome effect predictors add poly 5
  for(rdm in grep("^rdm\\d+$", colnames(res$data), value = TRUE) ){
    if(type=="poly"){
      over_mat <- as.matrix(poly(res$data[,rdm],nd),nrow=nrow(res$data), ncol=nd)
    }
    if(type=="rcs"){
      over_mat <- as.matrix(rms::rcs(res$data[,rdm],nd),nrow=nrow(res$data), ncol=nd)
    }
    if(type=="random"){
      over_mat <- matrix(rnorm(nrow(res$data)*2, 0, 1),nrow=nrow(res$data), ncol=2)
    }
    colnames(over_mat) <- paste0(rdm,c(1:ncol(over_mat)))
    res$data <- cbind(res$data, over_mat)
  }
 return(res)
}

detect_mal <- function(res_df, sim_condition){
  mal <- F
  
  best_df <- data.frame() 
  for(score in c("loodev", "nic", "aic", "bic", "nic")){
    best_size <- res_df$model_size[which(res_df[,score]==min(res_df[,score]))][1]
    best_score <- res_df[,score][which(res_df[,score]==min(res_df[,score]))][1]
    score_1se <- sd(res_df[,score])/sqrt(nrow(res_df))
    best_size_1se_min <- min(res_df$model_size[which(abs(res_df[,score]-min(res_df[,score]))<=score_1se)])
    best_size_1se_max <- max(res_df$model_size[which(abs(res_df[,score]-min(res_df[,score]))<=score_1se)])
    best_df <- bind_rows(best_df, data.frame(score,best_size, best_score, score_1se, best_size_1se_min,best_size_1se_max))
  }
  # library(tidyr)
  # library(dplyr)
  # library(ggplot2)
  # library(ggpubr)
  # res_df_long <- res_df %>%
  #   pivot_longer(
  #     cols = c(nic, aic, bic, dev, loodev),  # Specify columns to lengthen
  #     names_to = "score",  # New column for the names
  #     values_to = "value"  # New column for the values
  #   )
  # # ymin <- res_df_long$value[which(res_df_long$value == min(res_df_long$value[res_df_long$score=="dev"]) )][1]
  # # ymax <- res_df_long$value[res_df_long$score=="loodev"&res_df_long$model_size==1]
  # ggplot(res_df_long, aes(x = model_size, y = value, group = score, color = score)) +
  #   geom_line() +
  #   geom_line(data = filter(res_df_long, score == "loodev"), linetype = "dotted") +
  #   geom_vline(xintercept = sim_condition$n_ttl_betas, color = "grey") +
  #   scale_color_manual(values = c("nic" = "red", "aic" = "blue", "bic" = "orange", "dev" = "gray", "loodev" = "black")) +
  #   theme_minimal() +
  #   geom_errorbar(data = best_df, aes(x = best_size, xmin=best_size_1se_min, xmax=best_size_1se_max, y = best_score, color=score))+
  #   geom_point(data = best_df, aes(x = best_size, y = best_score, color=score), size = 3)+
  #   labs(title = paste0("iter = ",i), x = "Model Size", y = "Value", color = "Score")
  
  
  # sharp decrease in dev
  d <- c(0,diff(res_df$dev))
  for(di in c(sim_condition$n_ttl_betas:length(d)) ){
    e1 <- F
    e2 <- F
    if(d[di] < 5*mean(d[c( (di-3):(di-1) )]) & mean(d[c( (di-3):(di-1) )])<(-1) ){
      e1 <- T
    }
    # if drop in single step is larger than 50% cumulative drop
    if(abs(d[di]) > 0.5*(max(d[c(sim_condition$n_ttl_betas:di)]) - min(d[c(sim_condition$n_ttl_betas:di)])) ){
      e2 <- T
    }
    if(e1&e2){
      print(paste0("sharp decrease in dev at model size ",di))
      mal <- T
    }
  }
  
  
  # loodev best should not be too far away from data generating model
  e <- abs(best_df$best_size[best_df$score=="loodev"] - sim_condition$n_ttl_betas)
  if(e > max(res_df$model_size)/2){
    print("detected mal simulation: loodev too far from generating model")
    mal <- T
  }
  
  # loodev 1se not to wide
  e <- best_df$best_size_1se_max[best_df$score=="loodev"] - best_df$best_size_1se_min[best_df$score=="loodev"] 
  if(e >= 0.8* max(res_df$model_size)){
    print("detected mal simulation: loodev 1se too wide")
    mal <- T
  }
  return(mal)
}
