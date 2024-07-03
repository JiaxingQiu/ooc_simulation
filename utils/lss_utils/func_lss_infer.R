lss_infer <- function(
    lasso_obj = lasso_obj, # a grouped lasso regression model object
    lasso_by = lasso_by
){
  
  tune_trace_plot <- NULL
  coef_trace_plot <- NULL
  opt_model_df <- NULL
  
  if(lasso_by == "none"){
    
    # ---- lambda tuning trace plot ----
    plot.new() 
    plot(lasso_obj$tune_trace, main = "Lasso penalty\n\n", 
         xlab="Log Lambda",
         xlim=rev(c(min(log(lasso_obj$tune_trace$lambda)),max(log(lasso_obj$tune_trace$lambda)))))
    abline(v = log(lasso_obj$tune_trace$lambda.min), col = "red", lty = "dashed")
    abline(v = log(lasso_obj$tune_trace$lambda.1se), col = "blue", lty = "dashed")
    tune_trace_plot <- recordPlot()
    
    # ---- coefficient trace plot ----
    lasso_coef <- NULL
    for(s in seq(length(lasso_obj$coef_trace$lambda),1,-1) ){
      tmp_coeffs <- coef(lasso_obj$coef_trace, s=lasso_obj$coef_trace$lambda[s])
      coef_df <- data.frame(varname = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coef = tmp_coeffs@x)
      colnames(coef_df)[which(colnames(coef_df)=="coef")] <- paste0("s",s)
      if(is.null(lasso_coef)){
        lasso_coef <- coef_df
      }else{
        lasso_coef <- merge(lasso_coef, coef_df, by="varname", all.x=TRUE)
      }
    }
    rownames(lasso_coef) <- lasso_coef$varname
    lasso_coef <- lasso_coef[setdiff(rownames(lasso_coef),c("(Intercept)")),setdiff(colnames(lasso_coef),c("varname"))]
    lasso_coef <- as.data.frame( t(lasso_coef) )
    lasso_coef$s <- as.numeric( gsub("s","", rownames(lasso_coef)))
    lmd_s <- c()
    var_s <- c()
    for (var in setdiff(colnames(lasso_coef),"s") ){
      lmd_s <- c(lmd_s, min(lasso_coef$s[which(!is.na(lasso_coef[,var]))], na.rm=TRUE) )# max(lasso_coef$s[which(is.na(lasso_coef[,var]))], na.rm=TRUE)
      var_s <- c(var_s, var)
    }
    lambda_zero_coef <- data.frame(varname = var_s, 
                                   s = lmd_s,
                                   lambda = lasso_obj$coef_trace$lambda[lmd_s+1],
                                   log_lambda = log(lasso_obj$coef_trace$lambda[lmd_s+1]))
    lambda_zero_coef <- lambda_zero_coef[order(lambda_zero_coef$log_lambda),]
    rownames( lambda_zero_coef ) <- NULL
    plot.new() 
    plot(lasso_obj$coef_trace, main = "Lasso penalty\n\n",
         xvar = "lambda",
         xlim=rev(c(min(log(lasso_obj$tune_trace$lambda)),max(log(lasso_obj$tune_trace$lambda)))) )
    abline(v = log(lasso_obj$tune_trace$lambda.min), col = "red", lty = "dashed")
    abline(v = log(lasso_obj$tune_trace$lambda.1se), col = "blue", lty = "dashed")
    text(lambda_zero_coef[,'log_lambda'], 0, lambda_zero_coef[,'varname'], 
         cex=0.7, col="black",srt=90)
    coef_trace_plot <- recordPlot()
    
    # ---- selected variable table ----
    coef_df <- data.frame(coef=as.numeric(lasso_obj$opt_model$beta))
    coef_df$varname <- as.character(rownames(lasso_obj$opt_model$beta))
    opt_model_df <- coef_df[order(-abs(coef_df$coef)),c("varname","coef")]
  }
  if(lasso_by == "group"){
    
    # ---- lambda tuning trace plot ----
    plot.new() 
    plot(lasso_obj$tune_trace, main = "Lasso penalty\n\n", 
         xlab="Log Lambda",
         xlim=rev(c(min(log(lasso_obj$tune_trace$lambda)),max(log(lasso_obj$tune_trace$lambda))))) 
    abline(v = log(lasso_obj$tune_trace$lambda.min), col = "red", lty = "dashed")
    abline(v = log(lasso_obj$tune_trace$lambda.1se), col = "blue", lty = "dashed")
    tune_trace_plot <- recordPlot()
    
    # ---- coefficient trace plot ----
    
    lasso_coef <- as.data.frame( lasso_obj$coef_trace$beta)
    lasso_coef <- as.data.frame( t(lasso_coef) )
    lasso_coef$s <- as.numeric( gsub("s","", rownames(lasso_coef)))
    lmd_s <- c()
    var_s <- c()
    for (var in setdiff(colnames(lasso_coef),"s") ){
      lmd_s <- c(lmd_s, min(lasso_coef$s[which(!lasso_coef[,var]==0)], na.rm=TRUE) )
      var_s <- c(var_s, var)
    }
    lambda_zero_coef <- data.frame(varname = var_s, 
                                   s = lmd_s,
                                   lambda = lasso_obj$coef_trace$lambda[lmd_s+1],
                                   log_lambda = log(lasso_obj$coef_trace$lambda[lmd_s+1]))
    lambda_zero_coef <- lambda_zero_coef[order(lambda_zero_coef$log_lambda),]
    rownames( lambda_zero_coef ) <- NULL
    plot.new() 
    plot(lasso_obj$coef_trace, main = "Lasso penalty\n\n",
         xvar = "lambda",
         xlim=rev(c(min(log(lasso_obj$tune_trace$lambda)),max(log(lasso_obj$tune_trace$lambda)))) )
    abline(v = log(lasso_obj$tune_trace$lambda.min), col = "red", lty = "dashed")
    abline(v = log(lasso_obj$tune_trace$lambda.1se), col = "blue", lty = "dashed")
    text(lambda_zero_coef[,'log_lambda'], 0, lambda_zero_coef[,'varname'], 
         cex=0.7, col="black",srt=90)
    coef_trace_plot <- recordPlot()
    
    # ---- selected variable table ----
    coef_df <- data.frame(coef=as.numeric(lasso_obj$opt_model$beta))
    coef_df$varname <- as.character(rownames(lasso_obj$opt_model$beta))
    opt_model_df <- coef_df[order(-abs(coef_df$coef)),c("varname","coef")]
    
  }
  if(lasso_by == "cluster"){
    # ---- lambda tuning trace plot ----
    tune_df <- lasso_obj$tune_trace
    tune_by <- setdiff(colnames(tune_df),"lambda")
    tune_df$loss <- tune_df[,tune_by]
    tune_df <- tune_df[which(tune_df$lambda>0),]
    tune_df$panelty <- log(tune_df$lambda)#log(1/tune_df$lambda)
    tune_trace_plot <- ggplot(data = tune_df, aes(x=panelty, y=loss))+
      geom_line()+
      geom_vline(aes(xintercept=log(lasso_obj$opt_model$lambda.opt) ), color="red", linetype="dashed")+ #log(1/lasso_obj$opt_model$lambda.opt)
      geom_point(size=2,colour="red")+
      labs(x="Log Lambda", y=tune_by)+
      theme_minimal()
      
    # ---- coefficient trace plot ----
    coef_df_all <- data.frame()
    zero_df_all <- data.frame()
    for(var in setdiff(colnames(lasso_obj$coef_trace),c("lambda","(Intercept)")) ){
      coef_df <- lasso_obj$coef_trace[,c(var, "lambda")]
      coef_df <- coef_df[order(-abs(coef_df$lambda)),]
      coef_df$tag <- ifelse(coef_df[,var]==0,0,1)
      expand_df <- data.frame(var=0,lambda=max(lasso_obj$coef_trace$lambda,na.rm = TRUE)+10,tag=0)
      colnames(expand_df)[which(colnames(expand_df)=="var")] <- var
      coef_df <- bind_rows(expand_df,coef_df)
      coef_df$varname <- var
      coef_df$varvalue <- coef_df[,var]
      coef_df_all <- bind_rows(coef_df_all, coef_df[,c("lambda", "varname","varvalue")])
      zero_df <- coef_df[which(c(0,diff(coef_df$tag))==1)[1],c("varname", "lambda")]
      zero_df_all <- bind_rows(zero_df_all, zero_df)
    }
    coef_df_all <- coef_df_all[which(coef_df_all$lambda>0),]
    zero_df_all <- zero_df_all[which(zero_df_all$lambda>0),]
    coef_df_all$panelty <- log(coef_df_all$lambda)#log(1/coef_df_all$lambda)
    zero_df_all$panelty <- log(zero_df_all$lambda)#log(1/zero_df_all$lambda)
    coef_trace_plot <- ggplot(data = coef_df_all, aes(x=panelty, y=varvalue, color=varname))+
      geom_line()+
      geom_vline(aes(xintercept= log(lasso_obj$opt_model$lambda.opt) ), color="red", linetype="dashed")+
      geom_text(data=zero_df_all, aes(x=panelty,y=0, label=varname),hjust=0.5,vjust=-0.5,angle=90, colour="black")+
      labs(x="Log Lambda", y="Coefficient")+
      theme_minimal()+
      theme(legend.position = "none")
    lambda_zero_coef <- zero_df_all
    lasso_coef <- coef_df_all
    
    # ---- selected variable table ----
    coef_df <- as.data.frame(lasso_obj$opt_model$coefficients)
    colnames(coef_df) <- "coef"
    coef_df$varname <- rownames(coef_df)
    rownames( coef_df ) <- NULL
    opt_model_df <- coef_df[order(-abs(coef_df$coef)), ]
    
  }
  
  return(list(tune_trace_plot = tune_trace_plot,
              coef_trace_plot = coef_trace_plot,
              opt_model_df = opt_model_df,
              lambda_zero_coef = lambda_zero_coef,
              lasso_coef_trace = lasso_coef
              ))
}
