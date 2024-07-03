mdl_tradeoffs <- function(
  y_true, # vector of y_true values
  y_hat # vector of predicted y_hat, it can be probability, log odds, fold risk
  #rel_time = NULL # vector of continuous relative time from episodes
){
  y_hat_min <- min(y_hat, na.rm = TRUE)
  y_hat_max <- max(y_hat, na.rm = TRUE)
  step <- (y_hat_max-y_hat_min)/100
  thresholds <- seq(y_hat_min,y_hat_max,step)
  t_true <- mean(y_true, na.rm = TRUE)
  confusion_df_all <- data.frame()
  for (t in thresholds[2:length(thresholds)] ){
    # calculate TP, TN, FP, FN
    confusion_df <- data.frame(
      threshold = t,
      TP = sum(y_hat>=t&y_true>=t_true, na.rm = TRUE),
      TN = sum(y_hat<t&y_true<t_true, na.rm = TRUE),
      FP = sum(y_hat>=t&y_true<t_true, na.rm = TRUE),
      FN = sum(y_hat<t&y_true>=t_true, na.rm = TRUE)
    )
    confusion_df_all <- bind_rows(confusion_df_all,confusion_df)
  }
  # calculate x and y axis
  confusion_df_all <- confusion_df_all %>% mutate(recall = TP/(TP+FN),
                                                  precision = TP/(TP+FP),
                                                  specificity = TN/(FP+TN),
                                                  sensitivity = TP/(TP+FN),
                                                  alarm_rate = (TP+FP)/(TP+TN+FP+FN)) %>%
    as.data.frame()
  # add ending points for plotting
  confusion_df_all <- bind_rows(confusion_df_all,
                                data.frame(specificity=c(1,0), 
                                           sensitivity=c(0,1),
                                           precision=c(1,0),
                                           recall=c(0,1),
                                           alarm_rate=c(0,1)
                                           ))
  
  # ROC
  df_plot <- confusion_df_all%>% group_by(specificity) %>% summarise(sensitivity = max(sensitivity)) %>% as.data.frame()
  roc <- ggplot(data = df_plot, aes(x=1-specificity, y=sensitivity ) ) +
    geom_line() + xlim(0,1) + ylim(0,1) + theme_bw() +
    geom_abline(slope=1,intercept=0,linetype='dashed') +
    labs(x="1 - Specificity", y="Sensitivity", title="ROC")
  
  # PRC
  df_plot <- confusion_df_all%>% group_by(recall) %>% summarise(precision = max(precision)) %>% as.data.frame()
  prc <- ggplot(data = df_plot, aes(x=recall, y=precision ) ) +
    geom_line() + xlim(0,1) + ylim(0,1) + theme_bw() +
    geom_abline(slope=0,intercept=0,linetype='dashed') +
    labs(x="Recall", y="Precision", title = "PRC")
  
  # SAC (sensitivity VS alarm ratio)
  df_plot <- confusion_df_all%>% group_by(alarm_rate) %>% summarise(sensitivity = max(sensitivity)) %>% as.data.frame()
  sac <- ggplot(data = df_plot, aes(x=alarm_rate, y=sensitivity ) ) +
    geom_line() + xlim(0,1) + ylim(0,1) + theme_bw() +
    geom_abline(slope=1,intercept=0,linetype='dashed') +
    labs(x="Alarm Rate", y="Sensitivity", title = "SAC")
  
  return(list(roc = roc,
              prc = prc,
              sac = sac))
}
