p_ls <- list()
for(rn in c("lm","lr")){
  
  agg_df <- agg_df_ls[[rn]]
  agg_df <- agg_df %>% filter(sigma_rdm_fix_ratio == 1, ar1_phi == 0.3) %>% as.data.frame()
  
  plot_ls <- list()
  for(st in c("t","loopred","loodev") ){
    plot_df <- NULL
    for(en in c("_median","_q25", "_q75")){
      # time 
      tmpdf <- pivot_longer(agg_df[,c("id", paste0(st,c("0", "1", "2"),en) )], 
                            cols = starts_with(st)&ends_with(en), 
                            names_to = "model", 
                            values_to = paste0("model",en))
      tmpdf$model <- gsub(en,"", tmpdf$model)
      if(is.null(plot_df)){plot_df <- tmpdf}
      else{
        plot_df <- merge(plot_df, tmpdf)
      }
    }
    colnames(plot_df) <- gsub("_median","_m",colnames(plot_df))
    colnames(plot_df) <- gsub("_q25","_l",colnames(plot_df))
    colnames(plot_df) <- gsub("_q75","_u",colnames(plot_df))
      
    plot_df <- merge(plot_df,simulation_conditions, by="id", all.x=T) %>% as.data.frame()
    
    plot_df$n_obs_per_cluster <- factor(plot_df$n_obs_per_cluster, levels=sort(unique(plot_df$n_obs_per_cluster)) )
    levels(plot_df$n_obs_per_cluster) <- paste0(levels(plot_df$n_obs_per_cluster), " obs/cluster")
    plot_df <- plot_df %>%
      mutate(model_type = case_when(
        grepl("0$", model) ~ "RE",
        grepl("1$", model) ~ "GLM",
        grepl("2$", model) ~ "GEE",
        TRUE ~ NA_character_
      ))
    
    plot_df$model_type <- factor(plot_df$model_type, levels=c("GLM", "GEE","RE"))
    ylab <- NULL
    if(st=="t") ylab <- "Time (sec)"
    if(st=="loopred" & rn=="lr") ylab <- "looAUROC"
    if(st=="loopred" & rn=="lm") ylab <- "looMSE"
    if(st=="loodev") ylab <- "looDeviance /obs"
    
    plot_df <- plot_df %>%
      mutate(n_ttl_betas = case_when(
        model_type %in% c("GEE") ~ n_ttl_betas - 0.3,
        model_type %in% c("RE") ~ n_ttl_betas + 0.3,
        model_type %in% c("GLM") ~ n_ttl_betas + 0,
        TRUE ~ n_ttl_betas
      ))
    
    plot_ls[[st]] <- ggplot(data = plot_df, aes(x = n_ttl_betas, y = model_m, color = model_type)) + 
      geom_point(size=1) +
      geom_line(linewidth=0.3, linetype="dotted") + 
      geom_errorbar(aes(ymin = model_l, ymax = model_u),width=0.2) + 
      scale_x_continuous(limits = c(2, 8), breaks = c(3,5,7)) +
      # coord_trans(y = "sqrt") +
      facet_wrap(~ n_obs_per_cluster, ncol=3, nrow=1, scales="free_x") + 
      labs(x = NULL, 
           y = ylab,
           color = "Model") + 
      theme_bw()+
      scale_color_manual(values = scale_color) +
      theme(text = element_text(face = "bold"),
            axis.text = element_text(size=8),
            legend.title = element_text(size=10), 
            legend.text = element_text(size=8)) 
    
    
    
  }
  
  p <- ggarrange(plotlist = plot_ls, nrow=3, ncol=1, common.legend = T, legend = ifelse(rn=="lm","none","right"))
  
  p_ls[[rn]] <- annotate_figure(p, 
                                top=text_grob(ifelse(rn=="lm", "Gaussian", "Binomial"), size = 12, face = "bold"),
                                bottom = text_grob("Generating model size", size = 10, face = "bold"))
  
  
}


library(patchwork)
p <- p_ls[[1]] + p_ls[[2]] + 
  plot_layout(ncol = 2, widths = c(0.85, 1), guides = "collect") & 
  theme(plot.margin = unit(c(.3, .5, .2, .5), "cm")) #(top, right, bottom, left)

# p <- ggarrange(plotlist = p_ls, nrow=1, ncol=2, 
#                widths = c(0.9,1),
#                common.legend = T,legend = "right")
# 

library(cowplot)
p <- ggdraw(p) +
  draw_label("A. Cluster size", x = 0.01, y = 0.99, hjust = 0, vjust = 1, 
             fontface = 'bold', size = 13) 



