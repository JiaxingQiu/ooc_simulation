p_ls <- list()
for(rn in c("lm","lr")){
  pp_ls <- list()
  for(cs in c("raw", "unbalanced", "small", "rare")){
    plot_ls <- list()
    agg_df <- agg_df_ls[[rn]] %>% filter(case == cs)
    
      
    for(st in c("nc_rate","t","loopred","loodev") ){
      if(st =="nc_rate"){
        plot_df <- pivot_longer(agg_df[,c("id", paste0(st,c("0", "1", "2")) )], 
                                cols = starts_with(st), 
                                names_to = "model", 
                                values_to = "nc_rate")
      }else{
        plot_df <- NULL
        for(en in c("_median","_q25", "_q75","_q025", "_q975")){
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
        colnames(plot_df) <- gsub("_q25","_i",colnames(plot_df))
        colnames(plot_df) <- gsub("_q75","_r",colnames(plot_df))
        colnames(plot_df) <- gsub("_q025","_l",colnames(plot_df))
        colnames(plot_df) <- gsub("_q975","_u",colnames(plot_df))
      }
      
      plot_df <- merge(plot_df,simulation_conditions, by="id", all.x=T) %>% as.data.frame()
      
      
      plot_df <- plot_df %>%
        mutate(model_type = case_when(
          grepl("0$", model) ~ "RE",
          grepl("1$", model) ~ "GLM",
          grepl("2$", model) ~ "GEE",
          TRUE ~ NA_character_
        ))
      
      plot_df$model_type <- factor(plot_df$model_type, levels=c("GLM", "GEE","RE"))
      # fix ylab 
      ylab <- NULL
      if(cs =="raw"){
        if(st=="nc_rate"){
          ylab <- "Not converge (%)"
        } 
        if(st=="t"){
          ylab <- "Time (sec)"
        } 
        if(st=="loopred" & rn=="lr") {
          ylab <- "looAUC"
        }
        if(st=="loopred" & rn=="lm") {
          ylab <- "looMSE"
        }
        if(st=="loodev") {
          ylab <- "looDeviance /obs"
        }
      }
      # fix ylim
      if(st=="nc_rate") ylim <- c(0,0.7)
      if(st=="t"){
        if(rn=="lm") ylim <- c(0,110)
        if(rn=="lr") ylim <- c(0,3300)
      } 
      if(st=="loopred" & rn=="lr") {
        ylim <- c(0.4,1)
      }
      if(st=="loopred" & rn=="lm") {
        ylim <- c(150, 1000)
      }
      if(st=="loodev") {
        if(rn=="lm") ylim <- c(7.5, 10)
        if(rn=="lr") ylim <- c(0, 30.5)
      }
      
      plot_df <- plot_df %>%
        mutate(n_ttl_betas = case_when(
          model_type %in% c("GEE") ~ n_ttl_betas - 0.3,
          model_type %in% c("RE") ~ n_ttl_betas + 0.3,
          model_type %in% c("GLM") ~ n_ttl_betas + 0,
          TRUE ~ n_ttl_betas
        ))
      if(st=="nc_rate"){
        plot_ls[[st]] <- ggplot(data = plot_df, aes(x = n_ttl_betas, y = nc_rate, color = model_type)) + 
          geom_point(size=2, shape = 4) +
          geom_line(linewidth=0.3, linetype="dotted") + 
          scale_x_continuous(limits = c(2, 8), 
                             breaks = c(3,5,7)) +
          ylim(ylim) + 
          labs(subtitle = NULL,
               x = NULL, 
               y = ylab,
               color = "Model") + 
          theme_bw()+
          scale_color_manual(values = scale_color) +
          theme(plot.subtitle = element_text(size=10, face="bold"),
                text = element_text(face = "bold"),
                axis.text = element_text(size=8),
                axis.title = element_text(size=10),
                legend.title = element_text(size=10), 
                legend.text = element_text(size=8)) 
      }else{
        plot_ls[[st]] <- ggplot(data = plot_df, aes(x = n_ttl_betas, y = model_m, color = model_type)) + 
          # geom_point(size=1) +
          # geom_line(linewidth=0.3, linetype="dotted") + 
          # geom_errorbar(aes(ymin = model_l, ymax = model_u),width=0.2) + 
          geom_boxplot(
            aes(ymin = model_l,
                lower = model_i,
                middle = model_m,
                upper = model_r,
                ymax = model_u,
                group=n_ttl_betas), stat = "identity") +
          scale_x_continuous(limits = c(2, 8), 
                             breaks = c(3,5,7)) +
          ylim(ylim) + 
          labs(subtitle = NULL,
               x = NULL, 
               y = ylab,
               color = "Model") + 
          theme_bw()+
          scale_color_manual(values = scale_color) +
          theme(plot.subtitle = element_text(size=10, face="bold"),
                text = element_text(face = "bold"),
                axis.text = element_text(size=8),
                axis.title = element_text(size=10),
                legend.title = element_text(size=10), 
                legend.text = element_text(size=8)) 
      }
    }
    p <- ggarrange(plotlist = plot_ls, nrow=4, ncol=1, common.legend = T, legend = "none")
    stitle <- NULL
    if(cs=="raw") stitle <- "    Original"
    if(cs=="unbalanced") stitle <- "    Unbalanced clusters"
    if(cs=="small") stitle <- "    Small samples"
    if(cs=="rare") stitle <- "    Rare event"
    
    pp_ls[[cs]] <- annotate_figure(p, 
                                   top=text_grob(stitle, size = 10, face = "bold"),
                                   bottom = text_grob("Model size", size = 10, face = "bold"))
    
  }
  
  if(rn=="lm") pp_ls <- pp_ls[1:3]
  p <- ggarrange(plotlist = pp_ls, nrow=1, widths=c(1.1,1,1,1), legend.grob = get_legend(plot_ls[["t"]]), legend = "right")
  p_ls[[rn]] <- annotate_figure(p,top=text_grob(ifelse(rn=="lm", "Gaussian", "Binomial"), size = 13, face = "bold")) # , x=0, hjust=0
  
  
}


library(patchwork)
p <- p_ls[[1]] + p_ls[[2]] + 
  plot_layout(ncol = 2, widths = c(0.8,1), guides = "collect") & 
  theme(plot.margin = unit(c(.2, .2, .2, .2), "cm")) #(top, right, bottom, left)
library(cowplot)
p <- ggdraw(p) 

