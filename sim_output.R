setwd(dirname(rstudioapi::getSourceEditorContext()$path))
rm(list = ls())

library(dplyr)
library(rslurm)
library(ggplot2)
library(tidyr)
library(ggpubr)

lr_output_fname <- "output_ooc_run_lr.RDS"
lm_output_fname <- "output_ooc_run_lm.RDS"

agg_df_ls <- list()
for(rn in c("lr", "lm")){ 
  if(rn == "lr") output <- readRDS(paste0("./res/", lr_output_fname))
  if(rn == "lm") output <- readRDS(paste0("./res/", lm_output_fname))
  source("./sim_conditions.R")
  simulation_conditions <- simulation_conditions[,setdiff(colnames(simulation_conditions),"iter")]
  res_df <- as.data.frame(merge(output, simulation_conditions, by="id", all.x=T))
  res_df[res_df<0] <- NA 
  
  # deal with AUC < 0.5 for lr
  if(rn == "lr"){
    res_df[which(res_df$loopred0<0.5),c("loopred0", "loodev0")] <- NA
    res_df[which(res_df$loopred1<0.5),c("loopred1", "loodev1")] <- NA
    res_df[which(res_df$loopred2<0.5),c("loopred2", "loodev2")] <- NA
  }
  # add a column for nc_rate
  na_df <- res_df %>%
    group_by(id) %>% 
    summarise(nc_rate0 = mean(is.na(t0)),
              nc_rate1 = mean(is.na(t1)),
              nc_rate2 = mean(is.na(t2))) %>%
    as.data.frame()
  
  summa <- function(df){
    for(cl in c("t0", "t1", "t2", 
                "loopred0", "loopred1", "loopred2",
                "loodev0", "loodev1", "loodev2") ){
      if(grepl("dev",cl)) df[[cl]] <- df[[cl]] / df[["N"]]
      df[[paste0(cl,"_mean")]] <- mean(df[[cl]], na.rm=T)
      df[[paste0(cl,"_median")]] <- median(df[[cl]], na.rm=T)
      df[[paste0(cl,"_q25")]] <- quantile(df[[cl]],0.25, na.rm=T)
      df[[paste0(cl,"_q75")]] <- quantile(df[[cl]],0.75, na.rm=T)
      df[[paste0(cl,"_se")]] <- sd(df[[cl]], na.rm=T)/nrow(df)
    }
    return(df)
  }
  by_id <- group_by(res_df, id)
  agg_df <- do(by_id, summa(.))
  agg_df <- distinct(agg_df[,c("id",colnames(agg_df)[endsWith(colnames(agg_df), "_mean") |
                                                endsWith(colnames(agg_df), "_median") |
                                                endsWith(colnames(agg_df), "_q25") |
                                                endsWith(colnames(agg_df), "_q75") |
                                                endsWith(colnames(agg_df), "_se")]) ])
  agg_df <- as.data.frame(agg_df)
  agg_df <- merge(agg_df,na_df,by="id",all.x = T)
  agg_df_ls[[rn]] <- merge(agg_df, simulation_conditions, by="id", all.x=T)
}

scale_color <- c("GLM" = "red",  
                 "RE" = "blue", 
                 "GEE" = "darkorange")
source("./sim_plot_N.R")
p_n <- p
source("./sim_plot_phi.R")
p_p <- p
source("./sim_plot_sigma_ratio.R")
p_r <- p
source("./sim_plot_nc_rate.R")

p_all <- ggarrange(plotlist=list(p_n, p_p, p_r), ncol=1, nrow=3)
p_all <- annotate_figure(p_all, top = text_grob("Prediction performance & Computational efficiency\n of GLM, GEE (w/ AR1 correlation), and RE (generating model)", size = 14, face = "bold"))
p_all %>% ggsave(filename=paste0("./res/ooc_performance.png"), width = 10, height = 15, bg="white")


