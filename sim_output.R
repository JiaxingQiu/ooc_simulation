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
  colnames(res_df)
  
  summa <- function(df){
    for(cl in c("t0", "t1", "t2", 
                "loopred0", "loopred1", "loopred2",
                "loodev0", "loodev1", "loodev2") ){
      if(grepl("dev",cl)) df[[cl]] <- df[[cl]] / df[["N"]]
      df[[paste0(cl,"_mean")]] <- mean(df[[cl]], na.rm=T)
      df[[paste0(cl,"_median")]] <- median(df[[cl]], na.rm=T)
      df[[paste0(cl,"_q25")]] <- quantile(df[[cl]],0.025, na.rm=T)
      df[[paste0(cl,"_q75")]] <- quantile(df[[cl]],0.975, na.rm=T)
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

p_all <- ggarrange(plotlist=list(p_n, p_p, p_r), ncol=1, nrow=3)
p_all <- annotate_figure(p_all, top = text_grob("Leave-one-cluster-out prediction performance of\n GLM, GEE (w/ AR1 correlation), and RE (generating model)", size = 14, face = "bold"))
p_all %>% ggsave(filename=paste0("./res/ooc_performance.png"), width = 10, height = 16, bg="white")


