setwd(dirname(rstudioapi::getSourceEditorContext()$path))
rm(list = ls())

library(dplyr)
library(rslurm)
library(ggplot2)
library(tidyr)
library(ggpubr)


# lr
output <- readRDS("./res/output_ooc_run_lr.RDS")
source("./sim_conditions.R")
setdiff(unique(simulation_conditions$id),unique(output$id))

rest_id <- c(4, 7, 13, 16, 22, 25, 31, 34, 40, 43, 52, 58, 61, 67, 70, 76, 79)


# special lr
output <- readRDS("./res/output_ooc_run_lr_special.RDS")
source("./sim_conditions_special.R")
setdiff(unique(simulation_conditions$id),unique(output$id))
rest_id <- c(3)




