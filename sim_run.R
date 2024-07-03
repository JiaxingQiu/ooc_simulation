# setwd(dirname(rstudioapi::getSourceEditorContext()$path))

rm(list = ls())

library(dplyr)
library(rslurm)
library(Matrix)
library(lme4)
library(geepack)
library(MASS)
library(pROC)
library(tictoc)

list.of.packages <- c("dplyr",
                      "rslurm",
                      "MASS",
                      "lme4",
                      "Matrix",
                      "pROC",
                      "tictoc",
                      "geepack")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, lib = "/sfs/qumulo/qhome/jq2uw/R/goolf/4.3")



source("./sim_functions.R")
for(u in c("nic", "ass", "stp", "do")){
  path = paste0("./utils/",u,"_utils")
  flst = list.files( path)
  sapply(c(paste(path,flst,sep="/")), source, .GlobalEnv)
}
source("./sim_conditions.R")



# for test: 
sim_condition = simulation_conditions[which(simulation_conditions$id==11),] # 11


run_wrapper <- function(sim_condition, family) {
  results_list = list()
  for(i in 1:sim_condition$iter){
    tryCatch({
      res <- generate_data(sim_condition$n_cluster,
                           sim_condition$n_obs_per_cluster,
                           sim_condition$n_ttl_betas, 
                           sim_condition$fix_rdm_ratio,
                           sim_condition$sigma_fix,
                           sim_condition$sigma_rdm_fix_ratio,
                           sim_condition$ar1_phi,
                           sim_condition$na_rate,
                           family = family)
      
      # glme mixed effect model
      tic()
      m0 <- fit_eval_glmer(y = res$y,
                           c = res$c,
                           data = res$data,
                           family = family)
      t <- toc(log = TRUE)
      t0 <- as.numeric( t$toc - t$tic )
      
      
      # glm model 
      tic()
      m1 <- fit_eval_glm(y = res$y,
                         c = res$c,
                         data = res$data,
                         family = family)
      t <- toc(log = TRUE)
      t1 <- as.numeric( t$toc - t$tic )
      
      # gee model
      tic()
      m2 <- fit_eval_gee(y = res$y,
                         c = res$c,
                         data = res$data,
                         family = family)
      t <- toc(log = TRUE)
      t2 <- as.numeric( t$toc - t$tic )
      
      
      
      results_list[[i]] = list(id = sim_condition$id, 
                               iter = i, 
                               N = nrow(res$data),
                               t0 = t0,
                               t1 = t1,
                               t2 = t2,
                               loopred0 = m0$loopred,
                               loopred1 = m1$loopred,
                               loopred2 = m2$loopred,
                               loodev0 = m0$looDeviance,
                               loodev1 = m1$looDeviance,
                               loodev2 = m2$looDeviance)
      
    }, error = function(e){
      print(e)
      print(paste0("skip iteration ",i))
    })
  }
  results_list <- Filter(function(x) !is.null(x), results_list)
  toReturn = do.call("rbind", results_list)
  return(toReturn)
}



run_wrapper_lm <- function(sim_condition) run_wrapper(sim_condition, family="gaussian")
run_wrapper_lr <- function(sim_condition) run_wrapper(sim_condition, family="binomial")


sjob_lm = slurm_map(
  split(simulation_conditions, simulation_conditions$id),
  run_wrapper_lm,
  nodes=nrow(simulation_conditions),
  cpus_per_node = 1,
  jobname = "ooc_run_lm",
  submit = TRUE,
  preschedule_cores = F,
  slurm_options =
    c(account = "netlab", partition = "standard", time = "2-00:00:00"), # standard
  global_objects = lsf.str()
)
save(sjob_lm, file = "ooc_run_lm.Rdata")

sjob_lr = slurm_map(
  split(simulation_conditions, simulation_conditions$id),
  run_wrapper_lr,
  nodes=nrow(simulation_conditions),
  cpus_per_node = 1,
  jobname = "ooc_run_lr",
  submit = TRUE,
  preschedule_cores = F,
  slurm_options =
    c(account = "netlab", partition = "standard", time = "3-00:00:00"), # standard
  global_objects = lsf.str()
)
save(sjob_lr, file = "ooc_run_lr.Rdata")






