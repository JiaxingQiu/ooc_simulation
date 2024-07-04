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


# # # for test: 
# sim_condition = simulation_conditions[which(simulation_conditions$id==12),] 
# family =  "binomial" # "gaussian" #

# important:
# <simpleError: number of observations (=n_cluster * n_obs_per_cluster) <= number of random effects (= n_ttl_effect * n_cluster) for term (rdm1 + rdm2 + rdm3 + rdm4 + 1 | c); 
# the random-effects parameters and the residual variance (or scale parameter) are probably unidentifiable>

run_wrapper <- function(sim_condition, family) {
  results_list = list()
  for(i in 1:sim_condition$iter){
    res <- generate_data(sim_condition$n_cluster,
                         sim_condition$n_obs_per_cluster,
                         sim_condition$n_ttl_betas, 
                         sim_condition$fix_rdm_ratio,
                         sim_condition$sigma_fix,
                         sim_condition$sigma_rdm_fix_ratio,
                         sim_condition$ar1_phi,
                         sim_condition$na_rate,
                         family = family)
    
    # rare case scenario
    if("case" %in% colnames(sim_condition) & family == "binomial" ){
      if(sim_condition$case=="rare"){
        # control roughly the prevalence
        while(T){
          res <- generate_data(sim_condition$n_cluster,
                               sim_condition$n_obs_per_cluster,
                               sim_condition$n_ttl_betas, 
                               sim_condition$fix_rdm_ratio,
                               sim_condition$sigma_fix,
                               sim_condition$sigma_rdm_fix_ratio,
                               sim_condition$ar1_phi,
                               sim_condition$na_rate,
                               family = family)
          res$p <- res$p*0.015
          res$y <- rbinom(length(res$p), 1, res$p)
          prev <- length(unique(res$c[which(res$y==1)])) / length(unique(res$c))
          print(prev)
          if(prev > 0.1 & prev < 0.2) break
        }
      }
    }
    
    # glme mixed effect model
    m0 <- NULL
    tryCatch({
      tic()
      m0 <- fit_eval_glmer(y = res$y,
                           c = res$c,
                           data = res$data,
                           family = family)
      t <- toc(log = TRUE)
      t0 <- as.numeric( t$toc - t$tic )
    },error=function(e){
      print(e)
    })
    
    
    # glm model 
    m1 <- NULL
    tryCatch({
      tic()
      m1 <- fit_eval_glm(y = res$y,
                         c = res$c,
                         data = res$data,
                         family = family)
      t <- toc(log = TRUE)
      t1 <- as.numeric( t$toc - t$tic )
    },error=function(e){
      print(e)
    })
    
    # gee model
    m2 <- NULL
    tryCatch({
      tic()
      m2 <- fit_eval_gee(y = res$y,
                         c = res$c,
                         data = res$data,
                         family = family)
      t <- toc(log = TRUE)
      t2 <- as.numeric( t$toc - t$tic )
      
    },error=function(e){
      print(e)
    })
    
    
    results_list[[i]] = list(id = sim_condition$id, 
                             iter = i, 
                             N = nrow(res$data),
                             t0 = ifelse(is.null(m0),-1,t0),
                             t1 = ifelse(is.null(m1),-1,t1),
                             t2 = ifelse(is.null(m2),-1,t2),
                             loopred0 = ifelse(is.null(m0),-1,m0$loopred),
                             loopred1 = ifelse(is.null(m1),-1,m1$loopred),
                             loopred2 = ifelse(is.null(m2),-1,m2$loopred),
                             loodev0 = ifelse(is.null(m0),-1,m0$looDeviance),
                             loodev1 = ifelse(is.null(m1),-1,m1$looDeviance),
                             loodev2 = ifelse(is.null(m2),-1,m2$looDeviance) )
    
    
  }
  results_list <- Filter(function(x) !is.null(x), results_list)
  toReturn = do.call("rbind", results_list)
  return(toReturn)
}


run_wrapper_lm <- function(sim_condition) run_wrapper(sim_condition, family="gaussian")
run_wrapper_lr <- function(sim_condition) run_wrapper(sim_condition, family="binomial")


