setwd(dirname(rstudioapi::getSourceEditorContext()$path))

rm(list = ls())

library(dplyr)
library(rslurm)
library(Matrix)
library(lme4)
library(geepack)
library(MASS)
library(pROC)
library(tictoc)



source("./sim_functions.R")
for(u in c("nic", "ass", "stp", "do")){
  path = paste0("./utils/",u,"_utils")
  flst = list.files( path)
  sapply(c(paste(path,flst,sep="/")), source, .GlobalEnv)
}
source("./sim_conditions_rare_case.R")



# for test: 
sim_condition = simulation_conditions[which(simulation_conditions$id==11),] # 11
family="binomial"

#prevalence = c(0.1, 0.05)

run_wrapper <- function(sim_condition, family="binomial", prevalence=NULL) {
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
      
      # create rare case scenarios
      if(!is.null(prevalence)){
        # only keep prevalence % of 1
        idx1 <- sample(which(res$y==1), floor(length(res$y)*prevalence))
        idx0 <- which(res$y==0)
        rowidx <- sort(c(idx0, idx1))
        rownames(res$data) <- 1:nrow(res$data)
        
        res$y <- res$y[rowidx]
        res$c <- res$c[rowidx]
        res$data <- res$data[rowidx,]
        print(mean(res$y)) 
      }
      
      
        
      # glme mixed effect model
      m0 <- fit_eval_glmer(y = res$y,
                           c = res$c,
                           data = res$data,
                           family = family)
      
      
      # glm model 
      m1 <- fit_eval_glm(y = res$y,
                         c = res$c,
                         data = res$data,
                         family = family)
      t <- toc(log = TRUE)
      
      # gee model
      m2 <- fit_eval_gee(y = res$y,
                         c = res$c,
                         data = res$data,
                         family = family)
      t <- toc(log = TRUE)
      
      
      
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

