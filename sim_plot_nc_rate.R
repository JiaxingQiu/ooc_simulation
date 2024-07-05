

source("./sim_conditions.R")
agg_df <- agg_df_ls[["lm"]] %>% filter(nc_rate0>0)
simulation_conditions[which(simulation_conditions$id%in%unique(agg_df$id)),]
# important:
# <simpleError: number of observations (=n_cluster * n_obs_per_cluster) <= number of random effects (= n_ttl_effect * n_cluster) for term (rdm1 + rdm2 + rdm3 + rdm4 + 1 | c); 
# the random-effects parameters and the residual variance (or scale parameter) are probably unidentifiable>


source("./sim_conditions.R")
agg_df <- agg_df_ls[["lr"]] %>% filter(nc_rate0>0)
simulation_conditions <- merge(simulation_conditions, distinct(agg_df[,c("id","nc_rate0")]) )
simulation_conditions[which(simulation_conditions$id%in%unique(agg_df$id)),]
summary(simulation_conditions$nc_rate0)
# under general conditions with 100 iterations per each, the RE model had a non-convergence rate of 0.53 [0.34 - 0.65] 

source("./sim_conditions.R")
agg_df <- agg_df_ls[["lr"]] %>% filter(nc_rate1>0)
simulation_conditions <- merge(simulation_conditions, distinct(agg_df[,c("id","nc_rate1")]) )
simulation_conditions[which(simulation_conditions$id%in%unique(agg_df$id)),]
summary(simulation_conditions$nc_rate1)

source("./sim_conditions.R")
agg_df <- agg_df_ls[["lr"]] %>% filter(nc_rate2>0)
simulation_conditions <- merge(simulation_conditions, distinct(agg_df[,c("id","nc_rate2")]) )
simulation_conditions[which(simulation_conditions$id%in%unique(agg_df$id)),]
summary(simulation_conditions$nc_rate2)




