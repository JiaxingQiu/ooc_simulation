
# basis condition
n_cluster <- 50 # number of clusters
n_obs_per_cluster <- 25 # number of observations per cluster
n_ttl_betas <- c(3,5,7) # number of total effects
fix_rdm_ratio <- c(0.2) # proportion of fix effects
sigma_fix <- c(5) # fix effect beta variance 
sigma_rdm_fix_ratio <- 10 # random effect beta variance proportional to fixed effects
ar1_phi <- 0.6 # runif min, max = min + 0.2, 0-0.2 means low within-cluster correlation, 0.4-0.6 median correlation, 0.8-1 high correlation
na_rate <- 0
param_grid <- expand.grid(n_cluster = n_cluster,
                          n_obs_per_cluster = n_obs_per_cluster,
                          n_ttl_betas = n_ttl_betas,
                          fix_rdm_ratio = fix_rdm_ratio,
                          sigma_fix = sigma_fix,
                          sigma_rdm_fix_ratio = sigma_rdm_fix_ratio,
                          ar1_phi = ar1_phi,
                          na_rate = na_rate)
simulation_conditions <- as.data.frame(param_grid)


# special scenarios
# 0: basis condition (to compare with)
simulation_conditions0 <- simulation_conditions
simulation_conditions0$case <- "raw"
# 1: unbalanced cluster size 
simulation_conditions1 <- simulation_conditions
simulation_conditions1$na_rate <- 1
simulation_conditions1$case <- "unbalanced"
# 2. small sample
simulation_conditions2 <- simulation_conditions
simulation_conditions2$n_cluster <- 10
simulation_conditions2$case <- "small"
# 3. rare event
simulation_conditions3 <- simulation_conditions
simulation_conditions3$case <- "rare"
simulation_conditions <- rbind(simulation_conditions0, simulation_conditions1, simulation_conditions2, simulation_conditions3)
rm(simulation_conditions0, simulation_conditions1, simulation_conditions2, simulation_conditions3)


simulation_conditions$id <- seq(1:nrow(simulation_conditions))
simulation_conditions$iter <- 100


rest_id <- c(2,3,12)
simulation_conditions <- simulation_conditions[which(simulation_conditions$id %in% rest_id),]


