
# Parameters
n_cluster <- 100 # number of clusters
n_obs_per_cluster <- c(10, 50, 100) # number of observations per cluster
n_ttl_betas <- seq(5, 10) # number of total effects
fix_rdm_ratio <- c(0.2) # proportion of fix effects
sigma_fix <- c(5) # fix effect beta variance 
sigma_rdm_fix_ratio <- rev(c(0.5, 1, 10))  # random effect beta variance proportional to fixed effects
ar1_phi <- rev(c(0, 0.4, 0.8)) # runif min, max = min + 0.2, 0-0.2 means low within-cluster correlation, 0.4-0.6 median correlation, 0.8-1 high correlation
na_rate <- c(0) 

# Define the simulation conditions
param_grid <- expand.grid(n_cluster = n_cluster,
                          n_obs_per_cluster = n_obs_per_cluster,
                          n_ttl_betas = n_ttl_betas,
                          fix_rdm_ratio = fix_rdm_ratio,
                          sigma_fix = sigma_fix,
                          sigma_rdm_fix_ratio = sigma_rdm_fix_ratio,
                          ar1_phi = ar1_phi,
                          na_rate = na_rate)
simulation_conditions <- as.data.frame(param_grid)
simulation_conditions$cluster_strength <- NA
simulation_conditions$cluster_strength[which(simulation_conditions$n_obs_per_cluster==10 &
                                               simulation_conditions$sigma_rdm_fix_ratio==0.5 &
                                               simulation_conditions$ar1_phi==0)] <- "weak"

simulation_conditions$cluster_strength[which(simulation_conditions$n_obs_per_cluster==50 &
                                               simulation_conditions$sigma_rdm_fix_ratio==1 &
                                               simulation_conditions$ar1_phi==0.4 )] <- "moderate"

simulation_conditions$cluster_strength[which(simulation_conditions$n_obs_per_cluster==100 &
                                               simulation_conditions$sigma_rdm_fix_ratio==10 &
                                               simulation_conditions$ar1_phi==0.8 )] <- "strong"
simulation_conditions <- simulation_conditions[which(!is.na(simulation_conditions$cluster_strength)),]

simulation_conditions$id <- seq(1:nrow(simulation_conditions))
simulation_conditions$iter <- 100

