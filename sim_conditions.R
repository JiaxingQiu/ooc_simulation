
# Parameters
n_cluster <- rev(c(10, 50)) # number of clusters
n_obs_per_cluster <- rev(c(2, 5, 10, 25, 50, 150)) # number of observations per cluster
n_ttl_betas <- seq(5, 10) # number of total effects
fix_rdm_ratio <- c(0.2) # c(0.2, 0.5, 0.8) # proportion of fix effects
sigma_fix <- c(5) # fix effect beta variance # 5
sigma_rdm_fix_ratio <- rev(c(0.5, 1, 10))  # random effect beta variance proportional to fixed effects
ar1_phi <- rev(c(0, 0.4, 0.8)) # runif min, max = min + 0.2, 0-0.2 means low within-cluster correlation, 0.4-0.6 median correlation, 0.8-1 high correlation
na_rate <- c(0, 1) # 0.3, 0.7,  cluster wise missingness, rnorm centers, 0.3 means low level missingness, 0.7 means high level missingness, 0 means no missing
# 1 means unbalanced missingness, ranging from 0-1

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
simulation_conditions$id <- seq(1:nrow(param_grid))
simulation_conditions$iter <- 100




