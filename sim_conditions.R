
# Parameters
n_cluster <- 50 # number of clusters
n_obs_per_cluster <- rev(c(5, 25, 100)) # number of observations per cluster
n_ttl_betas <- c(3,5,7) # number of total effects
fix_rdm_ratio <- c(0.2) # c(0.2, 0.5, 0.8) # proportion of fix effects
sigma_fix <- c(5) # fix effect beta variance # 5
sigma_rdm_fix_ratio <- rev(c(0.5, 1, 10))  # random effect beta variance proportional to fixed effects
ar1_phi <- rev(c(0, 0.3, 0.6)) # runif min, max = min + 0.2, 0-0.2 means low within-cluster correlation, 0.3-0.5 median correlation, 0.6-0.8 high correlation
na_rate <- c(0) # 0.3, 0.7,  cluster wise missingness, rnorm centers, 0.3 means low level missingness, 0.7 means high level missingness, 0 means no missing
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
# # remove conditions where sigma_rdm_fix_ratio>0.5 and ar1_phi>0
# simulation_conditions <- simulation_conditions[which(!(simulation_conditions$sigma_rdm_fix_ratio>0.5&simulation_conditions$ar1_phi>0)),]
# # sanity check
# table(simulation_conditions[,c("sigma_rdm_fix_ratio", "ar1_phi")])
simulation_conditions$id <- seq(1:nrow(simulation_conditions))
simulation_conditions$iter <- 100

# special scenario can be analyzed under fewer clustering conditions
# small sample 10
# unbalanced 
# rare event

