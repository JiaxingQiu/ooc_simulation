source("./sim_run.R")
source("./sim_conditions_rest.R")


sjob_lr = slurm_map(
  split(simulation_conditions, simulation_conditions$id),
  run_wrapper_lr,
  nodes=nrow(simulation_conditions),
  cpus_per_node = 1,
  jobname = "ooc_run_lr_rest",
  submit = TRUE,
  preschedule_cores = F,
  slurm_options =
    c(account = "netlab", partition = "standard", time = "8-00:00:00"), 
  global_objects = lsf.str()
)
save(sjob_lr, file = "ooc_run_lr_rest.Rdata")

