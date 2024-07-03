source("./sim_run.R")

sjob_lr = slurm_map(
  split(simulation_conditions, simulation_conditions$id),
  run_wrapper_lr,
  nodes=nrow(simulation_conditions),
  cpus_per_node = 1,
  submit = TRUE,
  preschedule_cores = F,
  slurm_options =
    c(account = "netlab", partition = "standard", time = "3-00:00:00"), # standard
  global_objects = lsf.str()
)
save(sjob_lr, file = "ooc_run_lr.Rdata")



