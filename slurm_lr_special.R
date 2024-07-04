source("./sim_run.R")
source("./sim_conditions_special.R")

sjob = slurm_map(
  split(simulation_conditions, simulation_conditions$id),
  run_wrapper_lr,
  nodes=nrow(simulation_conditions),
  cpus_per_node = 1,
  jobname = "ooc_run_lr_special",
  submit = TRUE,
  preschedule_cores = F,
  slurm_options =
    c(account = "netlab", partition = "standard", time = "5-00:00:00"), 
  global_objects = lsf.str()
)
save(sjob, file = "ooc_run_lr_special.Rdata")
