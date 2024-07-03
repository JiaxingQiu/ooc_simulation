source("./sim_run.R")

sjob_lm = slurm_map(
  split(simulation_conditions, simulation_conditions$id),
  run_wrapper_lm,
  nodes=nrow(simulation_conditions),
  cpus_per_node = 1,
  jobname = "ooc_run_lm",
  submit = TRUE,
  preschedule_cores = F,
  slurm_options =
    c(account = "netlab", partition = "standard", time = "2-00:00:00"), # standard
  global_objects = lsf.str()
)
save(sjob_lm, file = "ooc_run_lm.Rdata")

