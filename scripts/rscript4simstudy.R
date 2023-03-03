
# callr -------------------------------------------------------------------
callr::rscript( # https://callr.r-lib.org/reference/rscript.html
  script = "scripts/purrr-simstudy.R",
  cmdargs = character(),
  libpath = .libPaths()[1],
  #repos = default_repos(),
  stdout = paste0(Sys.Date(), "_slurm.out"),  # Optionally a file name to send the standard output to
  stderr = paste0(Sys.Date(), "_error.out"),
  #poll_connection = TRUE,
  echo = FALSE,  # default=FALSE
  show = TRUE,
  #callback = NULL,
  #block_callback = NULL,
  #spinner = FALSE,
  #system_profile = FALSE,
  #user_profile = "project",
  #env = rcmd_safe_env(),
  timeout = as.difftime(12, units = "hours"),  # default=Inf
  #wd = "~/bestageing2022",  # defaults to the current working directory.
  fail_on_status = TRUE,  # default=FALSE
  color = TRUE,
  #...
)
