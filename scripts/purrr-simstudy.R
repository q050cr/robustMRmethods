

## ---------------------------
##
## Script name: purrr-simstudy.R
##
## Author: Christoph Reich
## Date Created: 2022-11-22
##
## Copyright (c) Christoph Reich, 2022
## Email: christoph.reich@med.uni-heidelberg.de
##
## ---------------------------
## Notes:
##   creates df `est_sim` with nsim* theta_estimates, nsim* standard errors of theta sims, nsim* pve for exposure
##      cols: "no_sim", "theta_sim_SAMPLE_SIZES", "theta_se_sim_SAMPLE_SIZES", "pve_X_sim_SAMPLE_SIZES"
##
##   implemented methods in simulation so far:
##          1. random IVW
##          2. weighted mode
##          3. CONMIX
##          4. MRMix
##
## ---------------------------

# set to TRUE if simulation should run when script is sourced
rm(list=ls())
runSimulation <- TRUE
nsim <- 5

library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(tibble)
library(MendelianRandomization)
library(TwoSampleMR)
library(mr.raps)
library(MRMix)
library(MRPRESSO)
library(penalized)
library(conflicted)
source("scripts/helper/MR_lasso.R") # Slob&Burgess 2020
source("scripts/helper/PVE_calculation_fns.R")
source("scripts/helper/simulation_helpers.R")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# GET DATA ----------------------------------------------------------------
# data
my_data_harm <- readRDS(file = dplyr::last(list.files("./output/Rdata/", pattern = "my_data_harm.rds", full.names = TRUE)))
# vars from MSc script
sim_vars <- readRDS(file = dplyr::last(list.files("./output/Rdata/", pattern = "_sim_vars.rds", full.names = TRUE)))
n <- sim_vars$n
colnames_pval <- sim_vars$colnames_pval
n.orig <- 334487

## PARAMS for simulation
set.seed(1234)
threshold <- 5 * 10^-6
n.orig_sqrt <- sqrt(n.orig)

# space to populate
est_sim_temp <- tibble(no_sim = 1:nsim)
est_sim <- tibble(no_sim = 1:nsim)

# create dataframes for each specified sample size ------------------------
gen_dfs <- map(1:length(colnames_pval), ~{
  threshold <- 5*10^-6
  numIV <- sum(my_data_harm[[colnames_pval[.x]]] < threshold)
  if (numIV<5) {
    return(NA)
  }
  sim_dat <- my_data_harm %>% 
    filter(!!sym(colnames_pval[.x]) < threshold) %>% 
    select(SNP, eaf.bmi, beta.bmi, se.bmi, p.value.bmi, eaf.gsd, beta.gsd, se.gsd, p.value.gsd) %>% 
    mutate(se_update_BMI = abs(beta.bmi / -qnorm(p=(p.value.bmi / 2), mean = 0) * sqrt(n.orig) / sqrt(n[.x])))
})

# run simulation ----------------------------------------------------------
T0 <- proc.time()[3]
for (dataset in 1:length(gen_dfs)) {
  nsim_list_dataset <- map(1:nsim, ~sim_function(dat = gen_dfs[[dataset]], index = dataset))  # returns nsim* results for MR simulations PER sample size 
  est_sim_temp <- populate_est_sim(nsim_list = nsim_list_dataset, nsim = nsim, no_cases = n[dataset])
  if(is.null(est_sim_temp)) next
  est_sim <- est_sim %>% cbind(est_sim_temp[-1])
  rm(est_sim_temp)
  #
  print(paste0("|||-----------------------Run finished for sample size: ", format(n[dataset], scientific = FALSE), " -----------------------|||"))
}
T1 <- proc.time()[3]
timediff_purrr <- T1 - T0

saveRDS(est_sim, paste0("./output/Rdata/", Sys.Date(), "_est_sim_purrr_NSIM_", nsim, ".rds"))

