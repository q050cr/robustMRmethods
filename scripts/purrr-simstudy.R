

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
##
## ---------------------------

rm(list=ls())
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
source("scripts/helper/PVE_calculation_fns.R")
source("scripts/helper/simulation_helpers.R")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# set to TRUE if simulation should run when script is sourced
runSimulation <- TRUE
subsetSampleSizes <- TRUE
runMRpresso <- TRUE  # runs foreeeever ;)
nsim <- 500
runGeneticArchitectures <- TRUE  # if set to TRUE, simulations run on 4 different genetic architectures (split: median maf, median beta)

# analysis set - create all combis
effect_size <- c("strong", "weak")
maf <- c("rare", "common")
all_quadrants <- tidyr::crossing(effect_size, maf)

# GET DATA ----------------------------------------------------------------
# data
my_data_harm <- readRDS(file = dplyr::last(list.files("./output/Rdata/harmonized-dat/", pattern = "my_data_harm.rds", full.names = TRUE)))
# vars from MSc script
sim_vars <- readRDS(file = dplyr::last(list.files("./output/Rdata/", pattern = "_sim_vars.rds", full.names = TRUE)))

# only simulate defined set of smaller sample sizes?
n <- sim_vars$n
colnames_pval <- sim_vars$colnames_pval
n.orig <- 334487

## subset?
subset.index <- c(8, 13, 18, 23, 33)
if (subsetSampleSizes == TRUE) {
  n <- as.integer(n[subset.index])
  colnames_pval <- colnames_pval[subset.index]
}


## PARAMS for simulation
set.seed(1234)
threshold <- 5 * 10^-6
n.orig_sqrt <- sqrt(n.orig)

# space to populate
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
    mutate(se_update_BMI = abs(beta.bmi / -qnorm(p=(p.value.bmi / 2), mean = 0) * sqrt(n.orig) / sqrt(n[.x]))) %>% 
    mutate(maf.bmi = ifelse(eaf.bmi >0.5, 1-eaf.bmi, eaf.bmi))
})

# run simulation ----------------------------------------------------------
T0 <- proc.time()[3]
if (runSimulation == TRUE) {
  for (dataset in 1:length(gen_dfs)) {
    # initialize in loop/ discard at the end
    est_sim_temp <- tibble(no_sim = 1:nsim)
    
    ## now purrr
    # returns nsim* results for MR simulations PER sample size
    nsim_list_dataset <- map(1:nsim, ~sim_function(dat = gen_dfs[[dataset]], index = dataset, runMRpresso = runMRpresso, subsetSampleSizes=subsetSampleSizes),      # RUN MR PRESSO???
                             .progress=list(clear=FALSE, 
                                            name=paste0("Progress of simulation of sample size: ", format(n[dataset], scientific = FALSE)))
    )
    est_sim_temp <- populate_est_sim(nsim_list = nsim_list_dataset, nsim = nsim, no_cases = n[dataset])
    if(is.null(est_sim_temp)) next
    est_sim <- est_sim %>% cbind(est_sim_temp[-1])
    rm(est_sim_temp, nsim_list_dataset)
    #
    print(paste0("|||-----------------------Run finished for sample size: ", format(n[dataset], scientific = FALSE), " -----------------------|||"))
  }
}
T1 <- proc.time()[3]
timediff_purrr <- T1 - T0
print(paste0("!!!!!--------------------------- ", round(timediff_purrr/60,2), " minutes passed for the simulation scenario! --------------------------------------!!!!!"))

saveRDS(est_sim, paste0("./output/Rdata/sim-results/", Sys.Date(), "_est_sim_purrr_NSIM_", nsim, ".rds"))

## with nsim=100 and the purrr script on the M1 Mac ∞ 5h (only 4 methods at that time)
#[1] "!!!!!--------------------------- 299.61 minutes passed for the simulation scenario! --------------------------------------!!!!!"


safe_filter <- safely(dplyr::filter, otherwise = NA, quiet = TRUE)  # needed because also NAs in `gen_dfs`
est_sim_architectures <- list()
T0_architecture <- proc.time()[3]
if (runGeneticArchitectures == TRUE) {
  for (architecture in 1:nrow(all_quadrants)) {
    effect_size <- all_quadrants[[architecture, 1]]
    maf <- all_quadrants[[architecture, 2]]
    
    ## filter for quadrant
    if (effect_size=="strong" & maf == "rare") {
      gen_dfs_filter <- map(.x = gen_dfs, .f = ~safe_filter(.x, abs(beta.bmi) > median(abs(beta.bmi)) &  maf.bmi < median(maf.bmi))$result)
    } else if (effect_size=="strong" & maf == "common") {
      gen_dfs_filter <- map(.x = gen_dfs, .f = ~safe_filter(.x, abs(beta.bmi) > median(abs(beta.bmi)) &  maf.bmi > median(maf.bmi))$result)
    } else if (effect_size=="weak" & maf == "common") {
      gen_dfs_filter <- map(.x = gen_dfs, .f = ~safe_filter(.x, abs(beta.bmi) < median(abs(beta.bmi)) &  maf.bmi > median(maf.bmi))$result)
    } else if (effect_size=="weak" & maf == "rare") {
      gen_dfs_filter <- map(.x = gen_dfs, .f = ~safe_filter(.x, abs(beta.bmi) < median(abs(beta.bmi)) &  maf.bmi < median(maf.bmi))$result)
    }
    ## initialize to store for later
    est_sim_naming <- paste0(toupper(effect_size), "_", toupper(maf),"_est_sim")
    est_sim_architectures[[est_sim_naming]] <- est_sim <- tibble(no_sim = 1:nsim)
    
    for (dataset in 1:length(gen_dfs_filter)) {
      # initialize in loop/ discard at the end
      est_sim_temp <- tibble(no_sim = 1:nsim)
      
      ## now purrr
      # returns nsim* results for MR simulations PER sample size
      nsim_list_dataset <- map(1:nsim, ~sim_function(dat = gen_dfs_filter[[dataset]], index = dataset, runMRpresso = runMRpresso, subsetSampleSizes=subsetSampleSizes),      # RUN MR PRESSO???
                               .progress=list(clear=FALSE, 
                                              name=paste0("Progress of simulation of sample size: ", format(n[dataset], scientific = FALSE)))
      )
      est_sim_temp <- populate_est_sim(nsim_list = nsim_list_dataset, nsim = nsim, no_cases = n[dataset])
      if(is.null(est_sim_temp)) next
      est_sim_architectures[[est_sim_naming]] <- est_sim_architectures[[est_sim_naming]] %>% cbind(est_sim_temp[-1])
      rm(est_sim_temp, nsim_list_dataset)
      #
      print(paste0("|||-----------------------Run finished for sample size: ", format(n[dataset], scientific = FALSE), " -----------------------|||"))
    }
  }
}

if (runGeneticArchitectures == TRUE) {
  T1_architecture <- proc.time()[3]
  timediff_purrr_archi <- T1_architecture - T0_architecture
  print(paste0("!!!!!--------------------------- ", round(timediff_purrr_archi/60,2), " minutes passed for the simulation scenario! --------------------------------------!!!!!"))
  
  saveRDS(est_sim_architectures, paste0("./output/Rdata/sim-results/", Sys.Date(), "_est_sim_architectures_purrr_NSIM_", nsim, ".rds"))
}

