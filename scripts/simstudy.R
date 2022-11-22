

## ---------------------------
##
## Script name: Simulation of betas and sample sizes
##
## Author: Christoph Reich
## Date Created: 2022-11-22
##
## Copyright (c) Christoph Reich, 2022
## Email: christoph.reich@med.uni-heidelberg.de
##
## ---------------------------
## Notes:
##   
## ---------------------------


# set to TRUE if simulation should run when script is sourced
runSimulation <- FALSE

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
source("scripts/MR_lasso.R")  # Slob&Burgess 2020
source("scripts/PVE_calculation_fns.R")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")



# GET DATA ----------------------------------------------------------------
# dataa
my_data_harm <- readRDS(file = "./output/Rdata/2022-11-22_my_data_harm.rds")
# MR estimates from full study
est <- readRDS(file = "./output/Rdata/2022-11-22_est.rds")
# vars from MSc script
sim_vars <- readRDS(file = "./output/Rdata/2022-11-22_sim_vars.rds")
n <- sim_vars$n
colnames_pval <- sim_vars$colnames_pval

###
# SIMULATION --------------------------------------------------------------
###
set.seed(1234)
## initialize
nsim <- 25
est_sim <- tibble(no_sim = 1:nsim)

## START SIM
if (runSimulation==TRUE) {
  T0 = proc.time()[3]
  for (i in 1:length(colnames_pval)) {
    threshold <- 5*10^-6
    se_update_BMI <- c()
    se_update_GBC <- c()
    
    numIV <- sum(my_data_harm[[colnames_pval[i]]] < threshold)
    
    if (numIV<5) { # start with 5 IVs only (sample size ~50,000)
      next
    }
    
    # FILTER DATA
    sim_dat <- my_data_harm %>% 
      filter(!!sym(colnames_pval[i]) < threshold) %>% 
      select(SNP, eaf.bmi, beta.bmi, p.value.bmi, eaf.gbc, beta.gbc, p.value.gbc)
    
    # get simulated se_{sample size}
    for (j in 1:nrow(sim_dat)) {
      # as done before in section "simulating on empirical data"
      se_update_BMI[j] <- abs(sim_dat[["beta.bmi"]][j]/-qnorm(p=(sim_dat[["p.value.bmi"]][i]/2),mean=0)*sqrt(n.orig)/sqrt(n[i]))
      se_update_GBC[j] <- abs(sim_dat[["beta.gbc"]][j]/-qnorm(p=(sim_dat[["p.value.gbc"]][i]/2),mean=0)*sqrt(n.orig)/sqrt(n[i]))
    }
    
    sim_dat <- sim_dat %>%
      bind_cols(tibble(se_update_BMI=se_update_BMI, se_update_GBC=se_update_GBC))
    
    # simulate beta values from normal_distribution # -----------------------------------------------------
    for (k in 1:nsim) {
      sim_dat <- sim_dat %>% 
        mutate(beta.bmi_sim_norm = purrr::map2_dbl(.x = beta.bmi, .y = se_update_BMI, .f = ~rnorm(n = 1, mean = .x, sd = .y))) %>% 
        mutate(beta.gbc_sim_norm = purrr::map2_dbl(.x = beta.gbc, .y = se_update_GBC, .f = ~rnorm(n = 1, mean = .x, sd = .y))) %>% 
        mutate(#explained_variance_X_sim = explained_variance_numeric(eaf = eaf.bmi, beta = beta.bmi), 
          explained_variance_X2_sim = explained_variance_numeric2(maf=eaf.bmi, 
                                                                  beta=beta.bmi_sim_norm, 
                                                                  se_beta=se_update_BMI, 
                                                                  samplesize=n[i]
          ),
          # explained_variance_GSD = explained_variance_binary(PA = eaf.gsd, 
          #                                                    RR1 = exp(beta.gsd), 
          #                                                    RR2 = exp(beta.gsd)^2, 
          #                                                    K = 0.1  # assumed prevalence in population
          #                                                    )$Vg,
          explained_variance_Y_sim = explained_variance_binary(PA = eaf.gbc, 
                                                               RR1 = exp(beta.gbc_sim_norm), 
                                                               RR2 = exp(beta.gbc_sim_norm)^2, 
                                                               K =1.2/100000  # assumed prevalence in population
          )$Vg
        )
      
      # PERFORM MR_analysis
      mr.obj = MendelianRandomization::mr_input(
        bx = sim_dat$beta.bmi_sim_norm, 
        bxse = sim_dat$se_update_BMI, 
        # outcome
        by = sim_dat$beta.gbc_sim_norm, 
        byse = sim_dat$se_update_GBC, 
        snps = sim_dat$SNP,
        exposure = "Body mass index",
        outcome = "Gallbladder Cancer"
      ) 
      # CONMIX (as reference)
      res = MendelianRandomization::mr_conmix(mr.obj)
      colname_theta <- paste0("theta_sim_", format(n[i], scientific=FALSE))
      colname_theta_se <- paste0("theta_se_sim_", format(n[i], scientific=FALSE))
      colname_pve_X <- paste0("pve_X_sim_", format(n[i], scientific=FALSE))
      colname_pve_Y <- paste0("pve_Y_sim_", format(n[i], scientific=FALSE))
      
      ## store results
      est_sim[k, colname_theta] <- res$Estimate
      CIlength = res$CIUpper-res$CILower
      if (length(CIlength)>1) print("conmix multimodal")
      est_sim[k, colname_theta_se] <- sum(CIlength)/1.96/2  ## Caution: this may be problematic 
      # add PVE cols
      est_sim[k, colname_pve_X] <- sum(sim_dat$explained_variance_X2_sim)
      est_sim[k, colname_pve_Y] <- sum(sim_dat$explained_variance_Y_sim)
      
      rm(res)
    }
    rm(se_update_BMI, se_update_GBC)
    print(paste0("|||-----------------------Run finished for sample size: ", format(n[i], scientific=FALSE), " -----------------------|||"))
  }  ## END OF SIMULATION
  
  
  T1 = proc.time()[3]
  timediff = T1-T0
  
  saveRDS(est_sim, paste0("../output/Rdata/", Sys.Date(), "_est_sim.rds"))
}


