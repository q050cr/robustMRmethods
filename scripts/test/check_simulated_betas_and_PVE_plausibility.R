


## simcheck CORRELATION simulated values




## ---------------------------
##
## Script name: simstudy.R
##
## Author: Christoph Reich
## Date Created: 2022-11-22
##
## Copyright (c) Christoph Reich, 2022
## Email: christoph.reich@med.uni-heidelberg.de
##
## ---------------------------
## Notes:
##   creates df `est_sim` with nsim* theta_estimates, nsim* standard errors of theta sims, nsim* pve for exposure and outcome
##      cols: "no_sim", "theta_sim_SAMPLE_SIZES", "theta_se_sim_SAMPLE_SIZES", "pve_X_sim_SAMPLE_SIZES", "pve_Y_sim_SAMPLE_SIZES"
## ---------------------------


# set to TRUE if simulation should run when script is sourced
runSimulation <- TRUE

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
source("scripts/MR_lasso.R") # Slob&Burgess 2020
source("scripts/PVE_calculation_fns.R")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")



# GET DATA ----------------------------------------------------------------
# data
my_data_harm <- readRDS(file = dplyr::last(list.files("./output/Rdata/", pattern = "my_data_harm.rds", full.names = TRUE)))
# MR estimates from full study
est <- readRDS(file = dplyr::last(list.files("./output/Rdata/", pattern = "_est.rds", full.names = TRUE)))
# vars from MSc script
sim_vars <- readRDS(file = dplyr::last(list.files("./output/Rdata/", pattern = "_sim_vars.rds", full.names = TRUE)))
n <- sim_vars$n
colnames_pval <- sim_vars$colnames_pval
n.orig <- 334487

###
# SIMULATION --------------------------------------------------------------
###
set.seed(1234)
## initialize
nsim <- 100
est_sim <- tibble(no_sim = 1:nsim)


## START SIM

sink(file = "scripts/test/output/simPVE_sd_observed_with_samplesize.txt")
if (runSimulation == TRUE) {
  T0 <- proc.time()[3]
  for (i in 1:length(colnames_pval)) {
    threshold <- 5 * 10^-6
    se_update_BMI <- c()
    se_update_GBC <- c()
    
    numIV <- sum(my_data_harm[[colnames_pval[i]]] < threshold)
    
    if (numIV < 5) { # start with 5 IVs only (sample size ~50,000)
      next
    }
    
    # FILTER DATA
    sim_dat <- my_data_harm %>%
      filter(!!sym(colnames_pval[i]) < threshold) %>%
      select(SNP, eaf.bmi, beta.bmi, se.bmi, p.value.bmi, eaf.gbc, beta.gbc, se.gbc, p.value.gbc)
    
    # get simulated se_{sample size}
    for (j in 1:nrow(sim_dat)) {
     # as done before in section "simulating on empirical data"
     se_update_BMI[j] <- abs(sim_dat[["beta.bmi"]][j] / -qnorm(p = (sim_dat[["p.value.bmi"]][i] / 2), mean = 0) * sqrt(n.orig) / sqrt(n[i]))
     se_update_GBC[j] <- abs(sim_dat[["beta.gbc"]][j] / -qnorm(p = (sim_dat[["p.value.gbc"]][i] / 2), mean = 0) * sqrt(n.orig) / sqrt(n[i]))
    }
    
    sim_dat <- sim_dat %>%
     bind_cols(tibble(se_update_BMI = se_update_BMI, se_update_GBC = se_update_GBC))
    
    # simulate beta values from normal_distribution # -----------------------------------------------------
    simulated_betas_df_X <- sim_dat %>%
      select(SNP, beta.bmi)
    simulated_betas_df_Y <- sim_dat %>%
      select(SNP, beta.gbc)
    
    for (k in 1:nsim) {
      sim_dat <- sim_dat %>% 
        mutate(beta.bmi_sim_norm = purrr::map2_dbl(.x = beta.bmi, .y = se_update_BMI, .f = ~rnorm(n = 1, mean = .x, sd = .y))) %>% 
        mutate(beta.gbc_sim_norm = purrr::map2_dbl(.x = beta.gbc, .y = se_update_GBC, .f = ~rnorm(n = 1, mean = .x, sd = .y))) %>% 
        mutate(#explained_variance_X_sim = explained_variance_numeric(eaf = eaf.bmi, beta = beta.bmi), 
          explained_variance_X2_sim = explained_variance_numeric2(maf=eaf.bmi, 
                                                                  beta=beta.bmi_sim_norm, 
                                                                  se_beta=se.bmi, 
                                                                  samplesize=n[length(n)]  # n[i]
          ),
          explained_variance_Y_sim = explained_variance_binary(PA = eaf.gbc, 
                                                               RR1 = exp(beta.gbc_sim_norm), 
                                                               RR2 = 1+(2*(exp(beta.gbc_sim_norm)-1)),  #exp(beta.gbc_sim_norm)^2, 
                                                               K =1.2/100000  # assumed prevalence in population
          )$Vg
        )
      
      # sum_PVE_Y_observed <- explained_variance_binary(
      #   PA = my_data_harm$eaf.gbc[my_data_harm[[colnames_pval[i]]] < 5*10^(-6)], 
      #   RR1 = exp(my_data_harm$beta.gbc[my_data_harm[[colnames_pval[i]]] < 5*10^(-6)]), 
      #   RR2 = 1+(2*(exp(my_data_harm$beta.gbc[my_data_harm[[colnames_pval[i]]] < 5*10^(-6)])-1)), 
      #   K =1.2/100000  # assumed prevalence in population
      # )$Vg %>% sum
      # 
      # PVE_Y <- explained_variance_binary(
      #   PA = sim_dat$eaf.gbc, 
      #   RR1 = exp(sim_dat$beta.gbc_sim_norm), 
      #   RR2 = 1+(2*(exp(sim_dat$beta.gbc_sim_norm)-1)), 
      #   K =1.2/100000  # assumed prevalence in population
      # )$Vg %>% sum
      # PVE_Y
      
      # store simulated beta values
      colname <- paste0("beta.bmi_sim", k)
      colnameY <- paste0("beta.gbc_sim", k)
      colname_pve_X <- paste0("pve_X_sim_", format(n[i], scientific=FALSE))
      colname_pve_Y <- paste0("pve_Y_sim_", format(n[i], scientific=FALSE))
      
      simulated_betas_df_X <- cbind(simulated_betas_df_X, sim_dat["beta.bmi_sim_norm"]) %>%
        rename(!!colname := beta.bmi_sim_norm) %>%
        as_tibble()
      
      simulated_betas_df_Y <- cbind(simulated_betas_df_Y, sim_dat["beta.gbc_sim_norm"]) %>%
        rename(!!colnameY := beta.gbc_sim_norm) %>%
        as_tibble()
      
      est_sim[k, colname_pve_X] <- sum(sim_dat$explained_variance_X2_sim, na.rm = TRUE)
      est_sim[k, colname_pve_Y] <- sum(sim_dat$explained_variance_Y_sim, na.rm = TRUE)
    }
    
    base::apply(simulated_betas_df_X %>% select(-c(SNP, beta.bmi)), 1, mean) -> mean_sim_betasX
    base::apply(simulated_betas_df_Y %>% select(-c(SNP, beta.gbc)), 1, mean) -> mean_sim_betasY
    
    simulated_betas_df_X <- simulated_betas_df_X %>% cbind(mean_sim_betasX)
    simulated_betas_df_Y <- simulated_betas_df_Y %>% cbind(mean_sim_betasY)
    
    simulated_betas_df_X %>% select(beta.bmi, mean_sim_betasX) -> checkX  # nearly perfect
    simulated_betas_df_Y %>% select(beta.gbc, mean_sim_betasY) -> checkY
    
    print(paste0("Correlation between observed beta_BMI & mean of ", nsim, " simulated beta_BMI: ", cor(checkX$beta.bmi,checkX$mean_sim_betasX)))
    print(paste0("Correlation between observed beta_GBC & mean of ", nsim, " simulated beta_GBC: ", cor(checkY$beta.gbc,checkY$mean_sim_betasY)))
    
    mean_PVE_X <- mean(est_sim[[colname_pve_X]], na.rm=TRUE)
    mean_PVE_Y <- mean(est_sim[[colname_pve_Y]], na.rm=TRUE)
    
    # vs observed PVE
    sum_PVE_X_observed <- explained_variance_numeric2(
      maf=my_data_harm$eaf.bmi[my_data_harm[[colnames_pval[i]]] < 5*10^(-6)], 
      beta=my_data_harm$beta.bmi[my_data_harm[[colnames_pval[i]]] < 5*10^(-6)], 
      se_beta=my_data_harm$se.bmi[my_data_harm[[colnames_pval[i]]] < 5*10^(-6)], 
      samplesize=n[length(n)]) %>% sum
    
    sum_PVE_Y_observed <- explained_variance_binary(
      PA = my_data_harm$eaf.gbc[my_data_harm[[colnames_pval[i]]] < 5*10^(-6)], 
      RR1 = exp(my_data_harm$beta.gbc[my_data_harm[[colnames_pval[i]]] < 5*10^(-6)]), 
      RR2 = 1+(2*(exp(my_data_harm$beta.gbc[my_data_harm[[colnames_pval[i]]] < 5*10^(-6)])-1)), 
      K =1.2/100000  # assumed prevalence in population
    )$Vg %>% sum
    
    print(paste0("Mean PVE_X with ", nsim, " simulated beta_BMIs: ", round(mean_PVE_X*100, 3), "% (observed: ", round(sum_PVE_X_observed*100, 3), "%)" ))
    print(paste0("Mean PVE_Y with ", nsim, " simulated beta_GBCs: ", round(mean_PVE_Y*100, 3), "% (observed: ", round(sum_PVE_Y_observed*100, 3), "%)" ))
    
    print(paste0("|||-----------------------Run finished for sample size: ", format(n[i], scientific=FALSE), " -----------------------|||"))
  } ## END OF SIMULATION
  
}
sink(file = NULL)

