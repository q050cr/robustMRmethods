
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
source("scripts/MR_lasso.R")  # Slob&Burgess 2020
source("scripts/PVE_calculation_fns.R")
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

###
# SIMULATION --------------------------------------------------------------
###
set.seed(1234)
## initialize
nsim <- 50
est_sim <- tibble(no_sim = 1:nsim)

## START SIM
if (runSimulation==TRUE) {
  T0 = proc.time()[3]
  for (i in 1:length(colnames_pval)) {
    threshold <- 5*10^-6
    se_update_BMI <- c()
    se_update_GSD <- c()
    se_update_GBC <- c()
    
    numIV <- sum(my_data_harm[[colnames_pval[i]]] < threshold)
    
    if (numIV<5) { # start with 5 IVs only (sample size ~50,000)
      next
    }
    
    # FILTER DATA
    sim_dat <- my_data_harm %>% 
      filter(!!sym(colnames_pval[i]) < threshold) %>% 
      select(SNP, 
             eaf.bmi, beta.bmi, se.bmi, p.value.bmi, 
             eaf.gsd, beta.gsd, se.gsd, p.value.gsd,
             eaf.gbc, beta.gbc, se.gbc ,p.value.gbc)
    
    # get simulated se_{sample size} --> this yields large standard errors for small sample sizes -> use original sd
    # for (j in 1:nrow(sim_dat)) {
    #   # as done before in section "simulating on empirical data"
    #   se_update_BMI[j] <- abs(sim_dat[["beta.bmi"]][j]/-qnorm(p=(sim_dat[["p.value.bmi"]][j]/2),mean=0)*sqrt(n.orig)/sqrt(n[i]))
    #   se_update_GSD[j] <- abs(sim_dat[["beta.gsd"]][j]/-qnorm(p=(sim_dat[["p.value.gsd"]][j]/2),mean=0)*sqrt(n.orig)/sqrt(n[i]))
    #   se_update_GBC[j] <- abs(sim_dat[["beta.gbc"]][j]/-qnorm(p=(sim_dat[["p.value.gbc"]][j]/2),mean=0)*sqrt(n.orig)/sqrt(n[i]))
    # }
    #
    # sim_dat <- sim_dat %>%
    #   bind_cols(tibble(se_update_BMI=se_update_BMI, se_update_GBC=se_update_GBC))
    
    # simulate beta values from normal_distribution # -----------------------------------------------------
    for (k in 1:nsim) {
      # scaling se_Y so that max(se_Y) = max(se_X)
      factor1 <- max(sim_dat$se.bmi)/max(sim_dat$se.gbc)
      factor2 <- max(sim_dat$se.bmi)/max(sim_dat$se.gsd)
      se.gbc_narrow <- sim_dat$se.gbc * factor1
      se.gsd_narrow <- sim_dat$se.gsd * factor2
      sim_dat <- sim_dat %>% 
        mutate(beta.bmi_sim_norm = purrr::map2_dbl(.x = beta.bmi, .y = se.bmi, .f = ~rnorm(n = 1, mean = .x, sd = .y))) %>% 
        # here we need to deal with standard deviation (standard errors for Y are  huge/ see mail 2022-12-16)
        mutate(beta.gsd_sim_norm = purrr::map2_dbl(.x = beta.gsd, .y = se.gsd_narrow, .f = ~rnorm(n = 1, mean = .x, sd = .y))) %>% 
        mutate(beta.gbc_sim_norm = purrr::map2_dbl(.x = beta.gbc, .y = se.gbc_narrow, .f = ~rnorm(n = 1, mean = .x, sd = .y))) %>% 
        mutate(#explained_variance_X_sim = explained_variance_numeric(eaf = eaf.bmi, beta = beta.bmi), 
          explained_variance_X2_sim = explained_variance_numeric2(maf=eaf.bmi, 
                                                                  beta=beta.bmi_sim_norm, 
                                                                  se_beta=se.bmi, 
                                                                  samplesize=n[length(n)]  # n[i]
          ),
          explained_variance_U_sim = explained_variance_binary(PA = eaf.gsd, 
                                                             RR1 = exp(beta.gsd_sim_norm), 
                                                             RR2 = 1+(2*(exp(beta.gsd_sim_norm)-1)), 
                                                             K = 0.1  # assumed prevalence in population
          )$Vg,
          explained_variance_Y_sim = explained_variance_binary(PA = eaf.gbc, 
                                                               RR1 = exp(beta.gbc_sim_norm), 
                                                               RR2 = 1+(2*(exp(beta.gbc_sim_norm)-1)),  #exp(beta.gbc_sim_norm)^2, 
                                                               K =1.2/100000  # assumed prevalence in population
          )$Vg
        )
      
      # PVE 
      colname_pve_X <- paste0("pve_X_sim_", format(n[i], scientific=FALSE))
      colname_pve_U <- paste0("pve_U_sim_", format(n[i], scientific=FALSE))
      colname_pve_Y <- paste0("pve_Y_sim_", format(n[i], scientific=FALSE))
      ## add PVE cols
      est_sim[k, colname_pve_X] <- sum(sim_dat$explained_variance_X2_sim)
      est_sim[k, colname_pve_U] <- sum(sim_dat$explained_variance_U_sim)
      est_sim[k, colname_pve_Y] <- sum(sim_dat$explained_variance_Y_sim)
      
      # PERFORM MR_analysis 
      ## A/ GSD  ---------------------------------------------------------------
      mr.obj.gsd = MendelianRandomization::mr_input(
        bx = sim_dat$beta.bmi_sim_norm, 
        bxse = sim_dat$se.bmi, 
        # outcome
        by = sim_dat$beta.gsd_sim_norm, 
        byse = sim_dat$se.gsd, 
        snps = sim_dat$SNP,
        exposure = "Body mass index",
        outcome = "Gallstone Disease"
      ) 
      
      # 1. random IVW
      ivw.res.gsd = MendelianRandomization::mr_ivw(mr.obj.gsd, model = "random")
      # colnames
      ivw.colname_theta_gsd <- paste0("IVW_theta_sim_BMI_GSD_", format(n[i], scientific=FALSE))
      ivw.colname_theta_se_gsd <- paste0("IVW_theta_se_sim_BMI_GSD_", format(n[i], scientific=FALSE))
      ivw.colname_pvalue_gsd <- paste0("IVW_mr_pval_BMI_GSD_", format(n[i], scientific=FALSE))
      # store results
      est_sim[k, ivw.colname_theta_gsd] <- ivw.res.gsd$Estimate
      est_sim[k, ivw.colname_theta_se_gsd] <- ivw.res.gsd$StdError
      est_sim[k, ivw.colname_pvalue_gsd] <- ivw.res.gsd$Pvalue
      
      # 2. weighted mode
      mode.res.gsd = MendelianRandomization::mr_mbe(mr.obj.gsd)
      # colnames
      mode.colname_theta_gsd <- paste0("MODE_theta_sim_BMI_GSD_", format(n[i], scientific=FALSE))
      mode.colname_theta_se_gsd <- paste0("MODE_theta_se_sim_BMI_GSD_", format(n[i], scientific=FALSE))
      mode.colname_pvalue_gsd <- paste0("MODE_mr_pval_BMI_GSD_", format(n[i], scientific=FALSE))
      # store results
      est_sim[k, mode.colname_theta_gsd] <- mode.res.gsd$Estimate
      est_sim[k, mode.colname_theta_se_gsd] <- mode.res.gsd$StdError
      est_sim[k, mode.colname_pvalue_gsd] <- mode.res.gsd$Pvalue
      
      # 3. CONMIX 
      conmix.res.gsd = MendelianRandomization::mr_conmix(mr.obj.gsd)
      # colnames
      conmix.colname_theta_gsd <- paste0("CONMIX_theta_sim_BMI_GSD_", format(n[i], scientific=FALSE))
      conmix.colname_theta_se_gsd <- paste0("CONMIX_theta_se_sim_BMI_GSD_", format(n[i], scientific=FALSE))
      conmix.colname_pvalue_gsd <- paste0("CONMIX_mr_pval_BMI_GSD_", format(n[i], scientific=FALSE))
      ## store results
      est_sim[k, conmix.colname_theta_gsd] <- conmix.res.gsd$Estimate
      CIlength = conmix.res.gsd$CIUpper-conmix.res.gsd$CILower
      if (length(CIlength)>1) print("conmix multimodal")
      est_sim[k, conmix.colname_theta_se_gsd] <- sum(CIlength)/1.96/2  ## Caution: this may be problematic 
      est_sim[k, conmix.colname_pvalue_gsd] <- conmix.res.gsd$Pvalue
      rm(conmix.res.gsd)  
      
      ## B/ GBC ---------------------------------------------------------------
      mr.obj.gbc = MendelianRandomization::mr_input(
        bx = sim_dat$beta.bmi_sim_norm, 
        bxse = sim_dat$se.bmi, 
        # outcome
        by = sim_dat$beta.gbc_sim_norm, 
        byse = sim_dat$se.gbc, 
        snps = sim_dat$SNP,
        exposure = "Body mass index",
        outcome = "Gallbladder Cancer"
      ) 
      
      # 1. random IVW
      ivw.res.gbc = MendelianRandomization::mr_ivw(mr.obj.gbc, model = "random")
      # colnames
      ivw.colname_theta_gbc <- paste0("IVW_theta_sim_BMI_GBC_", format(n[i], scientific=FALSE))
      ivw.colname_theta_se_gbc <- paste0("IVW_theta_se_sim_BMI_GBC_", format(n[i], scientific=FALSE))
      ivw.colname_pvalue_gbc <- paste0("IVW_mr_pval_BMI_GBC_", format(n[i], scientific=FALSE))
      # store results
      est_sim[k, ivw.colname_theta_gbc] <- ivw.res.gbc$Estimate
      est_sim[k, ivw.colname_theta_se_gbc] <- ivw.res.gbc$StdError
      est_sim[k, ivw.colname_pvalue_gbc] <- ivw.res.gbc$Pvalue
      
      # 2. weighted mode
      mode.res.gbc = MendelianRandomization::mr_mbe(mr.obj.gbc)
      # colnames
      mode.colname_theta_gbc <- paste0("MODE_theta_sim_BMI_GBC_", format(n[i], scientific=FALSE))
      mode.colname_theta_se_gbc <- paste0("MODE_theta_se_sim_BMI_GBC_", format(n[i], scientific=FALSE))
      mode.colname_pvalue_gbc <- paste0("MODE_mr_pval_BMI_GBC_", format(n[i], scientific=FALSE))
      # store results
      est_sim[k, mode.colname_theta_gbc] <- mode.res.gbc$Estimate
      est_sim[k, mode.colname_theta_se_gbc] <- mode.res.gbc$StdError
      est_sim[k, mode.colname_pvalue_gbc] <- mode.res.gbc$Pvalue
      
      # 3. CONMIX
      conmix.res.gbc <- tryCatch(
        expr = {
          MendelianRandomization::mr_conmix(mr.obj.gbc)
        },
        error = function(e){
          # the error we encounter is with "i=7" 
          #     Error in seq.default(from = CIMin, to = CIMax, by = CIStep) : 
          #         wrong sign in 'by' argument
          message(paste0("Caught an error with mr_conmix (BMI -> GBC)!\nSample size: ", n[i], 
                         "\nComputed CI's with predefined range of [-15;15]"))
          message("Below is the error message from R:")
          print(e) 
          return(MendelianRandomization::mr_conmix(mr.obj.gbc, CIMin = -15, CIMax=15))
        },
        warning = function(w){
          message('Caught an warning!')
          print(w)
        },
        finally = {
          #message('All done, quitting.')
        }
      )
      conmix.colname_theta_gbc <- paste0("CONMIX_theta_sim_BMI_GBC_", format(n[i], scientific=FALSE))
      conmix.colname_theta_se_gbc <- paste0("CONMIX_theta_se_sim_BMI_GBC_", format(n[i], scientific=FALSE))
      conmix.colname_pvalue_gbc <- paste0("CONMIX_mr_pval_BMI_GBC_", format(n[i], scientific=FALSE))
      
      ## store results
      est_sim[k, conmix.colname_theta_gbc] <- conmix.res.gbc$Estimate
      CIlength = conmix.res.gbc$CIUpper-conmix.res.gbc$CILower
      if (length(CIlength)>1) print("conmix multimodal")
      est_sim[k, conmix.colname_theta_se_gbc] <- sum(CIlength)/1.96/2  ## Caution: this may be problematic 
      est_sim[k, conmix.colname_pvalue_gbc] <- conmix.res.gbc$Pvalue
      
      rm(conmix.res.gbc, mode.res.gbc, mode.res.gsd, ivw.res.gbc, ivw.res.gsd)
    }
    
    #rm(se_update_BMI, se_update_GBC)
    print(paste0("|||-----------------------Run finished for sample size: ", format(n[i], scientific=FALSE), " -----------------------|||"))
  }  ## END OF SIMULATION
  
  T1 = proc.time()[3]
  timediff = T1-T0
  
  saveRDS(est_sim, paste0("./output/Rdata/", Sys.Date(), "_est_sim_NSIM_", nsim,".rds"))
}

