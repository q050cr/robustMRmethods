
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
source("scripts/helper/MR_lasso.R") # Slob&Burgess 2020
source("scripts/helper/PVE_calculation_fns.R")
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
threshold <- 5 * 10^-6
n.orig_sqrt <- sqrt(n.orig)

est_sim <- tibble(no_sim = 1:nsim)

## START SIM
if (runSimulation == TRUE) {
  T0 <- proc.time()[3]
  for (i in 1:length(colnames_pval)) {
    numIV <- sum(my_data_harm[[colnames_pval[i]]] < threshold)
    if (numIV < 5) next # start with 5 IVs only (sample size ~50,000)

    # FILTER DATA
    sim_dat <- my_data_harm %>%
      filter(!!sym(colnames_pval[i]) < threshold) %>%
      select(SNP, eaf.bmi, beta.bmi, se.bmi, p.value.bmi, eaf.gsd, beta.gsd, se.gsd, p.value.gsd) %>%
      mutate(se_update_BMI = abs(beta.bmi / -qnorm(p = (p.value.bmi / 2), mean = 0) * n.orig_sqrt / sqrt(n[i])))

    # simulate beta values from normal_distribution # -----------------------------------------------------
    for (k in 1:nsim) {
      sim_dat <- sim_dat %>%
        mutate(beta.bmi_sim_norm = purrr::map2_dbl(.x = beta.bmi, .y = se_update_BMI, .f = ~ rnorm(n = 1, mean = .x, sd = .y))) %>%
        mutate(explained_variance_X2_sim = explained_variance_numeric2(maf = eaf.bmi, beta = beta.bmi_sim_norm, se_beta = se_update_BMI, samplesize = n[i]))
      # PVE
      colname_pve_X <- paste0("pve_X_sim_", format(n[i], scientific = FALSE))
      ## add PVE cols
      est_sim[k, colname_pve_X] <- sum(sim_dat$explained_variance_X2_sim)

      # PERFORM MR_analysis
      ## A/ GSD  ---------------------------------------------------------------
      mr.obj.gsd <- MendelianRandomization::mr_input(
        bx = sim_dat$beta.bmi_sim_norm,
        bxse = sim_dat$se_update_BMI,
        # outcome
        by = sim_dat$beta.gsd,
        byse = sim_dat$se.gsd,
        snps = sim_dat$SNP,
        exposure = "Body mass index",
        outcome = "Gallstone Disease"
      )

      # 1. random IVW
      ivw.res.gsd <- MendelianRandomization::mr_ivw(mr.obj.gsd, model = "random")
      # colnames
      ivw.colname_theta_gsd <- paste0("IVW_theta_sim_BMI_GSD_", format(n[i], scientific = FALSE))
      ivw.colname_theta_se_gsd <- paste0("IVW_theta_se_sim_BMI_GSD_", format(n[i], scientific = FALSE))
      ivw.colname_pvalue_gsd <- paste0("IVW_mr_pval_BMI_GSD_", format(n[i], scientific = FALSE))
      # store results
      est_sim[k, ivw.colname_theta_gsd] <- ivw.res.gsd$Estimate
      est_sim[k, ivw.colname_theta_se_gsd] <- ivw.res.gsd$StdError
      est_sim[k, ivw.colname_pvalue_gsd] <- ivw.res.gsd$Pvalue

      # 2. weighted mode
      mode.res.gsd <- MendelianRandomization::mr_mbe(mr.obj.gsd)
      # colnames
      mode.colname_theta_gsd <- paste0("MODE_theta_sim_BMI_GSD_", format(n[i], scientific = FALSE))
      mode.colname_theta_se_gsd <- paste0("MODE_theta_se_sim_BMI_GSD_", format(n[i], scientific = FALSE))
      mode.colname_pvalue_gsd <- paste0("MODE_mr_pval_BMI_GSD_", format(n[i], scientific = FALSE))
      # store results
      est_sim[k, mode.colname_theta_gsd] <- mode.res.gsd$Estimate
      est_sim[k, mode.colname_theta_se_gsd] <- mode.res.gsd$StdError
      est_sim[k, mode.colname_pvalue_gsd] <- mode.res.gsd$Pvalue

      # 3. CONMIX
      conmix.res.gsd <- tryCatch(
        expr = {
          MendelianRandomization::mr_conmix(mr.obj.gsd)
        },
        error = function(e) {
          # the error we encounter is with "i=7"
          #     Error in seq.default(from = CIMin, to = CIMax, by = CIStep) :
          #         wrong sign in 'by' argument
          message(paste0(
            "Caught an error with mr_conmix (BMI -> GSD)!\nSample size: ", n[i],
            "\nComputed CI's with predefined range of [-15;15]"
          ))
          message("Below is the error message from R:")
          print(e)
          return(MendelianRandomization::mr_conmix(mr.obj.gbc, CIMin = -15, CIMax = 15))
        },
        warning = function(w) {
          message("Caught an warning!")
          print(w)
        },
        finally = {
          # message('All done, quitting.')
        }
      )
      # colnames
      conmix.colname_theta_gsd <- paste0("CONMIX_theta_sim_BMI_GSD_", format(n[i], scientific = FALSE))
      conmix.colname_theta_se_gsd <- paste0("CONMIX_theta_se_sim_BMI_GSD_", format(n[i], scientific = FALSE))
      conmix.colname_pvalue_gsd <- paste0("CONMIX_mr_pval_BMI_GSD_", format(n[i], scientific = FALSE))
      ## store results
      est_sim[k, conmix.colname_theta_gsd] <- conmix.res.gsd$Estimate
      CIlength <- conmix.res.gsd$CIUpper - conmix.res.gsd$CILower
      if (length(CIlength) > 1) print("conmix multimodal")
      est_sim[k, conmix.colname_theta_se_gsd] <- sum(CIlength) / 1.96 / 2 ## Caution: this may be problematic
      est_sim[k, conmix.colname_pvalue_gsd] <- conmix.res.gsd$Pvalue
      rm(conmix.res.gsd)

      # 4. MRMix  (takes time!!)
      # theta_temp_vec = seq(-0.5,0.5,by=0.01)
      mrmix.res.gsd <- MRMix::MRMix(
        betahat_x = sim_dat$beta.bmi_sim_norm,
        betahat_y = sim_dat$beta.gsd,
        sx = sim_dat$se_update_BMI,
        sy = sim_dat$se.gsd
      )
      mrmix.res_se.gsd <- MRMix::MRMix_se(
        betahat_x = sim_dat$beta.bmi_sim_norm,
        betahat_y = sim_dat$beta.gsd,
        sx = sim_dat$se_update_BMI,
        sy = sim_dat$se.gsd,
        theta = mrmix.res.gsd$theta, # estimate of causal effect, assuming the summary stats are standardized, theta represents increase in mean value of Y in s.d. unit of Y (for cont outcomes) or log-OR of Y (for binary outcomes)
        pi0 = mrmix.res.gsd$pi0,
        sigma2 = mrmix.res.gsd$sigma2
      )

      # colnames
      mrmix.colname_theta_gsd <- paste0("MRMIX_theta_sim_BMI_GSD_", format(n[i], scientific = FALSE))
      mrmix.colname_theta_se_gsd <- paste0("MRMIX_theta_se_sim_BMI_GSD_", format(n[i], scientific = FALSE))
      mrmix.colname_pvalue_gsd <- paste0("MRMIX_mr_pval_BMI_GSD_", format(n[i], scientific = FALSE))
      # store results
      est_sim[k, mrmix.colname_theta_gsd] <- mrmix.res.gsd$theta
      est_sim[k, mrmix.colname_theta_se_gsd] <- mrmix.res_se.gsd
      est_sim[k, mrmix.colname_pvalue_gsd] <- mrmix.res.gsd$pvalue_theta
      rm(mrmix.res.gsd, mrmix.res_se.gsd)
    }
    print(paste0("|||-----------------------Run finished for sample size: ", format(n[i], scientific = FALSE), " -----------------------|||"))
  } ## END OF SIMULATION

  T1 <- proc.time()[3]
  timediff <- T1 - T0

  saveRDS(est_sim, paste0("./output/Rdata/", Sys.Date(), "_est_sim_NSIM_", nsim, ".rds"))
  print(timediff)
}
