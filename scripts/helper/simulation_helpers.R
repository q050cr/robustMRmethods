
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
# source("scripts/helper/MR_lasso.R") # Slob&Burgess 2020, now implemented in `MendelianRandomization`
source("scripts/helper/PVE_calculation_fns.R")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

sim_function <- function(dat, index) {
  if(any(is.na(dat))) return(NULL)
  if(nrow(dat)<=3) return(NULL)  # cannot perform robust methods with one IV, "mr_presso method requires data on >3 variants."
  sim_vars <- readRDS(file = dplyr::last(list.files("./output/Rdata/", pattern = "_sim_vars.rds", full.names = TRUE)))
  n <- sim_vars$n
  
  return_vector <- c()
  dat <- dat %>% 
    mutate(beta.bmi_sim_norm = purrr::map2_dbl(.x = beta.bmi, .y = se_update_BMI, .f = ~rnorm(n = 1, mean = .x, sd = .y))) %>% 
    mutate(explained_variance_X2_sim = explained_variance_numeric2(maf = eaf.bmi, beta = beta.bmi_sim_norm, se_beta = se_update_BMI, samplesize = n[index]))
  # PVE 
  return_vector <- sum(dat$explained_variance_X2_sim)
  
  # PERFORM MR_analysis 
  ## A/ GSD  ---------------------------------------------------------------
  mr.obj.gsd = MendelianRandomization::mr_input(
    bx = dat$beta.bmi_sim_norm, 
    bxse = dat$se_update_BMI, 
    # outcome
    by = dat$beta.gsd, 
    byse = dat$se.gsd, 
    snps = dat$SNP,
    exposure = "Body mass index",
    outcome = "Gallstone Disease"
  ) 
  
  # 1. random IVW
  ivw.res.gsd = MendelianRandomization::mr_ivw(mr.obj.gsd, model = "random")
  # store results
  return_vector <- append(return_vector,ivw.res.gsd$Estimate)
  return_vector <- append(return_vector,ivw.res.gsd$StdError)
  return_vector <- append(return_vector, ivw.res.gsd$Pvalue)
  
  # 2. weighted mode
  mode.res.gsd = MendelianRandomization::mr_mbe(mr.obj.gsd)
  # store results
  return_vector <- append(return_vector, mode.res.gsd$Estimate)
  return_vector <- append(return_vector, mode.res.gsd$StdError)
  return_vector <- append(return_vector, mode.res.gsd$Pvalue)
  
  # 3. CONMIX 
  conmix.res.gsd <- tryCatch(
    expr = {
      MendelianRandomization::mr_conmix(mr.obj.gsd)
    },
    error = function(e){
      # the error we encounter is with "i=7" 
      #     Error in seq.default(from = CIMin, to = CIMax, by = CIStep) : 
      #         wrong sign in 'by' argument
      message(paste0("Caught an error with mr_conmix (BMI -> GSD)!\nSample size: ", n[index], 
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
  CIlength = conmix.res.gsd$CIUpper-conmix.res.gsd$CILower
  if (length(CIlength)>1) print("conmix multimodal")
  # store results
  return_vector <- append(return_vector, conmix.res.gsd$Estimate)
  return_vector <- append(return_vector, sum(CIlength)/1.96/2)
  return_vector <- append(return_vector, conmix.res.gsd$Pvalue)
  
  # 4. MRMix  (takes time!!)
  # theta_temp_vec = seq(-0.5,0.5,by=0.01)
  mrmix.res.gsd = MRMix::MRMix(betahat_x = dat$beta.bmi_sim_norm, 
                               betahat_y = dat$beta.gsd, 
                               sx = dat$se_update_BMI, 
                               sy = dat$se.gsd
  )
  mrmix.res_se.gsd = MRMix::MRMix_se(betahat_x = dat$beta.bmi_sim_norm,
                                     betahat_y = dat$beta.gsd, 
                                     sx = dat$se_update_BMI, 
                                     sy = dat$se.gsd,
                                     theta = mrmix.res.gsd$theta, # estimate of causal effect, assuming the summary stats are standardized, theta represents increase in mean value of Y in s.d. unit of Y (for cont outcomes) or log-OR of Y (for binary outcomes)
                                     pi0 = mrmix.res.gsd$pi0, 
                                     sigma2 = mrmix.res.gsd$sigma2
  )
  # store results
  return_vector <- append(return_vector, mrmix.res.gsd$theta)
  return_vector <- append(return_vector, mrmix.res_se.gsd)
  return_vector <- append(return_vector, mrmix.res.gsd$pvalue_theta)
  
  ## updated 2023-03-01
  # 5. weighted median
  median.res.gsd = MendelianRandomization::mr_median(mr.obj.gsd)
  return_vector <- append(return_vector, median.res.gsd$Estimate)
  return_vector <- append(return_vector, median.res.gsd$StdError)
  return_vector <- append(return_vector, median.res.gsd$Pvalue)
  
  # 6. Egger regression
  egger.res.gsd <- MendelianRandomization::mr_egger(mr.obj.gsd)
  return_vector <- append(return_vector,egger.res.gsd$Estimate)
  return_vector <- append(return_vector, egger.res.gsd$StdError.Est)
  return_vector <- append(return_vector, egger.res.gsd$Pvalue.Est)
  
  # 7. MR-robust
  robust.res.gsd <- MendelianRandomization::mr_ivw(mr.obj.gsd, "random", robust = TRUE) 
  return_vector <- append(return_vector,robust.res.gsd$Estimate)
  return_vector <- append(return_vector, robust.res.gsd$StdError)
  return_vector <- append(return_vector, robust.res.gsd$Pvalue)
  
  # 8. MR-PRESSO
  presso.df = data.frame(bx = dat$beta.bmi_sim_norm, by = dat$beta.gsd, bxse = dat$se_update_BMI, byse = dat$se.gsd)
  presso.res.gsd = MRPRESSO::mr_presso(BetaOutcome = "by", BetaExposure = "bx", SdOutcome = "byse", SdExposure = "bxse", 
                  OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = presso.df, NbDistribution = 500, SignifThreshold = 0.05)
  if (!is.na(presso.res.gsd$`Main MR results`[2,"Causal Estimate"]) & !is.na(presso.res.gsd$`Main MR results`[2,"Sd"])){
    return_vector <-  append(return_vector, presso.res.gsd$`Main MR results`[2,"Causal Estimate"])
    return_vector <- append(return_vector, presso.res.gsd$`Main MR results`[2,"Sd"])
    return_vector <- append(return_vector, presso.res.gsd$`Main MR results`[2,"P-value"])
  } else{  # if no outlier-corrected analysis was calculated
    return_vector <-  append(return_vector, presso.res.gsd$`Main MR results`[1,"Causal Estimate"])
    return_vector <- append(return_vector, presso.res.gsd$`Main MR results`[1,"Sd"])
    return_vector <- append(return_vector, presso.res.gsd$`Main MR results`[1,"P-value"])
  }

  # 9. MR-RAPS
  raps.res.gsd = mr.raps::mr.raps.overdispersed.robust(presso.df$bx, presso.df$by, presso.df$bxse, presso.df$byse,
                                     loss.function = "huber", k = 1.345, initialization = c("l2"), 
                                     suppress.warning = FALSE, diagnosis = FALSE, niter = 20, tol = .Machine$double.eps^0.5)
  return_vector <-  append(return_vector, raps.res.gsd$beta.hat)
  return_vector <-  append(return_vector, raps.res.gsd$beta.se)
  return_vector <-  append(return_vector, raps.res.gsd$beta.p.value)
  
  # 10. MR-Lasso
  lasso.res.gsd <- tryCatch(
    expr = {
      MendelianRandomization::mr_lasso(object = mr.obj.gsd)
    },
    error = function(e){
      # the error we encounter is with "i=7" 
      #     Error in seq.default(from = CIMin, to = CIMax, by = CIStep) : 
      #         wrong sign in 'by' argument
      message(paste0("Caught an error with mr_lasso (BMI -> GSD)!\nSample size: ", n[index]))
      message("Below is the error message from R:")
      print(e) 
      return(NA)
    },
    warning = function(w){
      message('Caught an warning!')
      print(w)
    },
    finally = {
      #message('All done, quitting.')
    }
  )
  if (is.na(lasso.res.gsd)) {
    return_vector <- append(return_vector, NA)
    return_vector <- append(return_vector, NA)
    return_vector <- append(return_vector, NA)
  } else{
    return_vector <- append(return_vector, lasso.res.gsd$Estimate)
    return_vector <- append(return_vector, lasso.res.gsd$StdError)
    return_vector <- append(return_vector, lasso.res.gsd$Pvalue)
  }
  ## return vector
  return(return_vector)
}


# POPULATE est_sim_temp --------------------------------------------------------
populate_est_sim <- function(nsim_list, nsim, no_cases) {
  if(is.null(unlist(nsim_list))) return(NULL)
  est_sim_temp <- tibble(no_sim = 1:nsim)
  for (nsims in 1:length(nsim_list)) {
    # cum_PVE_X
    est_sim_temp[nsims, paste0("pve_X_sim_", format(no_cases, scientific = FALSE))] <- nsim_list[[nsims]][1] 
    # 1. random IVW
    est_sim_temp[nsims, paste0("IVW_theta_sim_BMI_GSD_", format(no_cases, scientific = FALSE))] <- nsim_list[[nsims]][2] 
    est_sim_temp[nsims, paste0("IVW_theta_se_sim_BMI_GSD_", format(no_cases, scientific = FALSE))] <- nsim_list[[nsims]][3] 
    est_sim_temp[nsims, paste0("IVW_mr_pval_BMI_GSD_", format(no_cases, scientific = FALSE))] <- nsim_list[[nsims]][4] 
    # 2. weighted mode
    est_sim_temp[nsims, paste0("MODE_theta_sim_BMI_GSD_", format(no_cases, scientific = FALSE))] <- nsim_list[[nsims]][5] 
    est_sim_temp[nsims, paste0("MODE_theta_se_sim_BMI_GSD_", format(no_cases, scientific = FALSE))] <- nsim_list[[nsims]][6] 
    est_sim_temp[nsims, paste0("MODE_mr_pval_BMI_GSD_", format(no_cases, scientific = FALSE))] <- nsim_list[[nsims]][7] 
    # 3. CONMIX 
    est_sim_temp[nsims, paste0("CONMIX_theta_sim_BMI_GSD_", format(no_cases, scientific = FALSE))] <- nsim_list[[nsims]][8] 
    est_sim_temp[nsims, paste0("CONMIX_theta_se_sim_BMI_GSD_", format(no_cases, scientific = FALSE))] <- nsim_list[[nsims]][9] 
    est_sim_temp[nsims, paste0("CONMIX_mr_pval_BMI_GSD_", format(no_cases, scientific = FALSE))] <- nsim_list[[nsims]][10] 
    # 4. MRMix 
    est_sim_temp[nsims, paste0("MRMIX_theta_sim_BMI_GSD_", format(no_cases, scientific = FALSE))] <- nsim_list[[nsims]][11] 
    est_sim_temp[nsims, paste0("MRMIX_theta_se_sim_BMI_GSD_", format(no_cases, scientific = FALSE))] <- nsim_list[[nsims]][12] 
    est_sim_temp[nsims, paste0("MRMIX_mr_pval_BMI_GSD_", format(no_cases, scientific = FALSE))] <- nsim_list[[nsims]][13] 
    # 5. Median
    est_sim_temp[nsims, paste0("MEDIAN_theta_sim_BMI_GSD_", format(no_cases, scientific = FALSE))] <- nsim_list[[nsims]][14] 
    est_sim_temp[nsims, paste0("MEDIAN_theta_se_sim_BMI_GSD_", format(no_cases, scientific = FALSE))] <- nsim_list[[nsims]][15] 
    est_sim_temp[nsims, paste0("MEDIAN_mr_pval_BMI_GSD_", format(no_cases, scientific = FALSE))] <- nsim_list[[nsims]][16] 
    # 6. Egger 
    est_sim_temp[nsims, paste0("EGGER_theta_sim_BMI_GSD_", format(no_cases, scientific = FALSE))] <- nsim_list[[nsims]][17] 
    est_sim_temp[nsims, paste0("EGGER_theta_se_sim_BMI_GSD_", format(no_cases, scientific = FALSE))] <- nsim_list[[nsims]][18] 
    est_sim_temp[nsims, paste0("EGGER_mr_pval_BMI_GSD_", format(no_cases, scientific = FALSE))] <- nsim_list[[nsims]][19] 
    # 7. Robust 
    est_sim_temp[nsims, paste0("ROBUST_theta_sim_BMI_GSD_", format(no_cases, scientific = FALSE))] <- nsim_list[[nsims]][20] 
    est_sim_temp[nsims, paste0("ROBUST_theta_se_sim_BMI_GSD_", format(no_cases, scientific = FALSE))] <- nsim_list[[nsims]][21] 
    est_sim_temp[nsims, paste0("ROBUST_mr_pval_BMI_GSD_", format(no_cases, scientific = FALSE))] <- nsim_list[[nsims]][22] 
    # 8. Presso 
    est_sim_temp[nsims, paste0("PRESSO_theta_sim_BMI_GSD_", format(no_cases, scientific = FALSE))] <- nsim_list[[nsims]][23] 
    est_sim_temp[nsims, paste0("PRESSO_theta_se_sim_BMI_GSD_", format(no_cases, scientific = FALSE))] <- nsim_list[[nsims]][24] 
    est_sim_temp[nsims, paste0("PRESSO_mr_pval_BMI_GSD_", format(no_cases, scientific = FALSE))] <- nsim_list[[nsims]][25] 
    # 9. RAPS 
    est_sim_temp[nsims, paste0("RAPS_theta_sim_BMI_GSD_", format(no_cases, scientific = FALSE))] <- nsim_list[[nsims]][26] 
    est_sim_temp[nsims, paste0("RAPS_theta_se_sim_BMI_GSD_", format(no_cases, scientific = FALSE))] <- nsim_list[[nsims]][27] 
    est_sim_temp[nsims, paste0("RAPS_mr_pval_BMI_GSD_", format(no_cases, scientific = FALSE))] <- nsim_list[[nsims]][28] 
    # 10. lasso 
    est_sim_temp[nsims, paste0("LASSO_theta_sim_BMI_GSD_", format(no_cases, scientific = FALSE))] <- nsim_list[[nsims]][26] 
    est_sim_temp[nsims, paste0("LASSO_theta_se_sim_BMI_GSD_", format(no_cases, scientific = FALSE))] <- nsim_list[[nsims]][27] 
    est_sim_temp[nsims, paste0("LASSO_mr_pval_BMI_GSD_", format(no_cases, scientific = FALSE))] <- nsim_list[[nsims]][28] 
  }
  return(est_sim_temp)
}




