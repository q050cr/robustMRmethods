
## ---------------------------
##
## Script name: perform_mr.R
##
## Author: Christoph Reich
## Date Created: 2022-11-23
##
## Copyright (c) Christoph Reich, 2022
## Email: christoph.reich@med.uni-heidelberg.de
##
## ---------------------------
## Notes:
##   creates dataframe `est` with estimates from different robust MR methods: BMI -> GBC
##   creates dataframe `est.bmi.gsd` with estimates from different robust MR methods: BMI -> GSD
##        
## ---------------------------

library(tibble)
library(dplyr)
# mr methods
library(MendelianRandomization)
library(MRMix)
library(MRPRESSO)
library(mr.raps)
library(penalized)
source("./scripts/helper/MR_lasso.R")

## load data created in MSc.rmd
my_data_harm <- readRDS(file = dplyr::last(list.files("./output/Rdata/", pattern = "my_data_harm.rds", full.names = TRUE)))

# perform MR analysis -----------------------------------------------------
set.seed(123)

mr_methods = c("IVW", "median", "mode", "PRESSO", "Robust", "Lasso", "egger", "conmix", "MRMix", "RAPS")
est <- tibble(mr_methods=mr_methods, estimate=NA, se.estimate=NA, time=NA)
numIV <- nrow(my_data_harm)

if (numIV>2){
  mr.obj = MendelianRandomization::mr_input(
    # exposure
    bx = my_data_harm$beta.bmi, 
    bxse = my_data_harm$se.bmi, 
    # outcome
    by = my_data_harm$beta.gbc, 
    byse = my_data_harm$se.gbc, 
    snps = my_data_harm$SNP,
    exposure = "Body mass index",
    outcome = "Gallbladder Cancer"
  )
  # 1. IVW
  T0 = proc.time()[3]
  res = MendelianRandomization::mr_ivw(mr.obj, model = "random")
  T1 = proc.time()[3]
  est$estimate[est$mr_methods=="IVW"] = res$Estimate
  est$se.estimate[est$mr_methods=="IVW"] = res$StdError
  est$time[est$mr_methods=="IVW"] = T1-T0
  rm(res)
  
  # 2. median
  T0 = proc.time()[3]
  res = MendelianRandomization::mr_median(mr.obj)
  T1 = proc.time()[3]
  est$estimate[est$mr_methods=="median"] = res$Estimate
  est$se.estimate[est$mr_methods=="median"] = res$StdError
  est$time[est$mr_methods=="median"] = T1-T0
  rm(res)
  # 3. mode
  T0 = proc.time()[3]
  res = MendelianRandomization::mr_mbe(mr.obj)
  T1 = proc.time()[3]
  est$estimate[est$mr_methods=="mode"] = res$Estimate
  est$se.estimate[est$mr_methods=="mode"] = res$StdError
  est$time[est$mr_methods=="mode"] = T1-T0
  rm(res)
  # 4.MR-PRESSO    ## CAVE: 5 minutes per round, consider only calculating for defined thresholds
  T0 = proc.time()[3]
  res <- MRPRESSO:: mr_presso(BetaOutcome = "beta.gbc", 
                              BetaExposure = "beta.bmi", 
                              SdOutcome = "se.gbc", 
                              SdExposure = "se.bmi", 
                              OUTLIERtest = TRUE, 
                              DISTORTIONtest = TRUE, 
                              seed = 20221114,
                              data = as.data.frame(my_data_harm), ## must be a data.frame
                              NbDistribution = 1000,   # must be changed
                              SignifThreshold = 0.05)
  T1 = proc.time()[3]
  if (!is.na(res$`Main MR results`[2,"Causal Estimate"]) & !is.na(res$`Main MR results`[2,"Sd"])){
    est$estimate[est$mr_methods=="PRESSO"] =  res$`Main MR results`[2,"Causal Estimate"]
    est$se.estimate[est$mr_methods=="PRESSO"] = res$`Main MR results`[2,"Sd"]
  } else{
    est$estimate[est$mr_methods=="PRESSO"] =  res$`Main MR results`[1,"Causal Estimate"]
    est$se.estimate[est$mr_methods=="PRESSO"] = res$`Main MR results`[1,"Sd"]
  }
  est$time[est$mr_methods=="PRESSO"] = T1-T0
  rm(res)
  # 5. robust
  T0 = proc.time()[3]
  res = MendelianRandomization::mr_ivw(mr.obj,"random", robust = TRUE)
  T1 = proc.time()[3]
  est$estimate[est$mr_methods=="Robust"] = res$Estimate
  est$se.estimate[est$mr_methods=="Robust"] = res$StdError
  est$time[est$mr_methods=="Robust"] = T1-T0
  rm(res)
  # 6. MR-Lasso
  T0 = proc.time()[3]
  res = MR_lasso(betaYG = my_data_harm$beta.gbc, betaXG = my_data_harm$beta.bmi, sebetaYG = my_data_harm$se.gbc)
  T1 = proc.time()[3]
  est$estimate[est$mr_methods=="Lasso"] = res$ThetaEstimate
  est$se.estimate[est$mr_methods=="Lasso"] = res$ThetaSE
  est$time[est$mr_methods=="Lasso"] = T1-T0
  rm(res)
  # 7. Egger
  T0 = proc.time()[3]
  res = MendelianRandomization::mr_egger(mr.obj)
  T1 = proc.time()[3]
  est$estimate[est$mr_methods=="egger"] = res$Estimate
  est$se.estimate[est$mr_methods=="egger"] = res$StdError.Est
  est$time[est$mr_methods=="egger"] = T1-T0
  rm(res)
  # 8. contamination mixture
  T0 = proc.time()[3]
  res = MendelianRandomization::mr_conmix(mr.obj)
  T1 = proc.time()[3]
  est$estimate[est$mr_methods=="conmix"] = res$Estimate
  CIlength = res$CIUpper-res$CILower
  if (length(CIlength)>1) print("conmix multimodal")
  est$se.estimate[est$mr_methods=="conmix"] = sum(CIlength)/1.96/2 ## Caution: this may be problematic (Qi&Chatterjee2021)
  est$time[est$mr_methods=="conmix"] = T1-T0
  rm(res)
  # 9. MRMix
  # theta_temp_vec = seq(-0.5,0.5,by=0.01)
  T0 = proc.time()[3]
  res = MRMix::MRMix(betahat_x = my_data_harm$beta.bmi, 
                     betahat_y = my_data_harm$beta.gbc, 
                     sx =my_data_harm$se.bmi, 
                     sy = my_data_harm$se.gbc)
  res_se = MRMix::MRMix_se(betahat_x = my_data_harm$beta.bmi, 
                           betahat_y = my_data_harm$beta.gbc, 
                           sx =my_data_harm$se.bmi, 
                           sy = my_data_harm$se.gbc, 
                           theta = res$theta, # estimate of causal effect, assuming the summary stats are standardized, theta represents increase in mean value of Y in s.d. unit of Y (for cont outcomes) or log-OR of Y (for binary outcomes)
                           pi0 = res$pi0, 
                           sigma2 = res$sigma2
  )
  T1 = proc.time()[3]
  est$estimate[est$mr_methods=="MRMix"] = res$theta
  est$se.estimate[est$mr_methods=="MRMix"] = res_se
  est$time[est$mr_methods=="MRMix"] = T1-T0
  rm(res, res_se)
  # 10. MR-RAPS  
  ## Warning: The estimated overdispersion parameter is very small. Consider using the simple model without overdispersion.
  T0 = proc.time()[3]
  res = mr.raps::mr.raps.overdispersed.robust(b_exp = my_data_harm$beta.bmi, 
                                              b_out = my_data_harm$beta.gbc, 
                                              se_exp = my_data_harm$se.bmi, 
                                              se_out = my_data_harm$se.gbc, 
                                              loss.function = "huber", 
                                              k = 1.345, 
                                              initialization = c("l2"), 
                                              suppress.warning = FALSE, 
                                              diagnosis = FALSE, 
                                              niter = 20, 
                                              tol = .Machine$double.eps^0.5)
  T1 = proc.time()[3]
  est$estimate[est$mr_methods=="RAPS"] = res$beta.hat
  est$se.estimate[est$mr_methods=="RAPS"] = res$beta.se
  est$time[est$mr_methods=="RAPS"] = T1-T0
  rm(res)
}

## SAVE estimates
saveRDS(est, paste0("./output/Rdata/", Sys.Date(), "_est.rds"))


# BMI -> GSD --------------------------------------------------------------

mr_methods = c("IVW", "median", "mode", "PRESSO", "Robust", "Lasso", "egger", "conmix", "MRMix", "RAPS")
est.bmi.gsd <- tibble(mr_methods=mr_methods, estimate=NA, se.estimate=NA, time=NA)
numIV <- nrow(my_data_harm)

if (numIV>2){
  mr.obj = MendelianRandomization::mr_input(
    # exposure
    bx = my_data_harm$beta.bmi, 
    bxse = my_data_harm$se.bmi, 
    # outcome
    by = my_data_harm$beta.gsd, 
    byse = my_data_harm$se.gsd, 
    snps = my_data_harm$SNP,
    exposure = "Body mass index",
    outcome = "Gallstone disease"
  )
  # 1. IVW
  T0 = proc.time()[3]
  res = MendelianRandomization::mr_ivw(mr.obj, model = "random")
  T1 = proc.time()[3]
  est.bmi.gsd$estimate[est.bmi.gsd$mr_methods=="IVW"] = res$Estimate
  est.bmi.gsd$se.estimate[est.bmi.gsd$mr_methods=="IVW"] = res$StdError
  est.bmi.gsd$time[est.bmi.gsd$mr_methods=="IVW"] = T1-T0
  rm(res)
  
  # 2. median
  T0 = proc.time()[3]
  res = MendelianRandomization::mr_median(mr.obj)
  T1 = proc.time()[3]
  est.bmi.gsd$estimate[est.bmi.gsd$mr_methods=="median"] = res$Estimate
  est.bmi.gsd$se.estimate[est.bmi.gsd$mr_methods=="median"] = res$StdError
  est.bmi.gsd$time[est.bmi.gsd$mr_methods=="median"] = T1-T0
  rm(res)
  # 3. mode
  T0 = proc.time()[3]
  res = MendelianRandomization::mr_mbe(mr.obj)
  T1 = proc.time()[3]
  est.bmi.gsd$estimate[est.bmi.gsd$mr_methods=="mode"] = res$Estimate
  est.bmi.gsd$se.estimate[est.bmi.gsd$mr_methods=="mode"] = res$StdError
  est.bmi.gsd$time[est.bmi.gsd$mr_methods=="mode"] = T1-T0
  rm(res)
  # 4.MR-PRESSO    ## CAVE: 5 minutes per round, consider only calculating for defined thresholds
  T0 = proc.time()[3]
  res <- MRPRESSO:: mr_presso(BetaOutcome = "beta.gsd", 
                              BetaExposure = "beta.bmi", 
                              SdOutcome = "se.gsd", 
                              SdExposure = "se.bmi", 
                              OUTLIERtest = TRUE, 
                              DISTORTIONtest = TRUE, 
                              seed = 20221114,
                              data = as.data.frame(my_data_harm), ## must be a data.frame
                              NbDistribution = 1000,   # must be changed
                              SignifThreshold = 0.05)
  T1 = proc.time()[3]
  if (!is.na(res$`Main MR results`[2,"Causal Estimate"]) & !is.na(res$`Main MR results`[2,"Sd"])){
    est.bmi.gsd$estimate[est.bmi.gsd$mr_methods=="PRESSO"] =  res$`Main MR results`[2,"Causal Estimate"]
    est.bmi.gsd$se.estimate[est.bmi.gsd$mr_methods=="PRESSO"] = res$`Main MR results`[2,"Sd"]
  } else{
    est.bmi.gsd$estimate[est.bmi.gsd$mr_methods=="PRESSO"] =  res$`Main MR results`[1,"Causal Estimate"]
    est.bmi.gsd$se.estimate[est.bmi.gsd$mr_methods=="PRESSO"] = res$`Main MR results`[1,"Sd"]
  }
  est.bmi.gsd$time[est.bmi.gsd$mr_methods=="PRESSO"] = T1-T0
  rm(res)
  # 5. robust
  T0 = proc.time()[3]
  res = MendelianRandomization::mr_ivw(mr.obj,"random", robust = TRUE)
  T1 = proc.time()[3]
  est.bmi.gsd$estimate[est.bmi.gsd$mr_methods=="Robust"] = res$Estimate
  est.bmi.gsd$se.estimate[est.bmi.gsd$mr_methods=="Robust"] = res$StdError
  est.bmi.gsd$time[est.bmi.gsd$mr_methods=="Robust"] = T1-T0
  rm(res)
  # 6. MR-Lasso
  T0 = proc.time()[3]
  res = MR_lasso(betaYG = my_data_harm$beta.gsd, betaXG = my_data_harm$beta.bmi, sebetaYG = my_data_harm$se.gsd)
  T1 = proc.time()[3]
  est.bmi.gsd$estimate[est.bmi.gsd$mr_methods=="Lasso"] = res$ThetaEstimate
  est.bmi.gsd$se.estimate[est.bmi.gsd$mr_methods=="Lasso"] = res$ThetaSE
  est.bmi.gsd$time[est.bmi.gsd$mr_methods=="Lasso"] = T1-T0
  rm(res)
  # 7. Egger
  T0 = proc.time()[3]
  res = MendelianRandomization::mr_egger(mr.obj)
  T1 = proc.time()[3]
  est.bmi.gsd$estimate[est.bmi.gsd$mr_methods=="egger"] = res$Estimate
  est.bmi.gsd$se.estimate[est.bmi.gsd$mr_methods=="egger"] = res$StdError.Est
  est.bmi.gsd$time[est.bmi.gsd$mr_methods=="egger"] = T1-T0
  rm(res)
  # 8. contamination mixture
  T0 = proc.time()[3]
  res = MendelianRandomization::mr_conmix(mr.obj)
  T1 = proc.time()[3]
  est.bmi.gsd$estimate[est.bmi.gsd$mr_methods=="conmix"] = res$Estimate
  CIlength = res$CIUpper-res$CILower
  if (length(CIlength)>1) print("conmix multimodal")
  est.bmi.gsd$se.estimate[est.bmi.gsd$mr_methods=="conmix"] = sum(CIlength)/1.96/2 ## Caution: this may be problematic (Qi&Chatterjee2021)
  est.bmi.gsd$time[est.bmi.gsd$mr_methods=="conmix"] = T1-T0
  rm(res)
  # 9. MRMix
  # theta_temp_vec = seq(-0.5,0.5,by=0.01)
  T0 = proc.time()[3]
  res = MRMix::MRMix(betahat_x = my_data_harm$beta.bmi, 
                     betahat_y = my_data_harm$beta.gsd, 
                     sx =my_data_harm$se.bmi, 
                     sy = my_data_harm$se.gsd)
  res_se = MRMix::MRMix_se(betahat_x = my_data_harm$beta.bmi, 
                           betahat_y = my_data_harm$beta.gsd, 
                           sx =my_data_harm$se.bmi, 
                           sy = my_data_harm$se.gsd, 
                           theta = res$theta, # estimate of causal effect, assuming the summary stats are standardized, theta represents increase in mean value of Y in s.d. unit of Y (for cont outcomes) or log-OR of Y (for binary outcomes)
                           pi0 = res$pi0, 
                           sigma2 = res$sigma2
  )
  T1 = proc.time()[3]
  est.bmi.gsd$estimate[est.bmi.gsd$mr_methods=="MRMix"] = res$theta
  est.bmi.gsd$se.estimate[est.bmi.gsd$mr_methods=="MRMix"] = res_se
  est.bmi.gsd$time[est.bmi.gsd$mr_methods=="MRMix"] = T1-T0
  rm(res, res_se)
  # 10. MR-RAPS  
  ## Warning: The estimated overdispersion parameter is very small. Consider using the simple model without overdispersion.
  T0 = proc.time()[3]
  res = mr.raps::mr.raps.overdispersed.robust(b_exp = my_data_harm$beta.bmi, 
                                              b_out = my_data_harm$beta.gsd, 
                                              se_exp = my_data_harm$se.bmi, 
                                              se_out = my_data_harm$se.gsd, 
                                              loss.function = "huber", 
                                              k = 1.345, 
                                              initialization = c("l2"), 
                                              suppress.warning = FALSE, 
                                              diagnosis = FALSE, 
                                              niter = 20, 
                                              tol = .Machine$double.eps^0.5)
  T1 = proc.time()[3]
  est.bmi.gsd$estimate[est.bmi.gsd$mr_methods=="RAPS"] = res$beta.hat
  est.bmi.gsd$se.estimate[est.bmi.gsd$mr_methods=="RAPS"] = res$beta.se
  est.bmi.gsd$time[est.bmi.gsd$mr_methods=="RAPS"] = T1-T0
  rm(res)
}

## SAVE estimates
saveRDS(est.bmi.gsd, paste0("./output/Rdata/", Sys.Date(), "_est_BMI_GSD.rds"))
