

## ---------------------------
##
## Script name: PVE_calculation_fns.R
##
## Author: Christoph Reich
## Date Created: 2022-11-23
##
## Copyright (c) Christoph Reich, 2022
## Email: christoph.reich@med.uni-heidelberg.de
##
## ---------------------------
## Notes:
##   formulas to calculate the proportion of variance explained for numeric and 
##   binary traits      
## ---------------------------


# PVE FORMULAS ------------------------------------------------------------
explained_variance_numeric <- function(eaf, beta) {
  # eaf = allele frequency
  # beta = additive genetic effect
  explained_variance_numeric <- 2*eaf *(1-eaf) * beta^2
  return(explained_variance_numeric)
}

explained_variance_numeric2 <- function(maf, beta, se_beta, samplesize) {
  ### Teslovich Nature 2010
  # 1) maf = minor allele frequency
  # 2) beta = additive genetic effect
  # 3) se_beta = standard error of beta estimate
  # 4) samplesize = samplesize of study
  maf <- ifelse(maf>0.5, 1-maf, maf)
  
  explained_variance_numeric <- (2*beta^2*maf*(1-maf))/ 
    ( 2*beta^2*maf*(1-maf) + (se_beta^2)*2*samplesize*maf*(1-maf) )
  return(explained_variance_numeric)
}

explained_variance_binary <- function(PA,RR1,RR2,K) {
  # from So et al 2011
  # Calculates the variance in liability explained by a single biallelic loci.
  # PA: allele frequency of the risk allele (denoted A)
  # RR1: relative risk of Aa (one risk allele) compared to aa (no risk allele)
  # RR2: relative risk of AA (two risk alleles) compared to aa (no risk allele)
  # K:  overall probability of disease in population
  # Returns the variance explained (Vg) and the mean liability for each genotype (the overall liability is normalized to mean 0 and variance 1)
  Paa = (1-PA)^2
  PAa = 2*PA*(1-PA)
  PAA = PA^2
  muaa=0
  faa= K/(Paa + PAa*RR1 + PAA*RR2)
  fAa= RR1*faa
  fAA= RR2*faa 
  T = qnorm(1-faa) 
  muAa = T-qnorm(1-fAa)
  muAA = T-qnorm(1-fAA)
  mean.all= PAa*muAa+ PAA*muAA
  Vg= Paa*(muaa-mean.all)^2 + PAa*(muAa-mean.all)^2+ PAA*(muAA-mean.all)^2
  actual.Vg =  Vg/(1+Vg) 
  VR = 1-actual.Vg 
  actual.T = Paa*sqrt(VR)*qnorm(1-faa) + PAa*sqrt(VR)*qnorm(1-fAa) + PAA*sqrt(VR)*qnorm(1-fAA)
  actual.muaa = actual.T - sqrt(VR) * qnorm(1-faa)
  actual.muAa = actual.T - sqrt(VR) * qnorm(1-fAa)
  actual.muAA = actual.T - sqrt(VR) * qnorm(1-fAA)
  
  res <- list(Vg=actual.Vg,muaa=actual.muaa, muAa = actual.muAa, muAA=actual.muAA)
  res
} 
