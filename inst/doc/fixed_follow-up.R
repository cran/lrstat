## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(lrstat)

## -----------------------------------------------------------------------------
lrsamplesize(beta = 0.1, kMax = 1, criticalValues = 1.96, 
             allocationRatioPlanned = 3, accrualIntensity = 5, 
             lambda2 = 0.95/12, lambda1 = 0.3*0.95/12, 
             gamma1 = -log(1-0.1)/24, gamma2 = -log(1-0.1)/24, 
             accrualDuration = NA, followupTime = 26/4, 
             fixedFollowup = TRUE,
             typeOfComputation = "schoenfeld")

## -----------------------------------------------------------------------------
lrsamplesize(beta = 0.1, kMax = 1, criticalValues = 1.96, 
             allocationRatioPlanned = 3, accrualIntensity = 5, 
             lambda2 = 0.95/12, lambda1 = 0.3*0.95/12, 
             gamma1 = -log(1-0.1)/24, gamma2 = -log(1-0.1)/24, 
             accrualDuration = NA, followupTime = 26/4, 
             fixedFollowup = TRUE,
             typeOfComputation = "direct")

## -----------------------------------------------------------------------------
lrsim(kMax = 1, criticalValues = 1.96,  
      allocation1 = 3, allocation2 = 1,
      accrualIntensity = 5, 
      lambda2 = 0.95/12, lambda1 = 0.3*0.95/12, 
      gamma1 = -log(1-0.1)/24, gamma2 = -log(1-0.1)/24,
      n = 191, followupTime = 6.5, 
      fixedFollowup = TRUE,  
      plannedEvents = 39, 
      maxNumberOfIterations = 10000, seed = 12345)

lrsim(kMax = 1, criticalValues = 1.96,  
      allocation1 = 3, allocation2 = 1,
      accrualIntensity = 5, 
      lambda2 = 0.95/12, lambda1 = 0.3*0.95/12, 
      gamma1 = -log(1-0.1)/24, gamma2 = -log(1-0.1)/24,
      n = 127, followupTime = 6.5, 
      fixedFollowup = TRUE,  
      plannedEvents = 26, 
      maxNumberOfIterations = 10000, seed = 12345)

## -----------------------------------------------------------------------------
lrsim(kMax = 1, criticalValues = 1.96,  
      allocation1 = 3, allocation2 = 1,
      accrualIntensity = 5, 
      lambda2 = 0.95/12, lambda1 = 0.3*0.95/12, 
      gamma1 = -log(1-0.1)/24, gamma2 = -log(1-0.1)/24,
      n = 156, followupTime = 6.5, 
      fixedFollowup = TRUE,  
      plannedEvents = 32, 
      maxNumberOfIterations = 10000, seed = 12345)

