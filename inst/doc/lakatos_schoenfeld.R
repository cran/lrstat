## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(lrstat)
library(rpact)

## -----------------------------------------------------------------------------
caltime(nevents = 40, accrualDuration = 12, accrualIntensity = 200/12,
        lambda1 = -log(1-0.2)/12, lambda2 = -log(1-0.4)/12, 
        followupTime = 100)

## -----------------------------------------------------------------------------
lrpower(kMax = 1, criticalValues = 1.96, accrualDuration = 12, 
        accrualIntensity = 200/12, lambda1 = -log(1-0.2)/12, 
        lambda2 = -log(1-0.4)/12,  followupTime = 1.63)

## -----------------------------------------------------------------------------
getPowerSurvival(alpha = 0.025, directionUpper = FALSE, pi1 = 0.2, 
                 pi2 = 0.4, eventTime = 12, accrualTime = c(0, 12),
                 maxNumberOfSubjects = 200, maxNumberOfEvents = 40)

## -----------------------------------------------------------------------------
getSimulationSurvival(directionUpper = FALSE, pi1 = 0.2, pi2 = 0.4, 
                      eventTime = 12, accrualTime = 12, plannedEvents = 40,
                      maxNumberOfSubjects = 200, 
                      maxNumberOfIterations = 10000, seed = 12345)

