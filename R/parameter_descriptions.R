#' Parameter Description: accrualTime
#' @param accrualTime Accrual time intervals, must start with 0, e.g.,
#'   \code{c(0, 3)} breaks the time axis into 2 accrual intervals:
#'   [0, 3) and [3, Inf).
#' @name param_accrualTime
#' @keywords internal
NULL

#' Parameter Description: accrualIntensity
#' @param accrualIntensity A vector of accrual intensities, one for
#'   each accrual time interval.
#' @name param_accrualIntensity
#' @keywords internal
NULL

#' Parameter Description: accrualDuration
#' @param accrualDuration Duration of the enrollment period.
#' @name param_accrualDuration
#' @keywords internal
NULL

#' Parameter Description: piecewiseSurvivalTime
#' @param piecewiseSurvivalTime A vector that specifies the time intervals for
#'   the piecewise exponential survival distribution, must start with 0, e.g.,
#'   \code{c(0, 6)} breaks the time axis into 2 event intervals:
#'   [0, 6) and [6, Inf).
#'   Defaults to 0 for exponential distribution.
#' @name param_piecewiseSurvivalTime
#' @keywords internal
NULL

#' Parameter Description: allocationRatioPlanned
#' @param allocationRatioPlanned Allocation ratio for the active treatment
#'   versus control. Defaults to 1 for equal randomization.
#' @name param_allocationRatioPlanned
#' @keywords internal
NULL

#' Parameter Description: stratumFraction
#' @param stratumFraction A vector of stratum fractions that sum to 1.
#'   Defaults to 1 for no stratification.
#' @name param_stratumFraction
#' @keywords internal
NULL

#' Parameter Description: lambda
#' @param lambda A vector of hazard rates for the event, one for
#'   each analysis time interval.
#' @name param_lambda
#' @keywords internal
NULL

#' Parameter Description: lambda1
#' @param lambda1 A vector of hazard rates for the event for the
#'   active treatment group, one for each analysis time interval, by stratum.
#' @name param_lambda1
#' @keywords internal
NULL

#' Parameter Description: lambda2
#' @param lambda2 A vector of hazard rates for the event for the
#'   control group, one for each analysis time interval, by stratum.
#' @name param_lambda2
#' @keywords internal
NULL

#' Parameter Description: lambda1 stratified
#' @param lambda1 A vector of hazard rates for the event for the
#'   active treatment group, one for each analysis time interval, by stratum.
#' @name param_lambda1_stratified
#' @keywords internal
NULL

#' Parameter Description: lambda2_stratified
#' @param lambda2 A vector of hazard rates for the event for the
#'   control group, one for each analysis time interval, by stratum.
#' @name param_lambda2_stratified
#' @keywords internal
NULL

#' Parameter Description: gamma
#' @param gamma The hazard rate for exponential dropout or a vector of hazard
#'   rates for piecewise exponential dropout. Defaults to 0 for no dropout.
#' @name param_gamma
#' @keywords internal
NULL

#' Parameter Description: gamma1
#' @param gamma1 The hazard rate for exponential dropout or a vector of hazard
#'   rates for piecewise exponential dropout for the active treatment group.
#'   Defaults to 0 for no dropout.
#' @name param_gamma1
#' @keywords internal
NULL

#' Parameter Description: gamma2
#' @param gamma2 The hazard rate for exponential dropout or a vector of hazard
#'   rates for piecewise exponential dropout for the control group.
#'   Defaults to 0 for no dropout.
#' @name param_gamma2
#' @keywords internal
NULL

#' Parameter Description: followupTime
#' @param followupTime Follow-up time for the last enrolled subject.
#' @name param_followupTime
#' @keywords internal
NULL

#' Parameter Description: fixedFollowup
#' @param fixedFollowup Whether a fixed follow-up design is used.
#'   Defaults to 0 for variable follow-up.
#' @name param_fixedFollowup
#' @keywords internal
NULL

#' Parameter Description: minFollowupTime
#' @param minFollowupTime Follow-up time for the last enrolled subject.
#' @name param_minFollowupTime
#' @keywords internal
NULL

#' Parameter Description: maxFollowupTime
#' @param maxFollowupTime Follow-up time for the first enrolled subject.
#'   For fixed followup, \code{maxFollowupTime = minFollowupTime}.
#'   For variable followup,
#'   \code{maxFollowupTime = accrualDuration + minFollowupTime}.
#' @name param_maxFollowupTime
#' @keywords internal
NULL

#' Parameter Description: rho1
#' @param rho1 First parameter of the Fleming-Harrington family of weighted
#'   log-rank test. Defaults to 0 for conventional log-rank test.
#' @name param_rho1
#' @keywords internal
NULL

#' Parameter Description: rho2
#' @param rho2 Second parameter of the Fleming-Harrington family of weighted
#'   log-rank test. Defaults to 0 for conventional log-rank test.
#' @name param_rho2
#' @keywords internal
NULL

#' Parameter Description: numSubintervals
#' @param numSubintervals Number of sub-intervals to approximate the mean
#'   and variance of the weighted log-rank test score statistic.
#'   Defaults to 300. Specify a larger number for better approximation.
#' @name param_numSubintervals
#' @keywords internal
NULL

#' Parameter Description: kMax
#' @param kMax The maximum number of stages.
#' @name param_kMax
#' @keywords internal
NULL

#' Parameter Description: informationRates
#' @param informationRates The information rates fixed prior to the trial.
#'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
#' @name param_informationRates
#' @keywords internal
NULL

#' Parameter Description: criticalValues
#' @param criticalValues Upper boundaries on the z-test statistic scale
#'   for stopping for efficacy.
#' @name param_criticalValues
#' @keywords internal
NULL

#' Parameter Description: futilityBounds
#' @param futilityBounds Lower boundaries on the z-test statistic scale
#'   for stopping for futility at stages 1, ..., kMax-1. Defaults to
#'   \code{rep(-Inf, kMax-1)} if left unspecified.
#' @name param_futilityBounds
#' @keywords internal
NULL
