% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getDesignSurvivals.R
\name{lrschoenfeld}
\alias{lrschoenfeld}
\title{Schoenfeld Method for Log-Rank Test Sample Size Calculation}
\usage{
lrschoenfeld(
  beta = 0.2,
  kMax = 1L,
  informationRates = NA_real_,
  efficacyStopping = NA_integer_,
  futilityStopping = NA_integer_,
  criticalValues = NA_real_,
  alpha = 0.025,
  typeAlphaSpending = "sfOF",
  parameterAlphaSpending = NA_real_,
  userAlphaSpending = NA_real_,
  futilityBounds = NA_real_,
  typeBetaSpending = "none",
  parameterBetaSpending = NA_real_,
  userBetaSpending = NA_real_,
  hazardRatioH0 = 1,
  allocationRatioPlanned = 1,
  accrualTime = 0L,
  accrualIntensity = NA_real_,
  piecewiseSurvivalTime = 0L,
  stratumFraction = 1L,
  hazardRatio = NA_real_,
  lambda2 = NA_real_,
  gamma1 = 0L,
  gamma2 = 0L,
  followupTime = NA_real_,
  fixedFollowup = 0L,
  interval = as.numeric(c(0.001, 240)),
  spendingTime = NA_real_,
  rounding = 1L,
  calibrate = 1L,
  maxNumberOfIterations = 10000L,
  maxNumberOfRawDatasetsPerStage = 0L,
  seed = NA_integer_
)
}
\arguments{
\item{beta}{Type II error. Defaults to 0.2.}

\item{kMax}{The maximum number of stages.}

\item{informationRates}{The information rates in terms of number
of events for the conventional log-rank test and in terms of
the actual information for weighted log-rank tests.
Defaults to \code{(1:kMax) / kMax} if left unspecified.}

\item{efficacyStopping}{Indicators of whether efficacy stopping is allowed
at each stage. Defaults to true if left unspecified.}

\item{futilityStopping}{Indicators of whether futility stopping is allowed
at each stage. Defaults to true if left unspecified.}

\item{criticalValues}{Upper boundaries on the z-test statistic scale
for stopping for efficacy.}

\item{alpha}{The significance level. Defaults to 0.025.}

\item{typeAlphaSpending}{The type of alpha spending. One of the following:
"OF" for O'Brien-Fleming boundaries, "P" for Pocock boundaries,
"WT" for Wang & Tsiatis boundaries, "sfOF" for O'Brien-Fleming type
spending function, "sfP" for Pocock type spending function,
"sfKD" for Kim & DeMets spending function, "sfHSD" for Hwang,
Shi & DeCani spending function, "user" for user defined spending,
and "none" for no early efficacy stopping. Defaults to "sfOF".}

\item{parameterAlphaSpending}{The parameter value for the alpha spending.
Corresponds to Delta for "WT", rho for "sfKD", and gamma for "sfHSD".}

\item{userAlphaSpending}{The user defined alpha spending. Cumulative alpha
spent up to each stage.}

\item{futilityBounds}{Lower boundaries on the z-test statistic scale
for stopping for futility at stages 1, ..., \code{kMax-1}. Defaults to
\code{rep(-6, kMax-1)} if left unspecified. The futility bounds are
non-binding for the calculation of critical values.}

\item{typeBetaSpending}{The type of beta spending. One of the following:
"sfOF" for O'Brien-Fleming type spending function, "sfP" for Pocock type
spending function, "sfKD" for Kim & DeMets spending function,
"sfHSD" for Hwang, Shi & DeCani spending function, "user" for
user defined spending, and "none" for no early futility stopping.
Defaults to "none".}

\item{parameterBetaSpending}{The parameter value for the beta spending.
Corresponds to rho for "sfKD", and gamma for "sfHSD".}

\item{userBetaSpending}{The user defined beta spending. Cumulative beta
spent up to each stage.}

\item{hazardRatioH0}{Hazard ratio under the null hypothesis for the active
treatment versus control. Defaults to 1 for superiority test.}

\item{allocationRatioPlanned}{Allocation ratio for the active treatment
versus control. Defaults to 1 for equal randomization.}

\item{accrualTime}{A vector that specifies the starting time of
piecewise Poisson enrollment time intervals. Must start with 0, e.g.,
\code{c(0, 3)} breaks the time axis into 2 accrual intervals:
[0, 3) and [3, Inf).}

\item{accrualIntensity}{A vector of accrual intensities. One for
each accrual time interval.}

\item{piecewiseSurvivalTime}{A vector that specifies the starting time of
piecewise exponential survival time intervals. Must start with 0, e.g.,
\code{c(0, 6)} breaks the time axis into 2 event intervals:
[0, 6) and [6, Inf).
Defaults to 0 for exponential distribution.}

\item{stratumFraction}{A vector of stratum fractions that sum to 1.
Defaults to 1 for no stratification.}

\item{hazardRatio}{Hazard ratio under the alternative hypothesis for
the active treatment versus control.}

\item{lambda2}{A vector of hazard rates for the event in each analysis
time interval by stratum for the control group.}

\item{gamma1}{The hazard rate for exponential dropout, a vector of hazard
rates for piecewise exponential dropout applicable for all strata, or a
vector of hazard rates for dropout in each analysis time interval
by stratum for the active treatment group.}

\item{gamma2}{The hazard rate for exponential dropout, a vector of hazard
rates for piecewise exponential dropout applicable for all strata, or a
vector of hazard rates for dropout in each analysis time interval
by stratum for the control group.}

\item{followupTime}{Follow-up time for the last enrolled subject.}

\item{fixedFollowup}{Whether a fixed follow-up design is used.
Defaults to 0 for variable follow-up.}

\item{interval}{The interval to search for the solution of
followupTime. Defaults to \code{c(0.001, 240)}.}

\item{spendingTime}{A vector of length \code{kMax} for the error spending
time at each analysis. Defaults to missing, in which case, it is the
same as \code{informationRates}.}

\item{rounding}{Whether to round up sample size and events.
Defaults to 1 for sample size rounding.}

\item{calibrate}{Whether to use simulations to calibrate the number of
events calculated using the Schoenfeld formula.}

\item{maxNumberOfIterations}{The number of simulation iterations.
Defaults to 10000.}

\item{maxNumberOfRawDatasetsPerStage}{The number of raw datasets per
stage to extract.}

\item{seed}{The seed to reproduce the simulation results.
The seed from the environment will be used if left unspecified.}
}
\value{
A list of two components:
\itemize{
\item \code{analyticalResults}: An S3 class \code{lrpower} object for
the asymptotic power.
\item \code{simulationResults}: An S3 class \code{lrsim} object for
the empirical power.
}
}
\description{
Obtains the sample size and study duration by calibrating
the number of events calculated using the Schoenfeld formula
under the proportional hazards assumption.
}
\details{
This function calculates the sample size and study duration
by calibrating the number of events estimated using the
Schoenfeld formula under the proportional hazards assumption,
particularly when the hazard ratio is far away from one and/or
the allocation between groups is unequal.

For a fixed design, the Schoenfeld formula for the required
number of events is
\deqn{D = \frac{(\Phi^{-1}(1-\alpha) + \Phi^{-1}(1-\beta))^2}
{(\theta - \theta_0)^2 r(1-r)}}
where \eqn{D} is the total number of events required,
\eqn{\alpha} is the type I error rate,
\eqn{\beta} is the type II error rate,
\eqn{r} is the randomization probability for the active treatment group,
\eqn{\theta_0} and \eqn{\theta} are the log hazard ratios under
the null and alternative hypotheses, respectively.

The function first computes the number of events using the
Schoenfeld formula. If \code{calibrate} is set to 1, the
function uses simulations to calibrate the number of
events, accounting for scenarios where the Schoenfeld formula
may be inaccurate (e.g., when allocation is unequal or the hazard
ratio is extreme).

Let \eqn{D_{schoenfeld}} be the number of events calculated
by the Schoenfeld formula, and \eqn{D_{calibrated}}
be the calibrated number of events. The calibrated number of
events is calculated as
#' \deqn{D_{\text{calibrated}} =
\frac{\left\{\Phi^{-1}(1-\alpha) + \Phi^{-1}(1-\beta)\right\}^2}
{\left\{\Phi^{-1}(1-\alpha) +
\Phi^{-1}(1-\beta_{\text{schoenfeld}})\right\}^2}
D_{\text{schoenfeld}}}
where \eqn{\beta_{schoenfeld}} is the empirical type II error
estimated via simulation.

A second round of simulation is performed to obtain the
empirical power using the calibrated number of events.
}
\examples{

(lr1 <- lrschoenfeld(
  beta = 0.1, kMax = 2, alpha = 0.025,
  hazardRatioH0 = 1, allocationRatioPlanned = 1,
  accrualIntensity = 20, hazardRatio = 0.3,
  lambda2 = 1.9/12,
  gamma1 = -log(1-0.1)/24, gamma2 = -log(1-0.1)/24,
  fixedFollowup = 0, rounding = 1,
  calibrate = 0, maxNumberOfIterations = 1000,
  seed = 12345))

(lr2 <- lrschoenfeld(
  beta = 0.1, kMax = 2, alpha = 0.025,
  hazardRatioH0 = 1, allocationRatioPlanned = 1,
  accrualIntensity = 20, hazardRatio = 0.3,
  lambda2 = 1.9/12,
  gamma1 = -log(1-0.1)/24, gamma2 = -log(1-0.1)/24,
  fixedFollowup = 0, rounding = 1,
  calibrate = 1, maxNumberOfIterations = 1000,
  seed = 12345))

}
\author{
Kaifeng Lu, \email{kaifenglu@gmail.com}
}
