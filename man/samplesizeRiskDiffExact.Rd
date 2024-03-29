% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{samplesizeRiskDiffExact}
\alias{samplesizeRiskDiffExact}
\title{Sample size for exact unconditional test of risk difference}
\usage{
samplesizeRiskDiffExact(
  beta = NA_real_,
  riskDiffH0 = 0,
  pi1 = NA_real_,
  pi2 = NA_real_,
  allocationRatioPlanned = 1,
  alpha = 0.025
)
}
\arguments{
\item{beta}{The type II error.}

\item{riskDiffH0}{The risk difference under the null hypothesis.
Defaults to 0.}

\item{pi1}{The assumed probability for the active treatment group.}

\item{pi2}{The assumed probability for the control group.}

\item{allocationRatioPlanned}{Allocation ratio for the active treatment
versus control. Defaults to 1 for equal randomization.}

\item{alpha}{The one-sided significance level.}
}
\value{
A data frame with the following variables:
\itemize{
\item \code{alpha}: The specified one-sided significance level.
\item \code{attainedAlpha}: The attained one-sided significance level.
\item \code{power}: The power.
\item \code{n}: The sample size.
\item \code{riskDiffH0}: The risk difference under the null hypothesis.
\item \code{pi1}: The assumed probability for the active treatment group.
\item \code{pi2}: The assumed probability for the control group.
\item \code{allocationRatioPlanned}: Allocation ratio for the active
treatment versus control.
\item \code{zstatRiskDiffBound}: The critical value on the scale of
score test statistic for risk difference.
\item \code{pi2star}: The response probability in the control group
at which the critical value of the test statistic is attained.
}
}
\description{
Obtains the sample size given power for exact unconditional
test of risk difference.
}
\examples{

samplesizeRiskDiffExact(beta = 0.2, riskDiffH0 = -0.2,
                        pi1 = 0.8, pi2 = 0.8, alpha = 0.025)

}
\author{
Kaifeng Lu, \email{kaifenglu@gmail.com}
}
\keyword{internal}
