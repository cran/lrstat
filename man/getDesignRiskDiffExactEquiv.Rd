% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getDesignProportions.R
\name{getDesignRiskDiffExactEquiv}
\alias{getDesignRiskDiffExactEquiv}
\title{Power and sample size for exact unconditional test for
equivalence in risk difference}
\usage{
getDesignRiskDiffExactEquiv(
  beta = NA_real_,
  n = NA_real_,
  riskDiffLower = NA_real_,
  riskDiffUpper = NA_real_,
  pi1 = NA_real_,
  pi2 = NA_real_,
  allocationRatioPlanned = 1,
  alpha = 0.05
)
}
\arguments{
\item{beta}{The type II error.}

\item{n}{The total sample size.}

\item{riskDiffLower}{The lower equivalence limit of risk difference.}

\item{riskDiffUpper}{The upper equivalence limit of risk difference.}

\item{pi1}{The assumed probability for the active treatment group.}

\item{pi2}{The assumed probability for the control group.}

\item{allocationRatioPlanned}{Allocation ratio for the active treatment
versus control. Defaults to 1 for equal randomization.}

\item{alpha}{The significance level for each of the two one-sided
tests. Defaults to 0.05.}
}
\value{
A data frame with the following variables:
\itemize{
\item \code{alpha}: The specified significance level for each of the two
one-sided tests.
\item \code{attainedAlpha}: The attained significance level.
\item \code{power}: The power.
\item \code{n}: The sample size.
\item \code{riskDiffLower}: The lower equivalence limit of risk difference.
\item \code{riskDiffUpper}: The upper equivalence limit of risk difference.
\item \code{pi1}: The assumed probability for the active treatment group.
\item \code{pi2}: The assumed probability for the control group.
\item \code{riskDiff}: The risk difference.
\item \code{allocationRatioPlanned}: Allocation ratio for the active
treatment versus control.
\item \code{zstatRiskDiffLower}: The efficacy boundaries on the
z-test statistic scale for the one-sided null hypothesis on the
lower equivalence limit.
\item \code{zstatRiskDiffUpper}: The efficacy boundaries on the
z-test statistic scale for the one-sided null hypothesis on the
upper equivalence limit.
}
}
\description{
Obtains the power given sample size or obtains the sample
size given power for exact unconditional test of equivalence in risk
difference.
}
\examples{

getDesignRiskDiffExactEquiv(
  n = 200, riskDiffLower = -0.2, riskDiffUpper = 0.2,
  pi1 = 0.775, pi2 = 0.775, alpha = 0.05)

}
\author{
Kaifeng Lu, \email{kaifenglu@gmail.com}
}
