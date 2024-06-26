% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{powerRiskRatioExactEquiv}
\alias{powerRiskRatioExactEquiv}
\title{Power for exact unconditional test of equivalence in risk
ratio}
\usage{
powerRiskRatioExactEquiv(
  n = NA_integer_,
  riskRatioLower = NA_real_,
  riskRatioUpper = NA_real_,
  pi1 = NA_real_,
  pi2 = NA_real_,
  allocationRatioPlanned = 1,
  alpha = 0.05,
  calculateAttainedAlpha = 1L
)
}
\arguments{
\item{n}{The total sample size.}

\item{riskRatioLower}{The lower equivalence limit of risk ratio.}

\item{riskRatioUpper}{The upper equivalence limit of risk ratio.}

\item{pi1}{The assumed probability for the active treatment group.}

\item{pi2}{The assumed probability for the control group.}

\item{allocationRatioPlanned}{Allocation ratio for the active treatment
versus control. Defaults to 1 for equal randomization.}

\item{alpha}{The significance level for each of the two one-sided
tests. Defaults to 0.05.}

\item{calculateAttainedAlpha}{Whether to calculate the attained alpha.}
}
\value{
A data frame with the following variables:
\itemize{
\item \code{alpha}: The specified significance level for each of the two
one-sided tests.
\item \code{attainedAlphaH10}: The attained significance level under H10
if requested.
\item \code{attainedAlphaH20}: The attained significance level under H20
if requested.
\item \code{power}: The power.
\item \code{n}: The sample size.
\item \code{riskRatioLower}: The lower equivalence limit of risk ratio.
\item \code{riskRatioUpper}: The upper equivalence limit of risk ratio.
\item \code{pi1}: The assumed probability for the active treatment group.
\item \code{pi2}: The assumed probability for the control group.
\item \code{riskRatio}: The risk ratio.
\item \code{allocationRatioPlanned}: Allocation ratio for the active
treatment versus control.
\item \code{zstatRiskRatioLower}: The efficacy boundaries on the
z-test statistic scale for the one-sided null hypothesis on the
lower equivalence limit.
\item \code{zstatRiskRatioUpper}: The efficacy boundaries on the
z-test statistic scale for the one-sided null hypothesis on the
upper equivalence limit.
}
}
\description{
Obtains the power given sample size for exact unconditional
test of equivalence in risk ratio.
}
\examples{

powerRiskRatioExactEquiv(
  n = 200, riskRatioLower = 0.8, riskRatioUpper = 1.25,
  pi1 = 0.775, pi2 = 0.775, alpha = 0.05)

}
\author{
Kaifeng Lu, \email{kaifenglu@gmail.com}
}
\keyword{internal}
