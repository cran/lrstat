% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{zstatRiskRatio}
\alias{zstatRiskRatio}
\title{Miettinen-Nurminen Score Test Statistic for Two-Sample Risk Ratio}
\usage{
zstatRiskRatio(
  riskRatioH0 = 1,
  n1 = NA_real_,
  y1 = NA_real_,
  n2 = NA_real_,
  y2 = NA_real_
)
}
\arguments{
\item{riskRatioH0}{The risk ratio under the null hypothesis.
Defaults to 1.}

\item{n1}{The sample size for the active treatment group.}

\item{y1}{The number of responses for the active treatment group.}

\item{n2}{The sample size for the control group.}

\item{y2}{The number of responses for the control group.}
}
\value{
The value of the score test statistic.
}
\description{
Obtains the Miettinen-Nurminen score test statistic for
two-sample risk ratio possibly with stratification.
}
\details{
The Mantel-Haenszel sample size weights are used for stratified
samples.
}
\examples{

zstatRiskRatio(riskRatioH0 = 1, n1 = c(10,10), y1 = c(4,3),
               n2 = c(20,10), y2 = c(2,0))

}
\author{
Kaifeng Lu, \email{kaifenglu@gmail.com}
}
