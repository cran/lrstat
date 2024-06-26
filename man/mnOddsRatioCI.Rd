% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{mnOddsRatioCI}
\alias{mnOddsRatioCI}
\title{Miettinen-Nurminen score confidence interval for
two-sample odds ratio}
\usage{
mnOddsRatioCI(
  n1 = NA_real_,
  y1 = NA_real_,
  n2 = NA_real_,
  y2 = NA_real_,
  cilevel = 0.95
)
}
\arguments{
\item{n1}{The sample size for the active treatment group.}

\item{y1}{The number of responses for the active treatment group.}

\item{n2}{The sample size for the control group.}

\item{y2}{The number of responses for the control group.}

\item{cilevel}{The confidence interval level.}
}
\value{
A list with two components:
\itemize{
\item \code{data} A data frame containing the input sample size
and number of responses for each treatment group.
It has the following variables:
\itemize{
\item \code{n1}: The sample size for the active treatment group.
\item \code{y1}: The number of responses for the active treatment group.
\item \code{n2}: The sample size for the control group.
\item \code{y2}: The number of responses for the control group.
}
\item \code{estimates}: A data frame containing the point estimate
and confidence interval for odds ratio. It has the following
variables:
\itemize{
\item \code{scale}: The scale of treatment effect.
\item \code{estimate}: The point estimate.
\item \code{lower}: The lower limit of the confidence interval.
\item \code{upper}: The upper limit of the confidence interval.
\item \code{cilevel}: The confidence interval level.
}
}
}
\description{
Obtains the Miettinen-Nurminen score confidence
interval for two-sample odds ratio possibly with
stratification.
}
\details{
The Mantel-Haenszel sample size weights are used for stratified
samples.
}
\examples{

mnOddsRatioCI(n1 = c(10,10), y1 = c(4,3), n2 = c(20,10), y2 = c(2,0))

}
\author{
Kaifeng Lu, \email{kaifenglu@gmail.com}
}
