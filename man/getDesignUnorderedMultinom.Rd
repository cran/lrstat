% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getDesignProportions.R
\name{getDesignUnorderedMultinom}
\alias{getDesignUnorderedMultinom}
\title{Power and sample size for unordered multi-sample
multinomial response}
\usage{
getDesignUnorderedMultinom(
  beta = NA_real_,
  n = NA_real_,
  ngroups = NA_integer_,
  ncats = NA_integer_,
  pi = NA_real_,
  allocationRatioPlanned = NA_integer_,
  rounding = TRUE,
  alpha = 0.05
)
}
\arguments{
\item{beta}{The type II error.}

\item{n}{The total sample size.}

\item{ngroups}{The number of treatment groups.}

\item{ncats}{The number of categories of the multinomial response.}

\item{pi}{The matrix of response probabilities for the treatment groups.
It should have \code{ngroups} rows and \code{ncats-1} or \code{ncats}
columns.}

\item{allocationRatioPlanned}{Allocation ratio for the treatment groups.}

\item{rounding}{Whether to round up sample size. Defaults to 1 for
sample size rounding.}

\item{alpha}{The two-sided significance level. Defaults to 0.05.}
}
\value{
An S3 class \code{designUnorderedMultinom} object with the
following components:
\itemize{
\item \code{power}: The power to reject the null hypothesis.
\item \code{alpha}: The two-sided significance level.
\item \code{n}: The maximum number of subjects.
\item \code{ngroups}: The number of treatment groups.
\item \code{ncats}: The number of categories of the multinomial response.
\item \code{pi}: The response probabilities for the treatment groups.
\item \code{effectsize}: The effect size for the chi-square test.
\item \code{allocationRatioPlanned}: Allocation ratio for the treatment
groups.
\item \code{rounding}: Whether to round up sample size.
}
}
\description{
Obtains the power given sample size or obtains the sample
size given power for the chi-square test for unordered multi-sample
multinomial response.
}
\examples{

(design1 <- getDesignUnorderedMultinom(
  beta = 0.1, ngroups = 3, ncats = 4,
  pi = matrix(c(0.230, 0.320, 0.272,
                0.358, 0.442, 0.154,
                0.142, 0.036, 0.039),
              3, 3, byrow = TRUE),
  allocationRatioPlanned = c(2, 2, 1),
  alpha = 0.05))

}
\author{
Kaifeng Lu, \email{kaifenglu@gmail.com}
}
