% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getDesignProportions.R
\name{getDesignTwoOrdinal}
\alias{getDesignTwoOrdinal}
\title{Power and Sample Size for the Wilcoxon Test for Two-Sample
Ordinal Response}
\usage{
getDesignTwoOrdinal(
  beta = NA_real_,
  n = NA_real_,
  ncats = NA_integer_,
  pi1 = NA_real_,
  pi2 = NA_real_,
  allocationRatioPlanned = 1,
  rounding = TRUE,
  alpha = 0.05
)
}
\arguments{
\item{beta}{The type II error.}

\item{n}{The total sample size.}

\item{ncats}{The number of categories of the ordinal response.}

\item{pi1}{The prevalence of each category for the treatment group.
Only need to specify the valued for the first \code{ncats-1} categories.}

\item{pi2}{The prevalence of each category for the control group.
Only need to specify the valued for the first \code{ncats-1} categories.}

\item{allocationRatioPlanned}{Allocation ratio for the active treatment
versus control. Defaults to 1 for equal randomization.}

\item{rounding}{Whether to round up sample size. Defaults to 1 for
sample size rounding.}

\item{alpha}{The significance level. Defaults to 0.025.}
}
\value{
An S3 class \code{designTwoOrdinal} object with the following
components:
\itemize{
\item \code{power}: The power to reject the null hypothesis.
\item \code{alpha}: The two-sided significance level.
\item \code{n}: The maximum number of subjects.
\item \code{ncats}: The number of categories of the ordinal response.
\item \code{pi1}: The prevalence of each category for the treatment group.
\item \code{pi2}: The prevalence of each category for the control group.
\item \code{meanscore1}: The mean midrank score for the treatment group.
\item \code{meanscore2}: The mean midrank score for the control group.
\item \code{allocationRatioPlanned}: Allocation ratio for the active treatment
versus control.
\item \code{rounding}: Whether to round up sample size.
}
}
\description{
Obtains the power given sample size or obtains the sample
size given power for the Wilcoxon test for two-sample ordinal response.
}
\details{
A two-sample ordinal response design is used to test whether the
ordinal response distributions differ between two treatment arms.
Let \eqn{\pi_{gi}} denote the prevalence of category \eqn{i}
in group \eqn{g}, where \eqn{g=1} represents the treatment group and
\eqn{g=2} represents the control group.

The parameter of interest is
\deqn{\theta = \sum_{i=1}^{C} w_i (\pi_{1i} - \pi_{2i})}
where \eqn{w_i} is the midrank score for category \eqn{i}.
The Z-test statistic is given by
\deqn{Z = \hat{\theta}/\sqrt{Var(\hat{\theta})}}
where \eqn{\hat{\theta}} is the estimate of \eqn{\theta}.

The midrank score \eqn{w_i} for category \eqn{i} is calculated as:
\deqn{w_i = \sum_{j=1}^{i} \pi_j - 0.5\pi_i}
where \eqn{\pi_i = r\pi_{1i} + (1-r)\pi_{2i}} denotes the average
prevalence of category \eqn{i} across both groups,
and \eqn{r} is the randomization probability for the active treatment.

To understand the midrank score, consider \eqn{n\pi_i} subjects
in category \eqn{i}. The midrank score is the average rank of
these subjects:
\deqn{s_i = \frac{1}{n\pi_i} \sum_{j=1}^{n\pi_i} (
n\pi_1 + \cdots + n\pi_{i-1} + j)}
This simplifies to
\deqn{s_i = n\left(\sum_{j=1}^{i} \pi_j - 0.5\pi_i\right) + \frac{1}{2}}
By dividing by \eqn{n} and ignoring \eqn{\frac{1}{2n}}, we obtain
the midrank score \eqn{w_i}.

The variance of \eqn{\hat{\theta}} can be derived from the multinomial
distributions and is given by
\deqn{Var(\hat{\theta}) = \frac{1}{n}\sum_{g=1}^{2} \frac{1}{r_g}
\left\{\sum_{i=1}^{C} w_i^2\pi_{gi} - \left(\sum_{i=1}^{C} w_i\pi_{gi}
\right)^2\right\}}
where \eqn{r_g} is the randomization probability for group \eqn{g}.

The sample size is chosen such that the power to reject the null
hypothesis is at least \eqn{1-\beta} for a given
significance level \eqn{\alpha}.
}
\examples{

(design1 <- getDesignTwoOrdinal(
  beta = 0.1, ncats = 4, pi1 = c(0.55, 0.3, 0.1),
  pi2 = c(0.214, 0.344, 0.251), alpha = 0.025))

}
\author{
Kaifeng Lu, \email{kaifenglu@gmail.com}
}
