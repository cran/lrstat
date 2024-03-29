% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{nbstat}
\alias{nbstat}
\title{Negative binomial rate ratio}
\usage{
nbstat(
  time = NA_real_,
  rateRatioH0 = 1,
  allocationRatioPlanned = 1,
  accrualTime = 0L,
  accrualIntensity = NA_real_,
  piecewiseSurvivalTime = 0L,
  stratumFraction = 1L,
  kappa1 = NA_real_,
  kappa2 = NA_real_,
  lambda1 = NA_real_,
  lambda2 = NA_real_,
  gamma1 = 0L,
  gamma2 = 0L,
  accrualDuration = NA_real_,
  followupTime = NA_real_,
  fixedFollowup = 0L,
  nullVariance = 0L
)
}
\arguments{
\item{time}{A vector of calendar times for data cut.}

\item{rateRatioH0}{Rate ratio under the null hypothesis.}

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

\item{kappa1}{The dispersion parameter (reciprocal of the shape
parameter of the gamma mixing distribution) for the active treatment
group by stratum.}

\item{kappa2}{The dispersion parameter (reciprocal of the shape
parameter of the gamma mixing distribution) for the control group by
stratum.}

\item{lambda1}{The rate parameter of the negative binomial distribution
for the active treatment group by stratum.}

\item{lambda2}{The rate parameter of the negative binomial distribution
for the control group by stratum.}

\item{gamma1}{The hazard rate for exponential dropout, a vector of hazard
rates for piecewise exponential dropout applicable for all strata, or a
vector of hazard rates for dropout in each analysis time interval
by stratum for the active treatment group.}

\item{gamma2}{The hazard rate for exponential dropout, a vector of hazard
rates for piecewise exponential dropout applicable for all strata, or a
vector of hazard rates for dropout in each analysis time interval
by stratum for the control group.}

\item{accrualDuration}{Duration of the enrollment period.}

\item{followupTime}{Follow-up time for the last enrolled subject.}

\item{fixedFollowup}{Whether a fixed follow-up design is used.
Defaults to 0 for variable follow-up.}

\item{nullVariance}{Whether to calculate the variance for log rate ratio
under the null hypothesis.}
}
\value{
A list with two components:
\itemize{
\item \code{resultsUnderH1}: A data frame containing the following variables:
\itemize{
\item \code{time}: The analysis time since trial start.
\item \code{subjects}: The number of enrolled subjects.
\item \code{nevents}: The total number of events.
\item \code{nevents1}: The number of events in the active treatment
group.
\item \code{nevents2}: The number of events in the control group.
\item \code{ndropouts}: The total number of dropouts.
\item \code{ndropouts1}: The number of dropouts in the active treatment
group.
\item \code{ndropouts2}: The number of dropouts in the control group.
\item \code{nfmax}: The total number of subjects reaching maximum
follow-up.
\item \code{nfmax1}: The number of subjects reaching maximum follow-up
in the active treatment group.
\item \code{nfmax2}: The number of subjects reaching maximum follow-up
in the control group.
\item \code{exposure}: The total exposure time.
\item \code{exposure1}: The exposure time for the active treatment group.
\item \code{exposure2}: The exposure time for the control group.
\item \code{rateRatio}: The rate ratio of the active treatment group
versus the control group.
\item \code{vlogRate1}: The variance for the log rate
parameter for the active treatment group.
\item \code{vlogRate2}: The variance for the log rate
parameter for the control group.
\item \code{vlogRR}: The variance of log rate ratio.
\item \code{information}: The information of log rate ratio.
\item \code{zlogRR}: The Z-statistic for log rate ratio.
}
\item \code{resultsUnderH0} when \code{nullVariance = TRUE}: A data frame
with the following variables:
\itemize{
\item \code{time}: The analysis time since trial start.
\item \code{lambda1H0}: The restricted maximum likelihood estimate
of the event rate for the active treatment group.
\item \code{lambda2H0}: The restricted maximum likelihood estimate
of the event rate for the control group.
\item \code{rateRatioH0}: The rate ratio under H0.
\item \code{vlogRate1H0}: The variance for the log rate
parameter for the active treatment group under H0.
\item \code{vlogRate2H0}: The variance for the log rate
parameter for the control group under H0.
\item \code{vlogRRH0}: The variance of log rate ratio under H0.
\item \code{informationH0}: The information of log rate ratio under H0.
\item \code{zlogRRH0}: The Z-statistic for log rate ratio with variance
evaluated under H0.
\item \code{varianceRatio}: The ratio of the variance under H0 versus
the variance under H1.
\item \code{lambda1}: The true event rate for the active treatment group.
\item \code{lambda2}: The true event rate for the control group.
\item \code{rateRatio}: The true rate ratio.
}
\item \code{resultsUnderH0} when \code{nullVariance = FALSE}: A data frame
with the following variables:
\itemize{
\item \code{time}: The analysis time since trial start.
\item \code{rateRatioH0}: The rate ratio under H0.
\item \code{varianceRatio}: Equal to 1.
\item \code{lambda1}: The true event rate for the active treatment group.
\item \code{lambda2}: The true event rate for the control group.
\item \code{rateRatio}: The true rate ratio.
}
}
}
\description{
Obtains the number of subjects accrued, number of events,
number of dropouts, number of subjects reaching the maximum
follow-up, total exposure, and variance for log rate in each group,
rate ratio, variance, and Wald test statistic of
log rate ratio at given calendar times.
}
\details{
The probability mass function for a negative binomial distribution with
dispersion parameter \eqn{\kappa_i} and rate parameter \eqn{\lambda_i}
is given by
\deqn{P(Y_{ij} = y) = \frac{\Gamma(y+1/\kappa_i)}{\Gamma(1/\kappa_i) y!}
\left(\frac{1}{1 + \kappa_i \lambda_i t_{ij}}\right)^{1/\kappa_i}
\left(\frac{\kappa_i \lambda_i t_{ij}}
{1 + \kappa_i \lambda_i t_{ij}}\right)^{y},}
where \eqn{Y_{ij}} is the event count for subject \eqn{j} in
treatment group \eqn{i}, and \eqn{t_{ij}} is the exposure time for
the subject. If \eqn{\kappa_i=0}, the negative binomial distribution
reduces to the Poisson distribution.

For treatment group \eqn{i}, let \eqn{\beta_i = \log(\lambda_i)}.
The likelihood for \eqn{\{(\kappa_i, \beta_i):i=1,2\}} can be written as
\deqn{l = \sum_{i=1}^{2}\sum_{j=1}^{n_{i}}
\{\log \Gamma(y_{ij} + 1/\kappa_i) - \log \Gamma(1/\kappa_i) + y_{ij}
(\log(\kappa_i) + \beta_i) - (y_{ij} + 1/\kappa_i)
\log(1+ \kappa_i \exp(\beta_i) t_{ij})\}.}
It follows that
\deqn{\frac{\partial l}{\partial \beta_i} = \sum_{j=1}^{n_i}
\left\{y_{ij} - (y_{ij} + 1/\kappa_i)
\frac{\kappa_i \exp(\beta_i) t_{ij}}
{1 + \kappa_i \exp(\beta_i)t_{ij}}\right\},}
and
\deqn{-\frac{\partial^2 l}{\partial \beta_i^2} =
\sum_{j=1}^{n_i} (y_{ij} + 1/\kappa_i) \frac{\kappa_i \lambda_i t_{ij}}
{(1 + \kappa_i \lambda_i t_{ij})^2}.}
The Fisher information for \eqn{\beta_i} is
\deqn{E\left(-\frac{\partial^2 l}{\partial \beta_i^2}\right)
= n_i E\left(\frac{\lambda_i t_{ij}}
{1 + \kappa_i \lambda_i t_{ij}}\right).}
In addition, we can show that
\deqn{E\left(-\frac{\partial^2 l}
{\partial \beta_i \partial \kappa_i}\right) = 0.}
Therefore, the variance of \eqn{\hat{\beta}_i} is
\deqn{Var(\hat{\beta}_i) = \frac{1}{n_i} \left\{
E\left(\frac{\lambda_i t_{ij}}{1 + \kappa_i \lambda_i t_{ij}}\right)
\right\}^{-1}.}

To evaluate the integral, we need to obtain the distribution of the
exposure time,
\deqn{t_{ij} = \min(\tau - W_{ij}, C_{ij}, T_{fmax}),}
where \eqn{\tau} denotes the calendar time since trial start,
\eqn{W_{ij}} denotes the enrollment time for subject \eqn{j}
in treatment group \eqn{i}, \eqn{C_{ij}} denotes the time to dropout
after enrollment for subject \eqn{j} in treatment group \eqn{i}, and
\eqn{T_{fmax}} denotes the maximum follow-up time for
all subjects. Therefore,
\deqn{P(t_{ij} \geq t) = P(W_{ij} \leq \tau - t)P(C_{ij} \geq t)
I(t\leq T_{fmax}).}
Let \eqn{H} denote the distribution function of the enrollment time,
and \eqn{G_i} denote the survival function of the dropout time for
treatment group \eqn{i}. By the change of variables, we have
\deqn{E\left(\frac{\lambda_i t_{ij}}{1 + \kappa_i \lambda_i t_{ij}}
\right) = \int_{0}^{\tau \wedge T_{fmax}}
\frac{\lambda_i}{(1 + \kappa_i \lambda_i t)^2} H(\tau - t) G_i(t) dt.}
A numerical integration algorithm for a univariate function can be
used to evaluate the above integral.

For the restricted maximum likelihood (reml) estimate of
\eqn{(\beta_1,\beta_2)} subject to the
constraint that \eqn{\beta_1 - \beta_2 = \Delta}, we express the
log-likelihood in terms of \eqn{(\beta_2,\Delta,\kappa_1,\kappa_2)},
and takes the derivative of the log-likelihood function with respect
to \eqn{\beta_2}. The resulting score equation has asymptotic limit
\deqn{E\left(\frac{\partial l}{\partial \beta_2}\right) = s_1 + s_2,}
where
\deqn{s_1 = n r E\left\{\lambda1_1 t_{1j} - \left(\lambda_1t_{1j}
+ \frac{1}{\kappa_1}\right) \frac{\kappa_1 e^{\tilde{\beta}_2 +
\Delta}t_{1j}}{1 + \kappa_1 e^{\tilde{\beta}_2 +\Delta}t_{1j}}\right\},}
and
\deqn{s_2 = n (1-r) E\left\{\lambda_2 t_{2j} -
\left(\lambda_2 t_{2j} + \frac{1}{\kappa_2}\right)
\frac{\kappa_2 e^{\tilde{\beta}_2} t_{2j}}
{1 + \kappa_2 e^{\tilde{\beta}_2}t_{2j}}\right\}.}
Here \eqn{r} is the randomization probability for the active
treatment group. The asymptotic limit of the reml of \eqn{\beta_2}
is the solution \eqn{\tilde{\beta}_2} to
\eqn{E\left(\frac{\partial l}{\partial \beta_2}\right) = 0.}
}
\examples{
# Example 1: Variable follow-up design

nbstat(time = c(1, 1.25, 2, 3, 4),
       accrualIntensity = 1956/1.25,
       kappa1 = 5,
       kappa2 = 5,
       lambda1 = 0.7*0.125,
       lambda2 = 0.125,
       gamma1 = 0,
       gamma2 = 0,
       accrualDuration = 1.25,
       followupTime = 2.75)

# Example 2: Fixed follow-up design

nbstat(time = c(0.5, 1, 1.5, 2),
       accrualIntensity = 220/1.5,
       stratumFraction = c(0.2, 0.8),
       kappa1 = 3,
       kappa2 = 3,
       lambda1 = c(0.5*8.4, 0.6*10.5),
       lambda2 = c(8.4, 10.5),
       gamma1 = 0,
       gamma2 = 0,
       accrualDuration = 1.5,
       followupTime = 0.5,
       fixedFollowup = 1,
       nullVariance = 1)

}
\author{
Kaifeng Lu, \email{kaifenglu@gmail.com}
}
