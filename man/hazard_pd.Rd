% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrappers.R
\name{hazard_pd}
\alias{hazard_pd}
\title{Hazard Function for Progressive Disease (PD)}
\usage{
hazard_pd(
  piecewiseSurvivalTime = 0,
  hazard_pfs = 0.0578,
  hazard_os = 0.02,
  corr_pd_os = 0.5
)
}
\arguments{
\item{piecewiseSurvivalTime}{A vector that specifies the starting time of
piecewise exponential survival time intervals. Must start with 0, e.g.,
\code{c(0, 6)} breaks the time axis into 2 event intervals:
[0, 6) and [6, Inf).
Defaults to 0 for exponential distribution.}

\item{hazard_pfs}{A scalar or numeric vector specifying the
hazard(s) for PFS based on a pwexp distribution.}

\item{hazard_os}{A scalar or numeric vector specifying the
hazard(s) for overall survival (OS) based on a pwexp distribution.}

\item{corr_pd_os}{A numeric value specifying the correlation
between PD and OS times.}
}
\value{
A numeric vector representing the estimated hazard rates
for the pwexp distribution of PD.
}
\description{
Computes the hazard function of a piecewise exponential (pwexp)
distribution for progressive disease (PD), such that the
resulting hazard function for progression-free survival (PFS)
closely matches a given pwexp hazard for PFS.
}
\details{
This function determines the hazard vector \eqn{\lambda_{\text{pd}}}
for the pwexp distribution of PD, so that the implied survival
function for PFS time,
\eqn{T_{\text{pfs}} = \min(T_{\text{pd}}, T_{\text{os}})}, closely
matches the specified pwexp distribution for PFS with hazard vector
\eqn{\lambda_{\text{pfs}}}.

To achieve this, we simulate
\eqn{(Z_{\text{pd}}, Z_{\text{os}})} from
a standard bivariate normal distribution with correlation
\eqn{\rho}. Then, \eqn{U_{\text{pd}} = \Phi(Z_{\text{pd}})}
and \eqn{U_{\text{os}} = \Phi(Z_{\text{os}})} are generated, where
\eqn{\Phi} denotes the standard normal CDF.

The times to PD and OS are obtained via the inverse transform
method using quantile functions of the pwexp distribution:
\deqn{T_{\text{pd}} = \text{qpwexp}(U_{\text{pd}},u,\lambda_{\text{pd}})}
\deqn{T_{\text{os}} = \text{qpwexp}(U_{\text{os}},u,\lambda_{\text{os}})}
where \code{u = piecewiseSurvivalTime}.

The function solves for \eqn{\lambda_{\text{pd}}} such that
the survival function of \eqn{T_{\text{pfs}}} closely matches that
of a pwexp distribution with hazard \eqn{\lambda_{\text{pfs}}}:
\deqn{P(\min(T_{\text{pd}}, T_{\text{os}}) > t) = S_{\text{pfs}}(t)}
Since \deqn{Z_{\text{pd}} =
  \Phi^{-1}(\text{ppwexp}(T_\text{pd}, u, \lambda_{\text{pd}}))} and
\deqn{Z_{\text{os}} =
  \Phi^{-1}(\text{ppwexp}(T_\text{os}, u, \lambda_{\text{os}}))}
we have
\deqn{P(\min(T_{\text{pd}}, T_{\text{os}}) > t) =
P(Z_{\text{pd}} >
    \Phi^{-1}(\text{ppwexp}(t,u,\lambda_{\text{pd}})),
  Z_{\text{os}} >
    \Phi^{-1}(\text{ppwexp}(t,u,\lambda_{\text{os}})))}
while
\deqn{S_{\text{pfs}}(t) = 1 - \text{ppwexp}(t,u,\lambda_{\text{pfs}})}

Matching is performed sequentially at the internal cutpoints
\eqn{u_2, ..., u_J} and at the point
\eqn{u_J + \log(2)/\lambda_{\text{pfs},J}} for the final interval
to solve for \eqn{\lambda_{\text{pd},1}, \ldots,
\lambda_{\text{pd},J-1}} and \eqn{\lambda_{\text{pd},J}}, respectively.
}
\examples{
u <- c(0, 1, 3, 4)
lambda1 <- c(0.0151, 0.0403, 0.0501, 0.0558)
lambda2 <- 0.0145
rho <- 0.5
hazard_pd(u, lambda1, lambda2, rho)

}
\author{
Kaifeng Lu (\email{kaifenglu@gmail.com})
}
