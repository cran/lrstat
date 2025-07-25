% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getDesignMeans.R
\name{hedgesg}
\alias{hedgesg}
\title{Hedges' g Effect Size}
\usage{
hedgesg(tstat, m, ntilde, cilevel = 0.95)
}
\arguments{
\item{tstat}{The value of the t-test statistic for comparing two
treatment conditions.}

\item{m}{The degrees of freedom for the t-test.}

\item{ntilde}{The normalizing sample size to convert the
standardized treatment difference to the t-test statistic, i.e.,
\code{tstat = sqrt(ntilde)*meanDiff/stDev}.}

\item{cilevel}{The confidence interval level. Defaults to 0.95.}
}
\value{
A data frame with the following variables:
\itemize{
\item \code{tstat}: The value of the \code{t} test statistic.
\item \code{m}: The degrees of freedom for the t-test.
\item \code{ntilde}: The normalizing sample size to convert the
standardized treatment difference to the t-test statistic.
\item \code{g}: Hedges' \code{g} effect size estimate.
\item \code{varg}: Variance of \code{g}.
\item \code{lower}: The lower confidence limit for effect size.
\item \code{upper}: The upper confidence limit for effect size.
\item \code{cilevel}: The confidence interval level.
}
}
\description{
Obtains Hedges' g estimate and confidence interval of
effect size.
}
\details{
Hedges' \eqn{g} is an effect size measure commonly used in meta-analysis
to quantify the difference between two groups. It's an improvement
over Cohen's \eqn{d}, particularly when dealing with small sample sizes.

The formula for Hedges' \eqn{g} is \deqn{g = c(m) d} where \eqn{d}
is Cohen's \eqn{d} effect size estimate, and \eqn{c(m)} is the bias
correction factor, \deqn{d = (\hat{\mu}_1 - \hat{\mu}_2)/\hat{\sigma}}
\deqn{c(m) = 1 - \frac{3}{4m-1}.}
Since \eqn{c(m) < 1}, Cohen's \eqn{d} overestimates the true effect size,
\eqn{\delta = (\mu_1 - \mu_2)/\sigma.}
Since \deqn{t = \sqrt{\tilde{n}} d} we have
\deqn{g = \frac{c(m)}{\sqrt{\tilde{n}}} t} where \eqn{t}
has a noncentral \eqn{t} distribution with \eqn{m} degrees of freedom
and noncentrality parameter \eqn{\sqrt{\tilde{n}} \delta}.

The asymptotic variance of \eqn{g} can be approximated by
\deqn{Var(g) = \frac{1}{\tilde{n}} + \frac{g^2}{2m}.}
The confidence interval for \eqn{\delta}
can be constructed using normal approximation.

For two-sample mean difference with sample size \eqn{n_1} for the
treatment group and \eqn{n_2} for the control group, we have
\eqn{\tilde{n} = \frac{n_1n_2}{n_1+n_2}} and \eqn{m=n_1+n_2-2}
for pooled variance estimate.
}
\examples{

n1 = 7
n2 = 8
meanDiff = 0.444
stDev = 1.201
m = n1+n2-2
ntilde = n1*n2/(n1+n2)
tstat = sqrt(ntilde)*meanDiff/stDev

hedgesg(tstat, m, ntilde)

}
\references{
Larry V. Hedges. Distribution theory for Glass's estimator of
effect size and related estimators.
Journal of Educational Statistics 1981; 6:107-128.
}
\author{
Kaifeng Lu, \email{kaifenglu@gmail.com}
}
