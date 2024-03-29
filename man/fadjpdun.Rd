% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrappers.R
\name{fadjpdun}
\alias{fadjpdun}
\title{Adjusted p-values for Dunnett-based graphical approaches}
\usage{
fadjpdun(wgtmat, p, family = NULL, corr = NULL)
}
\arguments{
\item{wgtmat}{The weight matrix for intersection hypotheses.}

\item{p}{The raw p-values for elementary hypotheses.}

\item{family}{The matrix of family indicators for elementary hypotheses.}

\item{corr}{The correlation matrix that should be used for the parametric
test. Can contain NAs for unknown correlations between families.}
}
\value{
A matrix of adjusted p-values.
}
\description{
Obtains the adjusted p-values for graphical approaches
using weighted Dunnett tests.
}
\examples{

pvalues <- matrix(c(0.01,0.005,0.015,0.022, 0.02,0.015,0.010,0.023),
                  nrow=2, ncol=4, byrow=TRUE)
w <- c(0.5,0.5,0,0)
g <- matrix(c(0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0),
            nrow=4, ncol=4, byrow=TRUE)
wgtmat = fwgtmat(w,g)

family = matrix(c(1,1,0,0,0,0,1,1), nrow=2, ncol=4, byrow=TRUE)
corr = matrix(c(1,0.5,NA,NA, 0.5,1,NA,NA,
                NA,NA,1,0.5, NA,NA,0.5,1),
              nrow = 4, byrow = TRUE)
fadjpdun(wgtmat, pvalues, family, corr)

}
\references{
Frank Bretz, Martin Posch, Ekkehard Glimm, Florian Klinglmueller,
Willi Maurer, and Kornelius Rohmeyer. Graphical approach for multiple
comparison procedures using weighted Bonferroni, Simes, or
parameter tests. Biometrical Journal. 2011; 53:894-913.
}
\author{
Kaifeng Lu, \email{kaifenglu@gmail.com}
}
