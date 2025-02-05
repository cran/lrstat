% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{svdcpp}
\alias{svdcpp}
\title{Singular Value Decomposition of a Matrix}
\usage{
svdcpp(X, outtransform = 1L, decreasing = 1L)
}
\arguments{
\item{X}{A numeric matrix whose SVD decomposition is to be computed.}

\item{outtransform}{Whether the orthogonal matrices composing of the
left and right singular vectors are to be computed.}

\item{decreasing}{Whether the singular values should be sorted in
decreasing order and the corresponding singular vectors rearranged
accordingly.}
}
\value{
A list with the following components:
\itemize{
\item \code{d}: A vector containing the singular values of \eqn{X}.
\item \code{U}: A matrix whose columns contain the left singular vectors
of \eqn{X}.
\item \code{V}: A matrix whose columns contain the right singular vectors
of \eqn{X}.
}
}
\description{
Computes the singular-value decomposition of a
rectangular matrix.
}
\details{
Given \eqn{A \in R^{m\times n} (m \geq n)}, the following algorithm
overwrites \eqn{A} with \eqn{U^T A V = D}, where
\eqn{U\in R^{m\times m}} is orthogonal, \eqn{V \in R^{n\times n}} is
orthogonal, and \eqn{D \in R^{m\times n}} is diagonal.
}
\examples{

A <- matrix(c(1,0,0,0, 1,2,0,0, 0,1,3,0, 0,0,1,4), 4, 4)
svdcpp(A)

}
\references{
Gene N. Golub and Charles F. Van Loan.
Matrix Computations, second edition. Baltimore, Maryland:
The John Hopkins University Press, 1989, p.434.
}
\author{
Kaifeng Lu, \email{kaifenglu@gmail.com}
}
