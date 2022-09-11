#' @title Error spending functions
#' @description Obtains the error spent at the given information fractions
#' for the specified error spending function.
#'
#' @param t A vector of information fractions.
#' @param error Total error to spend.
#' @param sf Spending function. One of the following: "sfOF" for
#'  O'Brien-Fleming type spending function, "sfP" for Pocock type spending
#'  function, "sfKD" for Kim & DeMets spending function, and "sfHSD" for
#'  Hwang, Shi & DeCani spending function. Defaults to "sfOF".
#' @param sfpar Parameter for the spending function. Corresponds to rho for
#'  "sfKD" and gamma for "sfHSD".
#'
#' @return A vector of errors spent up to the interim look.
#'
#' @examples
#' errorSpent(t = 0.5, error = 0.025, sf = "sfOF")
#' errorSpent(t = c(0.5, 0.75, 1), error = 0.025, sf = "sfHSD", sfpar = -4)
#'
#' @export
errorSpent <- function(t, error, sf = "sfOF", sfpar = NA) {
  sapply(t, errorSpentcpp, error, sf, sfpar)
}



#' @title Quantile function of truncated piecewise exponential distribution
#' @description Obtains the quantile of a piecewise expoenential distribution
#' given that it exceeds a specified lower bound.
#'
#' @param probability The scalar probability corresponding to the quantile.
#' @inheritParams param_piecewiseSurvivalTime
#' @inheritParams param_lambda
#' @param lowerBound The left truncation time point for the survival time.
#' Defaults to 0 for no truncation.
#'
#' @return The quantile x such that
#' P(X > x | X > lowerBound) = 1 - probability.
#'
#' @examples
#' qtpwexp(probability = c(0.3, 0.5), piecewiseSurvivalTime = c(0, 6, 9, 15),
#'         lambda = c(0.025, 0.04, 0.015, 0.007))
#'
#' @export
qtpwexp <- function(probability, piecewiseSurvivalTime = 0, 
                    lambda = 0.0578, lowerBound = 0) {
  sapply(probability, qtpwexpcpp, piecewiseSurvivalTime, lambda, lowerBound)
}



#' @title Stagewise exit probabilities
#' @description Obtains the stagewise exit probabilities for both efficacy and
#' futility stopping.
#'
#' @param b Upper boundaries on the z-test statistic scale.
#' @param a Lower boundaries on the z-test statistic scale. Defaults to
#' \code{c(rep(-6.0, kMax-1), b[kMax])} if left unspecified, where
#' \code{kMax = length(b)}.
#' @param theta Stagewise parameter of interest, e.g., \code{-U/V} for
#'   weighted log-rank test, where \code{U} is the mean and \code{V} is
#'   the variance of the weighted log-rank test score statistic at each stage.
#'   For proportional hazards and conventional log-rank test, use the
#'   scalar input, \code{theta = -log(HR)}. Defaults to 0 corresponding to 
#'   the null hypothesis.
#' @param I Stagewise cumulative information, e.g., \code{V}, the variance
#'   of the weighted log-rank test score statistic at each stage. For
#'   conventional log-rank test, information can be approximated by
#'   \code{phi*(1-phi)*D}, where \code{phi} is the probability of being
#'   allocated to the active arm, and \code{D} is the total number of events 
#'   at each stage. Defaults to \code{seq(1, kMax)} if left unspecified. 
#'
#' @return A list of stagewise exit probabilities: one vector for efficacy
#' stopping probabilities, and the other vector for futility stopping
#' probabilities.
#'
#' @examples
#' exitprob(b = c(3.471, 2.454, 2.004), theta = -log(0.6), 
#'          I = c(50, 100, 150)/4)
#'          
#' exitprob(b = c(2.963, 2.359, 2.014), 
#'          a = c(-0.264, 0.599, 2.014), 
#'          theta = c(0.141, 0.204, 0.289), 
#'          I = c(81, 121, 160))
#'
#' @export
exitprob <- function(b, a = NA, theta = 0, I = NA) {
  exitprobcpp(b, a, theta, I)
}


#' @title Adjusted p-values for Bonferroni-based graphical approaches
#' @description Obtains the adjusted p-values for graphical approaches 
#' using weighted Bonferroni tests for fixed design.
#'
#' @param w The vector of initial weights for elementary hypotheses.
#' @param G The initial transition matrix.
#' @param p Raw p-values for elementary hypotheses.
#'
#' @return A matrix of adjusted p-values.
#'
#' @references
#' Frank Bretz, Willi Maurer, Werner Brannath and Martin Posch. A 
#' graphical approach to sequentially rejective multiple test 
#' procedures. Statistics in Medicine. 2009;28:586-604. 
#' 
#' @examples
#'
#' pvalues <- matrix(c(0.01,0.005,0.015,0.022, 0.02,0.015,0.010,0.023),
#'                   nrow=2, ncol=4, byrow=TRUE)
#' w <- c(0.5,0.5,0,0)
#' g <- matrix(c(0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0), 
#'             nrow=4, ncol=4, byrow=TRUE)
#' fadjpbon(w, g, pvalues)
#'
#' @export
fadjpbon <- function(w, G, p) {
  m = length(w)
  
  if (!is.matrix(p)) {
    p = matrix(p, ncol=m)
  }
  
  x = fadjpboncpp(w, G, p)
  if (nrow(x) == 1) {
    x = as.vector(x)
  }
  x
}

