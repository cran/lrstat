#' @title Print group sequential design
#' @description Prints the stopping boundaries and information inflation 
#' factor for group sequential design.
#'
#' @param x The design object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @export
print.design <- function(x, ...) {
  s = x$byStageResults;
  t = x$overallResults;
  k = length(s$informationRates)
  if (k>1) {
    df = t(data.frame(s$informationRates,
                      s$efficacyBounds,
                      s$futilityBounds,
                      s$cumulativeRejection,
                      s$cumulativeFutility,
                      s$cumulativeAlphaSpent,
                      s$efficacyP,
                      s$futilityP,
                      t$overallReject,
                      t$alpha,
                      t$drift,
                      t$inflationFactor))
    df[seq(9,12), -1] <- NA # only show overall      
    
    
    colnames(df) <- paste("stage", seq_len(ncol(df)), sep=" ")
  } else {
    
    df = t(data.frame(t$overallReject,
                      t$alpha,
                      s$efficacyBounds,
                      s$efficacyP,
                      t$drift))
    
    colnames(df) <- NA
  }
  rownames(df) <- sub("^[[:alpha:]][[:alnum:]]*.", "", rownames(df))
  print( round(df,3), ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print power and sample size results
#' @description Prints the summary statistics from power calculation.
#'
#' @param x The lrpower object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the summary statistics from power calculation.
#'
#' @keywords internal
#'
#' @export
print.lrpower <- function(x, ...) {
  s = x$byStageResults;
  t = x$overallResults;
  k = length(s$informationRates)
  if (k>1) {
    if (t$estimateHazardRatio) {
      df = t(data.frame(s$informationRates,
                        s$efficacyBounds,
                        s$futilityBounds,
                        s$cumulativeRejection,
                        s$cumulativeFutility,
                        s$cumulativeAlphaSpent,
                        s$numberOfEvents,
                        s$numberOfDropouts,
                        s$numberOfSubjects,
                        s$analysisTime,
                        s$efficacyHR,
                        s$futilityHR,
                        s$efficacyP,
                        s$futilityP,
                        s$information,
                        s$HR,
                        t$overallReject,
                        t$alpha,
                        t$numberOfEvents,
                        t$expectedNumberOfEvents,
                        t$numberOfDropouts,
                        t$expectedNumberOfDropouts,
                        t$numberOfSubjects,
                        t$expectedNumberOfSubjects,
                        t$studyDuration,
                        t$expectedStudyDuration,
                        t$accrualDuration,
                        t$followupTime,
                        t$fixedFollowup,
                        t$rho1,
                        t$rho2))
      df[seq(17,31), -1] <- NA # only show overall      
    } else {
      df = t(data.frame(s$informationRates,
                        s$efficacyBounds,
                        s$futilityBounds,
                        s$cumulativeRejection,
                        s$cumulativeFutility,
                        s$cumulativeAlphaSpent,
                        s$numberOfEvents,
                        s$numberOfDropouts,
                        s$numberOfSubjects,
                        s$analysisTime,
                        s$efficacyP,
                        s$futilityP,
                        s$information,
                        t$overallReject,
                        t$alpha,
                        t$numberOfEvents,
                        t$expectedNumberOfEvents,
                        t$numberOfDropouts,
                        t$expectedNumberOfDropouts,
                        t$numberOfSubjects,
                        t$expectedNumberOfSubjects,
                        t$studyDuration,
                        t$expectedStudyDuration,
                        t$accrualDuration,
                        t$followupTime,
                        t$fixedFollowup,
                        t$rho1,
                        t$rho2))
      df[seq(14,28), -1] <- NA # only show overall      
    }
    
    colnames(df) <- paste("stage", seq_len(ncol(df)), sep=" ")
  } else {
    if (t$estimateHazardRatio) {
      df = t(data.frame(t$overallReject,
                        t$alpha,
                        t$numberOfEvents,
                        t$numberOfDropouts,
                        t$numberOfSubjects,
                        t$studyDuration,
                        t$accrualDuration,
                        t$followupTime,
                        t$fixedFollowup,
                        t$rho1,
                        t$rho2,
                        s$efficacyBounds,
                        s$efficacyHR,
                        s$efficacyP,
                        s$information,
                        s$HR))
    } else {
      df = t(data.frame(t$overallReject,
                        t$alpha,
                        t$numberOfEvents,
                        t$numberOfDropouts,
                        t$numberOfSubjects,
                        t$studyDuration,
                        t$accrualDuration,
                        t$followupTime,
                        t$fixedFollowup,
                        t$rho1,
                        t$rho2,
                        s$efficacyBounds,
                        s$efficacyP,
                        s$information))
    }
    
    colnames(df) <- NA
  }
  rownames(df) <- sub("^[[:alpha:]][[:alnum:]]*.", "", rownames(df))
  print( round(df,3), ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print simulation results
#' @description Prints the summary statistics from simulation.
#'
#' @param x The lrsim object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the summary statistics from simulation runs.
#'
#' @keywords internal
#'
#' @export
print.lrsim <- function(x, ...) {
  s = x$overview
  k = length(s$numberOfEvents)
  if (k>1) {
    df = t(data.frame(s$cumulativeRejection,
                      s$cumulativeFutility,
                      s$numberOfEvents,
                      s$numberOfDropouts,
                      s$numberOfSubjects,
                      s$analysisTime,
                      s$overallReject,
                      s$expectedNumberOfEvents,
                      s$expectedNumberOfDropouts,
                      s$expectedNumberOfSubjects,
                      s$expectedStudyDuration))
    df[c(7,8,9,10,11), -1] <- NA # only show overall
    colnames(df) <- paste("stage", seq_len(ncol(df)), sep=" ")
  } else {
    df = t(data.frame(s$overallReject,
                      s$expectedNumberOfEvents,
                      s$expectedNumberOfDropouts,
                      s$expectedNumberOfSubjects,
                      s$expectedStudyDuration))
    colnames(df) <- NA
  }
  rownames(df) <- sub("^[[:alpha:]][[:alnum:]]*.", "", rownames(df))
  print( round(df,3), ..., na.print = "" , quote = FALSE )
  invisible(x)
}
