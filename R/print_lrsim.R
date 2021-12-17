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
  s = x$sumstat
  k = length(s$eventsPerStage)
  if (k>1) {
    df = t(data.frame(s$eventsPerStage,
                      s$expectedNumberOfEvents,
                      s$analysisTime,
                      s$expectedStudyDuration,
                      s$numberOfSubjects,
                      s$expectedNumberOfSubjects,
                      s$futilityPerStage,
                      s$rejectPerStage,
                      s$overallReject))
    df[c(2,4,6,9), -1] <- NA # only show overall
    colnames(df) <- paste("stage", seq_len(ncol(df)), sep=" ")
  } else {
    df = t(data.frame(s$expectedNumberOfEvents,
                      s$expectedStudyDuration,
                      s$expectedNumberOfSubjects,
                      s$overallReject))
    colnames(df) <- NA
  }
  rownames(df) <- sub("^[[:alpha:]][[:alnum:]]*.", "", rownames(df))
  print( df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}
