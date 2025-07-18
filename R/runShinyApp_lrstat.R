#' @title Run Shiny App
#' @description Runs the log-rank test power and sample size calculation
#' Shiny app.
#'
#' @return No return value, called for side effects.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
runShinyApp_lrstat <- function() {
  shiny::shinyAppDir(system.file("shinyApp", package = "lrstat"))
}
