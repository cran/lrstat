% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{rawdata}
\alias{rawdata}
\title{A simulated time-to-event data set with 10 replications}
\format{
An object of class \code{data.frame} with 4910 rows and 7 columns.
}
\usage{
rawdata
}
\description{
A simulated data set with stratification and delayed treatment effect:
\describe{
\item{\code{iterationNumber}}{The iteration number}
\item{\code{arrivalTime}}{The enrollment time for the subject}
\item{\code{stratum}}{The stratum for the subject}
\item{\code{treatmentGroup}}{The treatment group for the subject}
\item{\code{timeUnderObservation}}{The time under observation since
randomization}
\item{\code{event}}{Whether the subject experienced the event}
\item{\code{dropoutEvent}}{Whether the subject dropped out}
}
}
\keyword{datasets}
