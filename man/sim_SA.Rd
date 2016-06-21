\name{sim_SA}
\alias{sim_SA}
\docType{data}
\title{
Simulated data
}
\description{
Data simulated from the source attribution model. Intended for use demonstrating use of the \link{sourceR} package.
}
\usage{data("sim_SA")}
\format{
  The format is a list giving the prevalence data dataframe and source and human data dataframes:
\tabular{lll}{
  \emph{Parameter} \tab \emph{Format}\cr
 prev \tab Dataframe giving the prevalence values, source id and time identifiers for each source. \cr
 data \tab Dataframe giving the number of human cases for each type, time and location; and the \cr
 \tab number of positive source samples for each source, type, and time.\cr
 \tab ..$ Human   : number of human cases \cr
 \tab ..$ Type    : type id\cr
 \tab ..$ Time    : time id\cr
 \tab ..$ Location: location id\cr
 \tab ..$ Source1 : number of positive samples for each type (source 1)\cr
 \tab ..$ Source2 : number of positive samples for each type (source 2)\cr
 \tab ..$ Source3 : number of positive samples for each type (source 3)\cr
 \tab ..$ Source4 : number of positive samples for each type (source 4)\cr
 \tab ..$ Source5 : number of positive samples for each type (source 5)\cr
  }
}
\details{
Note that the source data differs between times, but is the same for each location within a time. This is necessary because the source prevalences only vary over time. If the source data is not the same for each location within a given time, the mcmc will return an error. 
}
\examples{
data(sim_SA)
}
\keyword{datasets}
