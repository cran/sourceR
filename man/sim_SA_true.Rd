\name{sim_SA_true}
\alias{sim_SA_true}
\docType{data}
\title{
True values of the parameters used to generate the simulated data set \link{sim_SA}.
}
\description{
True values of the parameters used to generate the simulated data set \link{sim_SA}.
}
\usage{data("sim_SA_true")}
\format{
  The format is:
\tabular{lll}{
  \emph{Parameter} \tab \emph{Format}\cr
  a \tab Dataframe with vectors giving the value, source, time and location identifiers for each source effect. \cr
  q \tab Dataframe with vectors giving the value and type identifiers for each source effect. \cr
  alpha \tab Gives the concentration parameter for the Dirichlet Process. \cr
  cluster \tab Dataframe with vectors giving the value and type identifiers for each source effect. \cr
  lambda_i \tab Dataframe with vectors giving the value, type, time and location identifiers for each lambda i (number of cases attributed to each type). \cr
  lambda_j \tab Dataframe with vectors giving the value, source, time and location identifiers for each lambda j (number of cases attributed to each source). \cr
}
}
\examples{
data(sim_SA_true)
}
\keyword{datasets}
