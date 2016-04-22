\name{campy}
\alias{campy}
\docType{data}
\title{
Campylobacter source attribution data set from the Manawatu, New Zealand.
}
\description{
Counts of human campylobacteriosis cases and positive source samples attributed to each bacterial subtype. For more information, see 
}
\usage{data("campy")}
\format{
  A data frame with 115 observations on the following 10 variables.
  \describe{
    \item{\code{Human}}{a numeric vector giving the number of human cases of campylobacteriosis attributed to each subtype \eqn{i}}
    \item{\code{ChickenA}}{a numeric vector giving the number of positive samples from Chicken supplier A for each subtype}
    \item{\code{ChickenB}}{a numeric vector giving the number of positive samples from Chicken supplier B for each subtype}
    \item{\code{ChickenC}}{a numeric vector giving the number of positive samples from Chicken supplier C for each subtype}
    \item{\code{Bovine}}{a numeric vector giving the number of positive samples from Bovine source samples for each subtype}
    \item{\code{Ovine}}{a numeric vector giving the number of positive samples from Ovine source samples for each subtype}
    \item{\code{Environment}}{a numeric vector giving the number of positive samples from environmental source samples for each subtype}
    \item{\code{Time}}{a factor with levels \code{1}}
    \item{\code{Location}}{a factor with levels \code{A}}
    \item{\code{Type}}{a factor with one level per subtype 
    }
  }
}
\references{
Mullner P, Jones G, Noble A, Spencer S, Hathaway S, French N (2009a). "Source Attribution of Food Borne Zoonoses in New Zealand: A Modified Hald Model." Risk Analysis, 29(7).
}
\examples{
data(campy)
}
\keyword{datasets}
