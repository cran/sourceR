\name{flatten}
\alias{flatten}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Flattens the McMC posterior from \link{saBayes} output into a data frame.
}
\description{
Flattens the McMC posterior from \link{saBayes} output into a data frame with one column per parameter.
}
\usage{
flatten(object)
}
\arguments{
  \item{object}{
an object that inherits from class \code{source_attribution}, such as a source_attribution object produced using \code{\link{saBayes}}.
}
}
\value{
Returns a data frame with one column per parameter in the object posterior list. The parameters will be named as follows:
\tabular{lll}{
  \emph{Parameter} \tab \emph{Name format}\cr
  Source effects (a): \tab \code{a_timeName.locationName.sourceName} \cr
  Type effects (q): \tab \code{q_typeName} \cr
  Theta values (theta): \tab \code{theta_typeName} \cr
  Cluster values (cluster): \tab \code{cluster_typeName} \cr
  Relative prevalences (r): \tab \code{r_timeName.typeName.sourceName} \cr
  Lambda i (li): \tab \code{li_timeName_locationName_typeName} \cr
  Lambda j (lj): \tab \code{lj_timeName_locationName_sourceName} \cr
  }
  where \code{timeName}, \code{locationName} and \code{typeName} are the values passed to the \code{time}, \code{location} and \code{type} arguments of \link{saBayes} prefixed with "time", "location", and "type" respectively. The \code{sourceName} values correspond to the source names passed to the \code{formula} argument of  \link{saBayes}.
}
\author{
Poppy Miller
}
\seealso{
\code{\link{saBayes}, \link{subset_posterior}, \link{summary}}
}
\examples{
##########################################################################
## Access simulated data set #############################################
##########################################################################
data(sim_SA)

##########################################################################
## Set priors ############################################################
##########################################################################

priors <- list(a = 1, r = 1, theta = c(0.01, 0.00001))

##########################################################################
## Run model #############################################################
##########################################################################

res <- saBayes(formula = Human~Source1+Source2+Source3+Source4+Source5, 
               time=~Time, location=~Location, type=~Type,
               data=sim_SA$data, priors = priors,
               alpha_conc = 1, prev = sim_SA$prev,
               likelihood_dist = "pois", n_iter = 20)
               
##########################################################################
#### Flatten #############################################################
##########################################################################

flat_res <- flatten(res)
names(flat_res)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
