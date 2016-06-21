\name{subset_posterior}
\alias{subset_posterior}

\title{
Subsets the posterior of the McMC posterior from \link{saBayes} output.
}
\description{
Subsets the posterior of the McMC posterior from \link{saBayes} output by parameter, time, location, type, source and iteration number.
}
\usage{
subset_posterior(object, params, t, l, i, j, iters)
}
\arguments{
  \item{object}{
an object that inherits from class \code{source_attribution}, such as a source_attribution object produced using \code{\link{saBayes}}.
}
  \item{params}{
A vector containing the names of the parameters to be returned. Default is all.
}
  \item{t}{
A vector containing the names of the times to be returned. Default is all.
}
  \item{l}{
A vector containing the names of the locations to be returned. Default is all.
}
  \item{i}{
A vector containing the names of the types to be returned. Default is all.
}
  \item{j}{
A vector containing the names of the sources to be returned. Default is all.
}
  \item{iters}{
A vector containing the iterations to be returned.
}
}
\value{
Returns the subsetted posterior in the same format as the input object. This can be passed to \link{summary.source_attribution} or \link{flatten} in the same way the full posterior can be.
}
\author{
Poppy Miller
}
\seealso{
\code{\link{saBayes}, \link{flatten}, \link{summary}}
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
#### Subset ##############################################################
##########################################################################

subset_res <- subset_posterior(res, params = c("a", "li", "lj"), 
                               t = "1", l = "B", j = c("Source2", "Source1"), 
                              i = c("47", "10"), iters = c(3:10))
flat_res <- flatten(subset_res)
str(subset_res)

}
