\name{summary.source_attribution}
\alias{summary.source_attribution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Summary function for \link{saBayes} output. 
}
\description{
This function generates a median and Chen Shao highest posterior density (HPD) region for each parameter in a source attribution model (\code{source_attribution} object). This algorithm requires unimodiality of the marginal posterior distribution for each parameter. If the marginal posterior distribution is not unimodel, the HPD region calculated may be invalid.
}
\usage{
\method{summary}{source_attribution}(object, alpha = 0.05, burn_in = 0, thin = 1, interval_type = "chen-shao", ...) %, interval_type = "chen-shao"
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  an object that inherits from class \code{source_attribution}, such as a source_attribution object produced using \code{\link{saBayes}}.
}
 \item{alpha}{
 a 100(1-alpha)\% interval has a probability content of 1-alpha.
}
\item{burn_in}{
 the number iterations to remove from the beginning of the MCMC chain prior to calculating the HPD intervals. A value of 0 results in no burn in being removed.
}
\item{thin}{
 the amount to thin the MCMC chain by prior to calculating the HPD intervals. A vlaue of 1 results in no thinning.
}
\item{interval_type}{
 The type of HPD to calculate. This defaults to Chen Shao, but can be set to SPIn
}
\item{...}{
Does nothing at the moment.
}
}
\value{
Returns a list with possible names \code{a}, \code{q}, \code{alpha}, \code{r}, \code{d}, \code{li}, \code{lj} and \code{lj_proportion}. The list contains the same nesting structure as \code{posterior}. 

  \tabular{lll}{
  \emph{Parameter} \tab \emph{Format and nesting} \tab \emph{Returns}\cr
  
  Source effects (\code{a}) \tab A nested list with a matrix for each time and location. eg \code{posterior$a[[t]][[l]][,j]} gives the median, lower and upper values of the Chen Shao HPD region for \code{a} \code{j} for time \code{t}, location \code{l}. \tab Always returned. \cr
  
  Type effects (\code{q}) \tab A matrix. eg \code{posterior$q[,i]} gives the the median, lower and upper values of the Chen Shao HPD region for \code{q} \code{i} for time \code{t}. \tab Returned if \code{names(posterior)} contains \code{q}. \cr
  
    \code{r} \tab A nested list containing an array for each time. eg \code{posterior$r[[t]][i,j,]} gives the  median, lower and upper values of the Chen Shao HPD region for \code{r} (type \code{i}, source \code{j}) at time \code{t}. \tab Returned if \code{names(posterior)} contains \code{r}. \cr
  
  Dispersion parameter (\code{d}) \tab An vector. \tab Returned if \code{names(posterior)} contains \code{d}. \cr
  
  Rate of human infection from type \code{i} (Lambda i: \code{li}) \tab A nested list with a matrix for each time and location. eg \code{posterior$li[[t]][[l]][,i]} gives the median, lower and upper values of the Chen Shao HPD region for \code{li} \code{i} for time \code{t}, location \code{l}. \tab Returned if \code{names(posterior)} contains \code{li}. \cr
  
  Rate of human infection from source \code{j} (Lambda j: \code{lj}) \tab A nested list with a matrix for each time and location. eg \code{posterior$lj[[t]][[l]][,j]} gives the median, lower and upper values of the Chen Shao HPD region for \code{lj} \code{j} for time \code{t}, location \code{l}. \tab Returned if \code{names(posterior)} contains \code{lj}. \cr
  
    Proportion of human infection from source \code{j} (\code{lj_proportion}) \tab A nested list with a matrix for each time and location. eg \code{posterior$lj_proportion[[t]][[l]][,j]} gives the median, lower and upper values of the Chen Shao HPD region for \code{lj_proportion} \code{j} for time \code{t}, location \code{l}. \tab Returned if \code{names(posterior)} contains \code{lj}. \cr
}

}
\references{
Chen, M.-H. and Shao, Q.-M. (1998). Monte Carlo estimation of Bayesian credible and HPD intervals, \emph{Journal of Computational and Graphical Statistics}, 7.
Liu Y, Gelman A, Zheng T (2015). "Simulation-efficient shortest probability intervals." Statistics and Computing.
}
\author{
Poppy Miller <poppy.p.miller@gmail.com>
}

\seealso{
\code{\link{saBayes}}
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
#### Summary #############################################################
##########################################################################

s <- summary(res, alpha = 0.05, burn_in = 1, thin = 1)
# print the median and HPD for the type effects
s$q            
# print the median and HPD for the lambda's for each source
# for time 2, location A
s$lj_proportion$time2$locationA 
## or
s$lj_proportion[[2]][[1]]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
