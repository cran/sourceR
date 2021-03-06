% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc.R
\docType{class}
\name{AdaptiveLogDirMRW}
\alias{AdaptiveLogDirMRW}
\title{AdaptiveLogDirMRW}
\format{
Object of \code{\link{R6Class}} with methods for updating a \link{DirichletNode} instance.
}
\value{
Object of \code{\link{AdaptiveLogDirMRW}}
}
\description{
This class implements an Adaptive Multi-site logarithmic Metropolis-Hastings random walk
algorithm, constrained so the parameter vector sums to 1.
}
\details{
An adaptive multivariate log-Gaussian proposal is used for $d-1$ elements of a $d$-dimensional parameter
vector contained in \code{node}, with the $d$th element updated to ensure that the vector sums to 1.
This makes the updater useful for Dirichlet distributed random variables, improving on \link{AdaptiveDirMRW} by
ensuring proposals do not go negative.

For details of the adaptive scheme, see Roberts and Rosenthal (2012) Examples of Adaptive MCMC. \emph{Journal of Computational
and Graphical Statistics}. \bold{18}:349--367.

Please note that no checks are performed as to the suitability of this
algorithm for a particular \link{StochasticNode}.  It is up to the user
to use the correct update algorithm for the appropriate nodes.
}
\section{Fields}{

\describe{
\item{\code{cov}}{the current covariance}

\item{\code{burnin}}{the number of updates to burn in}

\item{\code{tune}}{the current tuning matrix}

\item{\code{naccept}}{the number of accepted proposals}

\item{\code{ncalls}}{the number of times \code{update} has been called}

\item{\code{node}}{the node to which the updater is attached}
}}

\section{Methods}{

\describe{
  \item{\code{new(node, toupdate = function() 1:length(node$getData()), tune = rep(0.1, length(node$getData())), burning = 100)}}{constructor takes an instance of a \link{StochasticNode}
  node, function to choose the indices of the elements to update (by default all elements), initial tuning vector (diagonal of adaptive tuning matrix), and number of calls between
  adaptations.}
  \item{\code{update()}}{when called, updates \code{node}}
  \item{\code{acceptance()}}{return the acceptance rate}
}
}

\keyword{DAG}
\keyword{MCMC}
\keyword{node}
