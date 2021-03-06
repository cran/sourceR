% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc.R
\docType{class}
\name{AdaptiveMultiMRW}
\alias{AdaptiveMultiMRW}
\title{AdaptiveMultiMRW}
\format{
Object of \code{\link{R6Class}} with methods for updating a \link{Node} instance.
}
\value{
Object of \code{\link{AdaptiveMultiMRW}}
}
\description{
This class implements an Adaptive Multi-site Metropolis random walk
algorithm.
}
\details{
A multivariate Gaussian proposal is used, for which the proposal variance
is a scaled version of the evolving empirical posterior covariance matrix.
See Roberts and Rosenthal (2012) Examples of Adaptive MCMC. \emph{Journal of Computational
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
  \item{\code{new(node, tune = rep(0.1, length(node$getData())), burning = 100)}}{constructor takes an instance of a \link{StochasticNode}
  node, initial tuning vector (diagonal of adaptive tuning matrix), and number of burnin calls.}
  \item{\code{update()}}{when called, updates \code{node}}
  \item{\code{acceptance()}}{return the acceptance rate}
}
}

\keyword{DAG}
\keyword{MCMC}
\keyword{node}
