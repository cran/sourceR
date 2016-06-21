\name{heatmap_types}
\alias{heatmap_types}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Type effect dendrogram
}
\description{
Creates a dendrogram and heatmap for the type effects from an object of class source_attribution. This uses the heatmap.2 function from gplots.
}
\usage{
heatmap_types(object, iters, hclust_method = "complete", hm_colors = c("blue", "white"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object that inherits from class \code{source_attribution}, such as a source_attribution object produced using \code{\link{saBayes}
}}
  \item{iters}{
A vector containing the iterations to be used in constructing the graph. Default is all iteraitons in posterior.
}
  \item{hclust_method}{
See ?hclust for details. Note, the clustering is done on the group of each type effect (cluster column in the posterior) not the actual values of the type effects.
}
  \item{hm_colors}{
Gives the colours for completely dissimilar (dissimilarity value of 1), and identical (dissimilarity value of 0). All other values will be in between the two chosen colours. See ?colorRampPalette for more details.
}
}
\author{
Poppy Miller
}
\seealso{
\code{\link{saBayes}, \link{subset_posterior}, \link{summary}, \link{flatten}}
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
#### Plot heatmap ########################################################
##########################################################################
heatmap_types(object = res, hm_colors = c("purple", "white"))
}