\name{saBayes}
\alias{saBayes}
\title{
Fits a non-parametric Bayesian model to source attribution data
}
\description{
This function fits a non-parametric source attribution model for human cases of disease. It supports multiple types, sources, times and locations. The number of human cases for each type, time and location can follow a Poisson or Negative Binomial likelihood.
\deqn{y_{itl}\sim\textsf{Poisson}(\lambda_{itl})} 
or
\deqn{y_{itl}\sim\textsf{Negative Binomial}(\lambda_{itl},d)}
}
\usage{
saBayes(formula, time, location, type, data, priors, alpha_conc, prev, 
        likelihood_dist, n_iter, mcmc_params, initials, params_fix)
}

\arguments{

\item{formula}{
A formula object of the form \code{y ~ x1 + x2 + ... + xm}, where \code{y} is the name of the human cases column, and \code{x1, ..., xm} are the names of the source count columns in the data.
}

\item{time}{
A formula object of the form \code{~t}, where \code{t} is the name of the column containing the times in the data.
}

\item{location}{
A formula object of the form \code{~l}, where \code{l} is the name of the column containing the locations in the data.
}

\item{type}{
A formula object of the form \code{~s}, where \code{s} is the name of the column containing the (sub) types in the data.
}

  \item{data}{
  Data frame with columns for human cases, source counts (one column for each source), time, location, and type. The data for the human cases and source counts must be integers. The data for the time, location and type columns must be factors. The source counts are only allowed to vary over time, hence they must be repeated for each location within each time.
  
  \tabular{cccccc}{
\emph{Human Case Counts (y)} \tab \emph{Source 1 Counts (x1)} \tab \emph{Source 2 Counts (x2)}  \tab \emph{Time} \tab \emph{Location} \tab \emph{Type}\cr \cr
y(loc A, time 1, type 1) \tab x1(time 1, type 1) \tab x2(time 1, type 1)  \tab 1 \tab A \tab 1\cr
y(loc A, time 1, type 2) \tab x1(time 1, type 2) \tab x2(time 1, type 2)  \tab 1 \tab A \tab 2\cr
y(loc A, time 1, type 3) \tab x1(time 1, type 3) \tab x2(time 1, type 3)  \tab 1 \tab A \tab 3\cr \cr
y(loc B, time 1, type 1) \tab x1(time 1, type 1) \tab x2(time 1, type 1)  \tab 1 \tab B \tab 1\cr
y(loc B, time 1, type 2) \tab x1(time 1, type 2) \tab x2(time 1, type 2)  \tab 1 \tab B \tab 2\cr
y(loc B, time 1, type 3) \tab x1(time 1, type 3) \tab x2(time 1, type 3)  \tab 1 \tab B \tab 3\cr \cr
y(loc A, time 2, type 1) \tab x1(time 2, type 1) \tab x2(time 2, type 1)  \tab 2 \tab A \tab 1\cr
y(loc A, time 2, type 2) \tab x1(time 2, type 2) \tab x2(time 2, type 2)  \tab 2 \tab A \tab 2\cr
y(loc A, time 2, type 3) \tab x1(time 2, type 3) \tab x2(time 2, type 3)  \tab 2 \tab A \tab 3\cr \cr
y(loc B, time 2, type 1) \tab x1(time 2, type 1) \tab x2(time 2, type 1)  \tab 2 \tab B \tab 1\cr
y(loc B, time 2, type 2) \tab x1(time 2, type 2) \tab x2(time 2, type 2)  \tab 2 \tab B \tab 2\cr
y(loc B, time 2, type 3) \tab x1(time 2, type 3) \tab x2(time 2, type 3)  \tab 2 \tab B \tab 3\cr
}
}

\item{priors}{
A list giving the parameters of the priors for each parameter in the model with names \code{a}, \code{q}, \code{alpha}, \code{r}, and \code{d}. The priors for the source effect parameters and the relative prevalence parameters (\code{a} and \code{r}) can be given as a single number (to be used for the parameter at all times and locations) or a data frame. Please see the priors section of the vignette for further guidence on choosing priors for the parameters.

    \tabular{lllll}{
  \emph{Parameter} \tab \emph{Prior Distribution} \tab \emph{Prior Parameters} \tab \emph{Required}\cr
  
  \code{r} \tab Dirichlet(concentration) \tab A single positive number or a data frame with columns giving the prior values (named \code{value}), times (named the name of the time column in the data) and source and type id's (named \code{source_id} and the name given to the type id column in the data). If single number is supplied, it will be used for all times.\tab Always required. \cr
  
  Source effects (\code{a}) \tab Dirichlet(concentration) \tab A single positive number or a dataframe with columns giving the prior values (named \code{value}), times (named the name of the time column in the data), locations (named the name of the location column in the data) and the source id (named \code{source_id}).\tab Always required. \cr
  
  Type effects base distribution (\code{theta}) \tab DPM(Gamma(shape, rate), alpha) \tab Numerical vector of length 2 for the shape and rate of the Gamma base distribution. \tab Always required. \cr
  
  Dispersion (\code{d}) \tab Log Normal(meanlog, sdlog) \tab Numerical vector of length 2. \tab Required if \code{likelihood == "nbinom"}. \cr
}
}

\item{alpha_conc}{
A single positive number giving the concentration parameter for the base distribution of the type effects: DP(Gamma(shape, rate), alpha_conc). See the vignette for further guidence on how to choose the concentration parameter.
}

\item{prev}{
An optional data frame with columns giving the prevalences (named \code{value}) and the times (named the name of the time column in the data). Prevalences must be between 0 and 1. Prevalences are the proportion of samples that were positive for any type for a given source and time.
}

\item{likelihood_dist}{
Used to specify the likelihood distribution of the number of human cases from each type. Only Poisson (\code{"pois"}) and Negative Binomial (\code{"nbinom"}) are currently implemented. The model cannot have type effects if using a Negative Binomial distribution (ie the type effects will all be set to 1). This is done internally by setting \code{params_fix$type_effects} and \code{params_fix$alpha} to \code{TRUE}.
}

\item{n_iter}{
The total number of iterations if MCMC chain will run for.
}

\item{mcmc_params}{
An optional list containing parameters relating to the mcmc chain. 
    \tabular{lll}{
  \emph{Item} \tab \emph{Description} \tab \emph{Default} \cr
  
  \code{save_lambda} \tab If \code{TRUE} the lambda \code{i} and \code{j} parameters will be saved in the \code{posterior} list. These can be calculated from the other parameters after running the model. It is quicker and uses less memory to set \code{FALSE}. \tab \code{TRUE} \cr
  
  \code{burn_in} \tab The first k samples are discarded, where k is the value of \code{burn_in}. Discarding the burn in when running the model reduces storage requirements.\tab 0 \cr
  
  \code{thin} \tab Samples from every kth iteration will be stored, where k is the value of \code{thin}. Thinning whilst running the model reduces storage requirements. \tab 1 \cr
  
  \code{n_r} \tab The \code{r} matrices for each time are updated using single site updates for each \code{r_ijt}. \code{n_r} sets the number of \code{r_ijt}'s that are updated per iteration. \tab ceiling(0.4 * number of \code{r_ijt}'s) \cr
}
}

\item{initials}{
An optional list giving the starting values for the parameters.
    \tabular{lll}{
  \emph{Parameter} \tab \emph{Description} \tab \emph{Default} \cr
  
  \code{r} \tab A single positive number or a data frame with columns giving the prior values (named \code{value}), times (named the name of the time column in the data) and source and type id's (named \code{source_id} and the name given to the type id column in the data). If single number is supplied, it will be used for all times. \tab The default initial values are the maximum likelihood point estimates of \code{r} from the source matrix (i.e. \eqn{r_ij = x_ij / sum_i=1^n x_ij}). \cr
  
  Source effects (\code{a}) \tab A single positive number (which will be used as the starting value for all source effects) or a data frame with columns named \code{value} (containing the initial values), \code{source_id} (containing the source names which must be the same names as the source columns in the data) and columns giving the time and location for each parameter (named the same names as the time and location columns in the data). \tab The default initial values for the source effects are drawn the prior distribution (Dirichlet). \cr
  
  Type effects (\code{q}) \tab A data frame with columns giving the initial values (named \code{value}) and the type ids (named the same name as the type column in the data data frame). \tab Default is to initialise all type effects to be in a single group with a theta value drawn from \eqn{\theta_{1}\sim L(y|\theta_{1},a,p)\times Q_0(\theta_{1})}. \cr
  
  Dispersion parameter (\code{d}) \tab A single positive number. Only required if using a Negative Binomial likelihood. \tab Drawn from the prior (Log Normal with parameters specified in the \code{priors} list) for this parameter. \cr
}
}

\item{params_fix}{
An optional list specifying whether to fix the \code{r} matrices (\code{r}), by default it is set to \code{FALSE}.
}
}
\details{
This function fits a source attribution model for human cases of disease. It supports multiple types, sources, times and locations. The number of human cases for each type, time and location follows a Poisson or Negative Binomial likelihood.

\emph{Model}

\deqn{y_{itl}\sim\textsf{Poisson}(\lambda_{itl})} 
or
\deqn{y_{itl}\sim\textsf{Negative Binomial}(\lambda_{itl},d)}
where
\deqn{\lambda_{itl}=\sum_{j=1}^{m}\lambda_{ijtl}=q_{k(i)}\sum_{j=1}^{m}r_{ijt}\cdot \pi_{j}\cdot a_{jtl}}

The parameters are defined as follows: 
\deqn{a_{jtl}} is the unknown source effect for source \eqn{j}, time \eqn{t}, location \eqn{l}
\deqn{q_{k(i)}} is the unknown type effect for type \eqn{i} in group \eqn{k}.
\deqn{x_{ij}} is the known number of positive samples for each source \eqn{j} type\eqn{i} combination
\deqn{n_{ij}} is the known total number of samples for each source \eqn{j} type \eqn{i} combination
\deqn{\pi_{j}} is the fixed prevalence in source (i.e. the number of positive samples divided by the number of negative samples) \eqn{j}
\deqn{r_{ijt}}  is the unknown relative occurrence of type \eqn{i} on source \eqn{j}. 
A point estimate of \eqn{r_{ijt}} (which is used as the starting value for the MCMC chain) is
\deqn{r_{ijt}=x_{ij}/\sum_{i=1}^{n} x_{ij}}

\emph{Priors}
\deqn{x_{.j}\sim Multinomial(n_{j},p_{.j})}
\deqn{r_{.j}\sim Dirichlet(k_1,...,k_n)}
\deqn{a\sim Dirichlet(u_1,...,u_n)}
\deqn{q\sim DP(alpha, Gamma(alpha_{q},beta_{q}))}

\emph{Guidence on choosing priors and the concentration parameter for the DP}:

Dirichlet priors for the source effects and relative prevalences (r): symmetric Dirichlet priors (with a concentration parameter equal to 1) are often used when there is no prior knowledge favouring one component over another (it is uniform over all points in its support). 

Dirichlet process parameters: a DP is specified by a base distribution (Gamma in this model) and a positive real number alpha_conc called the concentration parameter. The base distribution is the expected value of the process; the concentration parameter specifies how strong the prior grouping is. In the limit alpha_conc tends to 0, all types will be assigned to one group, increasing alpha_conc makes a larger number of groups increasingly likely. A relatively uninformative base distribution can be specified with parameters c(0.01, 0.00001). 

See the vignette for further details about choosing priors.
}
\value{
Returns a list containing lists named \code{posterior} and \code{acceptance_rate}.
            
\item{posterior}{List with items named \code{a}, \code{q}, \code{cluster}, \code{theta}, \code{r}, \code{d}, \code{li} and \code{lj}. 

      \tabular{lll}{
  \emph{Parameter} \tab \emph{Format and nesting} \tab \emph{ }\cr
  
  Source effects (\code{a}) \tab A nested list with a matrix (dimensions (\code{n_iter}, \code{sources})) for each time and location. eg \code{posterior$a[[t]][[l]][n,j]} gives the value for \code{a} \code{j} for time \code{t}, location \code{l} at iteration \code{n}. \tab Always returned. \cr
  
  Type effects (\code{q}) \tab A matrix (dimensions (\code{n_iter}, \code{types})) containing the \code{theta} value for the cluster \code{g} that the type effect {qi} is assigned to in iteration \code{n}. eg \code{posterior$q[n,i]} gives the value for \code{q} \code{i} at iteration \code{n}. \tab Always returned. \cr
  
  Type effect cluster (\code{cluster}) \tab A matrix containing the group each type effect was assigned to at each iteration (dimensions (\code{n_iter}, \code{types})). eg \code{posterior$cluster[n,i]} gives the cluster for \code{q} \code{i} at iteration \code{n}. \tab Returned if \code{likelihood_dist == "pois"}. \cr
  
  Theta (\code{theta}) \tab A matrix containing the type effect parameters for each group (dimensions (\code{n_iter}, \code{gmax})). eg \code{posterior$theta[n,c]} gives the \code{theta} value for cluster \code{c} at iteration \code{n}. \tab Returned if \code{likelihood_dist == "pois"}. \cr
  
    \code{r} \tab A nested list containing an array (with dimensions(\code{types}, \code{sources}, \code{n_iter})) for each time. eg \code{posterior$r[[t]][i,j,n]} gives the value for \code{r} (type \code{i}, source \code{j}) at time \code{t}, iteration \code{n}. \tab Returned if \code{fix_r == FALSE}. \cr
  
  Dispersion parameter (\code{d}) \tab A vector. \tab Returned if \code{likelihood_dist == "nbinom"}. \cr
  
  Rate of human infection from type \code{i} (Lambda i: \code{li}) \tab A nested list with a matrix (dimensions (\code{n_iter}, \code{types})) for each time and location. eg \code{posterior$li[[t]][[l]][n,i]} gives the value for \code{li} \code{i} for time \code{t}, location \code{l} at iteration \code{n}. \tab Returned if \code{mcmc_params$save_lambda == TRUE}. \cr
  
  Rate of human infection from source \code{j} (Lambda j: \code{lj}) \tab A nested list with a matrix (dimensions (\code{n_iter}, \code{sources})) for each time and location. eg \code{posterior$lj[[t]][[l]][n,j]} gives the value for \code{lj} \code{j} for time \code{t}, location \code{l} at iteration \code{n}. \tab Returned if \code{mcmc_params$save_lambda == TRUE}. \cr
}
}
\item{acceptance_rate}{List with items named \code{a}, \code{r}, and \code{d}.

      \tabular{lll}{
  \emph{Parameter} \tab \emph{Format and nesting} \tab \emph{ }\cr
  
  Source effects (\code{a}) \tab A nested list with a positive numerical vector for each time and location. eg \code{acceptance_rate$a_centered[[t]][[l]][3]} gives the acceptance rate for \code{a} update for the third source effect at time \code{t}, location \code{l}. \tab Always returned. These should all be around 0.45. \cr
  
    \code{r} \tab A nested list containing a matrix (with dimensions(\code{types}, \code{sources})) for each time. eg \code{acceptance_rate$r[[t]][i,j]} gives the acceptance rate for \code{r} (type \code{i}, source \code{j}) at time \code{t}. \tab Only returned if \code{fix_r == FALSE}. These should all be around 0.45. \cr
  
  Dispersion parameter (\code{d}) \tab positive numerical vector of length 1. \tab Only returned if \code{likelihood_dist == "nbinom"} This should be around 0.45. \cr
}
}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Poppy Miller <p.miller@lancaster.ac.uk>
}

\seealso{
\code{\link{summary}}
}
\examples{
##########################################################################
## Access simulated data set #############################################
##########################################################################
#require(sourceR)
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
               
## or
## Run model using all parameters 
\dontrun{             
res <- saBayes(formula = Human~Source1+Source2+Source3+Source4+Source5, 
               time=~Time, location=~Location, type=~Type,
               data=sim_SA$data, priors = priors,
               alpha_conc = 1, prev = sim_SA$prev,
               likelihood_dist = "pois", n_iter = 5000,
               params_fix = list(type_effects = FALSE, r = FALSE),
               mcmc_params = list(save_lambda = TRUE, burn_in = 1, 
               thin = 1, n_r = 200))
}
               
##########################################################################
#### Results #############################################################
##########################################################################

## Example: Look at data #################################################
\dontrun{
res$posterior$cluster
}

## Example: Convert rates into proportions, example lambda j #############

lambda_j_proportion_t1_l1 <- apply(res$posterior$lj[[1]][[1]], 1, 
                                    function(x) x/sum(x))

## Example: Trace plot ###################################################
# trace plot for lambda (type = 3, time = 1, location = 2)
plot(res$posterior$li$time1$locationB[, "3"], type="l")
## or
plot(res$posterior$li[[1]][[2]][, 3], type="l")

## Example: Density plot #################################################
# density plot for type effect (type = 3)
plot(density(res$posterior$q[,"3"]))

## Example: Residual plot ################################################

summary_li <- summary(res)$li
median_li <- c(summary_li$time1$locationA[, "median"], 
               summary_li$time1$locationB[, "median"],
               summary_li$time2$locationA[, "median"],
               summary_li$time2$locationB[, "median"])

plot(median_li, sim_SA$data[, "Human"], 
     xlab="lambda i", ylab="Human cases")
     abline(0,1)
  
## Example: Box plots showing the type effect grouping ###################
  require(ggplot2)
  library(cluster)
  \dontrun{
  # subset your data by removing a burn in (of 1000) from the total number 
  # of iterations (5000) and thinning to have a maximum of 10000 values as 
  # this can be very slow
  sub_hclust <- seq(1000, 5000, by=ceiling(5000/10000))
  groups <- as.data.frame(t(res$posterior$cluster)) 
  groups <- as.data.frame(apply(groups, 2, function(x) as.factor(x)))
  
  # compute dissimilarity matrix for the type effect clusters
  disim_clust_g <- daisy(groups[,sub_hclust])
  clu <- hclust(disim_clust_g, method="complete")
  # hclust_groups contains the group id for each type effect, 
  # if they were split into 3 groups using the tree generated by the
  # hclust command above
  hclust_groups <- cutree(clu, k=3) 
  
  # Put the data in the format ggplot2 requires it to be in
  # i.e. a column for the values for q at each iteration (thinned)
  # a column for the factor name (type i) so that it will make boxplots
  # and a column for the group
  q_full <- data.frame(qi=c(as.vector(res$posterior$q[sub_hclust,])), 
    i=as.factor(rep(1:50, each=dim(res$posterior$q[sub_hclust,])[1])), 
    group=as.factor(rep(as.vector(hclust_groups), 
    each=dim(res$posterior$q[sub_hclust,])[1])))
                                         
  g <- ggplot(q_full, aes(x=i,y=qi, fill=group)) + 
    geom_boxplot(outlier.size=0.5) +
    theme_classic()+ ylab("density") + xlab("Type effect i")
  print(g)
  }
  
## Example: Dendrogram of type effect grouping ###########################
\dontrun{
  library(gplots)
  library(cluster)
  
  # subset your data by removing a burn in (of 1000) and thinning to 
  # have a maximum of 10000 values as this can be very slow
  sub_hclust <- seq(1000, 5000, by=ceiling(5000/10000))
  groups <- as.data.frame(t(res$posterior$cluster)) 
  groups <- as.data.frame(apply(groups, 2, function(x) as.factor(x)))
  rownames(groups) <- 1:50
  
  # compute dissimilarity matrix for the type effect clusters
  disim_clust_g <- daisy(groups[,sub_hclust])
  clu <- hclust(disim_clust_g, method="complete") 
  dend <- as.dendrogram(clu) 
  
  # OPTIONAL: change the colour of the heatmap. The lighter the colour,
  # the less likely two type effects are in the same group
  hmcols <- colorRampPalette(c("blue","white"))(299)
  # OPTIONAL: change the layout of the heatmap and dendrogram
  lmat <- rbind(c(0,3),c(2,1),c(0,4))
  lwid <- c(0.1,4)
  lhei <- c(1.5,4,0.4)
  layout(mat=lmat, widths=lwid, heights=lhei)

  heatmap_data <- as.matrix(disim_clust_g)

  rownames(heatmap_data) <- 1:50
  colnames(heatmap_data) <- 1:50
  
  heatmap.2(heatmap_data,
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            col=hmcols,           # use on color palette defined earlier 
            dendrogram="col",     # only draw a row dendrogram
            Colv=dend, 
            Rowv=dend,
            lmat=lmat, lwid=lwid, lhei=lhei
  )   
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Dirichlet Process }
\keyword{ Source Attribution }% __ONLY ONE__ keyword per line
