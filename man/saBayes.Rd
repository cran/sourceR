% rm sourceAttribution.pdf; R CMD Rd2pdf sourceAttribution creates the pdf of the help files

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
  \emph{Parameter} \tab \emph{Prior Distribution} \tab \emph{Prior Parameters}\cr
  
  \code{r} \tab Dirichlet(concentration) \tab A single positive number or a data \cr
  \tab \tab frame with columns giving the prior values \cr
  \tab \tab (named \code{value}), times (named the name \cr
  \tab \tab of the time column in the data) and source and\cr 
  \tab \tab type id's (named \code{source_id} and the name\cr 
  \tab \tab given to the type id column in the data). If a\cr
  \tab \tab single number is supplied, it will be used for\cr
  \tab \tab all times. \cr
  
  Source effects (\code{a}) \tab Dirichlet(concentration) \tab A single positive number or a dataframe \cr
  \tab \tab with columns giving the prior values (named  \cr
  \tab \tab \code{value}), times (named the name of the \cr
  \tab \tab time column in the data), locations (named the \cr
  \tab \tab name of the location column in the data) and \cr
  \tab \tab the source id (named \code{source_id}). \cr
  
  Type effects base \tab DPM(Gamma(shape, rate), \tab Numerical vector of length 2 for the shape and \cr
  distribution (\code{theta}) \tab alpha)\tab rate of the Gamma base distribution.\cr
  
  Dispersion (\code{d}) \tab Log Normal(meanlog, sdlog) \tab Numerical vector of length 2. Only required if \cr
  \tab \tab using a Negative Binomial likelihood.\cr
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
  
  \code{save_params} 
  \tab Specifies the parameters that should be returned in the \tab Save all parameters \cr
  \tab \code{posterior} list. The lambda parameters can be calculated\tab \code{c("a", "q", "r", "li", "d",} \cr
  \tab  from the other parameters post-hoc. It is quicker and uses\tab \code{"lj", "theta", "cluster"})\cr
  \tab less memory to save the values for fewer parameters. \tab \cr
  
  \code{burn_in} 
  \tab The first k samples are discarded, where k is the value of \tab 0 \cr
  \tab \code{burn_in}. Discarding the burn in when running the model \tab \cr
  \tab reduces storage requirements. \tab \cr
  
  \code{thin} 
  \tab Samples from every kth iteration will be stored, where k is \tab 1 \cr
  \tab the value of \code{thin}. Thinning whilst running the model \tab \cr
  \tab reduces storage requirements. \tab \cr
  
  \code{n_r} 
  \tab The \code{r} matrices for each time are updated using single \tab 40\% of the\cr
  \tab site updates for each \code{r_ijt}. \code{n_r} sets the number of \code{r_ijt}'s \tab r parameters.\cr
  \tab that are updated per iteration.\tab \cr
}
}

\item{initials}{
An optional list giving the starting values for the parameters.
    \tabular{lll}{
  \emph{Parameter} \tab \emph{Description} \cr
  
  \code{r} 
  \tab A single positive number or a data frame with columns giving the prior\cr
  \tab values (named \code{value}), times (named the name of the time column in\cr
  \tab the data) and source and type id's (named \code{source_id} and the name\cr
  \tab given to the type id column in the data). If single number is supplied,\cr
  \tab it will be used for all times. DEFAULT: the default initial values are the \cr
  \tab maximum likelihood point estimates of \code{r} from the source matrix \cr
  \tab (i.e. \eqn{r_ij = x_ij / sum_i=1^n x_ij}).\cr
  
  Source effects (\code{a}) 
  \tab A single positive number (which will be used as the starting value for \cr
  \tab all source effects) or a data frame with columns named \code{value} \cr
  \tab (containing the initial values), \code{source_id} (containing the source \cr
  \tab names which must be the same names as the source columns in the data) and \cr
  \tab columns giving the time and location for each parameter (named the same \cr
  \tab names as the time and location columns in the data). DEFAULT: The default \cr
  \tab initial values for the source effects are drawn the prior distribution \cr
  \tab (Dirichlet). \cr
  
  Type effects (\code{q}) 
  \tab A data frame with columns giving the initial values (named \code{value})\cr
  \tab and the type ids (named the same name as the type column in the data data \cr
  \tab frame). DEFAULT: initialise all type effects to be in a single group with \cr
  \tab a theta value drawn from \eqn{\theta_{1}\sim L(y|\theta_{1},a,p)\times Q_0(\theta_{1})}. \cr
  
  Dispersion parameter (\code{d}) 
  \tab A single positive number. Only required if using a Negative Binomial \cr
  \tab likelihood. DEFAULT: Drawn from the prior (Log Normal with parameters \cr
  \tab specified in the \code{priors} list) for this parameter. \cr
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
  \emph{Parameter} \tab \emph{Format and nesting} \tab \emph{}\cr
  
  Source effects (\code{a}) 
  \tab A nested list with a matrix (dimensions (\code{n_iter}, \tab Always returned. \cr
  \tab \code{sources})) for each time and location. \tab \cr
  \tab eg \code{posterior$a[[t]][[l]][n,j]} gives the value \tab \cr
  \tab for \code{a} \code{j} for time \code{t}, location \code{l} \tab \cr
  \tab at iteration \code{n}. \tab \cr
  
  Type effects (\code{q}) 
  \tab A matrix (dimensions (\code{n_iter}, \code{types})) \tab Always returned. \cr
  \tab containing the \code{theta} value for the cluster\tab \cr
  \tab \code{g} that the type effect {qi} is assigned to \tab \cr
  \tab in iteration \code{n}. eg \code{posterior$q[n,i]} \tab \cr
  \tab gives the value for \code{q} \code{i} at iteration \tab \cr
  \tab \code{n}.\tab \cr
  
  Type effect cluster (\code{cluster}) 
  \tab A matrix containing the group each type effect was \tab Returned if \cr
  \tab assigned to at each iteration (dimensions (\code{n_iter},  \tab \code{likelihood_dist}\cr
  \tab \code{types})). eg \code{posterior$cluster[n,i]}  \tab \code{== "pois"}\cr
  \tab gives the cluster for \code{q} \code{i} at iteration  \tab \cr
  \tab \code{n}.  \tab \cr
  
  Theta (\code{theta}) 
  \tab A matrix containing the type effect parameters for \tab Returned if  \cr
  \tab each group (dimensions (\code{n_iter}, \code{gmax})). \tab \code{likelihood_dist} \cr
  \tab eg \code{posterior$theta[n,c]} gives the \code{theta} \tab \code{ == "pois"}\cr
  \tab value for cluster \code{c} at iteration \code{n}. \tab \cr
  
    \code{r} 
    \tab A nested list containing an array (with dimensions \tab Returned if \cr
    \tab (\code{types}, \code{sources}, \code{n_iter})) for \tab \code{fix_r == FALSE}. \cr
    \tab each time. eg \code{posterior$r[[t]][i,j,n]} gives \tab \cr
    \tab the value for \code{r} (type \code{i}, source \code{j}) \tab \cr
    \tab at time \code{t}, iteration \code{n}. \tab \cr
  
  Dispersion parameter (\code{d}) \tab A vector. \tab Returned if \cr
  \tab \tab \code{likelihood_dist} \cr
  \tab \tab \code{ == "nbinom"}. \cr
  
  Rate of human infection from 
  \tab A nested list with a matrix (dimensions (\code{n_iter}, \tab Returned if \cr
  type \code{i} (Lambda i: \code{li}) \tab \code{types})) for each time and location. \tab  \code{mcmc_params$save_lambda} \cr
  \tab eg \code{posterior$li[[t]][[l]][n,i]} gives the value \tab \code{ == TRUE}\cr
  \tab for \code{li} \code{i} for time \code{t}, location \code{l}  \tab \cr
  \tab at iteration \code{n}. \tab \cr
  
  Rate of human infection from 
  \tab A nested list with a matrix (dimensions (\code{n_iter}, \tab Returned if \cr
  source \code{j} (Lambda j: \code{lj}) 
  \tab \code{sources})) for each time and location. \tab \code{mcmc_params$save_lambda}\cr
  \tab eg \code{posterior$lj[[t]][[l]][n,j]} gives the value \tab \code{ == TRUE}\cr
  \tab for \code{lj} \code{j} for time \code{t}, location \code{l} \tab \cr
  \tab at iteration \code{n}. \tab \cr
}
}
\item{acceptance_rate}{List with items named \code{a}, \code{r}, and \code{d}.

      \tabular{lll}{
  \emph{Parameter} \tab \emph{Format and nesting} \tab \emph{ }\cr
  
  Source effects (\code{a}) 
  \tab A nested list with a positive numerical vector \tab Always returned. \cr
  \tab for each time and location. \tab All should all be \cr
  \tab eg \code{acceptance_rate$a_centered[[t]][[l]][3]}\tab  around 0.45 \cr
  \tab gives the acceptance rate for \code{a} update \tab \cr
  \tab for the third source effect at time \code{t},  \tab \cr 
  \tab location \code{l}. \tab \cr
  
  \code{r} 
  \tab A nested list containing a matrix (with \tab Only returned if \cr
  \tab dimensions(\code{types}, \code{sources})) for \tab \code{fix_r == FALSE}. \cr
  \tab each time. eg \code{acceptance_rate$r[[t]][i,j]} \tab All should all be \cr
  \tab gives the acceptance rate for \code{r} (type \tab around 0.45 \cr
  \tab \code{i}, source \code{j}) at time \code{t}. \cr
  
  Dispersion parameter (\code{d}) \tab positive numerical vector of length 1. \tab Only returned if \cr
  \tab \tab \code{likelihood_dist == "nbinom"} \cr
  \tab \tab This should be \cr
  \tab \tab around 0.45 \cr
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
\code{\link{subset_posterior}, \link{flatten}, \link{summary}}
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
plot(res$posterior$li$time1$locationB[, "type3"], type="l")
## or
plot(res$posterior$li[[1]][[2]][, 3], type="l")

## Example: Density plot #################################################
# density plot for type effect (type = 3)
plot(density(res$posterior$q[,"type3"]))

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
  # if they were split into 5 groups using the tree generated by the
  # hclust command above
  hclust_groups <- cutree(clu, k=5) 
  
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
  lmat <- rbind(c(0,4),c(0,1),c(3,2),c(5,5))
  lwid <- c(0.1,4)
  lhei <- c(1,0.1,4,1)
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
            lmat=lmat, lwid=lwid, lhei=lhei,
            key = F
  )   
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Dirichlet Process }
\keyword{ Source Attribution }% __ONLY ONE__ keyword per line
