#####################################################
# Name: interface.R                                 #
# Author: Poppy Miller <p.miller@lancaster.ac.uk>   #
# Created: 2016-06-15                               #
# Copyright: Poppy Miller 2016                      #
# Purpose: Source attribution model interface       #
#####################################################

#' Runs the HaldDP source attribution model
#'
#' @docType class
#' @name HaldDP
#' @importFrom R6 R6Class
#' @importFrom grDevices col2rgb colorRampPalette
#' @importFrom stats median
#' @import dplyr
#' @export
#'
#' @return Object of \code{\link{HaldDP}} with methods for creating a HaldDP model,
#' running the model, and accessing and plotting the results.
#' @format Object of \code{\link{R6Class}} with methods for creating a HaldDP model,
#' running the model, and accessing and plotting the results.
#' @section Description:
#' This function fits a non-parametric Poisson source attribution model for human cases of
#' disease. It supports multiple types, sources, times and locations. The number of
#' human cases for each type, time and location follow a Poisson likelihood.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(data, k, priors, a_q, inits = NULL)}}{Constructor takes
#'
#'   \code{data} dataframe (with columns containing the number of human cases named
#'   \code{Human}, columns containing the number of positive source samples
#'   (one column per source), a column with the time id's named \code{Time},
#'   a column with the type id's named \code{Type}, and a column with the source
#'   location id's \code{Location}). The data for the human cases and source
#'   counts must be integers. The data for the time, location and type columns
#'   must be factors. The source counts are currently only allowed to vary over time,
#'   hence they must be repeated for each location within each time.
#'
#'   \code{k} prevalence dataframe (with columns named \code{Value, Time,
#'   Location and Source}). Prevalences must be between 0 and 1 as they are the
#'   proportion of samples that were positive for any type for a given source and time.
#'
#'   \code{priors} list with elements named \code{a_r}, \code{a_alpha}, \code{a_theta} and \code{b_theta},
#'   corresponding to the prior parameters for the \code{r}, \code{alpha}, and base
#'   distribution for the DP parameters respectively.
##   The \code{r} parameters have
##   a Dirichlet prior for each type/ time combination, whilst the \code{alpha}
##   parameters have a Dirichet prior for each time/ location combination. Therefore,
##   the prior values \code{a_alpha} can be either a single positive number (to be
##   used for all \code{alphas}) or a data frame with a single positive number for
##   each \code{alpha} parameter (in a column called \code{Value}) and columns named
##   \code{Time}, \code{Location} and \code{Source} containing the time, location and
##   source names for each \code{alpha} prior parameter. Similarly, the prior values \code{a_r}
##   can either be a single positive number or a data frame with columns named \code{Value},
##   \code{Time}, \code{Type} and \code{Source}.
##   The base distribution is Gamma distributed with shape and rate parameters given by \code{a_theta}
##   and \code{b_theta} respectively.
#'
#'   \tabular{lllll}{
#'   \emph{Parameter} \tab \emph{Prior Distribution} \tab \emph{Prior Parameters}\cr
#'   \code{a_r} \tab Dirichlet(concentration) \tab A single positive number or a data \cr
#'   \tab \tab frame with columns giving the prior values \cr
#'   \tab \tab (named \code{Value}), times (named Time) \cr
#'   \tab \tab and source ids (named \code{Source}. If a\cr
#'   \tab \tab single number is supplied, it will be used for\cr
#'   \tab \tab all times, sources an locations. \cr
#'
#'   \code{a_alpha} \tab Dirichlet(concentration) \tab A single positive number or a dataframe \cr
#'   \tab \tab with columns giving the prior values (named  \cr
#'   \tab \tab \code{value}), times (named the name of the \cr
#'   \tab \tab time column in the data), locations (named the \cr
#'   \tab \tab name of the location column in the data) and \cr
#'   \tab \tab the source id (named \code{source_id}). \cr
#'
#'   Type effects base \tab DPM(Gamma(shape, rate), \tab Numerical vector of length 2 for the shape and \cr
#'   distribution (\code{theta}) \tab alpha)\tab rate of the Gamma base distribution.\cr
#'   }
#'
#'   \code{a_q} the Dirichlet Process concentration parameter.
#'
#'   \code{inits} (optional) initial values for the mcmc algorithm. This is a list
#'   that may contain any of the following items: \code{alpha} (a data frame with
#'   columns named \code{Value} contining the initial values, \code{Time},
#'   \code{Location}, \code{Source}), \code{q} (a data frame with columns names
#'   \code{Value} contining the initial values and \code{Type}), and \code{r} (a data
#'   frame a column with the initial r values named \code{Value} (note these must
#'   sum to 1 for each source-time combination), a column with the source id's
#'   named \code{Source}, a column with the time id's named \code{Time}, a column
#'   with the type id's named \code{Type}.)
#'   An optional list giving the starting values for the parameters.
#'   \tabular{lll}{
#'   \emph{Parameter} \tab \emph{Description} \cr
#'   \code{r}
#'   \tab A data frame with columns giving the initial values (named \code{Value}),\cr
#'   \tab times (named Time) and source and type id's (named \code{Source} and Type. \cr
#'   \tab DEFAULT: the default initial values are the maximum likelihood point \cr
#'   \tab estimates of \code{r} from the source matrix (i.e. \eqn{r_ij = x_ij / sum_i=1^n x_ij}).\cr
#'   Source effects (\code{alpha})
#'   \tab A data frame with columns named \code{Value} (containing the initial values), \cr
#'   \tab \code{Source} (containing the source names) and columns giving the time and \cr
#'   \tab location for each parameter (named Location). DEFAULT: The default initial values\cr
#'   \tab for the source effects are drawn the prior distribution (Dirichlet). \cr
#'   Type effects (\code{q})
#'   \tab A data frame with columns giving the initial values (named \code{Value})\cr
#'   \tab and the type ids (named Type). DEFAULT: initialise all type effects to be in \cr
#'   \tab a single group with a theta value calculated as \cr
#'   \tab \eqn{\theta = sum(Human_itl) / sum_l=1^L(sum_t=1^T(sum_i=1^n(sum_j=1^m(alpha_jtl * r_ijt * k_jt))))}. \cr
#'   \tab i.e. \eqn{theta = sum(Human_itl) / sum(lambda_ijtl / theta)}}
#'   }
#'
#'   \item{\code{fit_params(n_iter = 1000, burn_in = 0, thin = 1,
#'   n_r = ceiling(private$nTypes * 0.2), params_fix = NULL)}}{when called, sets the mcmc
#'   parameters.
#'
#'   \code{n_iter} sets the number of iterations returned (after removing
#'   \code{burn_in} and thinning results by \code{thin} i.e. a total of
#'   (n_iter * thin) + burn_in iterations are run)
#'
#'   \code{n_r} is a positive
#'   integer that sets the total number of \code{r_{ijtl}} parameters to be updated
#'   at each time-location-source combination (the default is 20 percent updated
#'   per iteration)
#'
#'   \code{params_fix} is a list with a logical value for any of the model parameters
#'   (any of \code{"alpha", "r", or "q"}). For each parameter, if set to \code{TRUE},
#'   the parameters will be fixed at their initial values and not updated.}
#'
#'   \item{\code{update(n_iter, append = TRUE)}}{when called, updates the \code{HaldDP}
#'   model by running \code{n_iter} iterations.
#'
#'   If missing \code{n_iter}, the \code{n_iter} last set using \code{fit_params()}
#'   or \code{update()} is used.
#'
#'   \code{append}
#'   is a logical value which determines whether the next \code{n_iter} iterations
#'   are appended to any previous iterations, or overwrites them. When
#'   \code{append = TRUE}, the starting values are the last iteration and no
#'   \code{burn_in} is removed. Running the model for the first time, or changing any
#'   model or fitting parameters will set \code{append = FALSE}. }
#'
#'   \item{\code{print_data}}{returns a list containing the human data \code{y}
#'   (an array y[types, times, locations]), the source data \code{X} (an array X[types, sources, times]),
#'   the prevalence data (an array k[sources, times]), the type names, source names,
#'   time names, location names and number of different types, sources, times and locations.
#'   }
#'
#'   \item{\code{print_priors}}{returns a list containing the DP concentration
#'   parameter \code{a_q}, and the priors (R6 class with members named \code{a_alpha}
#'   (members are array \code{a_alpha[sources, times, locations]}), \code{a_r} (an array
#'   \code{a_r[types, sources, times]}), \code{a_theta} and \code{b_theta}).}
#'
#'   \item{\code{print_inits}}{returns an R6 class holding the initial values
#'   (members are \code{alpha} (an array \code{alpha[sources, times, locations]}),
#'   \code{theta} (an array \code{theta[types, iters]}), \code{s} (an array
#'   \code{s[types, iters]}), and \code{r} (an array \code{r[types, sources, times]})).}
#'
#'   \item{\code{print_fit_params}}{returns a list of fitting parameters (\code{n_iter},
#'   \code{append}, \code{burn_in}, \code{thin}, \code{params_fix} (R6 class with members
#'   \code{alpha}, \code{q}, \code{r})).}
#'
#'   \item{\code{print_acceptance}}{returns an R6 class containing the acceptance
#'   rates for each parameter (members are \code{alpha} (an array \code{alpha[sources, times, locations]}),
#'   and \code{r} (an array \code{r[types, sources, times]})).}
#'
#'   \item{\code{extract(params = c("alpha", "q", "s", "r", "lambda_i", "lambda_j", "lambda_j_prop"),
#'   times = NULL, locations = NULL, sources = NULL, types = NULL, iters = NULL,
#'   flatten = FALSE, drop = TRUE)}}{returns a list contining a subset of the parameters
#'   (determined by the \code{params} vector, \code{times}, \code{locations}, \code{sources}, \code{types} and \code{iters}).
#'
#'   If \code{flatten} is set to \code{TRUE}, it returns a dataframe with 1 column per
#'   parameter, otherwise it returns a list containing \code{params} containing a
#'   subset of the following arrays: \code{alpha[Sources, Times, Locations, iters]}, \code{q[Types, iters]},
#'   \code{s[Types, iters]}, \code{r[Types, Sources, Times, iters]},
#'   \code{lambda_i[Types, Times, Locations, iters]},
#'   \code{lambda_j[Sources, Times, Locations, iters]}.
#'
#'   \code{drop}
#'   determines whether to delete the dimensions of an array which have only one
#'   level when \code{flatten = FALSE}.}
#'
#'   \item{\code{summary(alpha = 0.05, params = c("alpha", "q", "s", "r", "lambda_i",
#'   "lambda_j" ,"lambda_j_prop"), times = NULL, locations = NULL, sources = NULL,
#'   types = NULL, iters = NULL, flatten = FALSE, drop = TRUE, CI_type = "chen-shao")}}{
#'   returns a list contining the
#'   median and credible intervals for a subset of the parameters. The default credible
#'   interval type are Chen-Shao (\code{"chen-shao"}) highest posterior density intervals (alternatives
#'   are \code{"percentiles"} and \code{"spin"}).
#'   See \code{extract} for details on the subsetting. \code{lambda_j_prop} returns the
#'   proportion of cases attributed to each source \code{j} and is calculated by dividing
#'   each iteration of \code{lambda_{jtl}} values by their sum within each time \code{t}
#'   and location \code{l}.}
#'
#'   \item{\code{plot_heatmap(iters, cols = c("blue","white"), hclust_method = "complete")}}{
#'   Creates a dendrogram and heatmap for the type effect groupings (\code{s} parameter
#'   in the model). This uses the heatmap.2 function from gplots.
#'
#'   \code{iters} is a vector containing the iterations to be used in constructing
#'   the graph. Default is all iterations in posterior.
#'
#'   \code{hclust_method} allows the user to select the method used by \code{stats::hclust} to
#'   cluster the type effect groupings \code{s}.
#'
#'   \code{cols} gives the colours for completely dissimilar (dissimilarity value
#'   of 1), and identical (dissimilarity value of 0). All other values will be in
#'   between the two chosen colours. See ?colorRampPalette for more details..}
#' }
#'
#' @section Details:
#' \describe{
#' This function fits a source attribution model for human cases of disease.
#' It supports multiple types, sources, times and locations. The number of human cases
#' for each type, time and location follows a Poisson or Negative Binomial likelihood.
#' \emph{Model}
#' \deqn{y_{itl}\sim\textsf{Poisson}(\lambda_{itl})}
#' where
#' \deqn{\lambda_{itl}=\sum_{j=1}^{m}\lambda_{ijtl}=q_{k(i)}\sum_{j=1}^{m}(r_{ijt}\cdot k_{j}\cdot alpha_{jtl})}
#'
#' The parameters are defined as follows:
#' \deqn{a_{jtl}} is the unknown source effect for source \eqn{j}, time \eqn{t}, location \eqn{l}
#' \deqn{q_{s(i)}} is the unknown type effect for type \eqn{i} in group \eqn{s}.
#' \deqn{x_{ij}} is the known number of positive samples for each source \eqn{j} type\eqn{i} combination
#' \deqn{n_{ij}} is the known total number of samples for each source \eqn{j} type \eqn{i} combination
#' \deqn{k_{j}} is the fixed prevalence in source (i.e. the number of positive samples
#' divided by the number of negative samples) \eqn{j}
#' \deqn{r_{ijt}}  is the unknown relative occurrence of type \eqn{i} on source \eqn{j}.
#'
#' \emph{Priors}
#' \deqn{r_{.jt}\sim Dirichlet(a\_r_{1jt},..., a\_r_{njt})}
#' \deqn{a_{tl}\sim Dirichlet(a\_alpha_{1tl},..., a\_alpha_{mtl})}
#' \deqn{q\sim DP(a_q, Gamma(a_{theta},b_{theta}))}
#' }
#'
#' @references Chen, M.-H. and Shao, Q.-M. (1998). Monte Carlo estimation of Bayesian
#' credible and HPD intervals, \emph{Journal of Computational and Graphical Statistics}, 7.
#' @references Liu Y, Gelman A, Zheng T (2015). "Simulation-efficient shortest probability
#' intervals." Statistics and Computing.
#' @author Chris Jewell and Poppy Miller \email{p.miller at lancaster.ac.uk}
#'
#' @examples
#' data(campy)
#' zero_rows <- which(apply(campy[,c(2:7)], 1, sum) == 0)
#' campy <- campy[-zero_rows,]
#' prevs <- data.frame(Value = c(181/ 239, 113/196, 109/127,
#'                               97/595, 165/552, 86/524),
#'                     Source = colnames(campy[, 2:7]),
#'                     Time = rep(1, 6),
#'                     Location = rep("A", 6))
#' priors <- list(a_alpha = 1, a_r = 0.1, a_theta = 0.01, b_theta = 0.00001)
#' res <- HaldDP$new(data = campy, k = prevs, priors = priors, a_q = 0.1)
#' res$fit_params(n_iter = 100, burn_in = 10, thin = 1)
#' res$update()
#'
#' dat <- res$print_data()
#' init <- res$print_inits()
#' prior <- res$print_priors()
#' acceptance <- res$print_acceptance()
#' fit_params <- res$print_fit_params()
#'
#' res$plot_heatmap(iters = 10:100, hclust_method = "complete")
#'
#' summarys <- res$summary(params = c("alpha", "q", "lambda_i"),
#'             times = "1", sources = c("ChickenA", "Bovine"),
#'             iters = 10:100, flatten = TRUE, CI_type = "chen-shao")
#'
#' posteriors <- res$extract(params = c("alpha", "r", "q", "lambda_j"),
#'             sources = c("ChickenB", "Ovine"),
#'             types = c("474", "52"),
#'             iters = 50:100, drop = FALSE, flatten = FALSE)

HaldDP <- R6::R6Class(
  "HaldDP",
  private = list(
    ## DATA
    y = NULL,
    X = NULL,
    ## 3D array [sources, type, time] giving the number of positive samples
    k = NULL,

    nTypes = NULL,
    nSources = NULL,
    nTimes = NULL,
    nLocations = NULL,

    namesTypes = NULL,
    namesSources = NULL,
    namesTimes = NULL,
    namesLocations = NULL,

    priors = NULL,
    inits = NULL, # Chain starting values

    ## MCMC PARAMETERS
    n_iter = NULL,
    n_iter_old = NULL,
    append = NULL,
    burn_in = NULL,
    thin = NULL,
    n_r = NULL,

    params_fix_class = R6::R6Class(
      "params_fix_class",
      public = list(
        q = NULL,
        alpha = NULL,
        r = NULL,
        initialize = function(params_fix)
        {
          self$q = params_fix$q
          self$alpha = params_fix$alpha
          self$r = params_fix$r
        }
      )
    ),
    params_fix = NULL,

    ## MODEL and RESULTS
    a_q = NULL,
    DPModel_impl = NULL,
    posterior = NULL,
    acceptance = NULL,

    ## MCMC updaters for each parameter
    r_updaters = NULL,
    q_updaters = NULL,
    alpha_updaters = NULL,
    updaters = NULL,
    save_chain_state = function(iter)
    {
      for (time in 1:private$nTimes) {
        for (location in 1:private$nLocations) {
          private$posterior$alpha[, time, location, iter] <-
            private$DPModel_impl$alphaNodes[[time]][[location]]$getData()
        }
        for (sources in 1:private$nSources) {
          # TODO: change to an lapply?
          private$posterior$r[, sources, time, iter] <-
            private$DPModel_impl$rNodes[[time]][[sources]]$getData()
        }
      }
      private$posterior$q[, iter] <-
        private$DPModel_impl$qNodes$getData()
      private$posterior$s[, iter] <- private$DPModel_impl$qNodes$s
    },
    create_posterior = function()
    {
      if (private$append == FALSE |
          isTRUE(all.equal(private$n_iter_old, 0)))
      {
        ## create posterior if no posterior exists or if append = FALSE
        private$posterior <-
          Posterior_HaldDP$new(
            nSources = private$nSources,
            nTimes = private$nTimes,
            nLocations = private$nLocations,
            nTypes = private$nTypes,
            n_iter = private$n_iter,
            namesSources = private$namesSources,
            namesTimes = private$namesTimes,
            namesLocations = private$namesLocations,
            namesTypes = private$namesTypes,
            namesIters = 1:private$n_iter
          )
      } else {
        new_post <- Posterior_HaldDP$new(
          nSources = private$nSources,
          nTimes = private$nTimes,
          nLocations = private$nLocations,
          nTypes = private$nTypes,
          n_iter = private$n_iter - private$n_iter_old,
          namesSources = private$namesSources,
          namesTimes = private$namesTimes,
          namesLocations = private$namesLocations,
          namesTypes = private$namesTypes,
          namesIters = (private$n_iter_old + 1):private$n_iter
        )
        private$posterior$alpha <-
          abind::abind(private$posterior$alpha, new_post$alpha, along = 4)
        names(dimnames(private$posterior$alpha)) <-
          names(dimnames(new_post$alpha))
        new_post$alpha <- NULL
        private$posterior$q <-
          abind::abind(private$posterior$q, new_post$q, along = 2)
        names(dimnames(private$posterior$q)) <-
          names(dimnames(new_post$q))
        new_post$q <- NULL
        private$posterior$s <-
          abind::abind(private$posterior$s, new_post$s, along = 2)
        names(dimnames(private$posterior$s)) <-
          names(dimnames(new_post$s))
        new_post$s <- NULL
        private$posterior$r <-
          abind::abind(private$posterior$r, new_post$r, along = 4)
        names(dimnames(private$posterior$r)) <-
          names(dimnames(new_post$r))
        new_post <- NULL
      }
    },
    assign_updaters = function()
    {
      r_updaters <- list()
      alpha_updaters <- list()
      q_updaters <- list()
      for (time in 1:private$nTimes) {
        alpha_updaters[[time]] <- list()
        for (location in 1:private$nLocations) {
          if (isTRUE(private$params_fix$alpha)) {
            alpha_updaters[[time]][[location]] <- NULL
          } else {
            alpha_updaters[[time]][[location]] <-
              AdaptiveLogDirMRW$new(
                private$DPModel_impl$alphaNodes[[time]][[location]],
                tune = rep(0.4, private$nSources),
                name = paste0(
                  "alpha: ",
                  "time ",
                  private$namesTimes[time],
                  ", location ",
                  private$namesLocations[location]
                )
              ) # TODO: User specificed starting variance?
          }
        }
        r_updaters[[time]] <- list()
        for (sources in 1:private$nSources) {
          # TODO: change to an lapply?
          if (isTRUE(private$params_fix$r)) {
            r_updaters[[time]][[sources]] <- NULL
          } else {
            # method of moments: beta to choose tuning values
            alpha <-
              private$priors$a_r[, sources, time] + private$X[, sources, time]
            alpha_0 <- sum(alpha)
            var_alphas <-
              (alpha * (alpha_0 - alpha)) / ((alpha_0 ^ 2) * (alpha_0 + 1))
            tune_val <- 100 * sqrt(var_alphas)
            r_updaters[[time]][[sources]] <-
              AdaptiveLogDirMRW$new(
                private$DPModel_impl$rNodes[[time]][[sources]],
                toupdate = function()
                  sample(private$nTypes, private$n_r),
                tune = tune_val,
                batchsize = 10,
                # adaptive batch size
                name = paste0(
                  "r: ",
                  "time ",
                  private$namesTimes[time],
                  ", source ",
                  private$namesSources[sources]
                )
              )
          }
        }
      }
      if (isTRUE(private$params_fix$q)) {
        q_updaters <- NULL
      } else {
        q_updaters <- PoisGammaDPUpdate$new(private$DPModel_impl$qNodes)
      }
      private$updaters <-
        unlist(c(q_updaters, r_updaters, alpha_updaters),
               recursive = TRUE)
    },
    set_data = function(data)
    {
      if (!is.data.frame(data))
        stop("data must be a data frame.")
      if (!all(c("Human", "Type", "Time", "Location") %in% colnames(data)))
        stop("The data must have columns names Human, Type, Time, and Location.")
      source_names <-
        colnames(data)[!(colnames(data) %in% c("Human", "Type", "Time", "Location"))]

      data <- na.omit(data)

      ## check all data values are positive, numeric, and whole numbers
      check_data <-
        all(sapply(data[, c("Human", source_names)], function(x)
          return(
            all(isFiniteInteger(x)) & all(x >= 0)
            )
          ))
      if (!check_data)
        stop("All human and source data values must be positive integers.")
      ## check that source data is the same for each location within time
      for (times in 1:length(unique(data$Time))) {
        if (length(unique(data$Location)) > 1) {
          tmp <- list()
          for (locations in 1:length(unique(data$Location))) {
            tmp[[locations]] <- subset(data, Time == unique(data$Time)[times] &
                                         Location == unique(data$Location)[locations])[, !(names(data) %in% c("Human", "Location", "Time"))]
          }
          combs <- combn(length(unique(data$Location)), 2)
          if (!all(sapply(1:ncol(combs) , function(x)
            nrow(setdiff(tmp[[combs[1, x]]], tmp[[combs[2, x]]])) == 0)))
            stop("Source data must be identical for each location within a time.")
        }
      }
      ## check all types have at least 1 source case over all times
      check_non0 <- matrix(NA, nrow = length(unique(data$Type)), ncol = length(unique(data$Time)))
      for (times in 1:length(unique(data$Time))) {
          tmp <- subset(data, Time == unique(data$Time)[times] & Location == unique(data$Location)[1])[, c(source_names)]
          tmp <- apply(tmp, 1, function(x) sum(x) <= 0)
          names(tmp) <- NULL
          if (times > 1) tmp <- tmp[match(subset(data, Time == unique(data$Time)[1] &
                                                   Location == unique(data$Location)[1])$Type,
                                          subset(data, Time == unique(data$Time)[times] &
                                                   Location == unique(data$Location)[1])$Type)] # same order of types
          check_non0[, times] <- tmp
      }
      ## if all source data for a particular time, type combo are 0, then stop
      if (!all(!apply(check_non0, 1, function(x) sum(x) >= length(unique(data$Time)))))
        stop("One or more type-time combinations has 0 source cases for all sources.")

      data$Time <- as.factor(data$Time)
      data$Location <- as.factor(data$Location)
      data$Type <- as.factor(data$Type)

      data$Type <- as.factor(data$Type)
      data$Time <- as.factor(data$Time)
      data$Location <- as.factor(data$Location)

      private$nTypes <- length(unique(data$Type))
      private$nTimes <- length(unique(data$Time))
      private$nLocations <- length(unique(data$Location))
      private$nSources <- length(source_names)

      if (dim(unique(data[, c("Type", "Time", "Location")]))[1] != (private$nTypes * private$nTimes * private$nLocations))
        stop(
          "data must have a non-NA count for the human and source sample columns for each time, type and location combination."
        )

      private$namesTypes <-
        as.factor(paste(gtools::mixedsort(unique(data$Type))))
      private$namesTimes <-
        as.factor(paste(gtools::mixedsort(unique(data$Time))))
      private$namesLocations <-
        as.factor(paste(gtools::mixedsort(unique(data$Location))))
      private$namesSources <-
        as.factor(paste(gtools::mixedsort(source_names)))

      private$y <- array(
        NA,
        dim = c(private$nTypes, private$nTimes, private$nLocations),
        dimnames = list(
          type = private$namesTypes,
          time = private$namesTimes,
          location = private$namesLocations
        )
      )
      private$X <- array(
        NA,
        dim = c(private$nTypes, private$nSources, private$nTimes),
        dimnames = list(
          type = private$namesTypes,
          source = private$namesSources,
          time = private$namesTimes
        )
      )

      ## Extract values from data frame and put into arrays. Surely this is very slow!! TODO: find better way!
      for (time in private$namesTimes) {
        for (type in private$namesTypes) {
          for (sources in private$namesSources) {
            tmp_X <-
              unique(data[which((data$Type == type) &
                                  (data$Time == time)),
                          sources])
            if (length(tmp_X) != 1L)
              stop("The data must have a non-NA row for each time and type for each source.")
            private$X[which(private$namesTypes == type),
                      which(private$namesSources == sources),
                      which(private$namesTimes == time)] <- tmp_X
          }
          for (location in private$namesLocations) {
            tmp_y <-
              data$Human[which((data$Type == type) &
                                 (data$Time == time) &
                                 (data$Location == location))]
            if (length(tmp_y) != 1L)
              stop(
                "The data must have a single non-NA row for each time, type, and location for each source."
              )
            private$y[which(private$namesTypes == type),
                      which(private$namesTimes == time),
                      which(private$namesLocations == location)] <-
              tmp_y
          }
        }
      }
    },
    set_k = function(k)
    {
      if (!is.data.frame(k))
        stop("k must be a data frame.")
      k <- na.omit(k)

      if (!all(c("Value", "Source", "Time") %in% colnames(k)))
        stop("k must have columns named Value, Source and Time.")

      k$Source <- as.factor(k$Source)
      k$Time <- as.factor(k$Time)

      if (all(paste(gtools::mixedsort(unique(k$Time))) != private$namesTimes))
        stop("The times in the Time column of k must be the same as those in the data.")
      if (all(paste(gtools::mixedsort(unique(k$Source))) != private$namesSources))
        stop("The sources in the Source column of k must be the same as those in the data.")

      if (dim(unique(k[, c("Time", "Source")]))[1] != (private$nTimes * private$nSources))
        stop("k must have a single number for each time, source combination.")

      ## check all data values are positive, numeric, and whole numbers
      if (!all(is.finite(k$Value)) |
          !all(k$Value >= 0) | !all(k$Value <= 1))
        stop("All prevalence values must be numbers between 0 and 1.")

      private$k <- array(
        NA,
        dim = c(private$nSources, private$nTimes),
        dimnames = list(
          source = private$namesSources,
          time = private$namesTimes
        )
      )

      ## Extract values from data frame and put into arrays. Surely this is very slow!! TODO: find better way!
      for (time in private$namesTimes) {
        for (sources in private$namesSources) {
          tmp_k <- k$Value[which(k$Time == time & k$Source == sources)]
          if (length(tmp_k) == 0L)
            stop("k must have a non-NA row for each time and source.")
          private$k[which(private$namesSources == sources), which(private$namesTimes == time)] <-
            tmp_k
        }
      }
    },
    set_a_q = function(a_q)
    {
      if (length(a_q) != 1 |
          !isFinitePositive(a_q) |
          a_q <= 0)
        stop("a_q should be a single positive number")
      private$a_q <- a_q
    },
    set_priors = function(priors)
    {
      if (!is.list(priors) |
          !all(c("a_alpha", "a_r", "a_theta", "b_theta") %in% names(priors)))
        stop("The priors must be a list with names a_alpha, a_r, a_theta, and b_theta.")

      ## theta priors
      if (length(priors$a_theta) != 1 |
          !isFinitePositive(priors$a_theta))
        stop("priors$a_theta must be a single positive number.")
      if (length(priors$b_theta) != 1 |
          !isFinitePositive(priors$b_theta))
        stop("priors$b_theta must be a single positive number.")

      ## r priors
      r_df <- is.data.frame(priors$a_r)
      if (isTRUE(r_df)) {
        ## check data frame has the right colnames
        if (!all(c("Value", "Type", "Time", "Source") %in% colnames(priors$a_r)))
          stop("priors$a_r must have columns names Value, Type, Time, and Source.")
        priors$a_r$Time <- as.factor(priors$a_r$Time)
        priors$a_r$Type <- as.factor(priors$a_r$Type)
        priors$a_r$Source <- as.factor(priors$a_r$Source)
        ## check the data frame has the right values in the columns
        if (!all(paste(gtools::mixedsort(unique(
          priors$a_r$Type
        ))) == private$namesTypes))
          stop("The levels of priors$a_r$Type must be the same as those in data$Type.")
        if (!all(paste(gtools::mixedsort(unique(
          priors$a_r$Time
        ))) == private$namesTimes))
          stop("The levels of priors$a_r$Time must be the same as those in data$Time.")
        if (!all(paste(gtools::mixedsort(unique(
          priors$a_r$Source
        ))) == private$namesSources))
          stop("The levels of priors$a_r$Source must be the same as those in data$Source.")
        priors$a_r <- na.omit(priors$a_r)
        ## check the dimensions are correct
        if (dim(unique(priors$a_r[, c("Type", "Time", "Source")]))[1] != (private$nTypes * private$nTimes * private$nSources))
          stop("priors$a_r must have a single number for each time, type and source combination.")

        # !all(private$test_integer(priors$a_r$Value)) |
        if (!all(priors$a_r$Value >= 0))
          stop("priors$a_r$Value must contain only positive numbers")
      } else {
        if (!isTRUE(length(priors$a_r) == 1) |
            !isFinitePositive(priors$a_r))
          stop("priors$a_r must be a data frame or a single number.")
      }

      a_r <- array(
        NA,
        dim = c(private$nTypes, private$nSources, private$nTimes),
        dimnames = list(
          type = private$namesTypes,
          source = private$namesSources,
          time = private$namesTimes
        )
      )

      ## alpha priors
      alpha_df <- is.data.frame(priors$a_alpha)
      if (isTRUE(alpha_df)) {
        ## check data frame has the right colnames
        if (!all(c("Value", "Source", "Time", "Location") %in% colnames(priors$a_alpha)))
          stop("a_alpha must have columns names Value, Source, Time, and Location.")

        priors$a_alpha$Time <- as.factor(priors$a_alpha$Time)
        priors$a_alpha$Location <-
          as.factor(priors$a_alpha$Location)
        priors$a_alpha$Source <- as.factor(priors$a_alpha$Source)

        ## check the data frame has the right values in the columns
        if (!all(paste(gtools::mixedsort(unique(
          priors$ar$Location
        ))) == private$namesLocations))
          stop("The levels of priors$a_alpha$Location must be the same as those in data$Location.")
        if (!all(paste(gtools::mixedsort(unique(priors$ar$Time))) == private$namesTimes))
          stop("The levels of priors$a_alpha$Time must be the same as those in data$Time.")
        if (!all(paste(gtools::mixedsort(unique(
          priors$ar$Source
        ))) == private$namesSources))
          stop("The levels of priors$a_alpha$Source must be the same as those in data$Source.")
        priors$a_alpha <- na.omit(priors$a_alpha)
        ## check the dimensions are correct
        if (dim(unique(priors$a_alpha[, c("Location", "Time", "Source")]))[1] != (private$nLocations * private$nTimes * private$nSources))
          stop(
            "priors$a_alpha must have a single number for each time, location and source combination."
          )

        if (!all(isFiniteInteger(priors$a_alpha$Value)) |
            !all(priors$a_alpha$Value >= 0))
          stop("priors$a_alpha$Value must contain only positive numbers")
      } else {
        if (!isTRUE(length(priors$a_alpha) == 1) |
            !isFinitePositive(priors$a_alpha))
          stop("priors$a_alpha must be a data frame or a single number.")
      }
      a_alpha <- array(
        NA,
        dim = c(private$nSources, private$nTimes, private$nLocations),
        dimnames = list(
          source = private$namesSources,
          time = private$namesTimes,
          location = private$namesLocations
        )
      )

      ## Extract values from priors alpha and r data frames for alpha and r and put into arrays.
      ## Surely this is very slow!! TODO: find better way!
      for (time in private$namesTimes) {
        for (sources in private$namesSources) {
          for (type in private$namesTypes) {
            if (isTRUE(r_df)) {
              tmp_a_r <-
                priors$a_r$Value[which((priors$a_r$Type == type) &
                                         (priors$a_r$Time == time) &
                                         (priors$a_r$Source == sources)
                )]

              if (length(tmp_a_r) != 1L)
                stop("priors$r must have a single non-NA row for each time, type and source.")
            } else {
              tmp_a_r <- priors$a_r
            }
            a_r[which(private$namesTypes == type),
                which(private$namesSources == sources),
                which(private$namesTimes == time)] <- tmp_a_r
          }
          for (location in private$namesLocations) {
            if (isTRUE(alpha_df)) {
              tmp_a_alpha <-
                priors$a_alpha$Value[which((priors$a_alpha$Source == sources) &
                                             (priors$a_alpha$Time == time) &
                                             (priors$a_alpha$Location == location)
                )]
              if (length(tmp_a_alpha) != 1L)
                stop("priors$alpha must have a single non-NA row for each time, location and source.")
            } else {
              tmp_a_alpha <- priors$a_alpha
            }
            a_alpha[which(private$namesSources == sources),
                    which(private$namesTimes == time),
                    which(private$namesLocations == location)] <-
              tmp_a_alpha
          }
        }
      }
      priors$a_alpha <- a_alpha
      priors$a_r <- a_r
      private$priors <- HaldDP_priors$new(priors)
    },
    set_inits = function(inits)
    {
      ## r values
      if (!("r" %in% names(inits))) {
        ## default is the source data matrix
        inits$r <-
          private$X + 0.000001  # Added for numeric stability -- only affects starting values for r
        for (times in 1:private$nTimes) {
          inits$r[, , times] <-
            apply(inits$r[, , times], 2, function(x)
              x / sum(x))
        }
      } else {
        if (!is.data.frame(inits$r))
          stop("inits$r must be a data frame.")
        inits$r <- na.omit(inits$r)
        if (!all(c("Type", "Source", "Time", "Value") %in% colnames(inits$r)))
          stop("inits$r must be a data frame with columns called Type, Source, Time and Value")

        inits$r$Type <- as.factor(inits$r$Type)
        inits$r$Source <- as.factor(inits$r$Source)
        inits$r$Time <- as.factor(inits$r$Time)

        if (!all(paste(gtools::mixedsort(unique(inits$r$Type))) == private$namesTypes) |
            !all(paste(gtools::mixedsort(unique(inits$r$Time))) == private$namesTimes) |
            !all(paste(gtools::mixedsort(unique(inits$r$Source))) == private$namesSources) |
            dim(inits$r)[1] != (private$nTypes * private$nSources * private$nTimes))
          stop(
            "inits$r must be a data frame with columns called Type, Source, Time and Value with one row per combination of type, source and time."
          )

        if (!all(isFiniteInteger(inits$r$Value)) |
            !all(inits$r$Value > 0) | !all(inits$r$Value < 1))
          stop("inits$r$Value must contain only numbers between 0 and 1.")

        inits_r <- array(
          NA,
          dim = c(private$nTypes, private$nSources, private$nTimes),
          dimnames = list(
            type = private$namesTypes,
            source = private$namesSources,
            time = private$namesTimes
          )
        )

        ## Extract values from inits r data frame and put into arrays.
        ## Surely this is very slow!! TODO: find better way!
        for (time in private$namesTimes) {
          for (sources in private$namesSources) {
            for (type in private$namesTypes) {
              tmp_init_r <-
                inits$r$Value[which((inits$r$Type == type) &
                                      (inits$r$Time == time) &
                                      (inits$r$Source == sources))]

              if (length(tmp_init_r) == 0L)
                stop("inits$r must have a non-NA row for each time, type and source.")

              inits_r[which(private$namesTypes == type),
                      which(private$namesSources == sources),
                      which(private$namesTimes == time)] <-
                tmp_init_r
            }
            ## check the initial values sum to 1 within each time and location
            if (!isTRUE(all.equal(sum(inits_r[, private$namesSources == sources, private$namesTimes == time]), 1, tol = 0.000001)))
              stop("inits$r must sum to 1 within each time and source.")
          }
        }
        inits$r <- inits_r
      }

      ## alpha inits
      inits_alpha <- array(
        NA,
        dim = c(private$nSources, private$nTimes, private$nLocations),
        dimnames = list(
          source = private$namesSources,
          time = private$namesTimes,
          location = private$namesLocations
        )
      )
      if (!("alpha" %in% names(inits))) {
        ## Surely this is very slow!! TODO: find better way!
        for (time in private$namesTimes) {
          for (location in private$namesLocations) {
            inits_alpha[, which(private$namesTimes == time),
                        which(private$namesLocations == location)] <-
              c(gtools::rdirichlet(1, private$priors$a_alpha[, time, location]))
          }
        }
      } else {
        if (!is.data.frame(inits$alpha))
          stop("inits$alpha must be a data frame.")
        inits$alpha <- na.omit(inits$alpha)
        if (!all(c("Location", "Source", "Time", "Value") %in% colnames(inits$alpha)))
          stop("inits$alpha must be a data frame with columns called Location, Source, Time and Value")

        inits$alpha$Location <- as.factor(inits$alpha$Location)
        inits$alpha$Source <- as.factor(inits$alpha$Source)
        inits$alpha$Time <- as.factor(inits$alpha$Time)

        if (!all(paste(gtools::mixedsort(unique(
          inits$alpha$Location
        ))) == private$namesLocations) |
        !all(paste(gtools::mixedsort(unique(
          inits$alpha$Time
        ))) == private$namesTimes) |
        !all(paste(gtools::mixedsort(unique(
          inits$alpha$Source
        ))) == private$namesSources) |
        dim(inits$alpha)[1] != (private$nLocations * private$nSources * private$nTimes))
          stop(
            "inits$alpha must be a data frame with columns called Location, Source, Time and Value with one row per combination of location, source and time."
          )
        if (!all(isFinitePositive(inits$alpha$Value)))
          stop("inits$alpha$Value must contain only positive numbers")

        ## Extract values from inits alpha data frame and put into arrays.
        ## Surely this is very slow!! TODO: find better way!
        for (time in private$namesTimes) {
          for (location in private$namesLocations) {
            for (sources in private$namesSources) {
              tmp_init_alpha <-
                inits$alpha$Value[which((inits$alpha$Location == location) &
                                          (inits$alpha$Time == time) &
                                          (inits$alpha$Source == sources)
                )]

              if (length(tmp_init_alpha) == 0L)
                stop("inits$alpha must have a non-NA row for each time, location and source.")

              if (length(tmp_init_alpha) != 1 |
                  !isFinitePositive(tmp_init_alpha))
                stop("All elements of inits$alpha must be positive numbers.")
              inits_alpha[which(private$namesSources == sources),
                          which(private$namesTimes == time),
                          which(private$namesLocations == location)] <-
                tmp_init_alpha
            }
            ## check the initial values sum to 1 within each time and location
            if (!isTRUE(all.equal(sum(inits_alpha[, which(private$namesTimes == time), which(private$namesLocations == location)]), 1)))
              stop("inits$alpha must sum to 1 within each time and location.")
          }
        }
      }
      inits$alpha <- inits_alpha

      ## q inits
      if (!("q" %in% names(inits))) {
        ## default to all in 1 group with mean solved using the initial parameters for r and alpha
        ## sum_{i=1}^n lambda_i = q * sum_{t=1}^T sum_{l=1}^L (sum_{i=1}^n sum_{j=1}^m p_{ijt} a_{jtl})
        ## sum_{i=1}^n y_i = q * sum_{t=1}^T sum_{l=1}^L (sum_{i=1}^n sum_{j=1}^m p_{ijt} a_{jtl})
        ## q = sum_{i=1}^n y_i / sum_{t=1}^T sum_{l=1}^L (sum_{i=1}^n sum_{j=1}^m p_{ijt} a_{jtl})
        sum_alpha_jtl_r_ijt <- 0
        for (time in 1:private$nTimes) {
          for (location in 1:private$nLocations) {
            sum_alpha_jtl_r_ijt <-
              sum_alpha_jtl_r_ijt + sum(inits$r[, , time] %*% inits$alpha[, time, location])
          }
        }
        q_val <- sum(private$y) / sum_alpha_jtl_r_ijt
        inits$theta <- q_val
        inits$s <- rep(1, private$nTypes)
      } else {
        if (!is.data.frame(inits$q))
          stop("inits$q must be a data frame.")
        inits$q <- na.omit(inits$q)
        if (!all(c("Type", "Value") %in% colnames(inits$q)))
          stop("inits$q must be a data frame with columns called Type and Value")
        inits$q$Type <- as.factor(inits$q$Type)
        if (!all(paste(gtools::mixedsort(unique(inits$q$Type))) == private$namesTypes) |
            dim(inits$q)[1] != private$nTypes)
          stop("inits$q must be a data frame with columns called Type and Value with one row per type.")

        if (!all(isFinitePositive(inits$q$Value)))
          stop("inits$q$Value must contain only positive numbers.")

        inits$q <- inits$q[gtools::mixedsort(inits$q$Type),]
        inits$theta <- sort(unique(inits$q$Value))
        inits$s <- as.numeric(as.factor(inits$q$Value))
        names(inits$s) <- private$namesTypes
      }
      private$inits <- HaldDP_inits$new(inits)
    },

    set_niter = function(n_iter)
    {
      ## check n_iter
      if (isFiniteInteger(n_iter) &
          length(n_iter) == 1 & n_iter > 0)
      {
        private$n_iter_old <-
          ifelse(private$append == TRUE, private$n_iter, 0)
        private$n_iter <- n_iter + private$n_iter_old
      } else
      {
        stop("n_iter is not a positive integer.")
      }
    },
    set_append = function(append)
    {
      if (isFiniteLogical(append) & length(append) == 1)
      {
        private$append <- append
      } else
      {
        stop("append is not a logical value")
      }
    },
    set_burn_in = function(burn_in)
    {
      if (!isFiniteInteger(burn_in) | burn_in < 0) {
        stop("burn_in is not a positive integer.")
      } else {
        private$burn_in <- burn_in
      }
    },
    set_thin = function(thin)
    {
      if (!isFiniteInteger(thin) | thin <= 0) {
        stop("thin is not a positive integer.")
      } else {
        private$thin <- thin
      }
    },
    set_n_r = function(n_r)
    {
      if (!isFiniteInteger(n_r) | n_r <= 0 | n_r > private$nTypes) {
        stop("n_r is not a positive real number or it is larger than the number of types.")
      } else {
        private$n_r <- n_r
      }
    },
    set_params_fix = function(params_fix)
    {
      if (is.null(params_fix))
      {
        params_fix$q <- FALSE     ## Logical: T/ F
        params_fix$alpha <- FALSE ## Logical: T/ F
        params_fix$r <- FALSE     ## Logical: T/ F
      } else {
        if ("q" %in% names(params_fix)) {
          if (length(params_fix$q) != 1 |
              !isFiniteLogical(params_fix$q)) {
            stop("params_fix$q does not have length 1. It should be a single logical value: T/ F.")
          }
        } else {
          params_fix$q <- FALSE
        }

        if ("alpha" %in% names(params_fix)) {
          if (length(params_fix$alpha) != 1 |
              !isFiniteLogical(params_fix$alpha)) {
            stop(
              "params_fix$alpha does not have length 1. It should be a single logical value: T/ F."
            )
          }
        } else {
          params_fix$alpha <- FALSE
        }

        if ("r" %in% names(params_fix)) {
          if (length(params_fix$r) != 1 |
              !isFiniteLogical(params_fix$r)) {
            stop("params_fix$r does not have length 1. It should be a single logical value: T/ F.")
          }
        } else {
          params_fix$r <- FALSE
        }
      }

      ## if already have some params fixed, check if any have changed, if yes, wipe posterior
      if (!is.null(private$params_fix$q)) {
        alpha_same <- params_fix$alpha == private$params_fix$alpha
        q_same <- params_fix$q == private$params_fix$q
        r_same <- params_fix$q == private$params_fix$q
        same <-
          ifelse(alpha_same == F |
                   q_same == F | r_same == F, FALSE, TRUE)

        ## if trying to append, but changing parameters
        if (same == FALSE & private$append == TRUE)
        {
          stop("Append must be set to false if changing params_fix, model parameters, or data.")
        }
      }
      private$params_fix <- private$params_fix_class$new(params_fix)
    },

    flatten_alpha = function(object)
    {
      res <- NULL
      ncol_old = 0
      for (times in 1:dim(object)[2]) {
        for (locations in 1:dim(object)[3]) {
          res <- cbind(res, t(object[, times, locations,]))
          colnames(res)[(ncol_old + 1):dim(res)[2]] <-
            paste(
              "alpha",
              colnames(res[, (ncol_old + 1):dim(res)[2]]),
              dimnames(object)$time[times],
              dimnames(object)$location[locations],
              sep = "_"
            )
          ncol_old <- dim(res)[2]
        }
      }
      return(res)
    },
    flatten_q_s = function(object, names)
    {
      res <- t(object)
      colnames(res) <- paste(names, colnames(res), sep = "_")
      return(res)
    },
    flatten_r = function(object)
    {
      res <- NULL
      ncol_old = 0
      for (times in 1:dim(object)[3]) {
        for (types in 1:dim(object)[1]) {
          res <- cbind(res, t(object[types, , times,]))
          colnames(res)[(ncol_old + 1):dim(res)[2]] <-
            paste("r",
                  colnames(res[, (ncol_old + 1):dim(res)[2]]),
                  dimnames(object)$time[times],
                  dimnames(object)$type[types],
                  sep = "_")
          ncol_old <- dim(res)[2]
        }
      }
      return(res)
    },
    flatten_lambda_i = function(object)
    {
      res <- NULL
      ncol_old = 0
      for (times in 1:dim(object)[2]) {
        for (locations in 1:dim(object)[3]) {
          res <- cbind(res, t(object[, times, locations,]))
          colnames(res)[(ncol_old + 1):dim(res)[2]] <-
            paste(
              "lambda",
              colnames(res[, (ncol_old + 1):dim(res)[2]]),
              dimnames(object)$time[times],
              dimnames(object)$location[locations],
              sep = "_"
            )
          ncol_old <- dim(res)[2]
        }
      }
      return(res)
    },
    flatten_lambda_j = function(object)
    {
      res <- NULL
      ncol_old = 0
      for (times in 1:dim(object)[2]) {
        for (locations in 1:dim(object)[3]) {
          res <- cbind(res, t(object[, times, locations,]))
          colnames(res)[(ncol_old + 1):dim(res)[2]] <-
            paste(
              "lambda",
              colnames(res[, (ncol_old + 1):dim(res)[2]]),
              dimnames(object)$time[times],
              dimnames(object)$location[locations],
              sep = "_"
            )
          ncol_old <- dim(res)[2]
        }
      }
      return(res)
    },

    calc_CI = function(x, alpha, CI_type)
    {
      # x is the MCMC output
      if (!is.atomic(alpha) |
          !isFinitePositive(alpha) |
          alpha > 1)
        stop("alpha must be a single logical value between 0 and 1.")

      switch(
        CI_type,
        "chen-shao" = ci_chenShao(x, alpha),
        "percentiles" = ci_percentiles(x, alpha),
        "SPIn" = ci_spin(x, alpha),
        stop(
          "The type of interval specified must be one of: chen-shao, percentiles, or SPIn."
        )
      )
    },

    calc_CI_alpha = function(object, alpha, CI_type)
    {
      res <- array(
        dim = c(dim(object)[1],
                dim(object)[2],
                dim(object)[3],
                3),
        dimnames = list(
          source = dimnames(object)$source,
          time = dimnames(object)$time,
          location = dimnames(object)$location,
          CI = c("median", "lower", "upper")
        )
      )
      for (times in 1:dim(object)[2]) {
        for (locations in 1:dim(object)[3]) {
          for (sources in 1:dim(object)[1]) {
            res[sources, times, locations, ] <-
              private$calc_CI(object[sources, times, locations,],
                              alpha, CI_type)
          }
        }
      }
      return(res)
    },
    calc_CI_q = function(object, alpha, CI_type)
    {
      res <-
        sapply(setNames(1:dim(object)[1], dimnames(object)$type), function(x)
          private$calc_CI(object[x, ], alpha, CI_type))
      rownames(res) <- c("median", "lower", "upper")
      return(res)
    },
    calc_CI_r = function(object, alpha, CI_type)
    {
      res <- array(
        dim = c(dim(object)[1],
                dim(object)[2],
                dim(object)[3],
                3),
        dimnames = list(
          type = dimnames(object)$type,
          source = dimnames(object)$source,
          time = dimnames(object)$time,
          CI = c("median", "lower", "upper")
        )
      )
      for (times in 1:dim(object)[3]) {
        for (sources in 1:dim(object)[2]) {
          for (types in 1:dim(object)[1]) {
            res[types, sources, times, ] <-
              private$calc_CI(object[types, sources, times, ],
                              alpha, CI_type)
          }
        }
      }
      return(res)
    },
    calc_CI_lambda_i = function(object, alpha, CI_type)
    {
      res <- array(
        dim = c(dim(object)[1],
                dim(object)[2],
                dim(object)[3],
                3),
        dimnames = list(
          type = dimnames(object)$type,
          time = dimnames(object)$time,
          location = dimnames(object)$location,
          CI = c("median", "lower", "upper")
        )
      )
      for (times in 1:dim(object)[2]) {
        for (locations in 1:dim(object)[3]) {
          for (types in 1:dim(object)[1])
            res[types, times, locations, ] <-
              private$calc_CI(object[types, times, locations, ],
                              alpha, CI_type)
        }
      }
      return(res)
    },
    calc_CI_lambda_j = function(object, alpha, CI_type)
    {
      res <- array(
        dim = c(dim(object)[1],
                dim(object)[2],
                dim(object)[3],
                3),
        dimnames = list(
          source = dimnames(object)$source,
          time = dimnames(object)$time,
          location = dimnames(object)$location,
          CI = c("median", "lower", "upper")
        )
      )
      for (times in 1:dim(object)[2]) {
        for (locations in 1:dim(object)[3]) {
          for (sources in 1:dim(object)[1])
            res[sources, times, locations, ] <-
              private$calc_CI(object[sources, times, locations, ],
                              alpha, CI_type)
        }
      }
      return(res)
    },

    check_extract_summary_params = function(params,
                                            times,
                                            locations,
                                            sources,
                                            types,
                                            iters,
                                            flatten)
    {
      if (!mode(params) %in% c("character") |
          length(params) > 7 |
          length(unique(params)) != length(params) |
          !all(
            unique(params) %in% c(
              "alpha",
              "q",
              "s",
              "r",
              "lambda_i",
              "lambda_j",
              "lambda_j_prop"
            )
          )) {
        stop(
          paste(
            "params must be a vector where elements are a subset of ",
            paste(
              c(
                "alpha",
                "q",
                "s",
                "r",
                "lambda_i",
                "lambda_j",
                "lambda_j_prop"
              ),
              collapse = ", "
            ),
            " with no repeated numbers."
          )
        )
      }

      if (is.null(iters)) {
        ## return all iters
        iters <- 1:private$n_iter
      } else {
        if (!mode(iters) %in% c("numeric") |
            !all(isFinitePositive(iters)) |
            !all(isFiniteInteger(iters)) |
            !all(iters <= private$n_iter) |
            !length(unique(iters)) == length(iters)) {
          stop(
            "iters must be a vector where all elements are integers between 0 and n_iter, with no repeated numbers."
          )
        }
      }

      if (is.null(times)) {
        ## return all times
        times <- private$namesTimes
      } else {
        if (!mode(times) %in% c("numeric", "character") |
            length(times) > private$nTimes |
            length(unique(times)) != length(times) |
            !all(unique(times) %in% private$namesTimes)) {
          stop(
            paste(
              "times must be a vector where elements can only be ",
              paste(private$namesTimes, collapse = ", "),
              " with no repeated values."
            )
          )
        }
      }

      if (is.null(locations)) {
        ## return all locations
        locations <- private$namesLocations
      } else {
        if (!mode(locations) %in% c("numeric", "character") |
            length(locations) > private$nLocations |
            length(unique(locations)) != length(locations) |
            !all(unique(locations) %in% private$namesLocations)) {
          stop(
            paste(
              "locations must be a vector where elements can only be ",
              paste(private$namesLocations, collapse = ", "),
              " with no repeated values."
            )
          )
        }
      }

      if (is.null(sources)) {
        ## return all sources
        sources <- private$namesSources
      } else {
        if (!mode(sources) %in% c("numeric", "character") |
            length(sources) > private$nSources |
            length(unique(sources)) != length(sources) |
            !all(unique(sources) %in% private$namesSources)) {
          stop(
            paste(
              "sources must be a vector where elements can only be ",
              paste(private$namesSources, collapse = ", "),
              " with no repeated values."
            )
          )
        }
      }

      if (is.null(types)) {
        ## return all sources
        types <- private$namesTypes
      } else {
        if (!mode(types) %in% c("numeric", "character") |
            length(types) > private$nTypes |
            length(unique(types)) != length(types) |
            !all(unique(types) %in% private$namesTypes)) {
          stop(
            paste(
              "types must be a vector where elements can only be ",
              paste(private$namesTypes, collapse = ", "),
              " with no repeated values."
            )
          )
        }
      }

      if (length(flatten) != 1 |
          !isFiniteLogical(flatten))
        stop("flatten must be a single logical value.")

      ## calculate the lambda's
      if ("lambda_i" %in% params) {
        private$posterior$calc_lambda_i(
          private$n_iter,
          private$nTimes,
          private$nLocations,
          private$nTypes,
          private$namesTimes,
          private$namesLocations,
          private$namesTypes,
          1:private$n_iter,
          private$k
        )
      }
      if ("lambda_j" %in% params) {
        private$posterior$calc_lambda_j(
          private$n_iter,
          private$nSources,
          private$nTimes,
          private$nLocations,
          private$namesSources,
          private$namesTimes,
          private$namesLocations,
          1:private$n_iter,
          private$k
        )
      }

      if ("lambda_j_prop" %in% params) {
        private$posterior$calc_lambda_j_prop(
          private$n_iter,
          private$nSources,
          private$nTimes,
          private$nLocations,
          private$namesSources,
          private$namesTimes,
          private$namesLocations,
          1:private$n_iter,
          private$k
        )
      }

      return(
        list(
          params = params,
          times = times,
          locations = locations,
          sources = sources,
          types = types,
          iters = iters,
          flatten = flatten
        )
      )
    }
  ),
  public = list(
    #TODO : change r init to an x/colsums(x) and remove from model!
    initialize = function(data, k, priors, a_q, inits = NULL)
    {
      ## read in data
      private$set_data(data) # sets y, X, names and number of sources, types, times and locations
      private$set_k(k) # sets prevalences

      ## read in priors
      private$set_priors(priors)

      private$set_a_q(a_q)

      ## read in initial values (must go after priors as uses priors to generate initial values)
      private$set_inits(inits)

      private$DPModel_impl <-
        DPModel_impl$new(
          y = private$y,
          X = private$X,
          R = private$inits$r,
          Time = private$namesTimes,
          Location = private$namesLocations,
          Sources = private$namesSources,
          Type = private$namesTypes,
          prev = private$k,
          a_q = private$a_q,
          a_theta = private$priors$a_theta,
          b_theta = private$priors$b_theta,
          a_r = private$priors$a_r,
          a_alpha = private$priors$a_alpha,
          s = private$inits$s,
          theta = private$inits$theta,
          alpha = private$inits$alpha
        )
      ## if the parameters are changed, any previous results are wiped to keep it consistant
      private$append <- NULL
      if (!is.null(private$posterior))
        private$create_posterior()
    },
    fit_params = function(n_iter = 1000,
                          burn_in = 0,
                          thin = 1,
                          n_r = ceiling(private$nTypes * 0.2),
                          params_fix = NULL)
    {
      private$set_append(FALSE)
      private$set_niter(n_iter)
      private$set_burn_in(burn_in)
      private$set_thin(thin)
      private$set_n_r(n_r)
      private$set_params_fix(params_fix)
      private$create_posterior()
    },
    update = function(n_iter, append = TRUE)
    {
      if (!missing(append)) {
        if (all(is.na(private$posterior$r))) {
          ## if first time running, set append to false
          private$set_append(FALSE)
        } else {
          private$set_append(append)
        }
      }

      ## when appending, we no longer need a burn in
      if (isTRUE(private$append)) {
        private$burn_in = 0
      }

      if (!missing(n_iter)) {
        # when run the first time n_iter is set using fit_params not update
        if (all(is.na(private$posterior$r)) &
            !n_iter == private$n_iter) {
          private$set_niter(n_iter)
          private$n_iter_old = 0
        } else {
          private$set_niter(n_iter)
        }
      } else {
        if (!all(is.na(private$posterior$r))) {
          # if no n_iter given in subsequent runs, run for another private$n_iter updates
          private$set_niter(private$n_iter)
        }
      }

      ## make posterior if append = F, otherwise extend its size
      private$create_posterior()

      ## assign updaters if running for the first time
      if (private$append == FALSE |
          isTRUE(all.equal(private$n_iter_old, 0))) {
        private$assign_updaters()
      }

      ## run mcmc
      ## user feedback and main MCMC loop
      total_iters <-
        (private$n_iter * private$thin) + private$burn_in
      iter_vec <- (private$n_iter_old + 1):total_iters

      saved_num <-
        private$n_iter_old + 1 # number of saved iterations, keeps track of where to put results in posterior

      pb <- txtProgressBar(max = total_iters, style = 3)
      for (i in iter_vec) {
        ## save results if iteration i is > burn in and a multiple of thin
        if (i >= (private$burn_in + iter_vec[1]) &
            (i - private$n_iter_old) %% private$thin == 0 &
            saved_num <= private$n_iter)
        {
          private$save_chain_state(saved_num)
          saved_num <- saved_num + 1
        }

        # Do the updates
        lapply(private$updaters, function(x)
          x$update())
        setTxtProgressBar(pb, i)
      }

      ## mcmc is finished, save and print acceptance rate summary for r and alpha
      private$acceptance <-
        Acceptance$new(
          nSources = private$nSources,
          nTimes = private$nTimes,
          nLocations = private$nLocations,
          nTypes = private$nTypes,
          namesSources = private$namesSources,
          namesTimes = private$namesTimes,
          namesLocations = private$namesLocations,
          namesTypes = private$namesTypes,
          paramsFix = private$params_fix
        )

      sapply(1:length(private$updaters), function(i) {
        tryCatch({
          acceptances <- private$updaters[[i]]$acceptance()
          names <- private$updaters[[i]]$name

          tmp <- strsplit(names, split = " ")[[1]]
          if (tmp[1] == "r:") {
            t_name <- substr(tmp[3], 1, nchar(tmp[3]) - 1)
            s_name <- tmp[5]
            private$acceptance$r[,
                                 which(s_name == dimnames(private$acceptance$r)$source),
                                 which(t_name == dimnames(private$acceptance$r)$time)] <-
              acceptances
          } else if (tmp[1] == "alpha:") {
            t_name <- substr(tmp[3], 1, nchar(tmp[3]) - 1)
            l_name <- tmp[5]
            private$acceptance$alpha[,
                                     which(t_name == dimnames(private$acceptance$alpha)$time),
                                     which(l_name == dimnames(private$acceptance$alpha)$location)] <-
              acceptances
          }

        },
        error = function(e) {
          NULL
        })
      })
      if (private$params_fix$alpha == FALSE) {
        cat("\nalpha acceptance: \n")
        print(private$acceptance$alpha)
      }
      if (private$params_fix$r == FALSE) {
        cat("\nr acceptance: \n")
        print(private$acceptance$r)
      }

    },

    ## Functions to access the data and results
    print_data = function()
    {
      return(
        list(
          y = private$y,
          X = private$X,
          k = private$k,

          namesType = private$namesTypes,
          namesSource = private$namesSources,
          namesTime = private$namesTimes,
          namesLocation = private$namesLocations,

          nTypes = private$nTypes,
          nSources = private$nSources,
          nTimes = private$nTimes,
          nLocations = private$nLocations
        )
      )
    },
    print_priors = function()
    {
      return(list(priors = private$priors,
                  a_q = private$a_q))
    },
    print_inits = function()
    {
      return(private$inits)
    },
    print_fit_params = function()
    {
      return(
        list(
          n_iter = private$n_iter,
          append = private$append,
          burn_in = private$burn_in,
          thin = private$thin,
          params_fix = private$params_fix
        )
      )
    },
    print_acceptance = function()
    {
      return(private$acceptance)
    },

    extract = function(params = c("alpha", "q", "s", "r", "lambda_i", "lambda_j", "lambda_j_prop"),
                       times = NULL,
                       locations = NULL,
                       sources = NULL,
                       types = NULL,
                       iters = NULL,
                       flatten = FALSE,
                       drop = TRUE)

    {
      params_checked <-
        private$check_extract_summary_params(params, times, locations,
                                             sources, types, iters, flatten)
      params <- params_checked$params
      times <- params_checked$times
      locations <- params_checked$locations
      sources <- params_checked$sources
      types <- params_checked$types
      iters <- params_checked$iters
      flatten <- params_checked$flatten

      if (!is.atomic(drop) |
          !isFiniteLogical(drop))
        stop("drop must be a single logical value.")

      if (flatten == FALSE) {
        return(lapply(setNames(params, params), function(x) {
          switch(
            x,
            "alpha" = private$posterior$alpha[sources, times, locations, iters, drop = drop],
            "q"     = private$posterior$q[types, iters, drop = drop],
            "s"     = private$posterior$s[types, iters, drop = drop],
            "r"     = private$posterior$r[types, sources, times, iters, drop = drop],
            "lambda_i" = private$posterior$lambda_i[types, times, locations, iters, drop = drop],
            "lambda_j" = private$posterior$lambda_j[sources, times, locations, iters, drop = drop],
            "lambda_j_prop" = private$posterior$lambda_j_prop[sources, times, locations, iters, drop = drop],
            stop("Unrecognised model component")
          )
        }))
      } else {
        res <- do.call(cbind, lapply(params, function(x) {
          switch(
            x,
            "alpha" = private$flatten_alpha(private$posterior$alpha[sources, times, locations, iters, drop = FALSE]),
            "q"     = private$flatten_q_s(private$posterior$q[types, iters, drop = FALSE], names = "q"),
            "r"     = private$flatten_r(private$posterior$r[types, sources, times, iters, drop = FALSE]),
            "lambda_i" = private$flatten_lambda_i(private$posterior$lambda_i[types, times, locations, iters, drop = FALSE]),
            "lambda_j" = private$flatten_lambda_j(private$posterior$lambda_j[sources, times, locations, iters, drop = FALSE]),
            "lambda_j_prop" = private$flatten_lambda_j(private$posterior$lambda_j_prop[sources, times, locations, iters, drop = FALSE]),
            stop("Unrecognised model component")
          )
        }))
        res <- as.data.frame(res)
        if ("s" %in% params) {
          tmp <-
            as.data.frame(private$flatten_q_s(private$posterior$s[types, iters, drop = FALSE], names = "s"))
          res <- cbind(res, tmp)
        }
        return(res)
      }
    },
    plot_heatmap = function(iters,
                            cols = c("blue", "white"),
                            hclust_method = "complete")
    {
      if (isTRUE(private$params_fix$q)) {
        stop(
          "params_fix$q is TRUE, therefore all values of the marginal posterior for q and s are the same (and a heatmap can't be plotted)."
        )
      }
      if (missing(iters)) {
        ## return all iters
        iters <- 1:private$n_iter
      } else {
        if (!mode(iters) %in% c("numeric") |
            !all(isFinitePositive(iters)) |
            !all(isFiniteInteger(iters)) |
            !all(iters <= private$n_iter) |
            !length(unique(iters)) == length(iters)) {
          stop(
            "iters must be a vector where all elements are integers between 0 and n_iter, with no repeated numbers."
          )
        }
      }

      # Draw the heatmap
      groups <- as.data.frame(private$posterior$s[, iters])
      clusterHeatMap(groups, cols, rownames(private$posterior$q), hclust_method)
    },
    summary = function(alpha = 0.05,
                       params = c("alpha", "q", "s", "r", "lambda_i", "lambda_j", "lambda_j_prop"),
                       times = NULL,
                       locations = NULL,
                       sources = NULL,
                       types = NULL,
                       iters = NULL,
                       flatten = FALSE,
                       CI_type = "chen-shao")

    {
      object <-
        self$extract(params, times, locations, sources, types, iters, flatten, drop = F)
      if (flatten == TRUE) {
        ## do summary on each column of array
        ## remove cols starting with s_ as no summary for s objects
        object <- object[, !grepl("s_", colnames(object))]
        return(do.call(cbind, lapply(setNames(1:dim(object)[2], colnames(object)), function(x)
          private$calc_CI(object[, x], alpha, CI_type))))
      } else {
        ## do each item in list seperately
        return(lapply(setNames(params, params), function(x) {
          switch(
            x,
            "alpha" = private$calc_CI_alpha(object$alpha, alpha, CI_type),
            "q"     = private$calc_CI_q(object$q, alpha, CI_type),
            "r"     = private$calc_CI_r(object$r, alpha, CI_type),
            "lambda_i" = private$calc_CI_lambda_i(object$lambda_i, alpha, CI_type),
            "lambda_j" = private$calc_CI_lambda_j(object$lambda_j, alpha, CI_type),
            "lambda_j_prop" = private$calc_CI_lambda_j(object$lambda_j_prop, alpha, CI_type)
          )
        }))
      }
    }
  )
)
