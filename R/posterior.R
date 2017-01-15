Posterior_HaldDP = R6::R6Class(
  "Posterior_HaldDP",
  public = list(
    q = NULL,
    s = NULL,
    alpha = NULL,
    r = NULL,
    lambda_i = NULL,
    lambda_j = NULL,
    lambda_j_prop = NULL,
    initialize = function(nSources,
                          nTimes,
                          nLocations,
                          nTypes,
                          n_iter,
                          namesSources,
                          namesTimes,
                          namesLocations,
                          namesTypes,
                          namesIters) {
      self$alpha <- array(
        dim = c(nSources,
                nTimes,
                nLocations,
                n_iter),
        dimnames = list(
          source = namesSources,
          time = namesTimes,
          location = namesLocations,
          iter = namesIters
        )
      )
      self$q <- array(
        dim = c(
          nTypes,
          n_iter
        ),
        dimnames = list(
          type = namesTypes, iter = namesIters
        )
      )

      self$s <- array(
        dim = c(
          nTypes,
          n_iter
        ),
        dimnames = list(
          type = namesTypes, iter = namesIters
        )
      )
      self$r <- array(
        dim = c(nTypes,
                nSources,
                nTimes,
                n_iter),
        dimnames = list(
          type = namesTypes, source = namesSources, time = namesTimes, iter = namesIters
        )
      )
    },
    calc_lambda_i = function(n_iter,
                             nTimes,
                             nLocations,
                             nTypes,
                             namesTimes,
                             namesLocations,
                             namesTypes,
                             namesIters,
                             k) {
      self$lambda_i <- array(
        dim = c(nTypes,
                nTimes,
                nLocations,
                n_iter),
        dimnames = list(
          type = namesTypes,
          time = namesTimes,
          location = namesLocations,
          iter = namesIters
        )
      )
      if (!is.null(self$q)) { # don't try to calculate before any iterations have occured
        for (i in 1:n_iter) {
          for (times in 1:nTimes) {
            for (locations in 1:nLocations) {
              self$lambda_i[, times, locations, i] <- self$q[, i] * self$r[, , times, i] %*% (k[, times] * self$alpha[, times, locations, i])
            }
          }
        }
      } else {
        self$lambda_i <- NULL
      }
    },

    calc_lambda_j = function(n_iter,
                             nSources,
                             nTimes,
                             nLocations,
                             namesSources,
                             namesTimes,
                             namesLocations,
                             namesIters,
                             k) {
      self$lambda_j <- array(
        dim = c(nSources,
                nTimes,
                nLocations,
                n_iter),
        dimnames = list(
          source = namesSources,
          time = namesTimes,
          location = namesLocations,
          iter = namesIters
        )
      )

      if (!is.null(self$q)) { # don't try to calculate before any iterations have occured
        for (i in 1:n_iter) {
          for (times in 1:nTimes) {
            for (locations in 1:nLocations) {
              self$lambda_j[, times, locations, i] <- self$alpha[, times, locations, i] * colSums(self$r[, , times, i] * self$q[, i]) * k[, times]
            }
          }
        }
      } else {
        self$lambda_j <- NULL
      }
    },

    calc_lambda_j_prop = function(n_iter,
                                  nSources,
                                  nTimes,
                                  nLocations,
                                  namesSources,
                                  namesTimes,
                                  namesLocations,
                                  namesIters,
                                  k) {
      if (is.null(self$lambda_j)) self$calc_lambda_j(n_iter,
                                                     nSources,
                                                     nTimes,
                                                     nLocations,
                                                     namesSources,
                                                     namesTimes,
                                                     namesLocations,
                                                     namesIters,
                                                     k)
      self$lambda_j_prop <- self$lambda_j
      for (iter in 1:dim(self$lambda_j)[4]) {
        for (times in 1:nTimes) {
          for (locations in 1:nLocations) {
            self$lambda_j_prop[,times,locations,iter] <-
              self$lambda_j[,times,locations,iter] / sum(self$lambda_j[,times,locations,iter])
          }
        }
      }
    }
  )
)
