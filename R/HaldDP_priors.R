HaldDP_priors = R6::R6Class(
  "priors",
  public = list(
    a_theta = NULL,
    b_theta = NULL,
    a_r = NULL,
    a_alpha = NULL,
    initialize = function(priors) {
      self$a_theta <- priors$a_theta ## single number
      self$b_theta <- priors$b_theta ## single number
      self$a_r <-
        priors$a_r ## DP prior for each type, source, time, location
      self$a_alpha <-
        priors$a_alpha ## DP prior for each source, time, location
    }
  )
)
