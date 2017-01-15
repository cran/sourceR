HaldDP_inits = R6::R6Class(
  "inits_class",
  public = list(
    theta = NULL,
    s = NULL,
    alpha = NULL,
    r = NULL,
    initialize = function(inits) {
      self$theta <-
        inits$theta    ## vector of length number of groups (i.e. no of unique values in s)
      self$s <- inits$s
      self$alpha <-
        inits$alpha    ## 3D array [source, time, location]
      self$r <-
        inits$r            ## 3D array [source, type, time] giving the relative prevalences
    }
  )
)
