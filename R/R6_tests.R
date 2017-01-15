test <- R6::R6Class(
  "test",
  private = list(
    x1 = NULL,
    x2 = NULL,
    x3 = R6::R6Class(
      "x3",
      public = list(
        a = NULL,
        b = NULL,
        initialize = function(x3, x1 = private$x1) {
          private$c = x3$c
          private$d = x3$d
          self$a = x1
          self$b = x3$b
        }
      ),
      private = list(
        c = NULL,
        d = NULL
      )
    )
  ),
  public = list(
    initialize = function(x1 = NA, x2 = NA, x3 = list())
    {
      private$x1 <- x1
      private$x2 <- x2
      if (is.null(private$x3$a)) {
        private$x3 <- private$x3$new(x3, x1)
      } else {
        stop("something is already here!")
      }
    },
    extract = function()
    {
      return(list(x1 = private$x1, x2 = private$x2, x3 = private$x3))
    }
  )
)

res <- test$new(x1 = 1, x2 = 2, list(a=3, b=4, c=5, d=6))

res
res$extract()
