Acceptance = R6::R6Class("acceptance",
                         public = list(
                           alpha = NULL,
                           r = NULL,
                           initialize = function(nSources,
                                                 nTimes,
                                                 nLocations,
                                                 nTypes,
                                                 namesSources,
                                                 namesTimes,
                                                 namesLocations,
                                                 namesTypes,
                                                 paramsFix) {
                             if (paramsFix$alpha == FALSE) {
                               self$alpha <- array(
                                 dim = c(nSources,
                                         nTimes,
                                         nLocations),

                                 dimnames = list(
                                   source = namesSources,
                                   time = namesTimes,
                                   location = namesLocations
                                 )
                               )
                             }

                             if (paramsFix$r == FALSE) {
                               self$r <- array(
                                 dim = c(nTypes,
                                         nSources,
                                         nTimes),
                                 dimnames = list(
                                   type = namesTypes,
                                   source = namesSources,
                                   time = namesTimes
                                 )
                               )
                             }
                           }
                         ))
