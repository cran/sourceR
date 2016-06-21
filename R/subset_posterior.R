subset_posterior <- function(object, params, t, l, i, j, iters) {
  object_subset <- list() #lapply(params, function(x) object$posterior[[x]])
  object_subset$posterior <- NULL
  
  check_names <- function(p, full, error_msg) {
    if (!all(p %in% full)) {
      stop(error_msg)
    }
  }
  
  if (missing(iters)) {
    n <- get_max_iter(object)
    iters <- c(1 : n)
  }
#   if (missing(t)) {
#     time_names <- calc_t_names(object)
#   }
#   if (missing(l)) {
#     location_names <- calc_l_names(object)
#   }
#   if (missing(j)) {
#     j <- calc_source_names(object)
#   }
#   if (missing(i)) {
#     i <- calc_type_names(object)
  # }
  
  
  if (!is.atomic(params)) stop("params must be a vector.")
  all_vals <- c("a", "q", "r", "li", "lj", "cluster", "theta", "d")
  get_trues <- sapply(all_vals, function(x) is.element(x, params))
  if (length(params) != sum(get_trues) || length(params) == 0) stop("params must only contain the values \"a\", \"q\", \"r\", \"li\", \"lj\", \"cluster\", \"theta\", \"d\" (or a subset of these values).")
  
  ##################################### q ######################################
  if ("q" %in% params) {
    if (!missing(i)) { # if missing i return all 
      i_names <- paste("type", i, sep = "")
      check_names(i_names, colnames(object$posterior$q), "The type names specified do not match the names of the types.")
      object_subset$posterior$q <- object$posterior$q[iters, i_names]
    } else {
      object_subset$posterior$q <- object$posterior$q[iters, ]
    }
  }
  
  ################################## cluster ###################################
  if ("cluster" %in% params) {
    if (!missing(i)) { # if missing i return all 
      i_names <- paste("type", i, sep = "")
      check_names(i_names, colnames(object$posterior$cluster), "The type names specified do not match the names of the types.")
      object_subset$posterior$cluster <- object$posterior$cluster[iters, i_names]
    } else {
      object_subset$posterior$cluster <- object$posterior$cluster[iters, ]
    }
  }
  ################################### theta ####################################
  if ("theta" %in% params) {
      object_subset$posterior$theta <- object$posterior$theta[iters, ]
  }
  
  ##################################### a ######################################
  if ("a" %in% params) {
    if (missing(t)) {
      time_names <- names(object$posterior$a)
    } else {
      time_names <- paste("time", t, sep = "")
      check_names(time_names, names(object$posterior$a), 
                  paste("The time names t must only contain the values \n", 
                        paste(names(object$posterior$a), collapse = ", "), "\n"))
    }
    
    if (missing(l)) {
      location_names <- names(object$posterior$a[[1]])
    } else {
      location_names <- paste("location", l, sep = "")
      check_names(location_names, names(object$posterior$a[[1]]), 
                  paste("The location names l must only contain the values \n", 
                        paste(names(object$posterior$a[[1]]), collapse = ", "), "\n"))
    }
    
    if (missing(j)) {
      j <- colnames(object$posterior$a[[1]][[1]])
    } else { 
      check_names(j, colnames(object$posterior$a[[1]][[1]]), 
                  paste("The source names j must only contain the values \n", 
                        paste(colnames(object$posterior$a[[1]][[1]]), collapse = ", "), "\n"))
    }
    object_subset$posterior$a <- lapply(time_names, function(x) lapply(location_names, function(y) object$posterior$a[[x]][[y]][iters, j]))
    ## rename list
    names(object_subset$posterior$a) <- time_names
    for (t1 in 1 : length(time_names)) {
        names(object_subset$posterior$a[[t1]]) <- location_names
    }
    }
  
  ##################################### r ######################################
  if ("r" %in% params) {
    if (missing(t)) {
      time_names <- names(object$posterior$r)
    } else {
      time_names <- paste("time", t, sep = "")
      check_names(time_names, names(object$posterior$r), 
                  paste("The time names t must only contain the values \n", 
                        paste(names(object$posterior$r), collapse = ", "), "\n"))
    }
    
    if (missing(j)) {
      j <- dimnames(object$posterior$r[[1]])[[2]]
    } else {
      check_names(j, dimnames(object$posterior$r[[1]])[[2]], 
                  paste("The source names j must only contain the values \n", 
                        paste(dimnames(object$posterior$r[[1]])[[2]], collapse = ", "), "\n"))
    }
    
    if (missing(i)) {
      i_names <- dimnames(object$posterior$r[[1]])[[1]]
    } else {
      i_names <- paste("type", i, sep = "")
      check_names(i_names, dimnames(object$posterior$r[[1]])[[1]], 
                  paste("The type names i must only contain the values \n", 
                        paste(dimnames(object$posterior$r[[1]])[[1]], collapse = ", "), "\n"))
    }
    
    object_subset$posterior$r <- lapply(time_names, function(x) object$posterior$r[[x]][i_names, j, iters])
    ## rename list
    names(object_subset$posterior$r) <- time_names
  }
  
  ##################################### li #####################################
  if ("li" %in% params) {
    if (missing(t)) {
      time_names <- names(object$posterior$li)
    } else {
      time_names <- paste("time", t, sep = "")
      check_names(time_names, names(object$posterior$li), 
                  paste("The time names t must only contain the values \n", 
                        paste(names(object$posterior$li), collapse = ", "), "\n"))
    }
    
    if (missing(l)) {
      location_names <- names(object$posterior$li[[1]])
    } else {
      location_names <- paste("location", l, sep = "")
      check_names(location_names, names(object$posterior$li[[1]]), 
                  paste("The location names l must only contain the values \n", 
                        paste(names(object$posterior$li[[1]]), collapse = ", "), "\n"))
    }
    
    if (missing(i)) {
      i_names <- colnames(object$posterior$li[[1]][[1]])
    } else {
      i_names <- paste("type", i, sep = "")
      check_names(i_names, colnames(object$posterior$li[[1]][[1]]), 
                  paste("The type names i must only contain the values \n", 
                        paste(colnames(object$posterior$li[[1]][[1]]), collapse = ", "), "\n"))
    }
    
    object_subset$posterior$li <- lapply(time_names, function(x) lapply(location_names, function(y) object$posterior$li[[x]][[y]][iters, i_names]))
    ## rename list
    names(object_subset$posterior$li) <- time_names
    for (t1 in 1 : length(time_names)) {
      names(object_subset$posterior$li[[t1]]) <- location_names
    }
  }
  
  ##################################### lj #####################################
  if ("lj" %in% params) {
    if (missing(t)) {
      time_names <- names(object$posterior$lj)
    } else {
      time_names <- paste("time", t, sep = "")
      check_names(time_names, names(object$posterior$lj), 
                  paste("The time names t must only contain the values \n", 
                        paste(names(object$posterior$lj), collapse = ", "), "\n"))
    }
    
    if (missing(l)) {
      location_names <- names(object$posterior$lj[[1]])
    } else {
      location_names <- paste("location", l, sep = "")
      check_names(location_names, names(object$posterior$lj[[1]]), 
                  paste("The location names l must only contain the values \n", 
                        paste(names(object$posterior$lj[[1]]), collapse = ", "), "\n"))
    }
    
    if (missing(j)) {
      j <- colnames(object$posterior$lj[[1]][[1]])
    } else {
      check_names(j, colnames(object$posterior$lj[[1]][[1]]), 
                  paste("The source names j must only contain the values \n", 
                        paste(colnames(object$posterior$lj[[1]][[1]]), collapse = ", "), "\n"))
    }

    object_subset$posterior$lj <- lapply(time_names, function(x) lapply(location_names, function(y) object$posterior$lj[[x]][[y]][iters, j]))
    ## rename list
    names(object_subset$posterior$lj) <- time_names
    for (t1 in 1 : length(time_names)) {
      names(object_subset$posterior$lj[[t1]]) <- location_names
    }
  }
  
  return(object_subset)
}