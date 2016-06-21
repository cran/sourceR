flatten <- function(object) {
  # always have source effects, sometimes have the others

  # number of times
  n_t <- calc_n_times(object)
  # number of locations
  n_l <- calc_n_locations(object)
  # number of sources
  n_sources <- calc_n_sources(object)
  
  n_a <- ifelse("a" %in% names(object$posterior), n_t * n_l * n_sources, 0)
  
  n_q <- n_clust <- n_theta <- ifelse("q" %in% names(object$posterior), 
                                      dim(object$posterior$q)[2], 0)
  
  n_r <- ifelse("r" %in% names(object$posterior), 
                dim(object$posterior$r[[1]])[1] * n_t * n_sources, 0) 
  
  n_li <- ifelse("li" %in% names(object$posterior), 
                 dim(object$posterior$li[[1]][[1]])[2] * n_t * n_l, 0) 
  
  n_lj <- ifelse("lj" %in% names(object$posterior), 
                 dim(object$posterior$lj[[1]][[1]])[2] * n_t * n_l, 0) 

  n_d <- ifelse("d" %in% names(object$posterior), 1, 0)
#   
#   num_params <- n_a + n_q + n_clust + n_theta + n_r + n_li + n_lj + n_d
  
  
  flat_object <- NULL
  if ("a" %in% names(object$posterior)) {
    flat_object <- plyr::ldply(list(object$posterior$a), data.frame)
    names(flat_object) <- paste("a_", names(flat_object), sep = "")
    object$posterior$a <- NULL
    n_named <- n_a
  }
  
  if ("q" %in% names(object$posterior)) {
    flat_object <- cbind(flat_object, object$posterior$q)
    object$posterior$q <- NULL
    vals <- (n_named + 1) : (n_named + n_q)
    names(flat_object)[vals] <- paste("q_", names(flat_object)[vals], sep = "")
    n_named <- n_named + n_q
  }
  
  if ("theta" %in% names(object$posterior)) {
    flat_object <- cbind(flat_object, object$posterior$theta)
    object$posterior$theta <- NULL
    vals <- (n_named + 1) : (n_named + n_theta)
    names(flat_object)[vals] <- paste("theta_", names(flat_object)[vals], sep = "")
    n_named <- n_named + n_theta
  }

  if ("cluster" %in% names(object$posterior)) { 
    flat_object <- cbind(flat_object, object$posterior$cluster)
    object$posterior$cluster <- NULL
    vals <- (n_named + 1) : (n_named + n_clust)
    names(flat_object)[vals] <- paste("cluster_", names(flat_object)[vals], sep = "")
    n_named <- n_named + n_clust
  } 

  # apply(res$posterior$r[[1]], 1, function(y) list(t(y))) creates 1 list per type with n_iter rows and n_sources cols
  # plyr::ldply(list(apply(res$posterior$r[[1]], 1, function(y) list(t(y)))), data.frame) creates a data.frame 
  # with one column for each source-type combination
  
  if ("r" %in% names(object$posterior)) {
    flat_object <- cbind(flat_object, plyr::ldply(list(lapply(object$posterior$r, function(x) 
      plyr::ldply(list(apply(x, 1, function(y) 
        list(t(y)))), data.frame))), data.frame))
    object$posterior$r <- NULL
    vals <- (n_named + 1) : (n_named + n_r)
    names(flat_object)[vals] <- paste("r_", names(flat_object)[vals], sep = "")
    n_named <- n_named + n_r
  }

  if ("li" %in% names(object$posterior)) {
    flat_object <- cbind(flat_object, plyr::ldply(list(object$posterior$li), data.frame))
    object$posterior$li <- NULL
    vals <- (n_named + 1) : (n_named + n_li)
    names(flat_object)[vals] <- paste("li_", names(flat_object)[vals], sep = "")
    n_named <- n_named + n_li
  }
  
  if ("lj" %in% names(object$posterior)) {
    flat_object <- cbind(flat_object, plyr::ldply(list(object$posterior$lj), data.frame))
    object$posterior$lj <- NULL
    vals <- (n_named + 1) : (n_named + n_lj)
    names(flat_object)[vals] <- paste("lj_", names(flat_object)[vals], sep = "")
    n_named <- n_named + n_lj
  }
  
  if ("d" %in% names(object$posterior)) {
    flat_object <- cbind(flat_object, object$posterior$d)
    object$posterior$d <- NULL
    vals <- (n_named + 1) : (n_named + n_d)
    names(flat_object)[vals] <- paste("d_", names(flat_object)[vals], sep = "")
    n_named <- n_named + n_d
  }

  return(flat_object)
}
