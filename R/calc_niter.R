get_max_iter <- function(object) {
  if ("a" %in% names(object$posterior)) {
    n <- dim(object$posterior$a[[1]][[1]])[1]
  } else if ("lj" %in% names(object$posterior)) {
    n <- dim(object$posterior$lj[[1]][[1]])[1]
  } else if ("q" %in% names(object$posterior)) {
    n <- dim(object$posterior$q)[1]
  } else if ("cluster" %in% names(object$posterior)) {
    n <- dim(object$posterior$cluster)[1]
  } else if ("theta" %in% names(object$posterior)) {
    n <- dim(object$posterior$theta)[1]
  } else if ("li" %in% names(object$posterior)) {
    n <- dim(object$posterior$li[[1]][[1]])[1]
  } else if ("r" %in% names(object$posterior)) {
    n <- dim(object$posterior$r[[1]])[3]
  } else if ("d" %in% names(object$posterior)) {
    n <- length(object$posterior$d)
  } else {
    stop("The object must contain the values \"a\", \"q\", \"r\", \"li\", \"lj\", \"cluster\", \"theta\", \"d\" (or a subset of these values).")
  }
  return(n)
}

calc_n_times <- function(object) {
  if ("a" %in% names(object$posterior)) {
    n <- length(object$posterior$a)
  } else if ("lj" %in% names(object$posterior)) {
    n <- length(object$posterior$lj)
  } else if ("li" %in% names(object$posterior)) {
    n <- length(object$posterior$li)
  } else if ("r" %in% names(object$posterior)) {
    n <- length(object$posterior$r)
  } else {
    stop("The object must contain the values \"a\", \"r\", \"li\", \"lj\" (or a subset of these values).")
  }
  return(n)
}

calc_n_locations <- function(object) {
  if ("a" %in% names(object$posterior)) {
    n <- length(object$posterior$a[[1]])
  } else if ("lj" %in% names(object$posterior)) {
    n <- length(object$posterior$lj[[1]])
  } else if ("li" %in% names(object$posterior)) {
    n <- length(object$posterior$li[[1]])
  } else {
    stop("The object must contain the values \"a\", \"li\", \"lj\" (or a subset of these values).")
  }
  return(n)
}

calc_n_sources <- function(object) {
  if ("a" %in% names(object$posterior)) {
    n <- dim(object$posterior$a[[1]][[1]])[2]
  } else if ("lj" %in% names(object$posterior)) {
    n <- dim(object$posterior$lj[[1]][[1]])[2]
  } else if ("r" %in% names(object$posterior)) {
    n <- dim(object$posterior$r[[1]])[2]
  } else {
    stop("The object must contain the values \"a\", \"r\", \"lj\" (or a subset of these values).")
  }
  return(n)
}

calc_n_types <- function(object) {
  if ("q" %in% names(object$posterior)) {
    n <- dim(object$posterior$q)[2]
  } else if ("cluster" %in% names(object$posterior)) {
    n <- dim(object$posterior$cluster)[2]
  } else if ("theta" %in% names(object$posterior)) {
    n <- dim(object$posterior$theta)[2]
  } else if ("li" %in% names(object$posterior)) {
    n <- dim(object$posterior$li[[1]][[1]])[2]
  } else if ("r" %in% names(object$posterior)) {
    n <- dim(object$posterior$r[[1]])[1]
  } else {
    stop("The object must contain the values \"a\", \"q\", \"r\", \"li\", \"lj\", \"cluster\", \"theta\", \"d\" (or a subset of these values).")
  }
  return(n)
}




calc_t_names <- function (object) {
  if ("a" %in% names(object$posterior)) {
    n <- names(object$posterior$a)
  } else if ("lj" %in% names(object$posterior)) {
    n <- names(object$posterior$lj)
  } else if ("li" %in% names(object$posterior)) {
    n <- names(object$posterior$li)
  } else if ("r" %in% names(object$posterior)) {
    n <- names(object$posterior$r)
  } else {
    stop("The object must contain the values \"a\", \"r\", \"li\", \"lj\" (or a subset of these values).")
  }
  return(n)
}

calc_l_names <- function (object) {
  if ("a" %in% names(object$posterior)) {
    n <- names(object$posterior$a[[1]])
  } else if ("lj" %in% names(object$posterior)) {
    n <- names(object$posterior$lj[[1]])
  } else if ("li" %in% names(object$posterior)) {
    n <- names(object$posterior$li[[1]])
  } else {
    stop("The object must contain the values \"a\", \"li\", \"lj\" (or a subset of these values).")
  }
  return(n)
}

calc_source_names <- function (object) {
  if ("a" %in% names(object$posterior)) {
    n <- colnames(object$posterior$a[[1]][[1]])
  } else if ("lj" %in% names(object$posterior)) {
    n <- colnames(object$posterior$lj[[1]][[1]])
  } else if ("r" %in% names(object$posterior)) {
    n <- dimnames(object$posterior$r[[1]])[[2]]
  } else {
    stop("The object must contain the values \"a\", \"r\", \"lj\" (or a subset of these values).")
  }
  return(n)
}

calc_type_names <- function(object) {
    if ("q" %in% names(object$posterior)) {
      n <- colnames(object$posterior$q)
    } else if ("cluster" %in% names(object$posterior)) {
      n <- colnames(object$posterior$cluster)
    } else if ("theta" %in% names(object$posterior)) {
      n <- colnames(object$posterior$theta)
    } else if ("li" %in% names(object$posterior)) {
      n <- colnames(object$posterior$li[[1]][[1]])
    } else if ("r" %in% names(object$posterior)) {
      n <- dimnames(object$posterior$r[[1]])[[1]]
    } else {
      stop("The object must contain the values \"a\", \"q\", \"r\", \"li\", \"lj\", \"cluster\", \"theta\", \"d\" (or a subset of these values).")
    }
    return(n)
}