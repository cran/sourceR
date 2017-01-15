#' True values for the parameters generating the simulated data.
#'
#' A list containing the true values of the parameters used to simulate the sim_SA_data dataset.
#'
#' @format A list with 5 items:
#' \describe{
#'   \item{alpha}{A dataframe with 24 rows and 4 variables: Value contains the true alpha values,
#'   Time, Location and Source contain the time, location and source id's respectively.}
#'   \item{q}{A dataframe with 91 rows and 2 variables: Value contains the true q values, and
#'   Type contains the type id's.}
#'   \item{lambda_i}{A dataframe with 364 rows and 4 variables: Value contains the true lambda_i values,
#'   Time, Location and Type contain the time, location and type id's respectively.}
#'   \item{lambda_j}{A dataframe with 24 rows and 4 variables: Value contains the true lambda_j values,
#'   Time, Location and Source contain the time, location and source id's respectively.}
#'   \item{r}{A dataframe with 2184 rows and 5 variables: Value contains the true r values,
#'   Time, Type, Location and Source contain the time, type, location and source id's respectively.}
#' }
"sim_SA_true"
