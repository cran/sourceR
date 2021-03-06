#' Human cases of campylobacteriosis and numbers of source samples positive for \emph{Campylobacter}.
#'
#' A dataset containing the number of human cases of campylobacteriosis and numbers of source samples
#' positive for \emph{Campylobacter} for each bacterial subtype.
#'
#' @format A list containing the human cases (`cases'), source samples (`sources'0 and, prevalences (`prev').
#'
#' \strong{cases:} data frame with 91 rows and 2 variables:
#' \describe{
#'   \item{Human}{number of human cases of campylobacteriosis between 2005-2008 in the Manawatu
#'   region of New Zealand}
#'   \item{Type}{MLST type id for the samples}
#' }
#'
#' \strong{sources:} data frame with 690 rows and 3 variables
#' \describe{
#'   \item{Count}{number of source samples positive for campylobacteriosis}
#'   \item{Source}{Source id for the samples}
#'   \item{Type}{MLST type id for the samples}
#' }
#'
#' \strong{prev:} data frame with 6 rows and 4 variables
#' \describe{
#'   \item{Value}{Prevalence value (number of positive samples divided by total number of samples)}
#'   \item{Source}{Source id for the samples}
#'   \item{n_positive}{number of positive source samples for campylobacter (PCR)}
#'   \item{n_total}{total number of source samples tested for campylobacter}
#' }
#'
"campy"
