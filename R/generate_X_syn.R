#' Generate X from data "syn"
#'
#' Generate data from synthetic clinical trial data with dimension n
#' and 30 variables.
#'
#' The data are generated based on the data-set data_syn which is not
#' exported from the package namespace. This data-set was generated
#' from pooled set of clinical trial data using the synthpop function
#' from the synthpop R package. data_syn has 1934 rows and 30
#' variables. Variable names have been removed. Continuous variables
#' have been scaled to the interval [0,1], levels of the categorical
#' variable have been blinded. For information on the synthetic data
#' which synthpop generates from, see the vignette: Description of
#' included synthetic data.
#' 
#' @param n Number of observations
#' @return A data frame of predictors
#' @import synthpop
#' @export
#' @examples
#' X <- generate_X_syn(n=100)
#'
generate_X_syn <- function(n){
  X <- synthpop::syn(data_syn, k = n, print.flag = FALSE)$syn
  return(X)
}
