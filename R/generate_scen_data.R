#' Generate simulation data from pre-specified scenarios
#'
#' Generate simulation data from simulation scenarios used in XYZ. The
#' scenarios were derived so that the prognostic part of the model on
#' control as an R^2 of 0.32 for continuous data and an AUC of 0.66
#' for binary data. The coefficient b0 was calculated in each scenario so
#' that the overall test for a treatment effect has a power of 50% for
#' 2.5% one-sided. There are four main scenarios (defined by the
#' expressions for prognostic and predictive part) replicated across
#' continuous and binary endpoints. Within each scenario there are 5
#' sub-scenarios corresponding to different selections of b1. The
#' third scenario in each of the sub-scenarios correspond to the scenario
#' where the interaction test (under the true model) has 80% power
#' (for 20% two-sided, or 10% one-sided as only b1 > 0 is
#' considered). The other 4 scenarios correspond to 0, 0.5, 1.5, 2 times
#' the b1 value that provides 80% power of the interaction test. In
#' total there are hence 4x2x5=40 scenarios.
#'
#' @param scen A row from the scen_param data set or scen_param_TTE data set for time-to-event data
#' @param include_truth Boolean, will the true treatment effect be included in the outcome data-set?
#' @param type For type == "sample" (default) X is generated using R package synthpop (using function
#'  generate_X_syn). For type == "resample" data are resampled from a large saved data-set generated
#'   from generate_X_syn (this option is considerably faster).
#' @return A data frame
#' @export
#' @examples
#' data(scen_param) ## scenarios used in XYZ
#' dat <- generate_scen_data(scen = scen_param[1, ])
generate_scen_data <- function(scen, include_truth = TRUE, type = c("sample", "resample")) {
  type <- match.arg(type)
  if(type == "resample"){
    ind <- sample(1:50000, 500)
    X <- X_large_pop[ind, ]
  }else{
    X <- generate_X_syn(n = 500)
  }
  trt <- generate_trt(n = 500, p_trt = 0.5)
  ## for survival cases
  lambda0 <- 0.0002
  cens_time <- function(n) {
    p <- sample(0:1, n, replace = TRUE, prob = c(0.05, 0.95))
    r1 <- runif(n, 0, 1000)
    r2 <- 1000 + (2000 - 1000) * rbeta(n, 1, 1.5)
    (1 - p) * r1 + p * r2
  }
  generate_y(
    X = X, trt = trt,
    prog = scen$prog, pred = scen$pred,
    b0 = scen$b0, b1 = scen$b1,
    type = scen$type,
    sigma_error = 1,
    cens_time = cens_time,
    lambda0 = lambda0,
    include_truth = include_truth
  )
}
