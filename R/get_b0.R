#' Calculate b0 parameter to achieve a given power to detect an overall
#' effect and the interaction effect
#'
#' Calculate the model parameters b0 (for a given b1) that provides a
#' given target power for testing for an overall treatment effect.
#' The calculation is performed conditional on a given data set of
#' covariates and treatment indicator. The power of the overall test
#' depends on the actual observed covariates. To make calculations
#' independent of the actual observed covariates (to focus on the
#' super-population) it is suggested to simulate a large number of
#' covariates (e.g. 100 times larger sample size than one would be
#' interested in), based on this one can calculate determine the
#' covariance matrix of the estimates. To obtain a covariance matrix
#' estimate for the sample size one is interested in one can then
#' scale up the covariance matrix estimate (e.g. by 100-fold in the
#' example above).
#'
#' @param X Matrix with all covariates
#' @param scal Scaling parameter for the covariance matrix
#' @param prog Character variable giving expression for prognostic
#'   effects (defined in terms of names in the X matrix)
#' @param pred Character variable giving expression for predictive
#'   effects (defined in terms of names in the X matrix)
#' @param trt Treatment effect indicator (same length as number of rows in X).
#' @param b1 The parameter b1 (calculation of b0 is for a given b1)
#' @param type type of data "continuous", "binary", "count" and
#'   "survival" are allowed here, but calculations are only
#'   implemented for "continuous" and "binary" currently.
#' @param power Vector of length 2, specifying the target power for the
#'   overall and the interaction test
#' @param alpha Vector of length 2, specifying the desired type 1
#'   error for the overall and the interaction test. The overall test
#'   is performed as a one-sided test (see also argument
#'   sign_better). The interaction test is performed as a two-sided
#'   test.
#' @param sign_better 1 if larger response is better, -1 is smaller
#'   response is better (used in power calculation for the overall
#'   test only)
#' @param sigma_error Residual variance assumed for type = "continuous"
#' @param theta overdispersion paramter for type = count". It is
#'   consistent with MASS::glm.nb and the (dprq)nbinom functions in R
#'   with the "size" parameter equal to theta. In this parameterization
#'   the negative binomial distribution converges to the Poisson
#'   distribution for theta goes to infinity.
#' @param lambda0 Intercept of exponential regression (on non-log scale)
#' @param cens_time Function to generate the censoring time, only
#'   needed for data_type = "survival"
#' @param t_mile Time point for comparing survival probabilities for overall
#'   test for treatment effect
#' @param interval Interval to search for b0 handed over to uniroot
#' @return Vector of model parameter b0
#' @export
#' @examples
#'
#' scal <- 100
#' X <- generate_X_dist(n = 500 * scal, p = 30, rho = 0.5)
#' trt <- generate_trt(n = 500 * scal, p_trt = 0.5)
#' prog <- "0.5*((X1=='Y')+X11)"
#' pred <- "X11>0.5"
#' get_b0(X, scal, prog, pred, trt,
#'   b1 = 0.5, type = "continuous",
#'   power = 0.9, alpha = 0.025, interval = c(-2, 2), sigma_error = 1
#' )
#' get_b0(X, scal, prog, pred, trt,
#'   b1 = 0.5, type = "binary",
#'   power = 0.9, alpha = 0.025, interval = c(-2, 2)
#' )
get_b0 <- function(X, scal, prog, pred, trt, b1, type,
                   power = 0.9, alpha = 0.025,
                   sign_better = 1, sigma_error, theta,
                   lambda0, cens_time, t_mile, interval = c(-5, 5)) {
  stopifnot(nrow(X) == length(trt), is.character(prog), is.character(pred))
  stopifnot(type %in% c("continuous", "binary", "count", "survival"))
  if (power < alpha) {
    stop("Power needs to be > alpha")
  }

  ## generate prognostic and predictive covariates
  pred_vals <- with(X, eval(parse(text = pred)))
  prog_vals <- with(X, eval(parse(text = prog)))

  opt <- uniroot(
    function(b0, ...) {
      calc_power(c(b0, b1), ...)[1] - power
    },
    scal = scal, prog_vals = prog_vals, trt = trt, pred_vals = pred_vals, alpha = alpha,
    type = type, sign_better = sign_better, sigma_error = sigma_error, theta = theta,
    lambda0 = lambda0, cens_time = cens_time, t_mile = t_mile, interval = interval
  )
  out <- opt$root
  outpow <- calc_power(c(out, b1),
    scal = scal, prog_vals = prog_vals, trt = trt, pred_vals = pred_vals,
    alpha = alpha, type = type,
    sign_better = sign_better, sigma_error = sigma_error, theta = theta,
    lambda0 = lambda0, cens_time = cens_time, t_mile = t_mile
  )[1]
  attr(out, "power_results") <- outpow
  if (abs(outpow - power) > 0.01) {
    message("optimization failed, try different starting value or optim_method")
  }
  out
}
