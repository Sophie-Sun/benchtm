#' Calculate b parameter to achieve a given power to detect an overall
#' effect and the interaction effect
#'
#' Calculate the model parameters b0 and b1 that provide a given
#' target power for testing for an overall treatment effect and an
#' interaction effect (under the true model). The calculation is
#' performed conditional on a given data set of covariates and
#' treatment indicator. It is assumed that the treatment indicator is
#' independent of the columns in X). The power of the overall test and
#' the interaction test will depend on the actual observed
#' covariates. To make calculations independent of the actual observed
#' covariates (to focus on the super-population) it is suggested to
#' simulate a large number of covariates (e.g. 100 times larger sample
#' size than one would be interested in), based on this one can
#' calculate an estimate of the covariance matrix of the estimates. To
#' obtain a covariance matrix estimate for the sample size one is
#' interested in one can then scale up the covariance matrix estimate
#' (e.g. by 100-fold in the example above). Note that there may not be
#' any solution or a unique solution, in case of non-convergence
#' (which particularly happens for the binary endpoint) it sometimes
#' helps to modify start and optim_method arguments.
#'
#' @param X Matrix with all covariates
#' @param scal Scaling parameter for the covariance matrix
#' @param prog Character variable giving expression for prognostic
#'   effects (defined in terms of names in the X matrix)
#' @param pred Character variable giving expression for predictive
#'   effects (defined in terms of names in the X matrix)
#' @param trt Treatment effect indicator (same length as number of rows in X).
#' @param type type of data "continuous", "binary", "count" and
#'   "survival" are allowed here, but calculations are only
#'   implemented for "continuous" and "binary".
#' @param power Vector of length 2, specifying the target power for the
#'   overall and the interaction test
#' @param alpha Vector of length 2, specifying the desired type 1
#'   error for the overall and the interaction test. The overall test
#'   is performed as a one-sided test (see also argument
#'   sign_better). The interaction test is performed as a two-sided
#'   test.
#' @param start Vector of length 2: Starting values for the parameter b.
#' @param sign_better 1 if larger response is better, -1 is smaller
#'   response is better (used in power calculation for the overall
#'   test only)
#' @param sigma_error Residual variance assumed for type = "continuous"
#' @param optim_method method argument in the optimization see ?optim
#' @return Vector of two model parameters b0 and b1
#' @export
#' @examples
#' scal <- 100
#' X <- generate_X_dist(n=500*scal, p=30, rho=0.5)
#' trt <- generate_trt(n=500*scal, p_trt = 0.5)
#' prog <- "0.5*((X1=='Y')+X11)"
#' pred <- "X11>0.5"
#' get_b(X, scal, prog, pred, trt, type="continuous",
#'       power = c(0.9, 0.9), alpha = c(0.025, 0.1),
#'       start=c(0,0), sigma_error=1)
#' get_b(X, scal, prog, pred, trt, type = "binary",
#'       power = c(0.9, 0.9), alpha = c(0.025, 0.1),
#'       start=c(0,0))
get_b <- function(X, scal, prog, pred, trt, type,
                  power = c(0.9, 0.8), alpha = c(0.025, 0.1), start=c(0,0),
                  sign_better = 1, sigma_error, optim_method = "Nelder-Mead"){

  stopifnot(nrow(X) == length(trt), is.character(prog), is.character(pred))
  stopifnot(type %in% c("continuous", "binary", "count", "survival"))
  
  if(type %in% c("count", "survival"))
    stop(sprintf("%s not implemented yet", type))
  if(any(power < alpha))
    stop("All elements in power need to be > alpha")
  
  ## generate prognostic and predictive covariates
  pred_vals <- with(X, eval(parse(text = pred)))
  prog_vals <- with(X, eval(parse(text = prog)))

  opt <- optim(function(b, ...){
    sq <- (calc_power(b, ...) - power)^2
    sqrt(sum(sq))
  }, par=start, scal = scal, prog_vals=prog_vals, trt=trt, pred_vals=pred_vals, alpha = alpha, 
  type = type, sign_better = sign_better, sigma_error = sigma_error, method = optim_method)
  out <- opt$par
  outpow <- calc_power(out, scal=scal, prog_vals=prog_vals, trt=trt, pred_vals=pred_vals,
                      alpha = alpha, type = type,
                      sign_better = sign_better, sigma_error = sigma_error)
  attr(out, "power_results") <- outpow
  if(any(abs(outpow-power) > 0.01))
    message("optimization failed, try different starting value or optim_method")
  out
}

