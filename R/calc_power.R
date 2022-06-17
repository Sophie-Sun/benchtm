#' Power calculation to detect overall treatment effect and interaction effect
#'
#' Calculate power for overall treatment effect (not utilizing any
#' covariates) and the power of the interaction test (under the true
#' model, which is Y=f_prog+(b0+b1*f_pred)*trt). The calculation is
#' performed conditional on a given data set of covariates (prognostic
#' effects specified via prog_vals; predictive effects specified via
#' pred_vals) and a treatment indicator. It is assumed that the
#' treatment indicator is independent of the columns in X). The power
#' of the overall test and the interaction test will depend on the
#' actual observed covariates. To make calculations independent of the
#' actual observed covariates (to focus on the super-population) it is
#' suggested to simulate a large number of covariates (e.g. 100 times
#' larger sample size than one would be interested in), based on this
#' one can determine the covariance matrix of the estimates. To obtain
#' a covariance matrix estimate for the sample size one is interested
#' in one can then scale up the covariance matrix estimate (e.g. by
#' 100-fold in the example above).
#'
#' @param b Vector of length 2: Coefficient for treatment indicator
#'   and predictive term in the true model
#' @param scal Scaling parameter for the covariance matrix
#' @param prog_vals Vector of main effects (one value per patient)
#' @param trt Vector of treatment indicators
#' @param pred_vals Vector of predictive effects (one value per
#'   patient)
#' @param alpha Vector of length 2, specifying the desired type 1
#'   error for the overall and the interaction test. The overall test
#'   is performed as a one-sided test (see also argument
#'   sign_better). The interaction test is performed as a two-sided
#'   test.
#' @param type type of data "continuous", "binary", "count" and
#'   "survival" are allowed here, but calculations are only
#'   implemented for "continuous" and "binary".
#' @param sign_better 1 if larger response is better, -1 is smaller
#'   response is better (used in power calculation for the overall
#'   test only)
#' @param sigma_error Residual variance assumed for type = "continuous"
#' @return Vector of two power values
#' @export
#' @examples
#' ###### generate a matrix of covariates
#' scal <- 100
#' X <- generate_X_dist(n=500*scal, p=30, rho=0.5)
#' trt <- generate_trt(n=500*scal, p_trt = 0.5)
#' prog <- "0.5*((X1=='Y')+X11)"
#' pred <- "X11>0.5"
#' pred_vals <- with(X, eval(parse(text = pred)))
#' prog_vals <- with(X, eval(parse(text = prog)))
#' alpha <- c(0.025,0.1)
#' calc_power(b=c(0,0), scal, prog_vals, trt, pred_vals, alpha,
#'            type="continuous", sigma_error=1, sign_better=1)
#' calc_power(b=c(0.2,0.3), scal, prog_vals, trt, pred_vals, alpha,
#'            type="continuous", sigma_error=1, sign_better=1)
#' calc_power(b=c(0,0), scal, prog_vals, trt, pred_vals, alpha,
#'            type="binary", sign_better=1)
#' calc_power(b=c(0.5,0.8), scal, prog_vals, trt, pred_vals, alpha,
#'            type="binary", sign_better=1)
calc_power <- function(b, scal, prog_vals, trt, pred_vals, alpha, type, sign_better, sigma_error){
  ## build design matrix
  X <- cbind(1, prog_vals, trt, trt*pred_vals)
  mu <- as.numeric(X%*%c(0,1,b)) # mean on linear predictor scale
  n1 <- sum(trt==1)
  n0 <- sum(trt==0)
  power <- numeric(2)
  if(type == "continuous"){ # calculations specific for cont
    delta_ov <- (b[1]+b[2]*mean(pred_vals))
    s0_2 <- var(prog_vals)+sigma_error^2 # var control
    s1_2 <- var(prog_vals+b[1]+b[2]*pred_vals)+sigma_error^2 # var trt
    Vc <- sigma_error^2*solve(crossprod(X)/scal) # variance estimate for true model (see generate_y + intercept)
    se_ov <- sqrt(s1_2/n1+s0_2/n0)*sqrt(scal)
  }
  if(type == "binary"){ # specific for binary
    p1 <- 1/(1+exp(-(prog_vals+b[1]+b[2]*pred_vals)))
    Ep1 <- mean(p1)
    p0 <- 1/(1+exp(-prog_vals))
    Ep0 <- mean(p0)
    delta_ov <- (Ep1 - Ep0)
    p <- 1/(1+exp(-mu))
    Vc <- solve(crossprod(X*sqrt(p*(1-p)))/scal) #solve(t(X)%*%diag(p*(1-p))%*%X) # variance estimate for true model
    v1 <- Ep1*(1-Ep1)/n1
    v0 <- Ep0*(1-Ep0)/n0
    se_ov <- sqrt(v1 + v0)*sqrt(scal)
  }
  ## power of gauss-test (assuming variances known)
  delta_ov <- delta_ov*sign_better
  power[1] <- 1-pnorm(qnorm(1-alpha[1]), delta_ov/se_ov, 1) 
  ## power to detect interaction term (under true model)
  se <- sqrt(Vc[4,4])
  power[2] <- 1-(pnorm(qnorm(1-alpha[2]/2), b[2]/se, 1)-pnorm(qnorm(alpha[2]/2), b[2]/se, 1))
  power
}
