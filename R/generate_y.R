#' Generate outcome variable y and merge with predictor data matrix X.
#'
#' Generate simulated data according to the model f(X) = f_prog(X) +
#' trt*(b0 + b1*f_pred(X)) for different outcome distributions.  For
#' continuous data f(X) is the conditional mean, for binary data the
#' logit response probility, for count data the log-rate and for
#' survival data the log-hazard rate.
#'
#' @param X matrix/dataframe of predictor variables
#' @param trt Binary treatment indicator variable
#' @param prog Character variable giving expression for prognostic
#'   effects (defined in terms of names in the X matrix)
#' @param pred Character variable giving expression for predictive
#'   effects (defined in terms of names in the X matrix)
#' @param b0 Treatment (main effect)
#' @param b1 Coefficient of the predictive effects defined in pred
#' @param sd_te Standard deviation of the treatment effects defined
#'   via pred. If given b1 is ignored. For binary data this is assumed
#'   to be on the log-odds ratio scale. For survival and count data
#'   this is assumed to be on the log scale.
#' @param type Outcome data type to generate ("continuous", "binary",
#'   "count" and "survival" are allowed here)
#' @param sigma_error Residual error, only needed for type =
#'   "continuous"
#' @param theta Overdispersion parameter, only needed for data_type =
#'   "count" (variance of neg bin in this parameterization is mu + mu^2/theta)
#' @param cens_time Maximum follow-up time, only needed for data_type
#'   = "survival"
#' @param include_truth boolean, will the true prognostic and predictive
#'   effects be included in the outcome data-set?
#' @param sign_better whether larger response is better (used to determine
#'   whether b1 is negative or positive if not given)
#' @return Data set
#' @import stats
#' @export
#' @examples
#' X <- generate_X_dist(n=10000, p=10, rho=0.5)
#' ## observational data set
#' trt <- generate_trt(n=nrow(X), type = "random", X=X, prop= "X2")
#' dat <- generate_y(X, trt, prog = "0.5*((X1=='Y')+X3)",
#'                   pred = "X3>0", b0 = 0, b1 = 1,
#'                   type = "continuous", sigma_error = 3)
#' 
#' 
#' #### generate data from user specified covariate X
#' X <- sapply(1:10, function(ii) {rnorm(500)})
#' X <- as.data.frame(X)
#' colnames(X) <- paste0("Z", 1:10)
#' trt <- generate_trt(nrow(X), p_trt = 0.5)
#' dat <- generate_y(X, trt, prog = "0.5*Z1",
#'                   pred = "(Z5>0)", b0 = 0, b1 = 1,
#'                   type = "binary")
#' glm(Y~trt*., data=dat, family=binomial)

generate_y <- function(X, trt, prog, pred, b0, b1 = NULL, sd_te = NULL,
                       type = c("binary","continuous", "count", "survival"),
                       sigma_error = 1.0,  theta = 1.0, cens_time = 2.0,
                       include_truth = FALSE,
                       sign_better = 1){
  stopifnot(nrow(X) == length(trt), sort(unique(trt)) == c(0,1), is.element(type, c("binary","continuous", "count", "survival")))
  n <- nrow(X)
  ## data generation is based on framework of f(Y) = prog + Trt(b0+b1*pred) + N(0,sigma_error)

  ## extract predictive variable values
  pred_val <- with(X, eval(parse(text = pred)))
  prog_val <- with(X, eval(parse(text = prog)))

  ## if sd_te is given, value of b1 is ignored
  if(!is.null(sd_te)){
    sd_pred <- sd(pred_val)
    b1 <- as.numeric(sd_te/sd_pred)*sign_better
  }
  dat <- cbind(trt,X)
  y_link <- prog_val + trt*(b0 + b1*pred_val)
  y_0 <- prog_val
  y_1 <- prog_val + b0 + b1*pred_val
  y_effect <- y_1 - y_0

  if(type == "continuous"){
    dat$Y <- y_link + stats::rnorm(n, 0, sigma_error)
    form_tmp <- paste0("mu = ",prog," + Trt*","(",round(b0,4),"+",
                       round(b1,4),"*(",pred,")) + N(0, ",sigma_error,")")
  }
  if(type == "binary"){
    dat$Y <- stats::rbinom(n, size = 1, prob = 1.0/(1.0 + exp(-y_link)))
    form_tmp <- paste0("logit(p) = ", prog, " + Trt*", "(", round(b0, 4), "+",
                       round(b1, 4), "*(",pred,"))")
  }
  if(type == "count"){ ## use negative binomial distribution
    dat$duration <- stats::rlnorm(n, 0, 0.3) ## observation time per patient (need to adjust for using offset in model)
    re <- stats::rgamma(n, theta, theta) ## gamma random effect per patient (-> will induce over-dispersion in count data)
    lambda <- dat$duration*re*exp(y_link)
    dat$Y <- stats::rpois(n, lambda)
    form_tmp <- paste0("log(lambda) = ", prog, " + Trt*", "(",
                       round(b0, 4), "+", round(b1, 4), "*(", pred, "))")
  }
  if(type == "survival"){
    lambda <- exp(y_link)
    dat$Y <- stats::rexp(n, lambda)
    dat$Y <- pmin(dat$Y, cens_time) ## censor for maximum follow-up time
    dat$event <- as.numeric(dat$Y < cens_time)
    form_tmp <- paste0("log(lambda) = ", prog, " + Trt*", "(", round(b0, 4),
                       "+", round(b1, 4), "*(", pred, "))")
  }

  if(include_truth){
    dat$trt_effect <- y_effect
    if (type == "binary") {
      dat$prob_diff <- 1/(1 + exp(-y_1)) - 1/(1 + exp(-y_0))
    }
    else {
      dat$prob_diff <- NA
    }
  }
  attr(dat, "form") <- form_tmp
  attr(dat, "b0") <- b0
  attr(dat, "b1") <- b1
  attr(dat, "pred") <- pred
  attr(dat, "prog") <- prog
  dat
}
