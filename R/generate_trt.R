#' Generate binary treatment variable
#' 
#' The proportion of assignment to treatment is p_trt. Generation can be done
#' fixing the proportion of patients in the treatment group to be as
#' close as possible to p_trt and then permuting (type = "exact") or
#' random sampling from a binomial distribution with probability p_trt
#' (type = "random"). For type = "random" one can also specify a
#' formula for logit(P(trt|X)) via "prop". The probability of obtaining 
#' treatment is then binomial with P(trt|X) = 1/(1+exp(-prop)). The variables
#' names specified in "prop" need to be handed over in a data frame "X".
#' 
#'
#' @param n Total sample size
#' @param p_trt Proportion of treatment assignment (default = 0.5)
#' @param type Either "exact" (proportion of p_trt is tried to match exactly), or
#' "random" then treatment assignment is generated from a binomial distribution.
#' @param X matrix/dataframe of predictor variables
#' @param prop A formula used to generate logit(P(trt|X)). Probability of
#' obtaining treatment is then Bin(1, P(trt|X)) with 
#' P(trt|X) = 1/(1+exp(-prop)). If type = "random" and prop is not specified 
#' p_trt is used in the binomial distribution
#' @return Vector of treatment indicators
#' @import stats
#' @export
#' @examples
#' ## by default try to match target trt proportion in observed trt proportion
#' trt <- generate_trt(10) # p_trt = 0.5 is default
#' table(trt)
#' trt <- generate_trt(10, p_trt = 0.2)
#' trt <- generate_trt(11, p_trt = 0.5)
#' table(trt) # remaining patient randomly allocated
#' ## example for random allocation
#' trt <- generate_trt(10, p_trt = 0.5, type = "random") # samples from binomial with p_trt = 0.5
#' table(trt) # observed allocation may deviate from p_trt=0.5 for small samples
#'
#' ## example, where treatment allocation depends on propensity model
#' X <- generate_X_dist(n=10000, p=10, rho=0.5)
#' ## Use propensity score for treatment as 1/(1+exp(-X2))
#' trt <- generate_trt(n=nrow(X), type = "random", X=X, prop = "X2")
#' ## fit propensity model
#' dat <- cbind(trt, X)
#' fit <- glm(trt~., data=dat, family=binomial)
#' summary(fit)

generate_trt <- function(n, p_trt = 0.5, type = c("exact", "random"), X = NULL, prop = NULL){
  type <- match.arg(type)
  rnd <- function(x) ## alternative rounding
    ifelse(x-floor(x) <= 0.5, floor(x), ceiling(x))
  if(type == "exact"){ # try to obtain exact proportion of p_trt
    if(!is.null(X)|!is.null(prop))
      stop("Need to specify type = \"random\" to use X and prop arguments.")
    ntrt <- rnd(p_trt*n)
    nctr <- rnd((1-p_trt)*n)
    trt <- rep(c(0,1), c(nctr, ntrt))
    if(n > ntrt+nctr) # remainders are 0.5,0.5 -> one pat missing
      trt <- c(trt, sample(0:1, 1)) # sample with prob 0.5
    trt <- sample(trt)
  }
  if(type == "random"){ # sample trt from binomial distr
    if(is.null(prop)){
      prop_linear <- log(p_trt/(1-p_trt))
    }
    if(!is.null(prop)){
      stopifnot(n == nrow(X))
      prop_linear <- with(X, eval(parse(text = prop)))
    }
    probs <- 1/(1+exp(-prop_linear))
    trt <- stats::rbinom(n, size = 1, prob = probs)
  }
  trt
}
