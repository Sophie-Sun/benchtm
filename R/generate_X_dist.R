#' Generate X from pre-specified distribution.
#'
#' Generate X from pre-specified distribution. X1 is a
#' Bernoulli distributed rv with probability 0.5. X2 an exponential
#' distribution with lambda=1. X3,X4,Z are standard normal variates
#' with corelation rho and X5=(Z>1), X6 is from a multinomial
#' distribution with 4 equally probable categories. X7-Xp are
#' independent standard normal distributions.
#'
#' @param n Number of observations
#' @param p Number of predictors requested
#' @param rho Correlation between variable 3-5
#' @return Data frame of predictors
#' @importFrom MASS mvrnorm
#' @import stats
#' @export
#' @examples
#' X <- generate_X_dist(n=100, p=10, rho=0.5)
generate_X_dist <- function(n, p, rho = 0.5){
  ############################################
  # X1~ ber(0,0.5)
  # X2~ exp(1)
  # X3, X4, Z ~ N((0,0,0),rho), X5 = (Z>1)
  # X6 ~ multinom(4, (0.25,0.25,0.25,0.25)), values = c("M1","M2","M3","M4")
  # X7,...~N(0,1)
  ###########################################
  # binary
  x_bin <- stats::rbinom(n, size = 1, prob = 0.50)
  x_bin <- ifelse(x_bin, "Y", "N")
  ## exponential
  x_exp <- stats::rexp(n, rate = 1)
  # normal
  Sigma <- matrix(rho, nrow = 3, ncol = 3) + diag(3) * (1-rho)
  x_norm_cor <- as.data.frame(MASS::mvrnorm(n, mu = c(0, 0, 0), Sigma = Sigma))
  x_norm_cor[,3] <- (x_norm_cor[, 3] > 0)
  x_norm_cor[,3] <- factor(ifelse(x_norm_cor[, 3], "Y", "N"))
  ##### x multinom
  x_multinom <- paste("M", sample(1:4, n, replace = TRUE), sep = "")
  x_norm_ind <- sapply(1:(p-6), function(ii){
    rnorm(n, 0, 1)
  })
  X <- as.data.frame(cbind(x_bin, x_exp, x_norm_cor, x_multinom, x_norm_ind))
  colnames(X) <- paste0("X", 1:(dim(X)[2]))
  return(X)
}

