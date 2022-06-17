## check dimension of the generated X
test_that("check dimension of X generation", {
  expect_equal(dim(generate_X_syn(100)), c(100, 30))
  expect_equal(dim(generate_X_dist(100,20)), c(100, 20))
})


## check dimension of the generated y
test_that("check dimension of y generation", {
  X <- generate_X_dist(n=100, p=10, rho=0.5)
  trt <- generate_trt(n=100, p_trt = 2/3)
  dat <- generate_y(X, trt, prog = "0.5*((X1=='Y')+X3)", pred = "X3>0", b0 = 0,
                    sd_te = 0.5,  type = "survival", sigma_error = 3,
                    include_truth = FALSE,sign_better = -1,cens_time = 2)
  expect_equal(dim(dat), c(100, 13))
})


## check power estimation
test_that("check estimation of power (no prog and pred effects)", {
  scal <- 100
  ## scenarios with no prognostic and no predictive effect
  ## (then we can compare against standard functions)
  X <- generate_X_dist(n=500*scal, p=10, rho=0.5)
  Z <- generate_trt(n=500*scal)
  prog <- "0.01*((X1=='Y')+X3)" ## essentially no prognostic effect
  pred <- "X3 > 0"
  pred_vals <- with(X, eval(parse(text = pred)))
  prog_vals <- with(X, eval(parse(text = prog)))
  ## scenario with no predictive effect (b2=0)
  calc1 <- calc_power(c(0.3,0), scal = scal, prog_vals, Z, pred_vals,
                     alpha=c(0.025,0.1), type="continuous", sigma_error=1,
                     sign_better=1)[1]
  exp1 <- power.t.test(delta=0.3, sd=1, sig.level=0.025,
                       alternative = "one.sided", n=250)$power
  expect_equal(calc1, exp1, tolerance = 0.01)

  ## binary case
  X <- generate_X_syn(n = 500*scal)
  X <- generate_X_dist(n=500*scal, p=30, rho=0.5)
  trt <- generate_trt(n=500*scal)
  prog <- "0.001*(X14-(X5=='b'))"
  pred <- "X14"
  pred_vals <- with(X, eval(parse(text = pred)))
  prog_vals <- with(X, eval(parse(text = prog)))
  ## scenario with no predictive effect (b2=0)
  calc2 <- calc_power(c(0.3,0), scal = scal, prog_vals, Z, pred_vals,
                     alpha=c(0.025,0.1), type="binary",
                     sign_better=1)[1]
  h <- pwr::ES.h(p1 = 1/(1+exp(-0.3)), p2 = 1/(1+exp(-0)))
  exp2 <- pwr::pwr.2p.test(h, n=250, sig.level = 0.025, alternative="greater")$power
  expect_equal(calc2, exp2, tolerance = 0.01)
})


test_that("check estimation of power (precalculated power values)", {
  s <- 100
  X <- generate_X_dist(n=100*s, p=15, rho=0.5)
  trt <- generate_trt(100*s)
  prog <- "0.5*(X1=='Y')+X11"
  pred <- "X11>0.5"
  pred_vals <- with(X, eval(parse(text = pred)))
  prog_vals <- with(X, eval(parse(text = prog)))
  calc1 <- calc_power(c(0.6,0.6), scal = s, prog_vals, trt, pred_vals,
                      alpha=c(0.025,0.1), type="continuous", sigma_error=1, sign_better=1)
  exp1 <- c(0.7362, 0.4953) ## pre-calculated based on simulation
  ## binary case
  calc2 <- calc_power(c(1.7,-1.6), scal = s, prog_vals, trt, pred_vals, alpha=c(0.025,0.1),
                      type="binary", sign_better=1)
  exp2 <- c(0.7487, 0.5131) ## pre-calculated based on simulation
  expect_equal(calc1, exp1, tolerance = 0.02)
  expect_equal(calc2, exp2, tolerance = 0.02)
})


## check calculation of b
test_that("check get_b", {
  s <- 100
  X <- generate_X_dist(n=100*s, p=15, rho=0.5)
  trt <- generate_trt(100*s)
  prog <- "0.5*(X1=='Y')+X11"
  pred <- "X11>0.5"
  ## continuous case
  calc1 <- get_b(X, scal=s, prog, pred, trt, type = "continuous", power = c(0.7362, 0.4953), sigma_error=1)
  exp1 <- c(0.6,0.6)
  ## binary case
  calc2 <- get_b(X, scal=s, prog, pred, trt, type = "binary", power = c(0.7487, 0.5131), start=c(0.5,-0.5))
  exp2 <- c(1.7,-1.6)
  expect_equal(as.numeric(calc1), exp1, tolerance = 0.05)
  expect_equal(as.numeric(calc2), exp2, tolerance = 0.05)
})


test_that("check get_b0", {
  s <- 100
  X <- generate_X_dist(n=50*s, p=30, rho=0.5)
  trt <- generate_trt(50*s)
  prog <- "X11"
  pred <- "X14"
  calc1 <- get_b0(X, scal = s, b1=-0.5, prog, pred, trt, type="continuous",
                  power = 0.5301, alpha = 0.05, sigma_error=5, interval=c(-5,5))
  exp1 <- 2.5
  X <- generate_X_dist(n=500*s, p=30, rho=0.5)
  trt <- generate_trt(500*s)
  calc2 <- get_b0(X, scal = s, b1=0.15, prog, pred, trt, type = "binary",
                  power = 0.4384, alpha = 0.1, interval = c(-5,5))
  exp2 <- 0.25
  expect_equal(as.numeric(calc1), exp1, tolerance = 0.05)
  expect_equal(as.numeric(calc2), exp2, tolerance = 0.05)
})
