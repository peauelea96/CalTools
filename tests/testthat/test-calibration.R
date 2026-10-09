test_that("linear calibration agrees with the analytic formula and sampling", {
  set.seed(17)
  X <- cbind(N = 1, x = rnorm(40), z = runif(40))
  d <- runif(40, 1, 4); q <- runif(40, .5, 2)
  t <- as.numeric(crossprod(X, d)) + c(1, 2, -1)
  v <- d * q
  expected <- d + v * as.numeric(X %*% solve(crossprod(X, X * v), t - crossprod(X, d)))
  z <- calibrate_weights(X, d, t, q)
  expect_equal(z$weights, expected, tolerance = 1e-9)
  expect_equal(as.numeric(X %*% z$lambda), (z$weights - d) / v, tolerance = 1e-9)
  if (requireNamespace("sampling", quietly = TRUE))
    expect_equal(z$g, sampling::calib(X, d, t, q, method = "linear"), tolerance = 1e-8)
})

test_that("signed and zero weights require a positive explicit metric", {
  X <- cbind(N = 1, x = 1:4); d <- c(-1, 0, 2, 3)
  expect_error(calibrate_weights(X, d, c(8, 22)), "metric_scale")
  z <- calibrate_weights(X, d, c(8, 22), method = "truncated",
                         lower = 0, upper = 5, metric_scale = 1)
  expect_true(is.na(z$g[2]))
  expect_equal(as.numeric(crossprod(X, z$weights)), c(8, 22), tolerance = 1e-7)
  expect_true(all(z$weights >= 0 & z$weights <= 5))
  expect_error(calibrate_weights(X, d, c(8, 22), method = "raking"), "positive")
  expect_error(calibrate_weights(X, rep(1, 4), c(8, 22), q = 0), "positive")
})

test_that("permanent clipping counterexample is feasible and solved correctly", {
  X <- rbind(c(-2, 3), c(0, 0), c(3, 0), c(0, 1))
  t <- c(1.9, 7.75); d <- rep(1, 4)
  witness <- c(1.95, 1, 29/15, 1.9)
  expect_equal(as.numeric(crossprod(X, witness)), t)
  unbounded <- calibrate_weights(X, d, t)$weights
  expect_true(all(unbounded[c(1, 3)] > 2))
  # Fixing both to 2 forces the first totals to 2, contradicting 1.9.
  z <- calibrate_weights(X, d, t, method = "truncated", lower = 0, upper = 2)
  ref <- calibrate_weights(X, d, t, method = "truncated", lower = 0, upper = 2, solver = "qp")
  expect_equal(z$weights, ref$weights, tolerance = 2e-6)
  expect_true(is_feasible(X, t, lower = 0, upper = 2))
})

test_that("random dual solutions agree with independent reference QP and KKT", {
  set.seed(123)
  for (k in seq_len(35)) {
    n <- 12; p <- 3
    X <- matrix(rnorm(n * p), n)
    L <- runif(n, -2, 0); U <- runif(n, .1, 3)
    witness <- L + runif(n) * (U - L)
    t <- as.numeric(crossprod(X, witness))
    d <- rnorm(n, 1, 3); v <- runif(n, .4, 3)
    a <- calibrate_weights(X, d, t, method = "truncated", lower = L, upper = U,
      metric_scale = v, solver = "dual", max_iter = 1000, tol = 1e-7)
    b <- calibrate_weights(X, d, t, method = "truncated", lower = L, upper = U,
      metric_scale = v, solver = "qp", tol = 1e-7)
    expect_equal(a$weights, b$weights, tolerance = 1e-5)
    expect_equal(a$objective, b$objective, tolerance = 1e-6)
    expect_equal(a$weights, pmax(L, pmin(U, d + v * as.numeric(X %*% a$lambda))), tolerance = 1e-6)
  }
})

test_that("redundant constraints, zero totals and negative totals are handled", {
  X <- cbind(a = c(1, 0, 1), b = c(0, 1, 1), c = c(1, 1, 2), zero = 0)
  d <- c(1, 2, 3)
  t <- c(a = 0, b = -1, c = -1, zero = 0)
  z <- calibrate_weights(X, d, t)
  expect_equal(as.numeric(crossprod(X, z$weights)), unname(t), tolerance = 1e-8)
  expect_equal(z$rank, 2)
  expect_error(calibrate_weights(X, d, c(a = 0, b = -1, c = 2, zero = 0)), "Inconsistent")
  expect_error(calibrate_weights(X, d, t[c(2, 1, 3, 4)]), "Names and order")
  z <- calibrate_weights(matrix(1, 1, 1), 2, 3, method = "truncated", lower = 3, upper = 3)
  expect_equal(z$weights, 3)
  z <- calibrate_weights(matrix(numeric(), 3, 0), c(1,2,3), numeric(), method = "truncated", lower = 0, upper = 2)
  expect_equal(z$weights, c(1,2,2))
})

test_that("nonlinear calibration uses correct derivatives with varying q", {
  set.seed(80)
  X <- cbind(1, rnorm(30), runif(30)); d <- runif(30, 1, 3); q <- runif(30, .3, 2)
  lambda <- c(.08, -.1, .05)
  wr <- d * exp(q * as.numeric(X %*% lambda)); t <- as.numeric(crossprod(X, wr))
  z <- calibrate_weights(X, d, t, q, "raking")
  expect_equal(z$weights, wr, tolerance = 1e-6)
  expect_equal(as.numeric(z$lambda), lambda, tolerance = 1e-6)
  l <- runif(30, -.2, .3); u <- runif(30, 1.5, 3)
  A <- (u-l)/((1-l)*(u-1))
  gl <- l + (u-l)*plogis(qlogis((1-l)/(u-l)) + A*q*as.numeric(X %*% lambda))
  z <- calibrate_weights(X, d, as.numeric(crossprod(X, d*gl)), q,
                         "logit", lower = l*d, upper = u*d)
  expect_equal(z$weights, d*gl, tolerance = 1e-6)
  expect_equal(as.numeric(z$lambda), lambda, tolerance = 1e-6)
})

test_that("distance functions include endpoint limits and correct raking name", {
  expect_equal(calibration_distance_from_factors(c(1,1), c(2,3), c(1,2), "raking"), 0)
  expect_equal(calibration_distance(c(0,2), c(1,1), method = "raking"), 2*log(2))
  expect_true(is.finite(calibration_distance(c(0,2), c(1,1), method = "logit", lower = 0, upper = 2)))
  expect_equal(calibration_distance(c(-1,2), c(1,1), method = "raking"), Inf)
  expect_equal(calibration_distance(c(1,1), c(1,1), method = "logit", lower = 0, upper = 2), 0)
})

test_that("input and nonconvergence checks never return uncalibrated success", {
  expect_error(calibrate_weights(matrix(1, 3, 1), c(1,1), 3), "d")
  expect_error(calibrate_weights(matrix(NA_real_, 3, 1), rep(1,3), 3), "finite")
  expect_error(calibrate_weights(matrix(1, 3, 1), rep(1,3), 3, max_iter = 0), "max_iter")
  expect_error(calibrate_weights(matrix(1, 3, 1), rep(1,3), 20,
    method = "truncated", lower = 0, upper = 2), "infeasible")
  expect_error(calibrate_weights(matrix(1, 3, 1), rep(1,3), -1,
    method = "raking", max_iter = 5), "did not converge")
})
