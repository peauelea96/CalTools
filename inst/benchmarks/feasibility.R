# Timings include input checks and factorization on identical feasible problems.
# An iteration limit is inconclusive. This is not an infeasible-case benchmark.
library(CalTools)
set.seed(20261009)
results <- list()
for (n in c(100L, 500L, 2000L)) {
  p <- 8L
  X <- cbind(1, matrix(rnorm(n * (p - 1L)), n))
  lower <- runif(n, 0, .3)
  upper <- runif(n, 1.7, 2)
  witness <- pmax(lower, pmin(upper, 1 + as.numeric(X %*% seq(-.7, .7, length.out=p))))
  totals <- as.numeric(crossprod(X, witness))
  for (method in c("projection", "auto", "lp")) {
    elapsed <- numeric(3)
    for (replicate in seq_along(elapsed)) {
      elapsed[replicate] <- system.time(result <- calibration_feasibility(
        X, totals, lower, upper, method=method, max_iter=5000,
        fallback=FALSE, tol=1e-8))[["elapsed"]]
      if (isTRUE(result$feasible)) stopifnot(result$max_scaled_residual <= 1e-8,
        all(result$weights >= lower-1e-7), all(result$weights <= upper+1e-7))
    }
    results[[length(results)+1L]] <- data.frame(n=n, p=p, method=method,
      median_seconds=median(elapsed), feasible=result$feasible, status=result$status,
      solver=result$solver, max_scaled_residual=result$max_scaled_residual)
  }
}
result <- do.call(rbind,results)
print(result,row.names=FALSE)
print(sessionInfo())
