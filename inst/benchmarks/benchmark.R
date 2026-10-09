# Run after installing CalTools: Rscript inst/benchmarks/benchmark.R
# Synthetic dense examples; elapsed times are machine/BLAS dependent.
library(CalTools)
set.seed(20261008)
results <- list()
for (n in c(200L, 1000L, 10000L)) {
  p <- 12L
  X <- cbind(1, matrix(rnorm(n * (p - 1L)), n))
  d <- runif(n, 1, 3)
  L <- rep(.5, n); U <- rep(4, n)
  # Construct a bounded optimum from known KKT multipliers.
  lambda <- c(.2, seq(-.5, .5, length.out = p-1L))
  wstar <- pmax(L, pmin(U, d + d * as.numeric(X %*% lambda)))
  target <- as.numeric(crossprod(X, wstar))
  solvers <- if (n <= 1000L) c("dual", "qp") else "dual"
  for (solver in solvers) {
    elapsed <- numeric(3)
    for (i in seq_along(elapsed)) {
      tm <- system.time(fit <- calibrate_weights(X, d, target, method="truncated",
        lower=L,upper=U,solver=solver,tol=1e-7,max_iter=1000))
      elapsed[i] <- tm[["elapsed"]]
      stopifnot(fit$max_scaled_residual <= 1e-7, max(abs(fit$weights - wstar)) < 1e-4)
    }
    results[[length(results)+1L]] <- data.frame(n=n,p=p,solver=solver,
      median_seconds=median(elapsed),max_scaled_residual=fit$max_scaled_residual,
      max_error_known_optimum=max(abs(fit$weights-wstar)),iterations=fit$iterations)
  }
}
result <- do.call(rbind,results)
print(result,row.names=FALSE)
cat("\nSession:\n")
print(sessionInfo())
# To persist benchmark values:
# write.csv(result, "benchmark-results.csv", row.names=FALSE)
