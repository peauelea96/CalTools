test_that("negative totals work with positive starting and final weights", {
  X <- cbind(N=rep(1,5),signed=-1:-5,contrast=c(-2,-1,0,1,2))
  d <- rep(2,5); q <- c(.4,.8,1,1.3,2)
  witness <- c(1,2,3,2,1)
  totals <- setNames(as.numeric(crossprod(X,witness)),colnames(X))
  expect_lt(totals["signed"],0);expect_equal(unname(totals["contrast"]),0)
  lower <- c(.1,.2,.3,.4,.5); upper <- c(4,5,6,5,4)
  for (method in c("linear","truncated","raking","logit")) {
    bounds <- if (method %in% c("truncated","logit")) list(lower=lower,upper=upper) else list()
    fit <- do.call(calibrate_weights,c(list(X=X,d=d,totals=totals,q=q,method=method),bounds))
    expect_lte(fit$max_scaled_residual,1e-8)
    # Flipping the units of a margin must leave the solution unchanged.
    flipped <- X; flipped[,2] <- -flipped[,2]
    flipped_totals <- totals; flipped_totals[2] <- -flipped_totals[2]
    reference <- do.call(calibrate_weights,c(list(X=flipped,d=d,totals=flipped_totals,q=q,method=method),bounds))
    expect_equal(fit$weights,reference$weights,tolerance=1e-6)
  }
  for (method in c("auto","projection","lp"))
    expect_true(calibration_feasibility(X,totals,lower,upper,method=method)$feasible)
})

test_that("all-negative targets do not cause false convergence", {
  X <- cbind(a=c(-1,0,-1),b=c(0,-1,-1))
  d <- rep(1,3);totals<-c(a=-4,b=-5)
  for (solver in c("auto","dual","qp")) {
    fit <- calibrate_weights(X,d,totals,method="truncated",lower=c(0,.2,.1),upper=c(5,6,4),solver=solver)
    expect_equal(as.numeric(crossprod(X,fit$weights)),unname(totals),tolerance=1e-7)
    expect_lte(fit$max_scaled_residual,1e-8)
  }
  expect_false(calibration_feasibility(matrix(1,3,1),-1,0,3)$feasible)
  expect_error(calibrate_weights(matrix(1,3,1),d,-1,method="truncated",lower=0,upper=3),"infeasible")
})
