test_that("projection counterexample is feasible but never labelled optimal", {
  X <- matrix(c(3,0,1,1,1,0,1,1,0,3),5,2)
  expect_warning(a <- calibrate(X,rep(1,5),c(.8,.4),method="truncated",lower=0,upper=1.2,
                                solver="projection"),"not certified")
  expect_equal(a$weights,c(0,1/3,1/15,11/15,0),tolerance=1e-10)
  expect_identical(a$status,"feasible_uncertified")
  expect_equal(a$objective,127/75,tolerance=1e-10)
  z <- calibrate(X,rep(1,5),c(.8,.4),method="truncated",lower=0,upper=1.2,
                 solver="projection",projection_action="refine")
  expect_equal(z$weights,c(1/35,38/105,4/105,71/105,0),tolerance=1e-7)
  expect_equal(z$objective,71/42,tolerance=1e-9)
  expect_gte(a$certificate$duality_gap+1e-10,1/350)
  expect_gte(a$certificate$kkt_objective_bound+1e-10,1/350)
  expect_error(calibrate(X,rep(1,5),c(.8,.4),method="truncated",lower=0,upper=1.2,
    solver="projection",projection_action="error"),"not certified")
})

test_that("wrong faces never establish infeasibility", {
  X <- t(rbind(c(1,-1,0),c(-2,0,3)))
  expect_error(calibrate(X,rep(1,3),c(1,1),method="truncated",lower=0,upper=1.2,
    solver="projection"),"does not prove infeasibility")
  a <- calibrate(X,rep(1,3),c(1,1),method="truncated",lower=0,upper=1.2,
    solver="projection",projection_action="refine")
  expect_equal(a$weights,c(1.2,.2,17/15),tolerance=1e-7)
  X <- t(rbind(c(2,0,-1,1,0),c(1,-3,0,0,0),c(1,0,0,1,0),c(2,0,0,0,3)))
  tt <- c(7/4,1,7/4,2)
  a <- calibration_feasibility(X,tt,0,1,method="fixing",fallback=FALSE)
  expect_true(is.na(a$feasible))
  b <- calibration_feasibility(X,tt,0,1,method="pfc",max_iter=10000,fallback=FALSE,tol=1e-7)
  expect_true(b$feasible)
  expect_lte(b$max_scaled_residual,1e-7)
  expect_true(calibration_feasibility(X,tt,0,1,method="fixing")$feasible)
  c <- calibration_feasibility(X,tt,0,1,method="pfc",max_iter=1,fallback=FALSE)
  expect_true(is.na(c$feasible))
})

test_that("PFC produces an independently checked separation certificate", {
  X <- rbind(c(1,0),c(0,1),c(1,1)); tt <- c(2,0)
  f <- calibration_feasibility(X,tt,0,1,method="pfc",fallback=FALSE)
  expect_false(f$feasible); expect_identical(f$solver,"pfc")
  v <- f$certificate$direction
  score <- sum(v*tt)-sum(pmax(as.numeric(X%*%v),0))
  expect_gt(score,0)
  expect_equal(score,f$certificate$separation)
  expect_lte(f$distance_lower,f$distance_upper+1e-12)
})

test_that("certificates allow fixed coordinates and nonunique multipliers", {
  X <- matrix(1,2,1)
  expect_true(calibration_certificate(X,c(1,1),2,c(2,0),lower=c(2,0),upper=c(2,0))$optimal)
  # No free coordinates; lambda=0 fails, but lambda=-2 certifies optimum.
  z <- calibration_certificate(X,c(3,1),0,c(0,0),lower=0,upper=2)
  expect_true(z$optimal)
  expect_lte(z$max_scaled_stationarity,1e-8)
  e <- matrix(numeric(),3,0)
  a <- calibrate(e,rep(1,3),numeric(),method="truncated",lower=0,upper=.5,solver="projection")
  expect_equal(a$weights,rep(.5,3))
})

test_that("mixed constraints use one QP and enforce inequalities independently", {
  X <- cbind(N=c(1,1),cap=c(1,0),floor=c(0,1))
  a <- calibrate(X,c(1,1),c(N=4,cap=1,floor=2),method="truncated",lower=0,upper=5,
                 direction=c("=","<=",">="))
  expect_equal(a$weights,c(1,3),tolerance=1e-8)
  expect_true(a$certificate$optimal); expect_equal(a$solver,"qp")
  expect_true(calibration_feasibility(X,c(N=4,cap=1,floor=2),0,5,
                                     direction=c("=","<=",">="))$feasible)
  expect_false(calibration_feasibility(X,c(N=4,cap=1,floor=5),0,5,
                                      direction=c("=","<=",">="))$feasible)
  expect_error(calibrate(X,c(1,1),c(N=4,cap=1,floor=2),method="raking",direction="<="),
               "quadratic")
  expect_error(calibrate(X,c(1,1),c(N=4,cap=1,floor=2),direction="<=",solver="dual"),
               "single QP")
  expect_error(calibration_feasibility(X,c(N=4,cap=1,floor=2),0,5,
                                      direction="<=",method="pfc"),"lp")
  expect_error(calibrate(X,c(1,1),c(N=4,cap=1,floor=2),direction=c("=","<=")),"direction")
})

test_that("mixed QP handles redundant equalities, signed weights and zero constraints", {
  X <- cbind(N=c(1,1,1),N2=c(2,2,2),cap=c(1,0,0),zero=0)
  a <- calibrate(X,c(-1,0,1),c(N=3,N2=6,cap=0,zero=0),metric_scale=1,
    lower=c(-2,-2,-2),upper=c(0,5,5),method="truncated",direction=c("=","=","<=","<="))
  expect_equal(a$weights,c(0,1,2),tolerance=1e-8)
  expect_true(a$certificate$optimal)
  f <- calibration_feasibility(matrix(1,1,1),-3,-Inf,Inf,direction="<=")
  expect_true(f$feasible); expect_lte(f$weights,-3+1e-8)
  expect_false(calibration_feasibility(matrix(0,1,1),-1,-Inf,Inf,direction="<=")$feasible)
  z <- calibrate(cbind(c(1,1),c(1,0)),c(1,1),c(3,1),method="truncated",
     lower=c(1,0),upper=c(1,4),direction=c("=","<="))
  expect_equal(z$weights,c(1,2),tolerance=1e-8)
})

test_that("nonlinear solver selection is explicit and truthful", {
  X <- matrix(1,3,1)
  expect_message(a <- calibrate(X,rep(1,3),6,method="raking"),"Newton")
  expect_equal(a$solver,"newton")
  expect_message(b <- calibrate(X,rep(1,3),6,method="logit",lower=0,upper=3,
                                solver="newton"),"Newton")
  expect_error(calibrate(X,rep(1,3),6,method="raking",solver="qp"),"newton")
})

test_that("chosen feasibility backend reaches all workflow stages", {
  X <- diag(3); colnames(X) <- c("A","B","C"); tt <- c(A=0,B=0,C=6)
  for (m in c("lp","pfc","fixing","projection")) {
    a <- automatic_calibration(X,tt,rep(1,3),lower=1,upper=3,bounds=c(1,3),
      search_limits=c(1,3),margin_groups=list(c("A","B","C")),ordered=TRUE,
      feasibility_method=m,max_depth=2)
    expect_length(a,1)
    expect_equal(a[[1]]$weights,rep(2,3),tolerance=1e-6)
  }
  expect_error(automatic_calibration(X,tt,rep(1,3),lower=1,upper=3,
    feasibility_method="invalid"),"arg")
})

test_that("random feasible problems agree with independent QP solutions", {
  for(seed in 1:12) {
    a <- generate_calibration_problem(n=25,categorical_groups=c(3,2),n_binary=1,
      n_continuous=2,weight_cv=.2,difficulty=.45,seed=seed)
    fit <- calibrate(a$X,a$d,a$totals,a$q,"truncated",a$lower,a$upper,solver="projection",
      projection_action="refine")
    qp <- calibrate(a$X,a$d,a$totals,a$q,"truncated",a$lower,a$upper,solver="qp")
    expect_equal(fit$weights,qp$weights,tolerance=1e-5)
    for(m in c("pfc","lp")) expect_true(calibration_feasibility(a$X,a$totals,a$lower,a$upper,method=m)$feasible)
    dirs <- rep(c("=","<=",">="),length.out=ncol(a$X))
    mix <- calibrate(a$X,a$d,a$totals,a$q,"truncated",a$lower,a$upper,direction=dirs)
    expect_true(mix$certificate$optimal)
    expect_lte(mix$objective,qp$objective+1e-6)
  }
})
