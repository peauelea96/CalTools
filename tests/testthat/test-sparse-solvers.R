test_that("new solvers agree with the quadratic reference and KKT", {
  solvers <- c("dykstra","dual_newton")
  if(requireNamespace("osqp",quietly=TRUE)) solvers <- c("osqp",solvers)
  set.seed(12)
  for(k in 1:20) {
    n<-sample(10:40,1); X<-cbind(1,matrix(rnorm(n*3),n)); d<-runif(n,.8,1.2)
    w<-runif(n,.5,1.5); tt<-drop(crossprod(X,w)); dir<-c("=","<=",">=","=")
    ref<-calibrate_weights(X,d,tt,method="truncated",lower=.2,upper=2,direction=dir,solver="qp")
    for(s in solvers) {
      a<-calibrate_quadratic(X,d,tt,lower=.2,upper=2,direction=dir,solver=s)
      expect_equal(a$weights,ref$weights,tolerance=1e-5,info=s)
      expect_equal(a$objective,ref$objective,tolerance=1e-6,info=s)
      expect_equal(a$status,"optimal_numerical",info=s)
      expect_lte(max(unlist(a$kkt[1:6])),1e-7)
    }
  }
})

test_that("redundancy, fixed bounds, negative and zero totals work", {
  X<-cbind(N=rep(1,8),x=-4:3,duplicate=2*(-4:3),zero=rep(0,8))
  d<-rep(1,8); w<-c(.8,rep(1,7)); t<-drop(crossprod(X,w))
  for(s in c("quadprog","osqp","dykstra","dual_newton")) {
    if(s=="osqp" && !requireNamespace("osqp",quietly=TRUE)) next
    fit<-calibrate_quadratic(Matrix::Matrix(X,sparse=TRUE),d,t,
      lower=c(.8,rep(.1,7)),upper=c(.8,rep(3,7)),solver=s)
    expect_equal(fit$status,"optimal_numerical",info=s)
    expect_equal(drop(crossprod(X,fit$weights)),unname(t),tolerance=1e-6,ignore_attr=TRUE)
    expect_equal(fit$weights[1],.8,tolerance=1e-7)
    expect_lte(max(unlist(fit$kkt[1:6])),1e-7)
  }
})

test_that("intervals, inactive inequalities and sparse matrices are supported", {
  X<-cbind(rep(1,5),-2:2)
  for(s in c("quadprog","osqp","dykstra","dual_newton")) {
    if(s=="osqp" && !requireNamespace("osqp",quietly=TRUE)) next
    a<-calibrate_quadratic(Matrix::Matrix(X,sparse=TRUE),rep(1,5),
      total_lower=c(5,-1),total_upper=c(6,1),lower=0,upper=3,solver=s)
    expect_equal(a$weights,rep(1,5),tolerance=1e-7)
    expect_equal(a$status,"optimal_numerical",info=s)
    expect_equal(calibration_diagnostics(a)$totals$violation,c(0,0),tolerance=1e-7)
  }
})

test_that("infeasibility is distinguished from iteration limits", {
  for(s in c("dykstra","dual_newton")) {
    a<-calibrate_quadratic(matrix(1,4,1),rep(1,4),10,lower=0,upper=1,solver=s,max_iter=10)
    expect_identical(a$status,"not_converged")
    expect_false(a$converged)
  }
  if(requireNamespace("osqp",quietly=TRUE)) {
    a<-calibrate_quadratic(matrix(1,4,1),rep(1,4),10,lower=0,upper=1,solver="osqp")
    expect_identical(a$status,"infeasible_certified")
    expect_true(a$infeasibility_certificate$verified)
    expect_true(all(is.na(a$weights)))
    a<-calibrate_quadratic(matrix(1,4,1),rep(1,4),10,solver="osqp",max_iter=1)
    expect_false(a$converged)
    expect_identical(a$status,"not_converged")
  }
})

test_that("API dispatch, warm starts and objective preserve conventions", {
  X<-cbind(1,1:8); d<-rep(2,8); t<-c(18,90)
  a<-calibrate_weights(X,d,t,method="truncated",solver="dykstra",lower=0,upper=5)
  expect_equal(a$g,a$weights/d)
  expect_equal(a$objective,sum((a$weights-d)^2/(2*d)))
  expect_equal(calibrate(X,d,t,method="truncated",solver="dual_newton",lower=0,upper=5)$weights,a$weights,tolerance=1e-6)
  expect_error(calibrate_weights(X,d,t,method="logit",solver="dykstra"),"quadratic")
  if(requireNamespace("osqp",quietly=TRUE)) {
    b<-calibrate_quadratic(X,d,t,solver="osqp")
    c<-calibrate_quadratic(X,d,t,solver="osqp",warm_start=b)
    expect_equal(c$weights,b$weights,tolerance=1e-6)
  }
})

test_that("ill conditioning never produces an unchecked optimal state", {
  set.seed(44); x<-rnorm(15); X<-cbind(1,x,x+1e-7*rnorm(15))
  target<-drop(crossprod(X,runif(15,.8,1.2)))
  for(s in c("dykstra","dual_newton","osqp")) {
    if(s=="osqp" && !requireNamespace("osqp",quietly=TRUE)) next
    a<-calibrate_quadratic(X,rep(1,15),target,lower=.1,upper=3,solver=s,max_iter=500)
    expect_true(a$status %in% c("optimal_numerical","approximate","not_converged"))
    if(a$converged) expect_lte(max(unlist(a$kkt[1:6])),1e-7)
  }
})

test_that("active vector bounds match quadprog", {
  set.seed(63); X<-cbind(1,seq(-1,1,length.out=20))
  d<-runif(20,.5,1.5); L<-seq(.1,.3,length.out=20); U<-seq(1.2,1.8,length.out=20)
  truth<-pmin(U,pmax(L,d+.8*X[,2])); target<-drop(crossprod(X,truth))
  ref<-calibrate_quadratic(X,d,target,lower=L,upper=U,solver="quadprog")
  expect_true(ref$converged)
  for(s in c("dykstra","dual_newton","osqp")) {
    if(s=="osqp" && !requireNamespace("osqp",quietly=TRUE)) next
    a<-calibrate_quadratic(X,d,target,lower=L,upper=U,solver=s)
    expect_equal(a$weights,ref$weights,tolerance=1e-5,info=s)
    expect_true(a$converged,info=s)
    expect_lte(max(unlist(a$kkt[1:6])),1e-7)
  }
})

test_that("empty controls and diagnostics for legacy fits are supported", {
  for(s in c("dykstra","dual_newton")) {
    a<-calibrate_quadratic(matrix(numeric(),3,0),rep(1,3),numeric(),lower=0,upper=.8,solver=s)
    expect_equal(a$weights,rep(.8,3)); expect_true(a$converged)
  }
  a<-calibrate_weights(matrix(1,3,1),rep(1,3),4)
  expect_equal(calibration_diagnostics(a)$totals$after,4)
  expect_equal(a$status,"optimal")
})

test_that("nonconstant q uses the same objective and orientation", {
  d<-c(1,2,3,4); q<-c(.5,2,1,3); X<-matrix(1,4,1)
  expected<-d+d*q*(15-sum(d))/sum(d*q)
  for(s in c("quadprog","dykstra","dual_newton","osqp")) {
    if(s=="osqp" && !requireNamespace("osqp",quietly=TRUE)) next
    a<-calibrate_quadratic(X,d,15,q=q,solver=s)
    expect_equal(a$weights,expected,tolerance=1e-6)
    expect_equal(a$objective,sum((expected-d)^2/(2*d*q)),tolerance=1e-7)
    expect_true(a$converged)
  }
})
