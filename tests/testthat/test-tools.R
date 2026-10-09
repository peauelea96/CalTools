test_that("LP feasibility supports signed bounds and zero/negative totals", {
  X <- cbind(1, -1:1)
  for (meth in c("auto", "lp")) {
    f <- calibration_feasibility(X, c(-3,0), -2, 0, method = meth)
    expect_true(f$feasible)
    expect_equal(as.numeric(crossprod(X,f$weights)),c(-3,0),tolerance=1e-7)
    expect_false(calibration_feasibility(X,c(9,0),-2,2,method=meth)$feasible)
    # Individual ranges pass but joint equalities are impossible.
    expect_false(calibration_feasibility(cbind(c(1,1),c(1,1)),c(1,2),0,3,method=meth)$feasible)
  }
  expect_true(is_feasible_lp(X,c(0,0),lower=-1,upper=1))
  expect_true(is_feasible_fast(X,c(0,0),lower=-1,upper=1))
  expect_error(calibration_feasibility(X,c(0,0),1,-1), "Bounds")
})

test_that("random constructed feasible LPs return verified witnesses", {
  set.seed(903)
  for (k in 1:20) {
    X <- matrix(rnorm(45),15,3); L <- runif(15,-3,-1); U <- runif(15,1,3)
    w <- runif(15,L,U);t <- as.numeric(crossprod(X,w))
    f <- calibration_feasibility(X,t,L,U,method="lp",tol=1e-7)
    expect_true(f$feasible)
    expect_lt(f$max_scaled_residual,1e-7)
  }
})

test_that("minimal directional relaxation has an analytic value", {
  z <- relax_bounds(matrix(1,3,1),9,0,2,lower_rate=0)
  expect_equal(z$status,"optimal")
  expect_equal(z$rho,1,tolerance=1e-7)
  expect_equal(sum(z$weights),9,tolerance=1e-7)
  z <- relax_bounds(matrix(1,2,1),-4,0,1,lower_rate=c(1,3),upper_rate=0)
  expect_equal(z$rho,1,tolerance=1e-7)
  expect_false(calibration_feasibility(matrix(1,2,1),-4,-.99*c(1,3),1)$feasible)
})

test_that("bound grids distinguish absolute and multiplicative endpoints", {
  X <- matrix(1,2,1)
  z <- analyze_feasibility_bounds(X,1,c(10,20),type="w-weight",method="lower",
    bounds=c(1,2), search_limits=0, step=.25)
  expect_equal(z$bounds,c(.5,2))
  z <- analyze_feasibility_bounds(X,7,c(1,1),type="w-weight",method="upper",
    bounds=c(0,1),search_limits=3.5,step=1)
  expect_equal(z$bounds,c(0,3.5))
  z <- analyze_feasibility_bounds(diag(2),c(-2,2),c(-1,1),type="g-weight",method="upper",
    bounds=c(0,1),search_limits=2,step=.5)
  expect_equal(z$lower,c(-2,0));expect_equal(z$upper,c(0,2))
})

test_that("diagnostics and aggregation preserve dimensions and names", {
  X <- cbind(a=c(1,0),b=c(0,1),c=c(1,1),zero=0)
  t <- c(a=2,b=3,c=5,zero=0)
  expect_false(analyze_inconsistency(X,t))
  expect_true(analyze_inconsistency(X,t+c(0,0,1,0)))
  rel <- analyze_collinearity(X,simple=TRUE)
  expect_equal(rel$rank,2)
  expect_equal(unname(X %*% t(rel$relations)),matrix(0,2,2),tolerance=1e-7)
  expect_equal(analyze_zero_sample(X)$names,"zero")
  z <- drop_margins(X,t,c("a","b","c"))
  expect_equal(dim(z$X),c(2L,1L))
  z <- drop_margins(X,t,colnames(X))
  expect_equal(dim(z$X),c(2L,0L))
  z <- merge_margins(X[,1:2],t[1:2],c("a","b"),"N")
  expect_equal(dim(z$X),c(2L,1L));expect_equal(unname(z$totals),5)
  expect_error(merge_margins(X,t,c("a","a")),"distinct")
  expect_error(merge_margins(X,t,c("a","b"),"c"),"unique")
  expect_length(adjacent_elements("a",1),0)
})

test_that("BFS reaches depth two and keeps every original residual", {
  X <- diag(3);colnames(X)<-c("a","b","c");t<-c(a=0,b=0,c=6)
  out <- automatic_calibration(X,t,rep(1,3),lower=1,upper=3,
    bounds=c(1,3),search_limits=c(1,3),margin_groups=list(c("a","b","c")),
    ordered=TRUE,max_solutions=3,max_depth=2)
  expect_length(out,1)
  expect_length(out[[1]]$history,2)
  expect_equal(sum(out[[1]]$weights),6,tolerance=1e-7)
  expect_equal(out[[1]]$residuals_original,c(a=2,b=2,c=-4),tolerance=1e-7)
  expect_equal(attr(out,"search_status"),"exhausted")
  pick <- interactive_calibration(X,t,rep(1,3),lower=1,upper=3,
    margin_groups=list(c("a","b","c")),ordered=TRUE,
    bounds=c(1,3),search_limits=c(1,3), choose=function(z) 1L)
  expect_equal(pick$weights,out[[1]]$weights)
})

test_that("recommendations never join a variable to itself on tied scores", {
  X<-diag(3);colnames(X)<-letters[1:3];t<-setNames(rep(1,3),letters[1:3])
  for (m in c("sample","lambda")) {
    z<-recommend_margin_merges(X,rep(1,3),t,margin_groups=list(letters[1:3]),ordered=FALSE,criterion=m)
    expect_length(unique(z[[1]]),2)
  }
})

test_that("ordered actions use supplied bounds and empty lists are defined", {
  X<-diag(2);colnames(X)<-c("a","b");t<-c(a=0,b=4)
  z<-check_modified_feasibility(X,t,c(1,1),1,3,c("a","b"),"join",bounds=c(1,3),search_limits=3,step=1)
  expect_equal(z$status,0)
  z<-check_sequential_feasibility(X,t,c(1,1),1,3,list(c("a","b")),"join",list("N"),
    bounds=c(1,3),search_limits=3,step=1)
  expect_equal(z$status,0);expect_equal(z$iter,1)
  z<-check_sequential_feasibility(X,t,c(1,1),1,3,list(),character(),bounds=c(1,3),search_limits=3,step=1)
  expect_equal(z$iter,0)
})

test_that("inclusion probability generator terminates at both endpoints", {
  expect_equal(random_inclusion_probabilities(4,0),rep(0,4))
  expect_equal(random_inclusion_probabilities(4,4),rep(1,4))
  set.seed(10)
  for (n in c(.1,1,5,9.99)) {
    p<-random_inclusion_probabilities(10,n)
    expect_equal(sum(p),n,tolerance=1e-12);expect_true(all(p>=0 & p<=1))
  }
  expect_error(random_inclusion_probabilities(3,4),"exceed")
})
