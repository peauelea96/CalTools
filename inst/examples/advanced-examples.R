# Executable examples from the 0.3.0 README. One expected projection warning.
library(CalTools)

X <- cbind(N=c(1,1), cap=c(1,0), floor=c(0,1))
fit <- calibrate(X, d=c(1,1), totals=c(N=4,cap=1,floor=2),
  method="truncated", lower=c(0,0), upper=c(2,5),
  direction=c("=","<=",">="))
fit$weights                  # c(1,3)
fit$solver                   # "qp"
fit$certificate$optimal      # TRUE, within the declared numerical tolerance
crossprod(X, fit$weights)    # N=4; cap=1 <= 1; floor=3 >= 2


X <- cbind(N=rep(1,4), x_lo=1:4, x_hi=1:4)
fit <- calibrate(X, rep(2,4), c(N=8,x_lo=21,x_hi=23),
  method="truncated", lower=0, upper=c(4,5,5,6),
  direction=c("=",">=","<="))
stopifnot(abs(sum(fit$weights)-8)<1e-7,
          sum((1:4)*fit$weights)>=21-1e-7,
          sum((1:4)*fit$weights)<=23+1e-7)


X <- matrix(1,3,1)
r <- calibrate(X, rep(1,3), 6, method="raking", solver="newton")
l <- calibrate(X, rep(1,3), 6, method="logit", lower=0, upper=3)
# Each call announces the actual nonlinear solver.
# calibrate(X, rep(1,3), 6, method="raking", solver="qp") is an error.


X <- matrix(c(3,0,1,1,1, 0,1,1,0,3), nrow=5)
a <- calibrate(X, rep(1,5), c(.8,.4), method="truncated",
  lower=0, upper=1.2, solver="projection")
# Warning: feasible projection; optimality not certified.
a$weights   # c(0, 1/3, 1/15, 11/15, 0)
a$status    # "feasible_uncertified"
a$objective # 127/75

b <- calibrate(X, rep(1,5), c(.8,.4), method="truncated",
  lower=0, upper=1.2, solver="projection", projection_action="refine")
b$weights   # c(1/35, 38/105, 4/105, 71/105, 0)
b$objective # 71/42; improvement = 1/350
b$projection$certificate # retains the first attempt for audit


certificate <- calibration_certificate(X, rep(1,5), c(.8,.4), a$weights,
  lower=0, upper=1.2)
certificate[c("optimal", "primal_feasible", "max_scaled_stationarity",
              "duality_gap", "kkt_objective_bound")]


X <- cbind(N=c(1,1), cap=c(1,0), floor=c(0,1))
f <- calibration_feasibility(X, c(N=4,cap=1,floor=2), 0, 5,
  direction=c("=","<=",">="))
f$feasible # TRUE
f$weights  # some witness, not necessarily optimal calibration weights

bad <- calibration_feasibility(X, c(N=4,cap=1,floor=5), 0, 5,
  direction=c("=","<=",">="))
bad$feasible # FALSE


X <- t(rbind(c(2,0,-1,1,0), c(1,-3,0,0,0),
             c(1,0,0,1,0), c(2,0,0,0,3)))
totals <- c(7/4,1,7/4,2)
witness <- c(1,0,1,3/4,0)
stopifnot(max(abs(crossprod(X,witness)-totals))<1e-12)
calibration_feasibility(X, totals, 0, 1, method="fixing", fallback=FALSE)$feasible
# NA: the centre projection fixes a wrong face.
calibration_feasibility(X, totals, 0, 1, method="pfc", max_iter=10000,
                        tol=1e-7, fallback=FALSE)$feasible
# TRUE


X <- rbind(c(1,0),c(0,1),c(1,1))
f <- calibration_feasibility(X,c(2,0),0,1,method="pfc",fallback=FALSE)
v <- f$certificate$direction
sum(v*c(2,0)) - sum(pmax(as.numeric(X%*%v),0)) # strictly positive
f[c("distance_lower","distance_upper","iterations","status")]


p <- generate_calibration_problem(n=80, categorical_groups=c(5,4),
  n_binary=2, n_continuous=2, n_zero=1, n_derived=1,
  weight_cv=.25, difficulty=.45, seed=2026)
# Make the already vector-valued bounds still more heterogeneous without
# removing the constructive witness.
L <- p$lower-seq(0,.2,length.out=80)
U <- p$upper+seq(.2,0,length.out=80)
a <- calibrate(p$X,p$d,p$totals,p$q,"truncated",L,U,solver="dual")
b <- calibrate(p$X,p$d,p$totals,p$q,"truncated",L,U,solver="projection",
               projection_action="refine")
stopifnot(max(abs(a$weights-b$weights))<1e-4)

# Preserve some equalities, replace others by upper/lower controls.
dir <- rep(c("=","<=",">="),length.out=ncol(p$X))
target <- p$totals
slack <- .02*pmax(1,abs(target))
target[dir=="<="] <- target[dir=="<="]+slack[dir=="<="]
target[dir==">="] <- target[dir==">="]-slack[dir==">="]
mix <- calibrate(p$X,p$d,target,p$q,"truncated",L,U,direction=dir)
stopifnot(mix$certificate$optimal, mix$objective<=a$objective+1e-6)


p <- generate_calibration_problem(n=40,categorical_groups=c(4,3),
  n_continuous=1,difficulty=.8,seed=112)
solutions <- automatic_calibration(p$X,p$totals,p$d,q=p$q,
  lower=p$lower,upper=p$upper,bounds=p$bounds,method="upper",
  search_limits=3,step=.1,margin_groups=p$margin_groups,ordered=p$ordered,
  max_solutions=2,max_depth=1,feasibility_method="pfc",
  feasibility_max_iter=500,feasibility_fallback=TRUE,
  solver="projection",projection_action="refine")
