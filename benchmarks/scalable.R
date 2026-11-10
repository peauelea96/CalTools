# Run from the package root with an installed CalTools >= 0.4.0:
# Rscript benchmarks/scalable.R
# R allocation counts exclude native OSQP allocations; not a peak RSS measure.
library(CalTools)
library(Matrix)
set.seed(20261021)
sizes <- c(1000L,3000L,10000L,50000L,100000L)
profiles <- data.frame(p=c(6L,16L,16L),inequalities=c(2L,6L,6L),density=c(1,.05,.1),tight=c(FALSE,FALSE,TRUE))
solvers <- c("quadprog","osqp","dykstra","dual_newton")
rows <- list()
for(n in sizes) for(k in seq_len(nrow(profiles))) {
  pr <- profiles[k,]; prep <- proc.time()[[3L]]
  X <- rsparsematrix(n,pr$p-1L,density=pr$density)
  X <- cbind(Matrix(1,n,1,sparse=TRUE),X)
  d <- runif(n,.8,1.2); truth <- d*runif(n,.7,1.3)
  if(pr$tight) truth <- pmin(1.2,pmax(.8,d+.15*as.numeric(X[,2])))
  target <- as.numeric(crossprod(X,truth))
  direction <- c(rep("=",pr$p-pr$inequalities),rep("<=",pr$inequalities))
  # Some inequalities active, some inactive; a known feasible truth remains.
  target[direction=="<="] <- target[direction=="<="]+rep(c(0,.1),length.out=pr$inequalities)
  input_time <- proc.time()[[3L]]-prep
  for(repetition in 1:3) for(s in solvers) {
    record <- data.frame(n=n,p=pr$p,inequalities=pr$inequalities,density=pr$density,
       solver=s,repetition=repetition,tight=pr$tight,active_bounds=NA_integer_,input_seconds=input_time,preparation_seconds=NA_real_,solve_seconds=NA_real_,
       total_seconds=NA_real_,R_allocated_MB=NA_real_,matrix_MB=as.numeric(object.size(X))/2^20,
       objective=NA_real_,KKT=NA_real_,iterations=NA_integer_,status="skipped_dense_limit")
    if(s=="quadprog" && n>2000) {rows[[length(rows)+1L]]<-record; next}
    if(s=="osqp" && !requireNamespace("osqp",quietly=TRUE)) {
      record$status<-"dependency_unavailable"; rows[[length(rows)+1L]]<-record; next
    }
    mem<-tempfile(); gc(); Rprofmem(mem)
    fit<-tryCatch(calibrate_quadratic(X,d,target,lower=if(pr$tight) .8 else .1,upper=if(pr$tight) 1.2 else 3,
      direction=direction,solver=s,tol=1e-6,max_iter=1000L),error=identity)
    Rprofmem(NULL)
    allocation<-suppressWarnings(as.numeric(sub(" .*","",readLines(mem,warn=FALSE))))
    unlink(mem); record$R_allocated_MB<-sum(allocation,na.rm=TRUE)/2^20
    if(inherits(fit,"error")) record$status<-conditionMessage(fit) else {
      record$preparation_seconds<-fit$timing$preparation
      record$solve_seconds<-fit$timing$solve
      record$total_seconds<-fit$timing$total
      record$objective<-fit$objective; record$KKT<-max(unlist(fit$kkt[1:6]))
      record$active_bounds<-sum(abs(fit$weights-fit$lower)<1e-5 | abs(fit$weights-fit$upper)<1e-5); record$iterations<-fit$iterations; record$status<-fit$status
    }
    rows[[length(rows)+1L]]<-record
    print(record[,c("n","p","solver","total_seconds","status")],row.names=FALSE)
    write.csv(do.call(rbind,rows),"benchmarks/results.csv",row.names=FALSE)
  }
}
write.csv(do.call(rbind,rows),"benchmarks/results.csv",row.names=FALSE)
writeLines(capture.output(sessionInfo()),"benchmarks/sessionInfo.txt")
