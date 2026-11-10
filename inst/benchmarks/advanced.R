# Run from the package project root after installation: Rscript inst/benchmarks/advanced.R
# Timings include validation, rank reduction and certification, but exclude generation.
library(CalTools)
dir.create("validation",showWarnings=FALSE)
rows <- list(); k <- 0L
record <- function(case,n,p,kind,requested,fun,repetitions=5L,batch=3L) {
  # Warm up code paths. Keep failure and warning outcomes in the measured record.
  run <- function() tryCatch(suppressMessages(suppressWarnings(fun())),error=function(e)e)
  run()
  for(rep in seq_len(repetitions)) {
    elapsed <- system.time(for(j in seq_len(batch)) ans <- run())[["elapsed"]]/batch
    failed <- inherits(ans,"error")
    k <<- k+1L
    rows[[k]] <<- data.frame(case=case,n=n,p=p,kind=kind,requested=requested,repetition=rep,
      seconds=elapsed,actual=if(failed) "error" else if(is.null(ans$solver)) "workflow" else ans$solver,
      status=if(failed) conditionMessage(ans) else if(is.null(ans$status)) attr(ans,"search_status") else ans$status,
      feasible=if(failed) NA else if(is.null(ans$feasible)) NA else ans$feasible,
      objective=if(failed||is.null(ans$objective)) NA_real_ else ans$objective,
      residual=if(failed||is.null(ans$max_scaled_residual)) NA_real_ else ans$max_scaled_residual)
  }
}
for(n in c(50L,200L,600L)) for(difficulty in c(.25,.49)) {
  a <- generate_calibration_problem(n=n,categorical_groups=c(5,4),n_binary=2,
    n_continuous=3,weight_cv=.3,difficulty=difficulty,seed=1200+n)
  for(solver in c("dual","qp","projection","projection_refine")) {
    ss <- if(solver=="projection_refine") "projection" else solver
    action <- if(solver=="projection_refine") "refine" else "report"
    record(paste0("generated_",difficulty),n,ncol(a$X),"calibration",solver,
      function() calibrate(a$X,a$d,a$totals,a$q,"truncated",a$lower,a$upper,
                           solver=ss,projection_action=action))
  }
  for(method in c("pfc","lp","projection","fixing"))
    record(paste0("generated_",difficulty),n,ncol(a$X),"feasibility",method,
      function() calibration_feasibility(a$X,a$totals,a$lower,a$upper,method=method,fallback=FALSE))
  dirs <- rep(c("=","<=",">="),length.out=ncol(a$X))
  record(paste0("generated_",difficulty),n,ncol(a$X),"mixed","qp",
    function() calibrate(a$X,a$d,a$totals,a$q,"truncated",a$lower,a$upper,direction=dirs))
  record(paste0("generated_",difficulty),n,ncol(a$X),"mixed_feasibility","lp",
    function() calibration_feasibility(a$X,a$totals,a$lower,a$upper,direction=dirs))
}
# These cases are deliberately not caught by checking a single margin's range.
X <- rbind(c(1,0),c(0,1),c(1,1))
for(m in c("pfc","lp","projection","fixing"))
  record("joint_infeasible",3,2,"feasibility",m,
    function() calibration_feasibility(X,c(2,0),0,1,method=m,fallback=FALSE))
X <- t(rbind(c(2,0,-1,1,0),c(1,-3,0,0,0),c(1,0,0,1,0),c(2,0,0,0,3)))
for(m in c("pfc","lp","projection","fixing"))
  record("boundary_feasible",5,4,"feasibility",m,
    function() calibration_feasibility(X,c(7/4,1,7/4,2),0,1,method=m,fallback=FALSE,max_iter=1000))
# Irreversible projection must not be compared as though it solved the same optimum.
X <- matrix(c(3,0,1,1,1,0,1,1,0,3),5,2)
for(s in c("dual","qp","projection")) record("suboptimal_projection",5,2,"calibration",s,
  function() calibrate(X,rep(1,5),c(.8,.4),method="truncated",lower=0,upper=1.2,solver=s))
# Measure the complete automatic search with its selectable feasibility backend.
a <- generate_calibration_problem(n=60,categorical_groups=c(4,3),n_continuous=1,
                                  difficulty=.8,seed=751)
for(m in c("pfc","lp","projection","fixing"))
  record("automatic_search",60,ncol(a$X),"workflow",m,function()
    automatic_calibration(a$X,a$totals,a$d,q=a$q,lower=a$lower,upper=a$upper,
      bounds=a$bounds,method="upper",search_limits=3,step=.1,
      margin_groups=a$margin_groups,ordered=a$ordered,max_solutions=2,max_depth=1,
      feasibility_method=m))
raw <- do.call(rbind,rows)
write.csv(raw,"validation/advanced-benchmarks.csv",row.names=FALSE)
key <- interaction(raw$case,raw$n,raw$kind,raw$requested,drop=TRUE)
summary <- do.call(rbind,lapply(split(raw,key),function(x) data.frame(
  x[1,c("case","n","p","kind","requested")],median_seconds=median(x$seconds),
  min_seconds=min(x$seconds),max_seconds=max(x$seconds),
  actual=paste(unique(x$actual),collapse=";"),status=paste(unique(x$status),collapse=";"))))
rownames(summary) <- NULL
write.csv(summary,"validation/advanced-benchmark-summary.csv",row.names=FALSE)
capture.output(sessionInfo(),file="validation/session-info-0.3.0.txt")
print(summary,row.names=FALSE)

# Compact status-aware table for the README; no optional packages needed.
selected <- summary[(summary$n==200 & summary$case=="generated_0.49") |
  summary$case %in% c("joint_infeasible","boundary_feasible","suboptimal_projection","automatic_search"),]
lines <- c("| Case | n | Task | Requested | Median ms | Actual route | Status |",
           "|---|---:|---|---|---:|---|---|")
for(i in seq_len(nrow(selected))) {
  x <- selected[i,]
  lines <- c(lines,paste0("| ",paste(c(x$case,x$n,x$kind,x$requested,
    sprintf("%.3f",1000*x$median_seconds),x$actual,x$status),collapse=" | ")," |"))
}
writeLines(lines,"validation/benchmark-table.md")
