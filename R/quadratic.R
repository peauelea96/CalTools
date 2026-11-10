# Sparse quadratic infrastructure. Constraint matrices have controls in rows.
.ct_qproblem <- function(X, d, totals, q, lower, upper, direction,
                        metric_scale, total_lower, total_upper) {
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Install package Matrix")
  if (length(dim(X)) != 2L || nrow(X) < 1L) stop("X must be a matrix with rows")
  if (!is.numeric(X) && !inherits(X,"Matrix")) stop("X must be numeric")
  X <- methods::as(Matrix::Matrix(X, sparse = TRUE), "dgCMatrix")
  if (any(!is.finite(X@x))) stop("X must be finite")
  n <- nrow(X); p <- ncol(X)
  d <- .ct_vec(d, n, "d", scalar = FALSE); q <- .ct_vec(q, n, "q")
  h <- .ct_metric(d, q, metric_scale)
  bounds <- .ct_bounds(lower, upper, n)
  direction <- .ct_direction(direction, p)
  if (is.null(total_lower) != is.null(total_upper))
    stop("Provide both total_lower and total_upper")
  if (!is.null(total_lower)) {
    tb <- .ct_bounds(total_lower, total_upper, p); lo <- tb$L; up <- tb$U
  } else {
    totals <- .ct_vec(totals, p, "totals", scalar = FALSE)
    if (!is.null(names(totals)) && !is.null(colnames(X)) &&
        !identical(names(totals), colnames(X))) stop("Names and order of totals must match columns of X")
    lo <- ifelse(direction == "<=", -Inf, totals)
    up <- ifelse(direction == ">=", Inf, totals)
  }
  R <- Matrix::t(X)
  list(X = X, R = R, d = d, h = h, q = q, L = bounds$L, U = bounds$U,
       lo = lo, up = up, n = n, p = p)
}
.ct_max <- function(x) max(c(0, x))
.ct_qdiagnostic <- function(p, w, y, z) {
  rw <- as.numeric(p$R %*% w)
  scale <- pmax(1, ifelse(is.finite(p$lo), abs(p$lo), 0),
                ifelse(is.finite(p$up), abs(p$up), 0))
  grad <- (w - p$d) / p$h
  station <- grad + as.numeric(Matrix::t(p$R) %*% y) + z
  comp <- function(y, x, l, u) {
    a <- y > 0 & is.finite(u); b <- y < 0 & is.finite(l)
    c(y[a] * (u[a] - x[a]), -y[b] * (x[b] - l[b]))
  }
  list(equalities = .ct_max(abs(rw[p$lo == p$up] - p$lo[p$lo == p$up]) / scale[p$lo == p$up]),
       inequalities = .ct_max(pmax(p$lo-rw, rw-p$up, 0) / scale),
       bounds = .ct_max(pmax(p$L-w, w-p$U, 0) / pmax(1,abs(w))),
       stationarity = .ct_max(abs(station)) / max(1, .ct_max(abs(grad))),
       complementarity = .ct_max(abs(c(comp(y,rw,p$lo,p$up), comp(z,w,p$L,p$U)))),
       dual_feasibility = .ct_max(c(pmax(y[!is.finite(p$up)],0),
         pmax(-y[!is.finite(p$lo)],0),pmax(z[!is.finite(p$U)],0),
         pmax(-z[!is.finite(p$L)],0))),
       achieved = rw)
}
.ct_osqp <- function(p, tol, max_iter, control, warm_start) {
  if (!requireNamespace("osqp", quietly = TRUE)) stop("Install optional package osqp")
  backend_start <- proc.time()[[3L]]
  P <- Matrix::Diagonal(p$n, 1/p$h)
  B <- rbind(p$R, Matrix::Diagonal(p$n))
  args <- list(eps_abs = tol * 0.1, eps_rel = tol * 0.1,
               max_iter = max_iter, verbose = FALSE,
               eps_prim_inf = tol * 0.1, eps_dual_inf = tol * 0.1)
  modern <- utils::packageVersion("osqp") >= "1.0.0"
  args[[if (modern) "polishing" else "polish"]] <- TRUE
  args <- utils::modifyList(args, control)
  model <- osqp::osqp(P = P, q = -p$d/p$h, A = B,
                     l = c(p$lo,p$L), u = c(p$up,p$U), pars = args)
  if (!is.null(warm_start)) {
    x <- if (is.list(warm_start)) warm_start$weights else warm_start
    x <- .ct_vec(x,p$n,"warm_start",scalar=FALSE)
    if (modern) model@WarmStart(x=x) else model$WarmStart(x=x)
  }
  tm <- proc.time()[[3L]]
  out <- if (modern) model@Solve() else model$Solve()
  elapsed <- proc.time()[[3L]]-tm
  state <- tolower(out$info$status)
  status <- if (state == "solved") "candidate" else if (state == "solved inaccurate")
    "approximate" else "not_converged"
  cert <- out$prim_inf_cert
  verified <- FALSE
  if (grepl("primal infeasible",state) && length(cert) == nrow(B) && all(is.finite(cert))) {
    cert <- cert / max(1,max(abs(cert)))
    l <- c(p$lo,p$L); u <- c(p$up,p$U)
    support <- sum(cert[cert>0]*u[cert>0]) + sum(cert[cert<0]*l[cert<0])
    # Finite-precision Farkas certificate, with a strict negative support margin.
    verified <- is.finite(support) && support < -tol &&
      max(abs(as.numeric(Matrix::t(B) %*% cert))) <= tol
    if (verified) status <- "infeasible_certified"
  }
  list(w=out$x, y=out$y[seq_len(p$p)], z=out$y[p$p+seq_len(p$n)],
       iterations=out$info$iter, status=status, backend=out$info,
       infeasibility_certificate=list(vector=cert,verified=verified,type="numerical_Farkas",tolerance=tol), preparation_time=tm-backend_start, solve_time=elapsed)
}

#' Sparse quadratic calibration
#'
#' Minimize sum((w-d)^2/(2*d*q)) with simultaneous equalities, inequalities,
#' intervals and individual weight bounds. Explicit metric_scale preserves the
#' package's alternative positive metric for signed initial weights.
#' @inheritParams calibrate_weights
#' @param totals Totals; may be NULL when interval endpoints are supplied.
#' @param total_lower,total_upper Optional per-control interval endpoints;
#'   these replace totals and direction. Both must be supplied.
#' @param solver auto, quadprog, osqp, dykstra or experimental dual_newton.
#'   auto uses quadprog for small dense problems, otherwise OSQP if installed,
#'   or Dykstra.
#' @param max_iter Iteration budget; exhaustion does not imply infeasibility.
#' @param control OSQP settings, or regularization and fallback_osqp for Newton.
#' @param warm_start OSQP weights or fit; for Newton a fit with compatible
#'   duals. Unsupported for Dykstra and quadprog.
#' @return A caltools_fit with weights, factors, objective, termination state,
#'   residuals, KKT diagnostics, duals and timing. Nonconvergence returns a result.
#' @details Dykstra retains corrections in the H-inverse metric. Newton
#'   uses a generalized Hessian, regularization and projected Armijo steps;
#'   no global convergence guarantee is asserted. OSQP inaccurate states
#'   remain approximate. Infeasibility certificates are numerical at tol.
#' @examples
#' X <- cbind(rep(1, 5), -2:2)
#' f <- calibrate_quadratic(X, rep(1, 5), c(6, 1),
#'   lower = 0, upper = 2, solver = "dykstra")
#' calibration_diagnostics(f)$totals
#' @export
calibrate_quadratic <- function(X, d, totals=NULL, q=1, lower=-Inf, upper=Inf,
  direction="=", solver="auto", tol=1e-7, max_iter=20000L,
  metric_scale=NULL, total_lower=NULL, total_upper=NULL,
  control=list(), warm_start=NULL) {
  start <- proc.time()[[3L]]
  .ct_scalar(tol,"tol"); .ct_scalar(max_iter,"max_iter",integer=TRUE)
  p <- .ct_qproblem(X,d,totals,q,lower,upper,direction,metric_scale,total_lower,total_upper)
  solver <- match.arg(solver,c("auto","quadprog","osqp","dykstra","dual_newton"))
  requested <- solver
  reason <- "explicit selection"
  if (solver == "auto") {
    sparse <- inherits(X,"sparseMatrix")
    density <- length(p$X@x)/max(1,p$n*p$p)
    if (p$n <= 200L && p$p <= 20L && !sparse && density > .2) {
      solver <- "quadprog"; reason <- "small dense reference problem"
    } else if (requireNamespace("osqp",quietly=TRUE)) {
      solver <- "osqp"; reason <- "sparse or larger problem; OSQP available"
    } else {
      solver <- "dykstra"; reason <- "OSQP unavailable; Dykstra uses sparse projections"
    }
  }
  fallback <- NULL
  use_fallback <- identical(control$fallback_osqp,TRUE)
  if (!is.null(control$fallback_osqp)) {
    if (solver != "dual_newton") stop("fallback_osqp applies only to dual_newton")
    control$fallback_osqp <- NULL
  }
  prepared <- proc.time()[[3L]]-start
  ans <- switch(solver,quadprog=.ct_quadprog(p,tol,max_iter,control,warm_start),osqp=.ct_osqp(p,tol,max_iter,control,warm_start),
                dykstra=.ct_dykstra(p,tol,max_iter,control,warm_start),
                dual_newton=.ct_dual_newton(p,tol,max_iter,control,warm_start))
  if (solver == "dual_newton" && ans$status != "candidate" && use_fallback &&
      requireNamespace("osqp",quietly=TRUE)) {
    fallback <- list(from="dual_newton",status=ans$status,iterations=ans$iterations,
      solve_time=ans$solve_time,preparation_time=ans$preparation_time)
    ans <- .ct_osqp(p,tol,max_iter,list(),ans$w)
    ans$solve_time <- ans$solve_time+fallback$solve_time
    ans$preparation_time <- ans$preparation_time+fallback$preparation_time
    solver <- "osqp"
  }
  if (grepl("infeasible",tolower(ans$backend$status)) || any(!is.finite(ans$w))) {
    ans$w <- rep(NA_real_,p$n); ans$y <- rep(NA_real_,p$p); ans$z <- rep(NA_real_,p$n)
    kkt <- list(equalities=Inf,inequalities=Inf,bounds=Inf,stationarity=Inf,
                complementarity=Inf,dual_feasibility=Inf,achieved=rep(NA_real_,p$p))
  } else kkt <- .ct_qdiagnostic(p,ans$w,ans$y,ans$z)
  error <- max(unlist(kkt[1:6]))
  status <- ans$status
  if (status == "candidate") status <- if (is.finite(error) && error <= tol) "optimal_numerical" else "approximate"
  g <- ans$w / p$d; g[p$d == 0] <- NA_real_
  structure(list(g=g,lambda=-ans$y,weights=ans$w,
    objective=sum((ans$w-p$d)^2/(2*p$h)),status=status,
    converged=status == "optimal_numerical",iterations=ans$iterations,
    solver=solver,requested_solver=requested,selection_reason=reason,method="truncated",lower=p$L,upper=p$U,
    residuals=kkt$achieved-ifelse(is.finite(p$lo),p$lo,p$up),
    max_scaled_residual=max(unlist(kkt[1:3])),metric_scale=p$h/p$q,q=p$q,
    kkt=kkt,duals=list(controls=ans$y,bounds=ans$z),
    timing=list(total=proc.time()[[3L]]-start,preparation=prepared+if(is.null(ans$preparation_time)) 0 else ans$preparation_time,solve=ans$solve_time),
    backend=ans$backend,fallback=fallback,infeasibility_certificate=ans$infeasibility_certificate,
    problem=p),class="caltools_fit")
}

.ct_blocks <- function(p) {
  eq <- which(p$lo == p$up)
  hi <- which(p$lo != p$up & is.finite(p$up))
  low <- which(p$lo != p$up & is.finite(p$lo))
  E <- p$R[eq,,drop=FALSE]
  C <- rbind(p$R[hi,,drop=FALSE], -p$R[low,,drop=FALSE])
  list(E=E,C=C,t=p$lo[eq],b=c(p$up[hi],-p$lo[low]),
       eq=eq,index=c(hi,low),sign=c(rep(1,length(hi)),rep(-1,length(low))))
}
.ct_map_dual <- function(p, blocks, alpha, beta) {
  y <- numeric(p$p); y[blocks$eq] <- alpha
  for (j in seq_along(beta)) y[blocks$index[j]] <-
    y[blocks$index[j]] + blocks$sign[j]*beta[j]
  y
}
.ct_dykstra <- function(p, tol, max_iter, control, warm_start) {
  if (length(control)) stop("Dykstra has no backend-specific control options")
  if (!is.null(warm_start)) stop("Dykstra must start from d; warm_start is unsupported")
  start <- proc.time()[[3L]]
  bl <- .ct_blocks(p); E <- bl$E
  HEt <- Matrix::Diagonal(p$n,p$h) %*% Matrix::t(E)
  gram <- as.matrix(E %*% HEt)
  # Factor the small control-space matrix once, retaining redundant controls.
  ss <- if (nrow(gram)) eigen(gram,symmetric=TRUE) else NULL
  keep <- if (is.null(ss)) logical() else ss$values > max(ss$values,0)*1e-12
  Ginv <- if (any(keep)) tcrossprod(sweep(ss$vectors[,keep,drop=FALSE],2,
                            sqrt(ss$values[keep]),"/")) else matrix(0,nrow(gram),nrow(gram))
  Ct <- methods::as(Matrix::t(bl$C),"dgCMatrix")
  columns <- lapply(seq_len(ncol(Ct)),function(j) {
    k <- seq.int(Ct@p[j]+1L,Ct@p[j+1L])
    if (Ct@p[j] == Ct@p[j+1L]) k <- integer()
    i <- Ct@i[k]+1L; v <- Ct@x[k]
    list(i=i,v=v,hv=p$h[i]*v,den=sum(p$h[i]*v*v))
  })
  alpha <- numeric(nrow(E)); beta <- numeric(nrow(bl$C))
  box <- numeric(p$n); w <- p$d; status <- "not_converged"
  solve_start <- proc.time()[[3L]]
  for (it in seq_len(max_iter)) {
    previous <- w
    if (length(alpha)) {
      z <- w + as.numeric(HEt %*% alpha)
      delta <- as.numeric(Ginv %*% (bl$t-as.numeric(E %*% z)))
      w <- z+as.numeric(HEt %*% delta); alpha <- -delta
    }
    z <- w+box; w <- .ct_clip(z,p$L,p$U); box <- z-w
    for (j in seq_along(beta)) {
      c <- columns[[j]]; i <- c$i
      if (c$den == 0) next
      z <- w[i]+c$hv*beta[j]
      beta[j] <- max(0,(sum(c$v*z)-bl$b[j])/c$den)
      w[i] <- z-c$hv*beta[j]
    }
    y <- .ct_map_dual(p,bl,alpha,beta)
    zdual <- box/p$h
    kkt <- .ct_qdiagnostic(p,w,y,zdual)
    if (max(unlist(kkt[1:6])) <= tol &&
        max(abs(w-previous)/pmax(1,abs(w))) <= tol) {
      status <- "candidate"; break
    }
    if (any(!is.finite(w))) break
  }
  list(w=w,y=y,z=zdual,iterations=it,status=status,
       backend=list(status=status),preparation_time=solve_start-start,solve_time=proc.time()[[3L]]-solve_start)
}

# Experimental projected semismooth Newton on the convex negative dual.
.ct_dual_newton <- function(p, tol, max_iter, control, warm_start) {
  unknown <- setdiff(names(control),c("regularization"))
  if (length(unknown)) stop("Unknown dual_newton control: ",paste(unknown,collapse=", "))
  reg <- if (is.null(control$regularization)) 1e-10 else control$regularization
  .ct_scalar(reg,"regularization")
  start <- proc.time()[[3L]]
  bl <- .ct_blocks(p); B <- rbind(bl$E,bl$C); rhs <- c(bl$t,bl$b)
  ne <- nrow(bl$E); m <- nrow(B); ii <- ne+seq_len(nrow(bl$C))
  y <- numeric(m)
  if (!is.null(warm_start)) {
    if (!is.list(warm_start) || length(warm_start$duals$controls) != p$p)
      stop("dual_newton warm_start requires a fit with compatible control duals")
    orig <- warm_start$duals$controls
    y <- c(orig[bl$eq],pmax(0,orig[bl$index]*bl$sign))
  }
  project <- function(y) {y[ii] <- pmax(0,y[ii]); y}
  evaluate <- function(y) {
    s <- as.numeric(Matrix::t(B) %*% y)
    raw <- p$d-p$h*s; w <- .ct_clip(raw,p$L,p$U)
    list(w=w,raw=raw,g=rhs-as.numeric(B %*% w),
         value=-sum((w-p$d)^2/(2*p$h)+s*w)+sum(rhs*y),
         z=-(w-p$d)/p$h-s)
  }
  full <- as.matrix(B %*% Matrix::Diagonal(p$n,p$h) %*% Matrix::t(B))
  lipschitz <- max(1,if(m) max(rowSums(abs(full))) else 0)
  state <- evaluate(y); status <- "not_converged"; stalled <- 0L
  solve_start <- proc.time()[[3L]]
  for (it in seq_len(max_iter)) {
    orig <- .ct_map_dual(p,bl,y[seq_len(ne)],y[ii])
    kkt <- .ct_qdiagnostic(p,state$w,orig,state$z)
    if (max(unlist(kkt[1:6])) <= tol) {status <- "candidate"; break}
    active <- c(seq_len(ne),ii[y[ii] > 0 | state$g[ii] < 0])
    free <- state$raw > p$L & state$raw < p$U
    M <- B[active,free,drop=FALSE]
    J <- as.matrix(M %*% Matrix::Diagonal(sum(free),p$h[free]) %*% Matrix::t(M))
    step <- numeric(m)
    if (length(active)) {
      diag(J) <- diag(J)+reg*max(1,max(abs(diag(J))))
      step[active] <- -.ct_pinv(J,state$g[active],tol=1e-14)
    }
    accepted <- FALSE
    for (strategy in 1:2) {
      if (strategy == 2) step <- -state$g/lipschitz
      rate <- 1
      for (ls in 1:45) {
        candidate <- project(y+rate*step); delta <- candidate-y
        slope <- sum(state$g*delta)
        next_state <- evaluate(candidate)
        if (slope < 0 && next_state$value <= state$value+1e-4*slope) {
          accepted <- TRUE; break
        }
        rate <- rate/2
      }
      if (accepted) break
    }
    if (!accepted) break
    change <- max(c(0,abs(next_state$w-state$w)/pmax(1,abs(state$w))))
    stalled <- if (change < .Machine$double.eps*16) stalled+1L else 0L
    y <- candidate; state <- next_state
    if (stalled >= 10L) break
  }
  orig <- .ct_map_dual(p,bl,y[seq_len(ne)],y[ii])
  list(w=state$w,y=orig,z=state$z,iterations=it,status=status,
       backend=list(status=status,experimental=TRUE),preparation_time=solve_start-start,solve_time=proc.time()[[3L]]-solve_start)
}

.ct_quadprog <- function(p,tol,max_iter,control,warm_start) {
  if (p$n > 2000L) stop("Dense quadprog is limited to 2000 observations; use osqp")
  if(length(control) || !is.null(warm_start)) stop("quadprog has no control or warm_start options")
  start <- proc.time()[[3L]]; bl <- .ct_blocks(p)
  X <- as.matrix(Matrix::t(rbind(bl$E,bl$C))); colnames(X) <- NULL
  fit <- tryCatch(.ct_mixed_fit(X,p$d,c(bl$t,bl$b),p$q,"truncated",p$L,p$U,
    p$h/p$q,tol,1e-10,c(rep("=",length(bl$t)),rep("<=",length(bl$b)))),error=identity)
  if (inherits(fit,"error")) return(list(w=.ct_clip(p$d,p$L,p$U),y=numeric(p$p),z=numeric(p$n),
    iterations=0L,status="not_converged",backend=list(status=conditionMessage(fit)),solve_time=proc.time()[[3L]]-start))
  y <- .ct_map_dual(p,bl,-fit$lambda[seq_along(bl$t)],-fit$lambda[length(bl$t)+seq_along(bl$b)])
  z <- -(fit$weights-p$d)/p$h-as.numeric(Matrix::t(p$R)%*%y)
  list(w=fit$weights,y=y,z=z,iterations=fit$iterations,status="candidate",
       backend=list(status="solved"),solve_time=proc.time()[[3L]]-start)
}

#' Inspect calibration totals, weights and optimality
#' @param fit A fit from calibrate_quadratic or a new sparse solver route.
#' @return A list with per-control totals, input/output weights and factors,
#'   multiplier tables, KKT diagnostics and summary statistics.
#' @export
calibration_diagnostics <- function(fit) {
  p <- fit$problem
  if (is.null(p)) stop("This fit has no stored problem; use calibration_certificate for legacy fits")
  before <- as.numeric(p$R %*% p$d); after <- as.numeric(p$R %*% fit$weights)
  list(totals=data.frame(control=if(is.null(colnames(p$X))) seq_len(p$p) else colnames(p$X),
    lower=p$lo,upper=p$up,before=before,after=after,
    violation=pmax(p$lo-after,after-p$up,0),lambda=fit$lambda),
    weights=data.frame(initial=p$d,calibrated=fit$weights,g=fit$g,change=fit$weights-p$d,
      lower=p$L,upper=p$U,bound_multiplier=fit$duals$bounds),
    factors=if(all(is.na(fit$g))) rep(NA_real_,5) else stats::quantile(fit$g,na.rm=TRUE),
    objective=fit$objective,status=fit$status,kkt=fit$kkt,timing=fit$timing)
}


.ct_augment_legacy <- function(fit,X,d,totals,started,direction) {
  p <- .ct_qproblem(X,d,totals,fit$q,fit$lower,fit$upper,direction,
                   fit$metric_scale,NULL,NULL)
  fit$problem <- p
  fit$duals <- list(controls=-fit$lambda,bounds=NULL)
  if (fit$method %in% c("linear","truncated")) {
    fit$duals$bounds <- -(fit$weights-d)/p$h+as.numeric(p$X %*% fit$lambda)
    fit$kkt <- .ct_qdiagnostic(p,fit$weights,fit$duals$controls,fit$duals$bounds)
  } else {
    fit$duals$bounds <- rep(NA_real_,length(d))
    fit$kkt <- list(note="Use the nonlinear distance gradient; quadratic KKT not applicable")
  }
  fit$termination <- if (fit$status == "optimal") "optimal_numerical" else "approximate"
  fit$timing <- list(total=proc.time()[[3L]]-started,preparation=NA_real_,solve=NA_real_)
  fit
}
