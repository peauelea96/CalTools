# Constraint convention: X[,j]'w direction[j] totals[j].
.ct_direction <- function(direction, p) {
  if (!is.character(direction) || anyNA(direction) ||
      !(length(direction) == 1L || length(direction) == p) ||
      any(!direction %in% c("=", "<=", ">=")))
    stop("direction must contain '=', '<=' or '>=': scalar or one per column")
  rep_len(direction, p)
}
.ct_violation <- function(X, w, totals, direction) {
  r <- as.numeric(crossprod(X, w)) - totals
  z <- abs(r)
  z[direction == "<="] <- pmax(r[direction == "<="], 0)
  z[direction == ">="] <- pmax(-r[direction == ">="], 0)
  if (!length(z)) 0 else max(z / pmax(1, abs(totals)))
}

#' Check a quadratic calibration optimality certificate
#' @inheritParams calibrate_weights
#' @param w Candidate final weights.
#' @param lambda Optional multipliers in the convention gradient = X lambda.
#' @param recover Search for alternative equality multipliers by LP when needed.
#' @return Primal and stationarity residuals, multiplier signs, complementarity,
#'   dual value and raw duality gap. optimal means certified within tol, not an
#'   exact-arithmetic proof. Error bounds require exact primal feasibility.
#' @details For equalities the KKT identity is w=clip(d+metric_scale*q*X lambda).
#'   Multipliers of <= constraints are nonpositive; those of >= are nonnegative.
#'   The gap bounds objective error only for exactly primal-feasible weights.
#'   Fixed coordinates impose no sign condition. Failed certification is not
#'   itself a proof of suboptimality. LP recovery handles nonunique multipliers.
#' @export
#' @examples
#' X <- matrix(1, 3, 1)
#' calibration_certificate(X, rep(1, 3), 6, rep(2, 3), lower=0, upper=3)
calibration_certificate <- function(X, d, totals, w, q = 1, lower = -Inf,
  upper = Inf, metric_scale = NULL, lambda = NULL, direction = "=",
  tol = 1e-8, recover = TRUE) {
  if (!is.logical(recover) || length(recover)!=1L || is.na(recover)) stop("recover must be TRUE or FALSE")
  dat <- .ct_data(X, totals); X <- dat$X; totals <- dat$totals
  n <- nrow(X); p <- ncol(X)
  d <- .ct_vec(d, n, "d", scalar = FALSE); w <- .ct_vec(w, n, "w", scalar = FALSE)
  v <- .ct_metric(d, .ct_vec(q, n, "q"), metric_scale)
  bu <- .ct_bounds(lower, upper, n); L <- bu$L; U <- bu$U
  direction <- .ct_direction(direction, p); .ct_scalar(tol, "tol")
  grad <- (w - d) / v
  fixed <- L == U
  low <- w <= L; high <- w >= U
  free <- !fixed & !low & !high
  if (is.null(lambda)) lambda <- .ct_pinv(X[free, , drop = FALSE], grad[free])
  lambda <- .ct_vec(lambda, p, "lambda", scalar = FALSE)
  stat <- function(la) {
    z <- grad - as.numeric(X %*% la)
    z[low] <- pmin(z[low], 0); z[high] <- pmax(z[high], 0); z[fixed] <- 0
    z
  }
  # Multipliers may be nonunique even after fitting all free coordinates.
  # Solve their linear feasibility problem, splitting lambda=lambda+ - lambda-.
  recovery <- "not_needed"
  if (recover && all(direction == "=") && p && max(abs(stat(lambda)), 0) > tol) {
    take <- which(!fixed)
    dirs <- ifelse(low[take], "<=", ifelse(high[take], ">=", "="))
    if (length(take)) {
      z <- lpSolve::lp("min", rep(0, 2*p), cbind(X[take,,drop=FALSE], -X[take,,drop=FALSE]),
                       dirs, grad[take])
      recovery <- paste0("lp_status_", z$status)
      if (z$status == 0L) lambda <- z$solution[seq_len(p)] - z$solution[p+seq_len(p)]
    }
  }
  eta <- as.numeric(X %*% lambda)
  z <- .ct_clip(d + v * eta, L, U)
  objective <- sum((w-d)^2/(2*v))
  dual <- sum(lambda*totals) + sum((z-d)^2/(2*v) - eta*z)
  r <- as.numeric(crossprod(X,w))-totals
  primal <- max(.ct_violation(X,w,totals,direction),
    pmax(L-w,w-U,0)/pmax(1,abs(w)))
  e <- stat(lambda)
  stationarity <- max(abs((w-z)/v)/pmax(1,abs(grad),abs(eta)), 0)
  sign_error <- max(c(0, pmax(lambda[direction=="<="],0), pmax(-lambda[direction==">="],0)))
  comp <- max(abs(lambda[direction!="="]*r[direction!="="]),0) / max(1,abs(objective))
  list(optimal = primal <= tol && stationarity <= tol && sign_error <= tol && comp <= tol,
    primal_feasible = primal <= tol, max_scaled_primal_violation = primal,
    max_scaled_stationarity = stationarity, multiplier_sign_error = sign_error,
    complementarity = comp, lambda = lambda, objective = objective, dual_value = dual,
    duality_gap = objective-dual, kkt_objective_bound = if(all(direction=="=")) sum(v*e^2)/2 else NA_real_,
    weighted_weight_error_bound = if(all(direction=="=")) sqrt(sum(v*e^2)) else NA_real_,
    bounds_require_exact_feasibility = TRUE, recovery = recovery, tolerance = tol)
}

# Irreversible affine projection on each selected face, always centred at d.
.ct_fixing <- function(X, totals, d, v, L, U, tol, max_iter) {
  w <- .ct_clip(d,L,U); fixed <- L==U
  la <- rep(0,ncol(X))
  for (it in seq_len(max_iter)) {
    free <- !fixed
    if (any(free)) {
      XF <- X[free,,drop=FALSE]
      rhs <- totals - as.numeric(crossprod(X[fixed,,drop=FALSE],w[fixed])) -
        as.numeric(crossprod(XF,d[free]))
      la <- .ct_pinv(crossprod(XF,XF*v[free]),rhs)
      w[free] <- d[free]+v[free]*as.numeric(XF%*%la)
    }
    if (.ct_error(X,w,totals)>tol)
      return(list(w=.ct_clip(w,L,U),lambda=la,iterations=it,converged=FALSE,status="inconsistent_selected_face"))
    bad <- !fixed & (w<L | w>U)
    if (!any(bad)) return(list(w=w,lambda=la,iterations=it,converged=TRUE,status="feasible"))
    w <- .ct_clip(w,L,U); fixed <- fixed|bad
  }
  list(w=w,lambda=la,iterations=max_iter,converged=.ct_error(X,w,totals)<=tol,status="iteration_limit")
}
.ct_projection_fit <- function(X,d,totals,q,v,L,U,red,tol,max_iter,action) {
  a <- .ct_fixing(red$A,red$b,d,v,L,U,tol*0.1,max_iter)
  cert <- calibration_certificate(X,d,totals,a$w,q,L,U,v/q,
                                   as.numeric(red$map%*%a$lambda),tol=tol)
  if (!cert$optimal && action == "refine") {
    fit <- calibrate_weights(X,d,totals,q,"truncated",L,U,v/q,"auto",tol=tol,max_iter=max_iter)
    fit$projection <- list(status=a$status,weights=a$w,certificate=cert,iterations=a$iterations)
    fit$solver <- paste0("projection+",fit$solver)
    return(fit)
  }
  if (!cert$primal_feasible)
    stop("Projection failed on a selected face; this does not prove infeasibility. Use projection_action='refine'")
  if (!cert$optimal) {
    if (action=="error") stop("Projection optimality was not certified")
    warning("Feasible projection; optimality not certified. Inspect $certificate or use projection_action='refine'",call.=FALSE)
  }
  g <- rep(NA_real_,length(d)); g[d!=0] <- a$w[d!=0]/d[d!=0]
  structure(list(g=g,lambda=cert$lambda,weights=a$w,
    residuals=setNames(as.numeric(crossprod(X,a$w))-totals,colnames(X)),
    max_scaled_residual=.ct_error(X,a$w,totals),objective=cert$objective,converged=TRUE,
    status=if(cert$optimal) "optimal" else "feasible_uncertified",iterations=a$iterations,
    method="truncated",solver="projection",rank=red$rank,lower=L,upper=U,
    metric_scale=v/q,q=q,certificate=cert),class="caltools_fit")
}

# One mixed-constraint backend: the existing strictly convex QP solver.
.ct_mixed_fit <- function(X,d,totals,q,method,lower,upper,metric_scale,tol,rank_tol,direction) {
  dat <- .ct_data(X,totals); X <- dat$X; totals <- dat$totals; n <- nrow(X)
  .ct_scalar(tol,"tol"); .ct_scalar(rank_tol,"rank_tol")
  d <- .ct_vec(d,n,"d",scalar=FALSE); q <- .ct_vec(q,n,"q"); v <- .ct_metric(d,q,metric_scale)
  if (n>2000L) stop("Dense mixed QP is limited to 2000 observations")
  if(is.null(lower)) lower <- if(method=="truncated") 0 else -Inf
  if(is.null(upper)) upper <- if(method=="truncated") 10*d else Inf
  bu <- .ct_bounds(lower,upper,n); L <- bu$L; U <- bu$U
  eq <- direction=="="; fixed <- which(L==U); I <- diag(n)
  E <- cbind(X[,eq,drop=FALSE],I[,fixed,drop=FALSE])
  red <- .ct_reduce(E,c(totals[eq],L[fixed]),rank_tol,tol)
  if(!red$consistent) stop("Inconsistent equality constraints or fixed weights")
  ii <- which(!eq); sg <- ifelse(direction[ii]==">=",1,-1)
  il <- which(is.finite(L)&L!=U); iu <- which(is.finite(U)&L!=U)
  C <- cbind(red$A,sweep(X[,ii,drop=FALSE],2,sg,"*"),I[,il,drop=FALSE],-I[,iu,drop=FALSE])
  rhs <- c(red$b,sg*totals[ii],L[il],-U[iu])
  if(!ncol(C)) return(calibrate_weights(X,d,totals,q,"linear",metric_scale=metric_scale,tol=tol))
  z <- tryCatch(quadprog::solve.QP(diag(1/v,n),d/v,C,rhs,meq=red$rank),error=function(e)e)
  if(inherits(z,"error")) stop("Mixed QP failed: ",conditionMessage(z),
    ". Use calibration_feasibility(..., direction=...) for a separate feasibility diagnosis.")
  # quadprog can return equality multipliers after internally flipping their
  # columns. Recover their signed values from stationarity in our original basis.
  non_eq <- seq.int(red$rank+1L,ncol(C))
  rest <- if(red$rank<ncol(C)) as.numeric(C[,non_eq,drop=FALSE]%*%z$Lagrangian[non_eq]) else rep(0,n)
  eqla <- as.numeric(crossprod(red$A,(z$solution-d)/v-rest))
  la <- numeric(ncol(X)); mapped <- as.numeric(red$map%*%eqla)
  la[eq] <- mapped[seq_len(sum(eq))]
  la[ii] <- sg*z$Lagrangian[red$rank+seq_along(ii)]
  cert <- calibration_certificate(X,d,totals,z$solution,q,L,U,v/q,la,direction,tol,recover=FALSE)
  if(!cert$optimal) stop("Mixed QP result failed the independent primal/KKT checks at tol")
  g <- rep(NA_real_,n); g[d!=0] <- z$solution[d!=0]/d[d!=0]
  structure(list(g=g,lambda=la,weights=z$solution,
    residuals=setNames(as.numeric(crossprod(X,z$solution))-totals,colnames(X)),
    max_scaled_residual=cert$max_scaled_primal_violation,objective=cert$objective,
    converged=TRUE,status="optimal",iterations=z$iterations[1],method=method,
    solver="qp",rank=red$rank,lower=L,upper=U,metric_scale=v/q,q=q,
    direction=direction,certificate=cert),class="caltools_fit")
}

# Accelerated projected residual minimization; certificates in scaled control units.
.ct_pfc <- function(X,totals,L,U,tol,max_iter) {
  scale <- pmax(1,abs(totals)); Z <- sweep(X,2,scale,"/"); b <- totals/scale
  lip <- max(svd(Z,nu=0,nv=0)$d,0)^2
  w <- L/2+U/2; y <- w; alpha <- 1
  best_upper <- Inf; best_lower <- 0; best <- w; certificate <- NULL
  result <- function(feasible,status,it) list(feasible=feasible,status=status,weights=best,
    max_scaled_residual=.ct_error(X,best,totals),solver="pfc",iterations=it,
    certificate=certificate,distance_lower=best_lower,distance_upper=best_upper,
    distance_units="Euclidean norm of totals scaled by pmax(1,abs(totals))")
  if(lip==0) return(result(if(max(abs(b),0)<=tol) TRUE else FALSE,"zero_matrix",0L))
  for(it in seq_len(max_iter)) {
    wn <- .ct_clip(y-as.numeric(Z%*%(as.numeric(crossprod(Z,y))-b))/lip,L,U)
    residual <- b-as.numeric(crossprod(Z,wn)); normr <- sqrt(sum(residual^2))
    if(normr<best_upper) {best_upper <- normr; best <- wn}
    if(max(abs(residual),0)<=tol) return(result(TRUE,"feasible",it))
    # Normalize the direction: avoids small squared residuals at the boundary.
    vv <- residual/normr; av <- as.numeric(Z%*%vv)
    support <- sum(pmax(L*av,U*av)); value <- sum(vv*b)-support
    # Conservative floating-point guard; not interval-arithmetic certification.
    roundoff <- 64*.Machine$double.eps*(nrow(X)+ncol(X)+1)*
      max(1,sum(abs(vv*b))+sum(pmax(abs(L*av),abs(U*av))))
    best_lower <- max(best_lower,0,value-roundoff)
    if(value>roundoff + tol*sum(abs(vv))) {
      certificate <- list(direction=vv/scale,scaled_direction=vv,
        separation=value,roundoff_guard=roundoff,tolerance_margin=tol*sum(abs(vv)),
        target=sum(vv*b),support=support,verified_in_exact_arithmetic=FALSE)
      return(result(FALSE,"infeasible_certified",it))
    }
    an <- (1+sqrt(1+4*alpha^2))/2
    y <- wn+(alpha-1)/an*(wn-w); w <- wn; alpha <- an
  }
  result(NA,"iteration_limit",max_iter)
}

# Standalone mixed feasibility is a shifted LP with checked primal witness.
.ct_mixed_feasibility <- function(X,totals,lower,upper,direction,tol) {
  n <- nrow(X); p <- ncol(X); bu <- .ct_bounds(lower,upper,n)
  L <- bu$L; U <- bu$U
  # Triplets avoid an n-by-n identity. Split unrestricted weights only when needed.
  finite_shift <- all(is.finite(L))
  shift <- if(finite_shift) L else rep(0,n)
  scale <- pmax(1,abs(totals))
  Z <- sweep(X,2,scale,"/")
  ix <- which(Z!=0,arr.ind=TRUE)
  trip <- cbind(ix[,2L],ix[,1L],Z[ix])
  if(!finite_shift) trip <- rbind(trip,cbind(ix[,2L],n+ix[,1L],-Z[ix]))
  rhs <- totals/scale-as.numeric(crossprod(Z,shift)); dirs <- direction
  il <- if(finite_shift) integer() else which(is.finite(L)); iu <- which(is.finite(U))
  ids <- c(il,iu); brhs <- c(L[il],U[iu]-shift[iu])
  if(length(ids)) {
    bs <- pmax(1,abs(brhs)); rows <- p+seq_along(ids)
    trip <- rbind(trip,cbind(rows,ids,1/bs))
    if(!finite_shift) trip <- rbind(trip,cbind(rows,n+ids,-1/bs))
    rhs <- c(rhs,brhs/bs); dirs <- c(dirs,rep(">=",length(il)),rep("<=",length(iu)))
  }
  nv <- if(finite_shift) n else 2*n
  # lpSolve needs at least one coefficient; zero rows are decided explicitly.
  nonzero <- unique(trip[,1L]); empty <- setdiff(seq_along(rhs),nonzero)
  if(length(empty)) {
    ok <- (dirs[empty]=="=" & abs(rhs[empty])<=tol) |
          (dirs[empty]=="<=" & rhs[empty]>=-tol) |
          (dirs[empty]==">=" & rhs[empty]<=tol)
    if(!all(ok)) return(list(feasible=FALSE,status="infeasible",weights=NULL,solver="zero-row"))
  }
  if(!nrow(trip)) return(list(feasible=TRUE,status="feasible",weights=shift,solver="lp",max_scaled_residual=0))
  z <- lpSolve::lp("min",rep(0,nv),const.dir=dirs,const.rhs=rhs,dense.const=trip)
  if(z$status!=0L) return(list(feasible=if(z$status==2L) FALSE else NA,
    status=if(z$status==2L) "infeasible" else paste0("lp_status_",z$status),weights=NULL,solver="lp"))
  w <- z$solution[seq_len(n)]+shift
  if(!finite_shift) w <- w-z$solution[n+seq_len(n)]
  err <- max(.ct_violation(X,w,totals,direction),pmax(L-w,w-U,0)/pmax(1,abs(w)))
  list(feasible=if(err<=tol) TRUE else NA,status=if(err<=tol) "feasible" else "numerical_failure",
       weights=w,solver="lp",max_scaled_residual=err,direction=direction)
}
