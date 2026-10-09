# Internal validation and numerical primitives. No packages are attached.
.ct_scalar <- function(x, name, integer = FALSE, zero = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < 0 || (!zero && x == 0) || (integer && x != floor(x)))
    stop(name, " must be a ", if (zero) "nonnegative" else "positive",
         if (integer) " integer" else " finite number", call. = FALSE)
  x
}
.ct_vec <- function(x, n, name, finite = TRUE, scalar = TRUE) {
  if (!is.numeric(x) || !is.null(dim(x)) || anyNA(x) ||
      (finite && any(!is.finite(x))) ||
      !(length(x) == n || (scalar && length(x) == 1L)))
    stop(name, " must be numeric, of length ", n,
         if (scalar) " (or scalar)" else "", call. = FALSE)
  if (length(x) == 1L) rep(x, n) else x
}
.ct_data <- function(X, totals = NULL) {
  X <- as.matrix(X)
  if (!is.numeric(X) || length(dim(X)) != 2L || nrow(X) == 0L ||
      any(!is.finite(X))) stop("X must be a finite numeric matrix with rows")
  storage.mode(X) <- "double"
  nm <- colnames(X)
  if (!is.null(nm) && (anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm)))
    stop("X column names must be nonempty and unique")
  if (!is.null(totals)) {
    totals <- .ct_vec(totals, ncol(X), "totals", scalar = FALSE)
    if (!is.null(names(totals)) && !is.null(nm) && !identical(names(totals), nm))
      stop("Names and order of totals must match columns of X")
  }
  list(X = X, totals = totals)
}
.ct_bounds <- function(L, U, n) {
  L <- .ct_vec(L, n, "lower", finite = FALSE)
  U <- .ct_vec(U, n, "upper", finite = FALSE)
  if (any(L > U) || any(L == Inf) || any(U == -Inf))
    stop("Bounds must satisfy -Inf <= lower <= upper <= Inf")
  list(L = L, U = U)
}
.ct_metric <- function(d, q, metric_scale) {
  if (is.null(metric_scale)) {
    if (any(d <= 0)) stop("Nonpositive d requires an explicit positive metric_scale; e.g. rep(1, length(d))")
    metric_scale <- d
  }
  metric_scale <- .ct_vec(metric_scale, length(d), "metric_scale")
  if (any(q <= 0) || any(metric_scale <= 0) || any(!is.finite(q * metric_scale)))
    stop("q and metric_scale must be strictly positive with finite product")
  metric_scale * q
}
.ct_error <- function(X, w, totals) {
  r <- as.numeric(crossprod(X, w)) - totals
  if (!length(r)) 0 else max(abs(r) / pmax(1, abs(totals)))
}
.ct_pinv <- function(A, b, tol = 1e-12) {
  if (!ncol(A)) return(numeric())
  ss <- svd(A)
  keep <- ss$d > tol * max(ss$d, 0)
  if (!any(keep)) return(rep(0, ncol(A)))
  as.numeric(ss$v[, keep, drop = FALSE] %*%
    (crossprod(ss$u[, keep, drop = FALSE], b) / ss$d[keep]))
}
# Whiten the constraint space once; return the map to original multipliers.
.ct_reduce <- function(X, totals, rank_tol = 1e-10, tol = 1e-8) {
  p <- ncol(X)
  if (!p) return(list(A = X, b = numeric(), map = matrix(0, 0, 0), rank = 0L,
                      consistent = TRUE, projected_totals = numeric(),
                      column_scale = numeric(), right_basis = matrix(0, 0, 0)))
  scale <- apply(abs(X), 2L, max)
  scale[scale == 0] <- 1
  Z <- sweep(X, 2L, scale, "/")
  s <- svd(Z, nu = min(dim(Z)), nv = min(dim(Z)))
  keep <- s$d > rank_tol * max(s$d, 0)
  V <- s$v[, keep, drop = FALSE]
  tn <- totals / scale
  recon <- as.numeric(V %*% crossprod(V, tn)) * scale
  consistent <- max(abs(totals - recon) / pmax(1, abs(totals))) <= tol
  map <- sweep(V, 1L, scale, "/")
  if (any(keep)) map <- sweep(map, 2L, s$d[keep], "/")
  list(A = s$u[, keep, drop = FALSE], b = as.numeric(crossprod(map, totals)),
       map = map, rank = sum(keep), consistent = consistent,
       projected_totals = recon, column_scale = scale, right_basis = V)
}
.ct_clip <- function(w, L, U) pmax(L, pmin(U, w))

# Minimize the negative calibration dual. Bound membership is recomputed at
# every evaluation: unlike permanent clipping, active constraints can leave.
.ct_dual <- function(A, b, d, v, method, L, U, tol, max_iter) {
  r <- ncol(A)
  if (method == "logit") {
    width <- U - L
    prob <- (d - L) / width
    offset <- stats::qlogis(prob)
    a <- v / (width * prob * (1 - prob))
    softplus <- function(z) pmax(z, 0) + log1p(exp(-abs(z)))
  }
  evaluate <- function(lambda) {
    eta <- as.numeric(A %*% lambda)
    if (method == "truncated") {
      raw <- d + v * eta
      w <- .ct_clip(raw, L, U)
      value <- sum(eta * w - (w - d)^2 / (2 * v)) - sum(b * lambda)
      h <- v * (raw > L & raw < U)
    } else if (method == "raking") {
      z <- v / d * eta
      w <- d * exp(z)
      value <- sum(d^2 / v * expm1(z)) - sum(b * lambda)
      h <- v / d * w
    } else {
      z <- offset + a * eta
      pr <- stats::plogis(z)
      w <- L + width * pr
      value <- sum(L * eta + width / a * (softplus(z) - softplus(offset))) - sum(b * lambda)
      h <- width * a * pr * (1 - pr)
    }
    list(w = w, value = value, grad = as.numeric(crossprod(A, w)) - b, h = h)
  }
  lambda <- rep(0, r)
  e <- evaluate(lambda)
  if (!r) return(list(w = e$w, lambda = lambda, iterations = 0L, converged = TRUE))
  done <- function(e) all(is.finite(e$grad)) && max(abs(e$grad) / pmax(1, abs(b))) <= tol
  it <- 0L
  for (it in seq_len(max_iter)) {
    if (done(e)) break
    H <- crossprod(A, A * e$h)
    step <- .ct_pinv(H, e$grad)
    if (!all(is.finite(step)) || sum(step * e$grad) <= 1e-14 * sum(e$grad^2))
      step <- e$grad / max(1, max(abs(H)))
    accepted <- FALSE
    for (j in 0:45) {
      trial <- lambda - 2^(-j) * step
      en <- evaluate(trial)
      if (is.finite(en$value) && en$value <= e$value - 1e-4 * 2^(-j) * sum(step * e$grad) +
          1e-13 * (1 + abs(e$value))) {
        accepted <- TRUE
        break
      }
    }
    if (!accepted) break
    lambda <- trial
    e <- en
  }
  if (!done(e)) {
    opt <- tryCatch(stats::optim(lambda, function(z) evaluate(z)$value,
      function(z) evaluate(z)$grad, method = "BFGS",
      control = list(maxit = max_iter, reltol = min(tol^2, 1e-14))), error = function(e) NULL)
    if (!is.null(opt)) {
      lambda <- opt$par
      e <- evaluate(lambda)
      it <- it + unname(opt$counts[1L])
    }
  }
  list(w = e$w, lambda = lambda, iterations = it, converged = done(e))
}

.ct_qp <- function(A, b, d, v, L, U) {
  n <- length(d)
  if (n > 2000L) stop("Dense QP fallback is limited to 2000 rows; use solver='dual' and increase max_iter")
  il <- which(is.finite(L)); iu <- which(is.finite(U))
  I <- diag(n)
  C <- cbind(A, I[, il, drop = FALSE], -I[, iu, drop = FALSE])
  z <- quadprog::solve.QP(diag(1 / v, n), d / v, C,
                          c(b, L[il], -U[iu]), meq = ncol(A))
  list(w = z$solution, lambda = z$Lagrangian[seq_len(ncol(A))],
       iterations = z$iterations[1L], converged = TRUE)
}
