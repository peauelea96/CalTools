#' Calibrate weights with explicit numerical diagnostics
#'
#' Solve calibration equations with linear, bounded quadratic, raking or logit
#' distances. Bounds refer to final weights, never implicitly to g-weights.
#' @param X Numeric matrix, observations in rows and margins in columns.
#' @param d Finite initial weights, one per row.
#' @param totals Known totals in column order; named totals must match exactly.
#' @param q Positive scalar or vector of distance factors.
#' @param method Calibration distance.
#' @param lower,upper Scalar or vector bounds on final weights.
#' @param metric_scale Positive quadratic scale. Default is d, requiring d>0.
#'   For signed or zero d specify a positive scale explicitly.
#' @param solver Quadratic solver: auto tries the low-dimensional dual then a
#'   dense QP fallback (at most 2000 observations); qp forces the reference QP.
#' @param tol Relative/absolute calibration tolerance using max(1,abs(totals)).
#' @param rank_tol Relative SVD rank tolerance after column scaling.
#' @param max_iter Iteration budget for each dual optimization phase.
#' @return A caltools_fit list. The first two entries are g and lambda for legacy
#'   indexing. Also contains weights, residuals, objective, status, iterations,
#'   rank, solver and bounds. g is NA where d is zero; use weights instead.
#' @details The quadratic objective is sum((w-d)^2/(2*metric_scale*q)).
#'   Raking and logit require positive d and do not accept metric_scale.
#'   Nonconvergence is an error, not an infeasibility claim. See the mathematical
#'   notes distributed with the package for dual and KKT derivations.
#' @export
#' @examples
#' X <- cbind(N = rep(1, 4), x = 1:4)
#' fit <- calibrate_weights(X, rep(2, 4), c(N = 8, x = 22))
#' fit$weights
#' crossprod(X, fit$weights)
#' calibrate_weights(X, c(-1, 0, 2, 3), c(N = 8, x = 22),
#'   method = "truncated", lower = 0, upper = 5, metric_scale = 1)
calibrate_weights <- function(X, d, totals, q = 1,
  method = c("linear", "truncated", "raking", "logit"),
  lower = NULL, upper = NULL, metric_scale = NULL,
  solver = c("auto", "dual", "qp"), tol = 1e-8, rank_tol = 1e-10,
  max_iter = 500L) {
  method <- match.arg(method); solver <- match.arg(solver)
  dat <- .ct_data(X, totals); X <- dat$X; totals <- dat$totals
  n <- nrow(X)
  d <- .ct_vec(d, n, "d", scalar = FALSE)
  q <- .ct_vec(q, n, "q")
  .ct_scalar(tol, "tol"); .ct_scalar(rank_tol, "rank_tol")
  .ct_scalar(max_iter, "max_iter", integer = TRUE)
  if (method %in% c("linear", "raking") && (!is.null(lower) || !is.null(upper)))
    stop("Bounds are only used by truncated and logit methods")
  if (method %in% c("raking", "logit") && (!is.null(metric_scale) || any(d <= 0)))
    stop("Raking/logit require strictly positive d and metric_scale=NULL")
  v <- .ct_metric(d, q, metric_scale)
  if (method %in% c("truncated", "logit")) {
    if (is.null(lower)) lower <- 0
    if (is.null(upper)) {
      if (any(d <= 0)) stop("Specify upper explicitly when d has nonpositive entries")
      upper <- 10 * d
    }
    bu <- .ct_bounds(lower, upper, n)
  } else bu <- list(L = rep(-Inf, n), U = rep(Inf, n))
  L <- bu$L; U <- bu$U
  if (method == "logit" && (any(!is.finite(L)) || any(!is.finite(U)) ||
                           any(L >= d) || any(d >= U)))
    stop("Logit requires finite vector bounds lower < d < upper")
  red <- .ct_reduce(X, totals, rank_tol, tol)
  if (!red$consistent) stop("Inconsistent calibration totals: linear dependencies are violated")
  used <- method
  if (method == "linear") {
    la <- .ct_pinv(crossprod(red$A, red$A * v), red$b - as.numeric(crossprod(red$A, d)))
    ans <- list(w = d + v * as.numeric(red$A %*% la), lambda = la,
                iterations = 1L, converged = TRUE)
  } else if (method == "truncated" && solver == "qp" && red$rank > 0L) {
    used <- "qp"
    ans <- .ct_qp(red$A, red$b, d, v, L, U)
  } else {
    used <- "dual"
    ans <- .ct_dual(red$A, red$b, d, v, method, L, U, tol * 0.1, max_iter)
    good <- ans$converged && .ct_error(X, ans$w, totals) <= tol
    if (!good && method == "truncated" && solver == "auto" && n <= 2000L) {
      qa <- tryCatch(.ct_qp(red$A, red$b, d, v, L, U), error = function(e) NULL)
      if (!is.null(qa)) { ans <- qa; used <- "qp-fallback" }
    }
  }
  err <- .ct_error(X, ans$w, totals)
  bounded <- all(ans$w >= L - tol * pmax(1, abs(ans$w))) &&
             all(ans$w <= U + tol * pmax(1, abs(ans$w)))
  if (!ans$converged || !is.finite(err) || err > tol || !bounded) {
    if (method == "truncated" && all(is.finite(L)) && all(is.finite(U))) {
      feas <- calibration_feasibility(X, totals, L, U, tol = tol, method = "lp")
      if (identical(feas$feasible, FALSE)) stop("Calibration problem is infeasible within the specified bounds")
    }
    stop("Calibration did not converge to the requested tolerance; this does not prove infeasibility")
  }
  w <- ans$w
  g <- rep(NA_real_, n); g[d != 0] <- w[d != 0] / d[d != 0]
  lambda <- as.numeric(red$map %*% ans$lambda)
  names(lambda) <- colnames(X)
  objective <- if (method %in% c("linear", "truncated")) sum((w - d)^2 / (2 * v)) else
    calibration_distance(w, d, q, method, lower = L, upper = U)
  structure(list(g = g, lambda = lambda, weights = w,
    residuals = setNames(as.numeric(crossprod(X, w)) - totals, colnames(X)),
    max_scaled_residual = err, objective = objective, converged = TRUE,
    status = "optimal", iterations = ans$iterations, method = method,
    solver = used, rank = red$rank, lower = L, upper = U,
    metric_scale = v / q, q = q), class = "caltools_fit")
}

#' Evaluate the calibration distance in final weights
#' @inheritParams calibrate_weights
#' @param w Final weights.
#' @return Nonnegative objective value, or Inf outside the distance domain.
#' @export
#' @examples
#' calibration_distance(c(2, 3), c(1, 2))
#' calibration_distance(c(0, 2), c(1, 1), method = "logit", lower = 0, upper = 2)
calibration_distance <- function(w, d, q = 1,
  method = c("linear", "truncated", "raking", "logit"),
  lower = NULL, upper = NULL, metric_scale = NULL) {
  method <- match.arg(method)
  d <- .ct_vec(d, length(w), "d", scalar = FALSE)
  w <- .ct_vec(w, length(d), "w", scalar = FALSE)
  q <- .ct_vec(q, length(d), "q")
  if (any(q <= 0)) stop("q must be positive")
  if (method %in% c("truncated", "logit")) {
    if (is.null(lower) || is.null(upper)) stop("Supply both bounds")
    bu <- .ct_bounds(lower, upper, length(d))
    if (any(w < bu$L) || any(w > bu$U)) return(Inf)
  }
  if (method %in% c("linear", "truncated"))
    return(sum((w - d)^2 / (2 * .ct_metric(d, q, metric_scale))))
  if (any(d <= 0) || !is.null(metric_scale)) stop("Raking/logit require positive d and no metric_scale")
  g <- w / d
  xlog <- function(x, y) { z <- numeric(length(x)); k <- x > 0; z[k] <- x[k] * log(x[k] / y[k]); z }
  if (method == "raking") {
    if (any(w < 0)) return(Inf)
    return(sum(d / q * (xlog(g, rep(1, length(g))) - g + 1)))
  }
  l <- bu$L / d; u <- bu$U / d
  if (any(!is.finite(l)) || any(!is.finite(u)) || any(l >= 1) || any(u <= 1))
    stop("Logit requires finite lower < d < upper")
  a <- (u - l) / ((1 - l) * (u - 1))
  sum(d / (q * a) * (xlog(g - l, 1 - l) + xlog(u - g, u - 1)))
}

#' Short calibration interface
#' @inheritParams calibrate_weights
#' @param description Print a short summary; no graphical side effects.
#' @return A caltools_fit; [[1]] is g and [[2]] is lambda.
#' @export
#' @examples
#' X <- matrix(1, 3, 1)
#' fit <- calibrate(X, rep(1, 3), 6)
#' fit$weights
calibrate <- function(X, d, totals, q = rep(1, length(d)),
  method = c("linear", "raking", "truncated", "logit"),
  lower = NULL, upper = NULL, description = FALSE, max_iter = 500,
  metric_scale = NULL, tol = 1e-8, solver = "auto") {
  z <- calibrate_weights(X, d, totals, q, match.arg(method), lower,
    upper, metric_scale, solver, tol = tol, max_iter = max_iter)
  if (description) print(z)
  z
}

#' Evaluate the calibration distance from multiplicative factors
#' @inheritParams calibrate_weights
#' @param g Multiplicative factors (only meaningful for nonzero d).
#' @return The calibration distance.
#' @export
#' @examples
#' calibration_distance_from_factors(c(1, 2), c(2, 3), method = "raking")
calibration_distance_from_factors <- function(g, d, q = 1,
  method = c("linear", "raking", "truncated", "logit"),
  lower = NULL, upper = NULL, metric_scale = NULL) {
  g <- .ct_vec(g, length(d), "g", scalar = FALSE)
  if (any(d == 0)) stop("Use calibration_distance(w, d, ...) when d contains zeros")
  calibration_distance(g * d, d, q, match.arg(method), lower, upper, metric_scale)
}

#' @export
print.caltools_fit <- function(x, ...) {
  cat("CalTools:", x$method, "calibration (", x$solver, ")\n", sep = " ")
  cat("Status:", x$status, " | rank:", x$rank,
      " | max scaled residual:", format(x$max_scaled_residual, digits = 3), "\n")
  cat("Weight range:", paste(format(range(x$weights), digits = 5), collapse = " to "), "\n")
  invisible(x)
}
