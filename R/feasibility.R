#' Check feasibility of calibration equations and finite vector bounds
#' @inheritParams calibrate_weights
#' @param method auto uses the low-dimensional quadratic dual then LP;
#'   projection uses alternating affine/box projections; lp uses LP directly.
#' @param fallback If projection exhausts max_iter, use LP when TRUE; otherwise
#'   return feasible=NA and status="iteration_limit".
#' @return A list with feasible (TRUE, FALSE, or NA for solver failure), status,
#'   witness weights, residual, and solver. FALSE is reported only from a checked
#'   separating certificate or LP infeasibility, never an iteration limit.
#' @details The LP uses w=L+z with z>=0 and sparse constraint triplets. The
#'   equations are scaled by max(1,abs(totals)); the tolerance does not reverse
#'   inequalities for negative totals. All bounds must be finite here.
#' @export
#' @examples
#' X <- cbind(N = rep(1, 3), x = -1:1)
#' calibration_feasibility(X, c(N = 3, x = 0), 0, 2)
#' calibration_feasibility(X, c(N = 9, x = 0), 0, 2)
calibration_feasibility <- function(X, totals, lower, upper,
  tol = 1e-8, method = c("auto", "projection", "lp"), max_iter = 500L, fallback = TRUE, rank_tol = 1e-10) {
  method <- match.arg(method)
  .ct_scalar(rank_tol, "rank_tol")
  dat <- .ct_data(X, totals); X <- dat$X; totals <- dat$totals
  .ct_scalar(tol, "tol"); .ct_scalar(max_iter, "max_iter", integer = TRUE)
  n <- nrow(X); p <- ncol(X)
  bu <- .ct_bounds(lower, upper, n); L <- bu$L; U <- bu$U
  if (any(!is.finite(L)) || any(!is.finite(U))) stop("Feasibility requires finite bounds")
  finish <- function(feasible, status, w = NULL, solver, certificate = NULL) {
    list(feasible = feasible, status = status, weights = w,
      max_scaled_residual = if (is.null(w)) NA_real_ else .ct_error(X, w, totals),
      solver = solver, certificate = certificate)
  }
  if (!is.logical(fallback) || length(fallback) != 1L || is.na(fallback))
    stop("fallback must be TRUE or FALSE")
  reduction <- .ct_reduce(X, totals, rank_tol = rank_tol, tol = tol)
  if (!reduction$consistent) return(finish(FALSE, "inconsistent_totals", solver = "rank",
    certificate = analyze_inconsistency(X, totals, tol, details = TRUE, rank_tol = rank_tol)))
  center <- L / 2 + U / 2
  if (!p || .ct_error(X, center, totals) <= tol)
    return(finish(TRUE, "feasible", center, "midpoint"))
  # Exact support of the box gives a cheap, independently checkable certificate.
  Z <- sweep(X, 2L, pmax(1, abs(totals)), "/")
  tn <- totals / pmax(1, abs(totals))
  lower <- colSums(pmax(Z, 0) * L + pmin(Z, 0) * U)
  upper <- colSums(pmax(Z, 0) * U + pmin(Z, 0) * L)
  bad <- which(tn < lower - tol | tn > upper + tol)
  if (length(bad)) return(finish(FALSE, "infeasible", solver = "box-certificate",
    certificate = list(columns = bad, min = lower[bad], max = upper[bad], target = tn[bad])))
  if (method == "auto") {
    red <- reduction
    if (red$consistent) {
      ans <- .ct_dual(red$A, red$b, center, rep(1, n), "truncated", L, U,
                      tol * 0.1, max_iter)
      if (all(is.finite(ans$w)) && .ct_error(X, ans$w, totals) <= tol)
        return(finish(TRUE, "feasible", ans$w, "dual-witness"))
    }
  }
  if (method == "projection") {
    weights <- center
    for (iteration in seq_len(max_iter)) {
      weights <- weights + as.numeric(reduction$A %*%
        (reduction$b - as.numeric(crossprod(reduction$A, weights))))
      weights <- .ct_clip(weights, L, U)
      if (.ct_error(X, weights, totals) <= tol) {
        result <- finish(TRUE, "feasible", weights, "alternating-projections")
        result$iterations <- iteration
        return(result)
      }
    }
    if (!fallback) {
      result <- finish(NA, "iteration_limit", weights, "alternating-projections")
      result$iterations <- max_iter
      return(result)
    }
  }
  # Shift negative weights: lpSolve variables themselves are nonnegative.
  ix <- which(Z != 0, arr.ind = TRUE)
  trip <- cbind(ix[, 2L], ix[, 1L], Z[ix])
  trip <- rbind(trip, cbind(p + seq_len(n), seq_len(n), rep(1, n)))
  rhs <- c(tn - as.numeric(crossprod(Z, L)), U - L)
  z <- lpSolve::lp(direction = "min", objective.in = rep(0, n),
    const.dir = c(rep("=", p), rep("<=", n)), const.rhs = rhs,
    dense.const = trip)
  if (z$status == 2L) return(finish(FALSE, "infeasible", solver = "lp"))
  if (z$status != 0L) return(finish(NA, paste0("lp_status_", z$status), solver = "lp"))
  w <- z$solution + L
  good <- .ct_error(X, w, totals) <= tol &&
    all(w >= L - tol * pmax(1, abs(L))) && all(w <= U + tol * pmax(1, abs(U)))
  if (!good) return(finish(NA, "numerical_failure", w, "lp"))
  finish(TRUE, "feasible", w, "lp")
}

.ct_feasible_bool <- function(z) {
  if (is.na(z$feasible)) stop("Feasibility solver is inconclusive: ", z$status)
  z$feasible
}
#' Boolean feasibility interfaces
#' @inheritParams calibration_feasibility
#' @param tol Calibration tolerance (absolute for totals of magnitude <=1).
#' @return Logical feasibility; numerical uncertainty raises an error.
#' @details The permanent-clipping rule and dense dual Farkas LP from the
#'   prototype are not used. Fast and regular wrappers share the same validated
#'   low-dimensional solver; use method='projection' in calibration_feasibility
#'   to request ordinary alternating projections.
#' @export
#' @examples
#' X <- matrix(-1, 3, 1)
#' is_feasible_lp(X, -6, lower = 0, upper = 3)
#' is_feasible(X, -6, lower = 0, upper = 3)
#' is_feasible_fast(X, -6, lower = 0, upper = 3)
is_feasible_lp <- function(X, totals, tol = 1e-6, lower, upper) {
  .ct_feasible_bool(calibration_feasibility(X, totals, lower, upper, tol, "lp"))
}
#' @rdname is_feasible_lp
#' @export
is_feasible <- function(X, totals, tol = 1e-6, lower, upper) {
  .ct_feasible_bool(calibration_feasibility(X, totals, lower, upper, tol, "auto"))
}
#' @rdname is_feasible_lp
#' @export
is_feasible_fast <- function(X, totals, tol = 1e-6, lower, upper) {
  is_feasible(X, totals, tol, lower, upper)
}

#' Minimally expand finite vector bounds along a specified direction
#' @inheritParams calibration_feasibility
#' @param lower_rate,upper_rate Nonnegative rates: L(rho)=L-rho*lower_rate,
#'   U(rho)=U+rho*upper_rate. At least one rate must be positive.
#' @return A list with rho, expanded bounds, a witness and status. rho is the
#'   minimum common expansion along this direction, not a general nearest box.
#' @export
#' @examples
#' relax_bounds(matrix(1, 3, 1), 9, 0, 2, lower_rate = 0)
relax_bounds <- function(X, totals, lower, upper,
                         lower_rate = 1, upper_rate = 1, tol = 1e-8) {
  dat <- .ct_data(X, totals); X <- dat$X; totals <- dat$totals
  .ct_scalar(tol, "tol")
  n <- nrow(X); p <- ncol(X)
  bu <- .ct_bounds(lower, upper, n)
  if (any(!is.finite(c(bu$L, bu$U)))) stop("Use finite initial bounds")
  lr <- .ct_vec(lower_rate, n, "lower_rate"); ur <- .ct_vec(upper_rate, n, "upper_rate")
  if (any(lr < 0 | ur < 0) || !any(lr + ur > 0)) stop("Expansion rates must be nonnegative and not all zero")
  # z = w-L+rho*lr >= 0, z-rho*(lr+ur) <= U-L.
  Z <- sweep(X, 2L, pmax(1, abs(totals)), "/")
  ix <- which(Z != 0, arr.ind = TRUE)
  trip <- cbind(ix[, 2L], ix[, 1L], Z[ix])
  trip <- rbind(trip, cbind(seq_len(p), rep(n + 1L, p), -as.numeric(crossprod(Z, lr))),
    cbind(p + seq_len(n), seq_len(n), 1),
    cbind(p + seq_len(n), rep(n + 1L, n), -(lr + ur)))
  trip <- trip[trip[, 3L] != 0, , drop = FALSE]
  rhs <- c(totals / pmax(1, abs(totals)) - as.numeric(crossprod(Z, bu$L)), bu$U - bu$L)
  fit <- lpSolve::lp("min", c(rep(0, n), 1),
    const.dir = c(rep("=", p), rep("<=", n)), const.rhs = rhs, dense.const = trip)
  if (fit$status != 0L) return(list(status = if (fit$status == 2L) "infeasible" else "numerical_failure",
                                   rho = NA_real_, solver_status = fit$status))
  rho <- fit$solution[n + 1L]
  L <- bu$L - rho * lr; U <- bu$U + rho * ur
  w <- fit$solution[seq_len(n)] + L
  if (.ct_error(X, w, totals) > tol || any(w < L - tol * pmax(1, abs(L))) ||
      any(w > U + tol * pmax(1, abs(U))))
    return(list(status = "numerical_failure", rho = rho))
  list(status = "optimal", rho = rho, lower = L, upper = U, weights = w)
}
