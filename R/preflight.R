#' Diagnose a calibration problem before searching or fitting
#' @inheritParams calibrate_weights
#' @return A caltools_preflight list with status, consistency, feasibility and
#'   zero_sample. Status distinguishes inconsistent_totals (no weight vector can
#'   satisfy the equalities), bound_infeasible (equalities are consistent but the
#'   box is too restrictive), feasible, and numerical_failure.
#' @details Incompatible totals are diagnosed before bound optimization. No
#'   bound relaxation can repair inconsistent equalities. An explicit merge can
#'   change the equalities, but does not preserve the original separate totals.
#' @export
#' @examples
#' X <- cbind(A = c(1, 0), B = c(0, 1), N = c(1, 1))
#' calibration_preflight(X, c(A = 4, B = 6, N = 11), 0, 10)
#' calibration_preflight(X, c(A = 4, B = 6, N = 10), 0, 3)
#' calibration_preflight(X, c(A = 4, B = 6, N = 10), 0, 10)
calibration_preflight <- function(X, totals, lower, upper, tol = 1e-8,
                                   rank_tol = 1e-10) {
  data <- .ct_data(X, totals)
  bounds <- .ct_bounds(lower, upper, nrow(data$X))
  if (any(!is.finite(c(bounds$L, bounds$U)))) stop("Preflight requires finite bounds")
  consistency <- analyze_inconsistency(data$X, data$totals, tol,
                                       details = TRUE, rank_tol = rank_tol)
  zero_sample <- analyze_zero_sample(data$X, tol = 0)
  if (consistency$inconsistent) {
    feasibility <- list(feasible = FALSE, status = "inconsistent_totals",
      solver = "rank", certificate = consistency, weights = NULL)
    status <- "inconsistent_totals"
  } else {
    feasibility <- calibration_feasibility(data$X, data$totals,
                                            bounds$L, bounds$U, tol, rank_tol = rank_tol)
    status <- if (is.na(feasibility$feasible)) "numerical_failure" else
      if (feasibility$feasible) "feasible" else "bound_infeasible"
  }
  structure(list(status = status, consistency = consistency,
    feasibility = feasibility, zero_sample = zero_sample), class = "caltools_preflight")
}

#' @export
print.caltools_preflight <- function(x, ...) {
  cat("CalTools preflight:", x$status, "\n")
  print(x$consistency)
  if (x$status == "bound_infeasible")
    cat("Equalities are consistent, but no weights satisfy the supplied bounds.\n")
  if (x$status == "numerical_failure")
    cat("The bounded feasibility check was inconclusive.\n")
  invisible(x)
}

#' @export
print.caltools_search <- function(x, ...) {
  cat("CalTools search:", length(x), "candidate(s);", attr(x, "search_status"), "\n")
  print(attr(x, "preflight"))
  if (length(x)) for (i in seq_along(x)) {
    cat("Candidate", i, "| merges:", length(x[[i]]$history),
        "| objective:", format(x[[i]]$fit$objective, digits = 5), "\n")
  }
  invisible(x)
}
