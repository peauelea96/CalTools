#' Find the closest algebraically consistent calibration totals
#'
#' @param X Numeric sample matrix with one row per unit and one column per margin.
#' @param totals Numeric target vector, in exactly the order of the columns of X.
#' @param scales Positive scalar or vector of length ncol(X). The objective is
#'   sum(((adjusted - totals) / scales)^2) / 2. With scales=1 the projection
#'   is Euclidean in original total units. Larger scales allow larger changes.
#' @param details Return a structured report instead of just the new total vector.
#' @param rank_tol Relative numerical rank tolerance after the same column
#'   scaling used by analyze_inconsistency().
#' @return A named numeric vector of consistent totals by default. With
#'   details=TRUE, a caltools_total_projection object with totals, adjusted,
#'   adjustment, residuals, distance, objective, rank_X, and initially_consistent.
#' @details This projects the supplied totals onto the (numerically determined)
#'   row space of X. In particular, it does NOT impose positivity or bounds on
#'   calibration weights. Call calibration_feasibility() separately for those.
#'   Unlike analyze_inconsistency()$projected_totals, the default projection
#'   minimizes unscaled Euclidean distance. Algebraically, if S=diag(scales),
#'   the solution is S*Q*Q'*solve(S,totals), where Q spans S^{-1} row(X).
#'   Relative column scaling is used only to determine numerical rank; it
#'   does not change the distance objective.
#' @export
#' @examples
#' X <- cbind(A = c(1, 0), B = c(0, 1), N = c(1, 1))
#' totals <- c(A = 4, B = 6, N = 11)
#' closest_consistent_totals(X, totals)
#' closest_consistent_totals(X, totals, scales = c(1, 1, 10), details = TRUE)
closest_consistent_totals <- function(X, totals, scales = 1,
                                      details = FALSE, rank_tol = 1e-10) {
  dat <- .ct_data(X, totals)
  X <- dat$X
  totals <- dat$totals
  .ct_scalar(rank_tol, "rank_tol")
  if (!is.logical(details) || length(details) != 1L || is.na(details))
    stop("details must be TRUE or FALSE", call. = FALSE)
  p <- ncol(X)
  if (!p) {
    # A zero-column matrix has one and only one possible total vector: numeric().
    scales <- .ct_vec(scales, 0L, "scales")
    adjusted <- numeric()
    red <- .ct_reduce(X, totals, rank_tol = rank_tol)
  } else {
    scales <- .ct_vec(scales, p, "scales")
    if (any(scales <= 0)) stop("scales must be strictly positive", call. = FALSE)
    red <- .ct_reduce(X, totals, rank_tol = rank_tol)
    r <- red$rank
    if (r == p) {
      adjusted <- as.numeric(totals)
    } else if (r == 0L) {
      adjusted <- rep(0, p)
    } else {
      # diag(column_scale) * V spans row(X). Apply S^-1 and project
      # orthogonally in the standardized total coordinates.
      B <- sweep(red$right_basis, 1L, red$column_scale / scales, "*")
      Q <- qr.Q(qr(B, LAPACK = FALSE))
      standardized <- totals / scales
      adjusted <- as.numeric(Q %*% crossprod(Q, standardized)) * scales
    }
  }
  labels <- colnames(X)
  if (is.null(labels)) labels <- names(totals)
  if (is.null(labels)) labels <- paste0("margin_", seq_len(p))
  adjusted <- setNames(adjusted, labels)
  if (!details) return(adjusted)
  original <- setNames(as.numeric(totals), labels)
  change <- adjusted - original
  scaled_change <- if (p) change / scales else numeric()
  structure(list(totals = original, adjusted = adjusted,
                 adjustment = change, residuals = original - adjusted,
                 scales = setNames(scales, labels),
                 distance = sqrt(sum(scaled_change^2)),
                 objective = sum(scaled_change^2) / 2,
                 rank_X = red$rank, initially_consistent = red$consistent),
            class = "caltools_total_projection")
}

#' @export
print.caltools_total_projection <- function(x, ...) {
  cat("Closest consistent totals | rank:", x$rank_X,
      "| initially consistent:", x$initially_consistent, "\n")
  cat("Scaled Euclidean distance:", format(x$distance, digits = 7), "\n")
  print(data.frame(original = x$totals, adjusted = x$adjusted,
                   adjustment = x$adjustment), ...)
  invisible(x)
}
