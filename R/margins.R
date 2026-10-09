#' Diagnose and explain incompatible calibration totals
#' @inheritParams calibrate_weights
#' @param details Return a structured diagnostic rather than a logical flag.
#' @param tol Relative/absolute tolerance; NULL uses 1e-8.
#' @return With details=FALSE, TRUE if the totals are incompatible. Otherwise a
#'   caltools_inconsistency list with rank_X, rank_augmented, relations,
#'   coefficients, incompatibilities, projected_totals and residuals.
#'   Each coefficient row c describes X c = 0; consistency requires c' totals=0.
#'   The incompatibilities data frame reports the observed discrepancy, tolerance,
#'   involved margins and sample residual for each violated relation.
#' @details The rank comparison is performed by extending the row space of the
#'   column-scaled matrix with the totals row. rank_X uses rank_tol; the appended
#'   row increases the rank by one when its projection residual exceeds tol in
#'   original total units. This deliberately separates rank and total tolerances.
#'   Reported relations need not form a minimal or independent conflict set.
#'   They locate conflicting groups, not a uniquely responsible incorrect total.
#'   Projected totals are a diagnostic only: input totals are never overwritten.
#' @export
#' @examples
#' X <- cbind(group_A = c(1, 0), group_B = c(0, 1), N = c(1, 1))
#' totals <- c(group_A = 4, group_B = 6, N = 11)
#' analyze_inconsistency(X, totals)
#' report <- analyze_inconsistency(X, totals, details = TRUE)
#' report$incompatibilities
#' report$coefficients
analyze_inconsistency <- function(X, totals, tol = NULL, details = FALSE,
                                  rank_tol = 1e-10) {
  data <- .ct_data(X, totals); X <- data$X; totals <- data$totals
  if (is.null(tol)) tol <- 1e-8
  .ct_scalar(tol, "tol"); .ct_scalar(rank_tol, "rank_tol")
  if (!is.logical(details) || length(details) != 1L || is.na(details))
    stop("details must be TRUE or FALSE")
  reduction <- .ct_reduce(X, totals, rank_tol, tol)
  if (!details) return(!reduction$consistent)
  p <- ncol(X)
  margin_labels <- colnames(X)
  if (is.null(margin_labels)) margin_labels <- if (p) paste0("margin_", seq_len(p)) else character()
  coefficients <- matrix(numeric(), 0L, p, dimnames = list(NULL, margin_labels))
  relations <- data.frame(relation = character(), discrepancy = numeric(),
    tolerance = numeric(), sample_residual = numeric(), violated = logical(),
    stringsAsFactors = FALSE)
  relations$margins <- I(list())
  # Row j of the orthogonal null-space projector gives a relation whose total
  # discrepancy is precisely the jth component of the consistency residual.
  projector <- diag(p) - tcrossprod(reduction$right_basis)
  keys <- character()
  for (j in seq_len(p)) {
    if (sqrt(sum(projector[j, ]^2)) <= 100 * .Machine$double.eps) next
    raw <- projector[j, ] / reduction$column_scale
    magnitude <- max(abs(raw))
    if (magnitude == 0) next
    coefficient <- raw / magnitude
    nonzero <- abs(projector[j, ]) > 100 * .Machine$double.eps
    # Only use the mask for labels. Keep full numerical coefficients as evidence.
    # Fix orientation by the first substantial coefficient; choosing a maximum
    # is unstable when several coefficients tie up to floating-point error.
    orientation <- sign(coefficient[which(abs(coefficient) >
      100 * .Machine$double.eps)[1L]])
    coefficient <- coefficient * orientation
    discrepancy <- sum(coefficient * totals)
    threshold <- tol * max(1, abs(totals[j])) /
      (reduction$column_scale[j] * magnitude)
    key <- paste(signif(coefficient, 10), collapse = ":")
    previous <- match(key, keys)
    if (!is.na(previous)) {
      relations$tolerance[previous] <- min(relations$tolerance[previous], threshold)
      relations$violated[previous] <- abs(relations$discrepancy[previous]) > relations$tolerance[previous]
      next
    }
    selected <- which(nonzero)
    terms <- vapply(selected, function(k) {
      value <- signif(abs(coefficient[k]), 7)
      if (value == 1) margin_labels[k] else
        paste0(format(value, trim = TRUE), " * ", margin_labels[k])
    }, "")
    signs <- ifelse(coefficient[selected] < 0, " - ", " + ")
    signs[1L] <- if (coefficient[selected[1L]] < 0) "-" else ""
    label <- paste0(signs, terms, collapse = "")
    row <- data.frame(relation = paste0(label, " = 0"), discrepancy = discrepancy,
      tolerance = threshold, sample_residual = max(abs(as.numeric(X %*% coefficient))),
      violated = abs(discrepancy) > threshold, stringsAsFactors = FALSE)
    row$margins <- I(list(margin_labels[nonzero]))
    relations <- rbind(relations, row)
    coefficients <- rbind(coefficients, coefficient)
    keys <- c(keys, key)
  }
  rownames(coefficients) <- if (nrow(coefficients)) paste0("relation_", seq_len(nrow(coefficients))) else character()
  rownames(relations) <- rownames(coefficients)
  residuals <- setNames(totals - reduction$projected_totals, margin_labels)
  structure(list(inconsistent = !reduction$consistent, rank_X = reduction$rank,
    rank_augmented = reduction$rank + as.integer(!reduction$consistent),
    relations = relations, coefficients = coefficients,
    incompatibilities = relations[relations$violated, , drop = FALSE],
    projected_totals = setNames(reduction$projected_totals, margin_labels),
    residuals = residuals, max_scaled_residual = if (p) max(abs(residuals) / pmax(1, abs(totals))) else 0,
    tol = tol, rank_tol = rank_tol), class = "caltools_inconsistency")
}

#' @export
print.caltools_inconsistency <- function(x, ...) {
  cat("Calibration totals:", if (x$inconsistent) "incompatible" else "consistent", "\n")
  cat("Rank of X:", x$rank_X, " | rank after appending totals:", x$rank_augmented, "\n")
  if (nrow(x$incompatibilities))
    print(x$incompatibilities[, c("relation", "discrepancy", "tolerance"), drop = FALSE], row.names = FALSE)
  invisible(x)
}

#' Diagnose linear dependencies and zero columns
#' @inheritParams calibrate_weights
#' @param rank_tol Relative rank tolerance after scaling columns.
#' @param simple Rescale each dependency so its smallest nonzero coefficient is one in magnitude.
#' @return A list: multicollinear, relations (one row per dependence, in original
#'   units) and indices of redundant columns, plus rank. First three positions
#'   retain their legacy roles, but indices now mean column indices.
#' @export
#' @examples
#' X <- cbind(A = c(1, 0), B = c(0, 1), N = c(1, 1), empty = 0)
#' analyze_collinearity(X, simple = TRUE)
#' analyze_zero_sample(X)
analyze_collinearity <- function(X, rank_tol = 1e-7, simple = FALSE) {
  X <- .ct_data(X)$X; .ct_scalar(rank_tol, "rank_tol")
  p <- ncol(X)
  if (!p) return(list(multicollinear = FALSE, relations = NULL, indices = integer(), rank = 0L))
  sc <- apply(abs(X), 2, max); sc[sc == 0] <- 1
  Z <- sweep(X, 2, sc, "/")
  # QR supplies explicit relations without forming the squared-condition Gram matrix.
  qrZ <- qr(Z, tol = rank_tol, LAPACK = FALSE)
  rank <- qrZ$rank
  indep <- qrZ$pivot[seq_len(rank)]
  dep <- setdiff(seq_len(p), indep)
  if (!length(dep)) return(list(multicollinear = FALSE, relations = NULL, indices = integer(), rank = rank))
  rel <- matrix(0, length(dep), p, dimnames = list(NULL, colnames(X)))
  for (k in seq_along(dep)) {
    j <- dep[k]; rel[k, j] <- 1 / sc[j]
    if (rank) rel[k, indep] <- -.ct_pinv(Z[, indep, drop = FALSE], Z[, j], rank_tol) / sc[indep]
    if (simple) {
      nz <- abs(rel[k, ]) > rank_tol * max(abs(rel[k, ]))
      rel[k, !nz] <- 0
      rel[k, ] <- rel[k, ] / min(abs(rel[k, nz]))
    }
  }
  list(multicollinear = TRUE, relations = rel, indices = dep, rank = rank)
}

#' @rdname analyze_collinearity
#' @param tol Absolute threshold for identifying zero columns.
#' @export
analyze_zero_sample <- function(X, tol = 1e-6) {
  X <- .ct_data(X)$X; .ct_scalar(tol, "tol", zero = TRUE)
  ind <- which(colSums(abs(X) > tol) == 0)
  list(indices = ind, names = colnames(X)[ind])
}

.ct_named <- function(X, totals) {
  dat <- .ct_data(X, totals)
  if (is.null(colnames(dat$X)) || is.null(names(dat$totals)))
    stop("Name both X columns and totals")
  dat
}
.ct_names <- function(x, available) {
  if (!is.character(x) || anyNA(x) || anyDuplicated(x) || !all(x %in% available))
    stop("Select distinct existing column names")
}
#' Suppress or aggregate calibration margins
#' @inheritParams calibrate_weights
#' @param delete Names of columns to suppress.
#' @param join Distinct names of at least two columns to aggregate by summation.
#' @param merged_name Name for the resulting aggregate.
#' @return A named list X, totals, retaining positional access [[1]], [[2]].
#'   Zero-column matrices are supported. Aggregation weakens the original
#'   constraints: original individual totals are no longer guaranteed.
#' @export
#' @examples
#' X <- cbind(a = c(1, 1, 0), b = c(0, 0, 1))
#' merge_margins(X, c(a = 4, b = 2), c("a", "b"), "N")
#' drop_margins(X, c(a = 4, b = 2), "b")
drop_margins <- function(X, totals, delete) {
  dat <- .ct_named(X, totals); .ct_names(delete, colnames(dat$X))
  keep <- setdiff(colnames(dat$X), delete)
  list(X = dat$X[, keep, drop = FALSE], totals = dat$totals[keep])
}
#' @rdname drop_margins
#' @export
merge_margins <- function(X, totals, join, merged_name = NULL) {
  dat <- .ct_named(X, totals); .ct_names(join, colnames(dat$X))
  if (length(join) < 2L) stop("join must contain at least two distinct names")
  if (is.null(merged_name)) merged_name <- paste(join, collapse = "+")
  if (!is.character(merged_name) || length(merged_name) != 1L || is.na(merged_name) ||
      !nzchar(merged_name) || merged_name %in% setdiff(colnames(dat$X), join))
    stop("merged_name must be a unique nonempty name")
  z <- drop_margins(dat$X, dat$totals, join)
  X <- cbind(z$X, rowSums(dat$X[, join, drop = FALSE]))
  colnames(X)[ncol(X)] <- merged_name
  list(X = X, totals = c(z$totals, setNames(sum(dat$totals[join]), merged_name)))
}

#' Initial linear calibration multipliers
#' @inheritParams calibrate_weights
#' @return Named multiplier vector. Multipliers are not unique when margins are
#'   dependent and their magnitudes depend on column units.
#' @export
#' @examples
#' linear_multipliers(cbind(1, 1:4), rep(1, 4), c(4, 11))
linear_multipliers <- function(X, d, totals, q = rep(1, length(d)), metric_scale = NULL) {
  calibrate_weights(X, d, totals, q, "linear", metric_scale = metric_scale)$lambda
}

#' Immediate neighbors in an ordered vector
#' @param values Ordered vector.
#' @param position One-based position.
#' @return Its existing adjacent elements (empty for a singleton).
#' @export
#' @examples
#' adjacent_elements(c("young", "middle", "older"), 2)
adjacent_elements <- function(values, position) {
  .ct_scalar(position, "position", integer = TRUE)
  if (position > length(values)) stop("position is outside values")
  ix <- c(position - 1L, position + 1L)
  values[ix[ix >= 1L & ix <= length(values)]]
}

#' Generate random inclusion probabilities with prescribed sum
#' @param N Population size, a positive integer.
#' @param n Desired sample size, between zero and N; may be noninteger.
#' @return Vector in [0,1] summing to n. The distribution is not uniform on the
#'   fixed-sum polytope and differs from the old rejection-based generator.
#' @export
#' @examples
#' set.seed(1)
#' sum(random_inclusion_probabilities(100, 10))
random_inclusion_probabilities <- function(N, n) {
  .ct_scalar(N, "N", integer = TRUE); .ct_scalar(n, "n", zero = TRUE)
  if (n > N) stop("n must not exceed N")
  if (n == 0) return(rep(0, N))
  if (n == N) return(rep(1, N))
  size <- stats::rexp(N)
  left <- seq_len(N); out <- numeric(N); remain <- n
  repeat {
    v <- remain * size[left] / sum(size[left])
    certain <- v >= 1
    if (!any(certain)) { out[left] <- v; break }
    out[left[certain]] <- 1; remain <- remain - sum(certain)
    left <- left[!certain]
  }
  out
}
