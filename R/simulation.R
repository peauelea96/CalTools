#' Generate calibration stress-test problems with controlled difficulty
#'
#' @param n Number of sampled units (rows).
#' @param categorical_groups Numbers of categories in each one-hot family;
#'   e.g. c(5, 4) creates two complete partitions with 5 and 4 margins.
#' @param n_binary,n_continuous,n_positive,n_zero,n_derived Numbers of Bernoulli,
#'   signed normal, positive continuous, identically zero, and exactly linearly
#'   dependent columns, respectively. Derived columns are sums of base columns.
#' @param include_intercept Add a constant-one N column.
#' @param base_weight Mean starting weight (>=1).
#' @param weight_cv Approximate coefficient of variation of d (0 gives constant d).
#' @param bounds Pair of positive g-weight factors with lower <= 1 < upper.
#' @param difficulty Real number in [0,1]. Up to 0.5 a known weight vector
#'   satisfies all starting bounds; beyond 0.5 a selected margin has a total
#'   above its largest possible value under the starting bounds. The distance
#'   beyond that upper limit increases with difficulty.
#' @param stress_pattern "local" stresses the smallest nonempty category of the
#'   first one-hot family (useful for margin aggregation); "global" stresses
#'   every unit (useful for testing bound relaxation). Without categories,
#'   "local" falls back to "global".
#' @param max_excess Extra upper g-factor excursion at difficulty=1; positive.
#' @param category_skew Value in [0,1] determining imbalance in category
#'   probabilities (0 is uniform).
#' @param ordered_groups Logical scalar or vector, specifying whether each
#'   categorical group may be merged only between adjacent categories.
#' @param inconsistency Either "none", "linear", or "zero". The last two
#'   deliberately make the totals algebraically inconsistent after generating
#'   the underlying consistent system. Required extra columns are added.
#' @param inconsistency_size Relative size of the deliberately introduced
#'   total error (strictly greater than 1e-6).
#' @param seed Optional integer seed. A supplied seed does not modify the
#'   caller's random-number generator state.
#' @return List containing X, totals, d, q, lower, upper, bounds, margin_groups,
#'   ordered, true_weights, difficulty, stress_pattern, stress_margin,
#'   expected_status, and legacy-style aliases Xs, total, L_bounds,
#'   U_bounds, mar_2_join, hierarchy. The reference true_weights satisfy all
#'   equations before an inconsistency is deliberately introduced.
#' @details This is a constructive benchmark, NOT a random population model.
#'   difficulty <= 0.5 guarantees bounded feasibility when inconsistency="none".
#'   difficulty > 0.5 guarantees bounded infeasibility when
#'   inconsistency="none" by exceeding an explicit nonnegative margin's box
#'   support. The exact amount of bound relaxation required depends on X and d.
#'   An algebraically inconsistent problem can never be repaired by bounds alone.
#'   Input factor bounds are converted to final-weight bounds via d.
#' @export
#' @examples
#' a <- generate_calibration_problem(n = 30, categorical_groups = c(5, 4),
#'   n_binary = 2, n_continuous = 1, n_zero = 1,
#'   difficulty = 0.8, seed = 42)
#' analyze_inconsistency(a$X, a$totals)
#' calibration_preflight(a$X, a$totals, a$lower, a$upper)$status
#' analyze_feasibility_bounds(a$X, a$totals, a$d, bounds = a$bounds,
#'   type = "g-weight", method = "upper", search_limits = 3, step = 0.1)$status
#' b <- generate_calibration_problem(n = 12, categorical_groups = c(3, 2),
#'   inconsistency = "linear", seed = 13)
#' analyze_inconsistency(b$X, b$totals)
generate_calibration_problem <- function(n = 30L, categorical_groups = c(5L, 4L),
    n_binary = 0L, n_continuous = 1L, n_positive = 0L,
    n_zero = 0L, n_derived = 0L, include_intercept = TRUE,
    base_weight = 5, weight_cv = 0, bounds = c(0.9, 1.1),
    difficulty = 0, stress_pattern = c("local", "global"),
    max_excess = 1, category_skew = 0.3, ordered_groups = FALSE,
    inconsistency = c("none", "linear", "zero"),
    inconsistency_size = 0.1, seed = NULL) {
  .ct_scalar(n, "n", integer = TRUE)
  for (arg in c("n_binary", "n_continuous", "n_positive", "n_zero", "n_derived"))
    .ct_scalar(get(arg), arg, integer = TRUE, zero = TRUE)
  if (!is.numeric(categorical_groups) || anyNA(categorical_groups) ||
      any(!is.finite(categorical_groups)) ||
      any(categorical_groups < 1) || any(categorical_groups != floor(categorical_groups)))
    stop("categorical_groups must be a vector of positive integers", call. = FALSE)
  if (!is.logical(include_intercept) || length(include_intercept) != 1L ||
      is.na(include_intercept)) stop("include_intercept must be TRUE or FALSE", call. = FALSE)
  .ct_scalar(base_weight, "base_weight")
  if (base_weight < 1) stop("base_weight must be >= 1", call. = FALSE)
  .ct_scalar(weight_cv, "weight_cv", zero = TRUE)
  bounds <- .ct_vec(bounds, 2L, "bounds", scalar = FALSE)
  if (!(bounds[1L] > 0 && bounds[1L] <= 1 && bounds[2L] > 1))
    stop("bounds must satisfy 0 < lower <= 1 < upper", call. = FALSE)
  .ct_scalar(difficulty, "difficulty", zero = TRUE)
  if (difficulty > 1) stop("difficulty must be in [0,1]", call. = FALSE)
  .ct_scalar(max_excess, "max_excess")
  .ct_scalar(category_skew, "category_skew", zero = TRUE)
  if (category_skew > 1) stop("category_skew must be in [0,1]", call. = FALSE)
  if (!is.logical(ordered_groups) || anyNA(ordered_groups) ||
      !(length(ordered_groups) == 1L ||
        length(ordered_groups) == length(categorical_groups)))
    stop("ordered_groups must be a logical scalar or one value per categorical group",
         call. = FALSE)
  ordered_groups <- rep_len(ordered_groups, length(categorical_groups))
  .ct_scalar(inconsistency_size, "inconsistency_size")
  if (inconsistency_size <= 1e-6)
    stop("inconsistency_size must exceed 1e-6", call. = FALSE)
  stress_pattern <- match.arg(stress_pattern)
  inconsistency <- match.arg(inconsistency)
  if (is.null(seed)) {
    # Respect the calling session's RNG state.
  } else {
    .ct_scalar(seed, "seed", integer = TRUE, zero = TRUE)
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) previous_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) assign(".Random.seed", previous_seed, envir = .GlobalEnv) else
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
          rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
  }
  add_columns <- function(A, B, labels) {
    B <- as.matrix(B)
    if (nrow(B) != n || ncol(B) != length(labels)) stop("internal column-size error")
    out <- cbind(A, B)
    colnames(out) <- c(colnames(A), labels)
    out
  }
  X <- matrix(numeric(), nrow = n, ncol = 0L)
  if (include_intercept) X <- add_columns(X, matrix(1, n, 1), "N")
  margin_groups <- vector("list", length(categorical_groups))
  memberships <- vector("list", length(categorical_groups))
  for (j in seq_along(categorical_groups)) {
    k <- as.integer(categorical_groups[j])
    prob <- exp(-category_skew * log(10) * (seq_len(k) - 1L) / max(1, k - 1L))
    membership <- sample.int(k, n, replace = TRUE, prob = prob)
    labels <- paste0("CAT", j, "_", seq_len(k))
    indicators <- outer(membership, seq_len(k), "==") * 1
    X <- add_columns(X, indicators, labels)
    margin_groups[[j]] <- labels
    memberships[[j]] <- membership
  }
  for (j in seq_len(n_binary))
    X <- add_columns(X, matrix(stats::rbinom(n, 1, stats::runif(1, 0.15, 0.85)), ncol = 1L),
                     paste0("BIN", j))
  for (j in seq_len(n_continuous))
    X <- add_columns(X, matrix(stats::rnorm(n), ncol = 1L), paste0("SIGNED", j))
  for (j in seq_len(n_positive))
    X <- add_columns(X, matrix(stats::runif(n, 1, 10), ncol = 1L), paste0("POS", j))
  if (n_zero > 0L)
    X <- add_columns(X, matrix(0, n, n_zero), paste0("ZERO", seq_len(n_zero)))
  if (n_derived > 0L) {
    available <- which(colSums(abs(X)) > 0)
    if (!length(available)) stop("n_derived requires at least one nonzero base column")
    for (j in seq_len(n_derived)) {
      a <- sample(available, 1L)
      b <- sample(available, 1L)
      X <- add_columns(X, matrix(X[, a] + 0.5 * X[, b], ncol = 1L),
                       paste0("DERIVED", j))
    }
  }
  if (inconsistency == "linear") {
    available <- which(colSums(abs(X)) > 0)
    if (!length(available)) {
      X <- add_columns(X, matrix(1, n, 1), "CONFLICT_BASE")
      available <- ncol(X)
    }
    a <- available[1L]
    X <- add_columns(X, matrix(X[, a], ncol = 1L), "LINEAR_CONFLICT")
  }
  if (inconsistency == "zero" && !any(colSums(abs(X)) == 0))
    X <- add_columns(X, matrix(0, n, 1), "ZERO_CONFLICT")
  if (!ncol(X)) stop("Generate at least one calibration column", call. = FALSE)
  sigma <- sqrt(log(1 + weight_cv^2))
  d <- if (sigma == 0) rep(base_weight, n) else
    pmax(1, base_weight * exp(stats::rnorm(n, -sigma^2 / 2, sigma)))
  if (any(!is.finite(d)))
    stop("Generated initial weights are not finite; reduce weight_cv or base_weight",
         call. = FALSE)
  lower <- d * bounds[1L]
  upper <- d * bounds[2L]
  stress_margin <- NA_character_
  chosen <- rep(TRUE, n)
  if (stress_pattern == "local" && length(margin_groups)) {
    group <- margin_groups[[1L]]
    sizes <- colSums(X[, group, drop = FALSE] * d)
    eligible <- which(sizes > 0)
    idx <- eligible[which.min(sizes[eligible])]
    stress_margin <- group[idx]
    chosen <- X[, stress_margin] == 1
  } else {
    stress_pattern <- "global"
    # A nonzero nonnegative column certifies global box infeasibility when
    # every witness weight exceeds its upper bound. Signed-only columns do not.
    eligible <- which(colSums(X < 0) == 0 & colSums(X > 0) > 0)
    if (difficulty > 0.5 && !length(eligible))
      stop("Difficulty > 0.5 with global stress requires a nonzero nonnegative margin; add an intercept, categorical family or positive column",
           call. = FALSE)
    if (length(eligible)) stress_margin <- colnames(X)[eligible[1L]]
  }
  # For local stress, use a single nonempty category. For global stress,
  # every individual exceeds its own bound beyond the threshold.
  z <- stats::runif(n, 0.5, 1)
  g <- rep(1, n)
  if (difficulty <= 0.5) {
    g[chosen] <- 1 + 2 * difficulty * (bounds[2L] - 1) * z[chosen]
  } else {
    g[chosen] <- bounds[2L] + (2 * difficulty - 1) * max_excess * z[chosen]
  }
  true_weights <- d * g
  totals <- setNames(as.numeric(crossprod(X, true_weights)), colnames(X))
  if (any(!is.finite(c(lower, upper, true_weights, totals))))
    stop("Generated values overflowed; reduce max_excess, base_weight or weight_cv",
         call. = FALSE)
  if (inconsistency != "none") {
    conflict_name <- if (inconsistency == "linear") "LINEAR_CONFLICT" else
      colnames(X)[which(colSums(abs(X)) == 0)[1L]]
    totals[conflict_name] <- totals[conflict_name] +
      inconsistency_size * max(1, abs(totals[conflict_name]))
  }
  expected_status <- if (inconsistency != "none") "inconsistent_totals" else
    if (difficulty <= 0.5) "feasible" else "bound_infeasible"
  list(X = X, Xs = X, totals = totals, total = totals, d = d, q = rep(1, n),
       lower = lower, upper = upper, L_bounds = lower, U_bounds = upper,
       bounds = bounds, margin_groups = margin_groups,
       mar_2_join = margin_groups, ordered = ordered_groups,
       hierarchy = ordered_groups,
       true_weights = true_weights, difficulty = difficulty,
       stress_pattern = stress_pattern, stress_margin = stress_margin,
       inconsistency = inconsistency, expected_status = expected_status)
}
