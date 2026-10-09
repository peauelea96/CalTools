.ct_gbounds <- function(d, bounds, type) {
  if (type == "w-weight") return(list(L = rep(bounds[1L], length(d)), U = rep(bounds[2L], length(d))))
  if (any(d == 0)) stop("g-weight bounds are undefined for zero d; use w-weight bounds")
  list(L = pmin(d * bounds[1L], d * bounds[2L]), U = pmax(d * bounds[1L], d * bounds[2L]))
}
.ct_grid <- function(bounds, search_limits, step, method) {
  bounds <- .ct_vec(bounds, 2L, "bounds", scalar = FALSE)
  .ct_scalar(step, "step")
  if (bounds[1L] > bounds[2L]) stop("bounds must be ordered")
  limits <- if (method == "both") .ct_vec(search_limits, 2L, "search_limits", scalar = FALSE) else {
    search_limits <- .ct_vec(search_limits, 1L, "search_limits")
    if (method == "upper") c(bounds[1L], search_limits) else c(search_limits, bounds[2L])
  }
  if (limits[1L] > bounds[1L] || limits[2L] < bounds[2L]) stop("Search limits must expand the initial bounds")
  count <- max(ceiling((bounds[1L] - limits[1L]) / step), ceiling((limits[2L] - bounds[2L]) / step))
  if (count > 100000) stop("Grid exceeds 100000 steps; increase step or use relax_bounds")
  k <- 0:count
  cbind(pmax(limits[1L], bounds[1L] - k * step), pmin(limits[2L], bounds[2L] + k * step))
}

#' Search an expanding grid of common bounds
#' @inheritParams calibrate_weights
#' @param tol Feasibility tolerance.
#' @param type Whether common bounds multiply d (g-weight) or apply directly.
#' @param method Which endpoints to expand.
#' @param bounds Initial lower and upper common endpoints.
#' @param search_limits Outermost endpoint(s), scalar for upper/lower, pair for both.
#' @param step Positive grid spacing; the final endpoint is always included.
#' @return bounds, status (0 found, 1 exhausted), vector bounds and witness when
#'   successful. For both, searches a simultaneous expanding path, not all pairs.
#' @export
#' @examples
#' analyze_feasibility_bounds(matrix(1, 3, 1), 9, rep(1, 3),
#'   type = "w-weight", method = "upper", bounds = c(0, 2),
#'   search_limits = 4, step = 0.25)
analyze_feasibility_bounds <- function(X, totals, d, tol = 1e-6,
  type = c("g-weight", "w-weight"), method = c("upper", "lower", "both"),
  bounds, search_limits, step) {
  type <- match.arg(type); method <- match.arg(method)
  dat <- .ct_data(X, totals); d <- .ct_vec(d, nrow(dat$X), "d", scalar = FALSE)
  grid <- .ct_grid(bounds, search_limits, step, method)
  # Feasibility is monotone along this nested path: binary search grid index.
  check <- function(k) {
    bu <- .ct_gbounds(d, grid[k, ], type)
    f <- calibration_feasibility(dat$X, dat$totals, bu$L, bu$U, tol)
    list(ok = .ct_feasible_bool(f), bounds = grid[k, ], lower = bu$L,
         upper = bu$U, weights = f$weights)
  }
  hi <- nrow(grid); last <- check(hi)
  if (!last$ok) return(list(bounds = bounds, status = 1L))
  lo <- 1L
  while (lo < hi) {
    mid <- floor((lo + hi) / 2)
    if (check(mid)$ok) hi <- mid else lo <- mid + 1L
  }
  out <- check(lo); out$ok <- NULL; out$status <- 0L; out$grid_index <- lo
  out
}

#' Modify margins and search for feasible bounds
#' @inheritParams analyze_feasibility_bounds
#' @param lower,upper Actual starting vector bounds, checked first.
#' @param margin_names Selected margin names.
#' @param action Suppress or join the selected margins.
#' @param merged_name Optional new aggregate name.
#' @return List with X, totals, bounds, vector bounds and status.
#' @export
#' @examples
#' X <- cbind(A = c(1, 0), B = c(0, 1))
#' check_modified_feasibility(X, c(A = 0, B = 4), c(1, 1), 1, 3,
#'   c("A", "B"), action = "join", bounds = c(1, 3),
#'   search_limits = 3, step = 1)
check_modified_feasibility <- function(X, totals, d, lower, upper, margin_names,
  action = c("suppress", "join"), merged_name = NULL, tol = 1e-6,
  type = c("g-weight", "w-weight"), method = c("upper", "lower", "both"),
  bounds, search_limits, step) {
  action <- match.arg(action); type <- match.arg(type); method <- match.arg(method)
  if (action == "suppress" && !is.null(merged_name)) stop("merged_name only applies to joins")
  z <- if (action == "suppress") drop_margins(X, totals, margin_names) else
    merge_margins(X, totals, margin_names, merged_name)
  first <- calibration_feasibility(z$X, z$totals, lower, upper, tol)
  if (.ct_feasible_bool(first)) return(c(z, list(bounds = bounds, status = 0L,
    lower = rep_len(lower, nrow(z$X)), upper = rep_len(upper, nrow(z$X)), weights = first$weights)))
  out <- analyze_feasibility_bounds(z$X, z$totals, d, tol, type, method, bounds, search_limits, step)
  c(z, out)
}

#' Apply an ordered list of margin changes
#' @inheritParams check_modified_feasibility
#' @param margin_names_list List of margin selections, one per action.
#' @param actions Character vector of actions.
#' @param merged_names List of aggregate names or NULLs, one per action.
#' @return First feasible modified problem, with iter and action history; if
#'   unsuccessful, the last attempted problem with status=1 (not the original).
#' @export
#' @examples
#' X <- cbind(A = c(1, 0), B = c(0, 1))
#' check_sequential_feasibility(X, c(A = 0, B = 4), c(1, 1), 1, 3,
#'   list(c("A", "B")), "join", list("N"), bounds = c(1, 3),
#'   search_limits = 3, step = 1)
check_sequential_feasibility <- function(X, totals, d, lower, upper,
  margin_names_list, actions, merged_names = rep(list(NULL), length(margin_names_list)), tol = 1e-6,
  type = c("g-weight", "w-weight"), method = c("upper", "lower", "both"),
  bounds, search_limits, step) {
  type <- match.arg(type); method <- match.arg(method)
  if (!is.list(margin_names_list) || length(actions) != length(margin_names_list) ||
      length(merged_names) != length(margin_names_list)) stop("Action lists must have equal lengths")
  z <- list(X = X, totals = totals, bounds = bounds, status = 1L, iter = 0L)
  history <- list()
  for (i in seq_along(margin_names_list)) {
    z <- check_modified_feasibility(z$X, z$totals, d, lower, upper, margin_names_list[[i]],
      actions[i], merged_names[[i]], tol, type, method, bounds, search_limits, step)
    history[[i]] <- list(action = actions[i], variables = margin_names_list[[i]])
    z$iter <- i
    if (z$status == 0L) break
  }
  z$history <- history
  z
}

.ct_families <- function(families, ordered, nm) {
  if (!is.list(families) || !is.logical(ordered) || anyNA(ordered) ||
      length(ordered) != length(families)) stop("Supply a list of families and one logical ordered flag per family")
  alln <- unlist(families, use.names = FALSE)
  if (length(alln)) .ct_names(alln, nm)
  invisible(TRUE)
}
.ct_pairs <- function(families, ordered) {
  out <- list()
  for (i in seq_along(families)) {
    f <- families[[i]]
    if (length(f) < 2L) next
    pairs <- if (ordered[i]) lapply(seq_len(length(f) - 1L), function(k) f[c(k, k + 1L)]) else
      combn(f, 2L, simplify = FALSE)
    out <- c(out, pairs)
  }
  out
}
.ct_update_family <- function(families, join, newname) {
  lapply(families, function(f) {
    if (!any(f %in% join)) return(f)
    first <- min(which(f %in% join))
    f[first] <- newname
    f[!f %in% join]
  })
}

#' Suggest pairs of margins to aggregate
#' @inheritParams calibrate_weights
#' @param margin_groups Disjoint list of families of column names.
#' @param ordered One flag per family: TRUE permits only adjacent merges.
#' @param criterion Heuristic: scaled multipliers or sample support.
#' @return List of distinct-name pairs. Scores are heuristic, not optimality
#'   guarantees; lambda scores use a least-squares target projection if inconsistent.
#' @export
#' @examples
#' X <- cbind(A = c(1, 0), B = c(0, 1))
#' recommend_margin_merges(X, c(1, 1), c(A = 2, B = 2),
#'   margin_groups = list(c("A", "B")), ordered = FALSE)
recommend_margin_merges <- function(X, d, totals, q = 1, margin_groups, ordered,
                                  criterion = c("lambda", "sample"), metric_scale = NULL) {
  criterion <- match.arg(criterion)
  dat <- .ct_named(X, totals); X <- dat$X; totals <- dat$totals
  .ct_families(margin_groups, ordered, colnames(X))
  d <- .ct_vec(d, nrow(X), "d", scalar = FALSE)
  if (criterion == "sample") score <- colSums(X != 0) else {
    q <- .ct_vec(q, length(d), "q"); v <- .ct_metric(d, q, metric_scale)
    red <- .ct_reduce(X, totals)
    la <- .ct_pinv(crossprod(red$A, red$A * v), red$b - as.numeric(crossprod(red$A, d)))
    score <- abs(as.numeric(red$map %*% la)) * pmax(1e-15, sqrt(colSums(X^2)))
  }
  names(score) <- colnames(X)
  out <- list()
  for (i in seq_along(margin_groups)) {
    f <- margin_groups[[i]]
    if (length(f) < 2L) next
    worst <- if (criterion == "sample") f[which.min(score[f])] else f[which.max(score[f])]
    possible <- if (ordered[i]) adjacent_elements(f, match(worst, f)) else setdiff(f, worst)
    best <- if (criterion == "sample") possible[which.max(score[possible])] else possible[which.min(score[possible])]
    out[[length(out) + 1L]] <- c(worst, best)
  }
  out
}

#' Search calibration problems after permitted margin aggregations
#' @inheritParams recommend_margin_merges
#' @inheritParams analyze_feasibility_bounds
#' @param lower,upper Initial final-weight bounds.
#' @param max_solutions Maximum number of candidate solutions.
#' @param max_depth Maximum number of successive pair merges.
#' @param max_nodes Maximum evaluated nodes, guarding combinatorial growth.
#' @param on_inconsistency How to handle incompatible totals: stop returns an
#'   empty search with a diagnostic; merge explicitly allows searching problems
#'   with weakened constraints after permitted aggregations.
#' @return List of candidate lists with fit, weights, g, X, totals, bounds,
#'   history and residuals_original. Attributes searched_nodes and search_status
#'   identify whether the search exhausted the tree or stopped at a limit.
#' @details Breadth-first search explores all allowed pairs, not just heuristic
#'   recommendations. It never silently deletes a nonzero population totals.
#'   Every candidate reports its modified constraints and original residuals.
#' @export
#' @examples
#' X <- cbind(a = c(1, 0), b = c(0, 1))
#' automatic_calibration(X, c(a = 0, b = 4), c(1, 1), lower = 1,
#'   upper = 3, bounds = c(1, 3), search_limits = c(1, 3),
#'   margin_groups = list(c("a", "b")), ordered = TRUE)
automatic_calibration <- function(X, totals, d, q = 1, lower, upper,
  method = c("both", "upper", "lower"), bounds = c(0.1, 10),
  search_limits = c(0.1, 15), step = 1, margin_groups = list(),
  ordered = rep(FALSE, length(margin_groups)), tol = 1e-6,
  max_solutions = 3L, max_depth = 2L, type = c("g-weight", "w-weight"),
  metric_scale = NULL, max_nodes = 1000L, on_inconsistency = c("stop", "merge")) {
  method <- match.arg(method); type <- match.arg(type)
  dat <- .ct_named(X, totals); X <- dat$X; totals <- dat$totals
  .ct_families(margin_groups, ordered, colnames(X))
  .ct_scalar(max_solutions, "max_solutions", integer = TRUE); .ct_scalar(max_depth, "max_depth", TRUE, TRUE)
  .ct_scalar(max_nodes, "max_nodes", integer = TRUE)
  .ct_grid(bounds, search_limits, step, method)
  .ct_bounds(lower, upper, nrow(X))
  d <- .ct_vec(d, nrow(X), "d", scalar = FALSE)
  .ct_metric(d, .ct_vec(q, length(d), "q"), metric_scale)
  on_inconsistency <- match.arg(on_inconsistency)
  initial_check <- calibration_preflight(X, totals, lower, upper, tol)
  if ((initial_check$status == "inconsistent_totals" && on_inconsistency == "stop") ||
      initial_check$status == "numerical_failure") {
    return(structure(list(), class = c("caltools_search", "list"),
      preflight = initial_check, searched_nodes = 0L, search_status = initial_check$status))
  }
  queue <- list(list(X = X, totals = totals, families = margin_groups,
    history = list(), groups = setNames(lapply(colnames(X), identity), colnames(X))))
  seen <- character(); solutions <- list(); nodes <- 0L; head <- 1L
  status <- "exhausted"
  while (head <= length(queue)) {
    if (nodes >= max_nodes) { status <- "max_nodes"; break }
    node <- queue[[head]]; head <- head + 1L
    key <- paste(sort(vapply(node$groups, function(g) paste(sort(g), collapse = "\034"), "")), collapse = "\035")
    if (key %in% seen) next
    seen <- c(seen, key); nodes <- nodes + 1L
    # Incompatible equalities cannot be repaired by moving weight bounds.
    incompatible <- analyze_inconsistency(node$X, node$totals, tol)
    bu <- NULL
    if (!incompatible) {
      f <- if (!length(node$history)) initial_check$feasibility else
        calibration_feasibility(node$X, node$totals, lower, upper, tol)
      bu <- list(lower = lower, upper = upper)
    }
    if (!incompatible && !.ct_feasible_bool(f)) {
      found <- analyze_feasibility_bounds(node$X, node$totals, d, tol, type,
        method, bounds, search_limits, step)
      if (found$status == 0L) bu <- found else bu <- NULL
    }
    if (!is.null(bu)) {
      fit <- calibrate_weights(node$X, d, node$totals, q, "truncated", bu$lower,
        bu$upper, metric_scale = metric_scale, tol = tol)
      solutions[[length(solutions) + 1L]] <- list(fit = fit, weights = fit$weights,
        g = fit$g, X = node$X, totals = node$totals,
        lower = fit$lower, upper = fit$upper, history = node$history,
        residuals_original = setNames(as.numeric(crossprod(X, fit$weights)) - totals, names(totals)))
      if (length(solutions) >= max_solutions) { status <- "max_solutions"; break }
    }
    if (length(node$history) >= max_depth) next
    for (pair in .ct_pairs(node$families, ordered)) {
      merged <- merge_margins(node$X, node$totals, pair)
      newname <- tail(colnames(merged$X), 1L)
      groups <- node$groups[setdiff(names(node$groups), pair)]
      groups[[newname]] <- unlist(node$groups[pair], use.names = FALSE)
      queue[[length(queue) + 1L]] <- list(X = merged$X, totals = merged$totals,
        families = .ct_update_family(node$families, pair, newname),
        history = c(node$history, list(list(action = "join", variables = pair, name = newname))), groups = groups)
    }
  }
  class(solutions) <- c("caltools_search", "list")
  attr(solutions, "preflight") <- initial_check
  attr(solutions, "searched_nodes") <- nodes
  attr(solutions, "search_status") <- status
  solutions
}

#' Interactive selection of audited calibration candidates
#' @inheritParams automatic_calibration
#' @param choose Optional function accepting candidates and returning an index
#'   (or zero to cancel). Allows scripted use and testing without readline.
#' @param ... Additional arguments passed to automatic_calibration.
#' @return Selected candidate with preflight diagnostics; NULL on cancellation;
#'   an empty caltools_search with diagnostics if no candidates exist.
#' @details Replaces the unbounded legacy menus with selection from reproducible
#'   candidates. For explicit deletion use drop_margins first.
#' @export
#' @examples
#' X <- cbind(A = c(1, 0), B = c(0, 1))
#' interactive_calibration(X, c(A = 0, B = 4), c(1, 1), lower = 1, upper = 3,
#'   margin_groups = list(c("A", "B")), ordered = TRUE,
#'   bounds = c(1, 3), search_limits = c(1, 3), choose = function(candidates) 1L)
interactive_calibration <- function(X, totals, d, q = 1, lower, upper,
  margin_groups = list(), ordered = rep(FALSE, length(margin_groups)), tol = 1e-6,
  choose = NULL, ...) {
  if (is.null(choose) && !interactive()) stop("Use automatic_calibration or provide choose in noninteractive sessions")
  candidates <- automatic_calibration(X, totals, d, q, lower, upper,
    margin_groups = margin_groups, ordered = ordered, tol = tol, ...)
  if (is.null(choose)) print(attr(candidates, "preflight"))
  if (!length(candidates)) return(candidates)
  if (is.null(choose)) {
    for (i in seq_along(candidates)) {
      cat("Candidate", i, "\n"); print(candidates[[i]]$history)
      print(candidates[[i]]$fit)
      print(candidates[[i]]$residuals_original)
    }
    k <- suppressWarnings(as.integer(readline("Candidate number (0 to cancel): ")))
  } else k <- choose(candidates)
  if (length(k) != 1L || !is.numeric(k) || is.na(k) || k != floor(k) || k < 0 || k > length(candidates))
    stop("Invalid candidate index")
  if (k == 0) return(NULL)
  selected <- candidates[[k]]
  selected$preflight <- attr(candidates, "preflight")
  selected
}
