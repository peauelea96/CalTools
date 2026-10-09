test_that("Euclidean closest consistent totals have the exact analytic solution", {
  X <- cbind(A = c(1, 0), B = c(0, 1), N = c(1, 1))
  target <- c(A = 4, B = 6, N = 11)
  adjusted <- closest_consistent_totals(X, target)
  expect_equal(unname(adjusted), c(13/3, 19/3, 32/3), tolerance = 1e-12)
  expect_identical(names(adjusted), colnames(X))
  expect_false(analyze_inconsistency(X, adjusted))
  expect_true(analyze_inconsistency(X, target))
  report <- closest_consistent_totals(X, target, details = TRUE)
  expect_s3_class(report, "caltools_total_projection")
  expect_equal(report$adjusted, adjusted)
  expect_equal(report$adjustment, adjusted - target)
  expect_equal(report$distance, sqrt(1/3), tolerance = 1e-12)
  expect_output(print(report), "Closest consistent totals")
  # Orthogonality means every already-consistent total is at least as far.
  for (w in list(c(0, 0), c(4, 6), c(100, -20))) {
    consistent <- as.numeric(crossprod(X, w))
    expect_gte(sum((consistent - target)^2) + 1e-10,
               sum((adjusted - target)^2))
  }
})

test_that("Scales change the geometry, without changing algebraic consistency", {
  X <- cbind(A = c(1, 0), B = c(0, 1), N = c(1, 1))
  target <- c(A = 4, B = 6, N = 11)
  a <- closest_consistent_totals(X, target, scales = c(1, 1, 10))
  expect_equal(unname(a), c(4 + 1/102, 6 + 1/102, 11 - 100/102),
               tolerance = 1e-12)
  expect_false(analyze_inconsistency(X, a))
  expect_error(closest_consistent_totals(X, target, scales = c(1, 0, 2)),
               "strictly positive")
  expect_error(closest_consistent_totals(X, target, scales = 1:2), "length 3")
  expect_error(closest_consistent_totals(X, target, details = NA), "details")
})

test_that("Projection handles full rank, zero columns and no columns", {
  expect_equal(closest_consistent_totals(diag(3), 1:3),
               setNames(1:3, paste0("margin_", 1:3)))
  X <- cbind(A = c(1, -1), zero = c(0, 0))
  corrected <- closest_consistent_totals(X, c(A = -7, zero = 3))
  expect_equal(unname(corrected), c(-7, 0))
  expect_false(analyze_inconsistency(X, corrected))
  expect_equal(closest_consistent_totals(matrix(0, 2, 2), c(5, -2)),
               setNames(c(0, 0), paste0("margin_", 1:2)))
  empty <- closest_consistent_totals(matrix(numeric(), 3, 0), numeric())
  expect_length(empty, 0)
})

test_that("Generator shapes, aliases, reproducibility and RNG isolation work", {
  set.seed(14)
  state <- .Random.seed
  args <- list(n = 40, categorical_groups = c(5, 4), n_binary = 3,
               n_continuous = 2, n_positive = 1, n_zero = 2,
               n_derived = 2, weight_cv = .25, difficulty = .25,
               ordered_groups = c(TRUE, FALSE), seed = 84)
  a <- do.call(generate_calibration_problem, args)
  b <- do.call(generate_calibration_problem, args)
  expect_identical(.Random.seed, state)
  expect_identical(a$X, b$X)
  expect_identical(a$totals, b$totals)
  expect_equal(nrow(a$X), 40)
  expect_equal(ncol(a$X), 1 + 5 + 4 + 3 + 2 + 1 + 2 + 2)
  expect_identical(names(a$totals), colnames(a$X))
  expect_equal(a$X, a$Xs)
  expect_equal(a$totals, a$total)
  expect_equal(a$L_bounds, a$lower)
  expect_equal(a$U_bounds, a$upper)
  expect_equal(a$mar_2_join, a$margin_groups)
  expect_equal(a$hierarchy, c(TRUE, FALSE))
  expect_equal(a$ordered, c(TRUE, FALSE))
  expect_equal(length(a$margin_groups), 2)
  expect_true(all(a$true_weights >= a$lower - 1e-12))
  expect_true(all(a$true_weights <= a$upper + 1e-12))
  expect_equal(unname(a$totals), as.numeric(crossprod(a$X, a$true_weights)))
  expect_equal(a$expected_status, "feasible")
  expect_equal(calibration_preflight(a$X, a$totals, a$lower, a$upper)$status,
               "feasible")
})

test_that("Difficulty controls a provable initial bound obstruction", {
  for (pattern in c("local", "global")) {
    a <- generate_calibration_problem(n = 25, categorical_groups = c(3, 2),
           n_continuous = 0, difficulty = .9, stress_pattern = pattern, seed = 7)
    expect_false(analyze_inconsistency(a$X, a$totals))
    expect_equal(a$expected_status, "bound_infeasible")
    expect_equal(calibration_preflight(a$X, a$totals, a$lower, a$upper)$status,
                 "bound_infeasible")
    if (pattern == "local") {
      j <- match(a$stress_margin, colnames(a$X))
      expect_gt(a$totals[j], sum(a$X[, j] * a$upper))
    } else {
      expect_true(all(a$true_weights > a$upper))
    }
    # A sufficiently expanded upper box contains the constructive witness.
    expect_true(calibration_feasibility(a$X, a$totals, a$lower,
                                        pmax(a$upper, a$true_weights))$feasible)
  }
  a <- generate_calibration_problem(n = 10, categorical_groups = integer(),
          n_continuous = 1, difficulty = .8, stress_pattern = "local", seed = 8)
  expect_equal(a$stress_pattern, "global")
  expect_equal(a$expected_status, "bound_infeasible")
})

test_that("Generator can force inconsistent linear or empty margins", {
  for (kind in c("linear", "zero")) {
    a <- generate_calibration_problem(n = 15, categorical_groups = c(3, 2),
               n_zero = 0, n_continuous = 0, inconsistency = kind, seed = 99)
    expect_true(analyze_inconsistency(a$X, a$totals))
    expect_equal(a$expected_status, "inconsistent_totals")
    expect_equal(calibration_preflight(a$X, a$totals, a$lower, a$upper)$status,
                 "inconsistent_totals")
    adjusted <- closest_consistent_totals(a$X, a$totals)
    expect_false(analyze_inconsistency(a$X, adjusted))
    if (kind == "zero") expect_true(any(colSums(abs(a$X)) == 0))
  }
  expect_error(generate_calibration_problem(difficulty = 1.5), "difficulty")
  expect_error(generate_calibration_problem(bounds = c(1.1, 2)), "bounds")
  expect_error(generate_calibration_problem(categorical_groups = integer(),
     include_intercept = FALSE, n_continuous = 1,
     difficulty = .8, stress_pattern = "global", seed = 2),
     "nonzero nonnegative margin")
})

test_that("Generated benchmark inputs are accepted by automatic_calibration", {
  a <- generate_calibration_problem(n = 12, categorical_groups = 3,
         n_continuous = 0, difficulty = .75, seed = 123)
  out <- automatic_calibration(a$X, a$totals, a$d, a$q,
    lower = a$lower, upper = a$upper, bounds = a$bounds,
    method = "upper", search_limits = 3, step = .25,
    margin_groups = a$margin_groups, ordered = a$ordered,
    max_solutions = 1, max_depth = 1, max_nodes = 10)
  expect_s3_class(out, "caltools_search")
  expect_equal(attr(out, "preflight")$status, "bound_infeasible")
})
