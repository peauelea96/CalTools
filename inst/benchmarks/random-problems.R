# Run with an installed CalTools package; this file does not run on package load.
# Examples of repeated controlled stress tests for new public APIs.
library(CalTools)

benchmark_statuses <- function(repetitions = 5L, n = 30L,
                               difficulty_grid = c(0, .25, .5, .6, .8, 1),
                               stress_pattern = "local") {
  rows <- list()
  k <- 0L
  for (difficulty in difficulty_grid) for (rep in seq_len(repetitions)) {
    problem <- generate_calibration_problem(
      n = n, categorical_groups = c(5, 4), n_binary = 2,
      n_continuous = 2, n_zero = 1, n_derived = 1,
      difficulty = difficulty, stress_pattern = stress_pattern,
      seed = rep + 10000L * k)
    seconds <- system.time({
      screening <- calibration_preflight(problem$X, problem$totals,
                                         problem$lower, problem$upper)
    })[["elapsed"]]
    stopifnot(identical(screening$status, problem$expected_status))
    k <- k + 1L
    rows[[k]] <- data.frame(repetition = rep, n = n,
      margins = ncol(problem$X), difficulty = difficulty,
      stress_pattern = problem$stress_pattern,
      expected = problem$expected_status, observed = screening$status,
      elapsed_seconds = seconds)
  }
  do.call(rbind, rows)
}

# set.seed is not necessary: seed= is local to each generated problem.
results <- benchmark_statuses()
print(results)

# Separate algebraic contradictions from box failures.
inconsistent <- generate_calibration_problem(n = 20,
  categorical_groups = c(3, 2), inconsistency = "linear", seed = 990)
stopifnot(analyze_inconsistency(inconsistent$X, inconsistent$totals))
repaired <- closest_consistent_totals(inconsistent$X, inconsistent$totals)
stopifnot(!analyze_inconsistency(inconsistent$X, repaired))

# Sample a small difficult case for automated calibration (uses new API names).
problem <- generate_calibration_problem(n = 16,
  categorical_groups = c(4), n_continuous = 0,
  difficulty = .8, stress_pattern = "local", seed = 112)
solutions <- automatic_calibration(problem$X, problem$totals, problem$d,
  q = problem$q, lower = problem$lower, upper = problem$upper,
  bounds = problem$bounds, method = "upper", search_limits = 3,
  step = .1, margin_groups = problem$margin_groups,
  ordered = problem$ordered, max_solutions = 2, max_depth = 2)
print(solutions)
