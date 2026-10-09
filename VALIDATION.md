# CalTools 0.2.0 validation record

Date: 2026-10-09. R 4.3.3, Ubuntu 24.04.3, x86_64, system BLAS/LAPACK.
Synthetic data only. Windows, macOS and the minimum supported R version were
not executed locally. The supplied GitHub workflow can check those platforms.

## Checks executed

- R CMD check --no-manual: Status: OK, with no errors, warnings or notes.
- Source-package installation and namespace loading succeeded.
- 25 testthat blocks, 299 assertions: all passed, no failures, warnings or skips.
- All 24 public functions have roxygen source, generated help and example calls.
- README.Rmd rendered successfully to README.md with all enabled chunks executed.
- Analytic linear solution and comparison with sampling 2.10.
- 35 random bounded problems with signed starting weights compared to quadprog
  1.5-8 and independently constructed KKT optima.
- 20 constructed feasible LP problems, infeasible cases and conflicting totals,
  using lpSolve 5.6.20.
- Raking and vector-bound logit with varying q and known multiplier solutions.
- Negative totals with positive initial weights, vector bounds, zero totals,
  all-negative targets, and invariance under negating a column and its total.
- Inconsistent-total coefficient evidence, scaled columns, full rank, zero
  columns, empty constraints and separate rank/total tolerances.
- Preflight distinction between algebraic inconsistency and restrictive bounds;
  default workflow stop and explicitly authorized aggregation.
- Permanent-clipping counterexample, alternating-projection witnesses,
  inconclusive iteration limit and LP fallback.
- Search requiring two successive merges, family updates, original residuals,
  scripted interactive selection and probability-generator endpoints.

The source archive is checked in a fresh directory so documentation from the
previous version cannot contaminate its examples. The final check log is
included at validation/R-CMD-check.log. PDF manual generation is not part of
`R CMD check --no-manual`; generated function help and examples are checked.

## Dense bounded calibration benchmark

Median elapsed seconds over three complete calls, 12 margins, active bounds,
and a known optimum constructed from KKT multipliers:

| Rows | Reduced dual | Reference QP | Maximum dual weight error |
|---:|---:|---:|---:|
| 200 | 0.002 s | 0.015 s | 5.95e-14 |
| 1000 | 0.002 s | 1.649 s | 8.88e-14 |
| 10000 | 0.021 s | Not run | 4.55e-14 |

This illustrates the benefit over a dense sample-dimension QP in these cases.
It does not establish an advantage over all other implementations. In the 0.1.0
review, one separate easy 1000-row case was faster in the original R prototype
(0.001 s versus 0.002 s); that historical observation remains relevant. The
original C++ wrappers were not benchmarked as a complete calibration solver.

## Bounded feasibility benchmark

Identical constructed feasible problems per method, 8 margins, three complete
calls per method, tolerance 1e-8, projection limit 5000 with fallback disabled.
Timings include validation and constraint-space decomposition.

| Rows | Alternating projections | Reduced dual (auto) | Shifted primal LP |
|---:|---:|---:|---:|
| 100 | 0.007 s | 0.001 s | 0.002 s |
| 500 | 0.009 s | 0.002 s | 0.010 s |
| 2000 | 0.021 s | 0.004 s | 0.087 s |

All nine method/problem combinations returned checked feasible witnesses.
Projection residuals were below 1e-8; dual and LP residuals below 1.4e-12.
At 100 rows the LP was faster than ordinary projections. These measurements
do not support a universal ranking, especially on infeasible or nearly tangent
problems. Farkas' lemma itself is not an algorithm whose runtime can be measured.
The comparison is with the stated primal LP, not the old dense certificate LP.
Millisecond timings depend on timer resolution, hardware, BLAS and system load.

Reproducible scripts are in inst/benchmarks/. Raw measurements, test summaries,
and sessionInfo are in validation/. The test RDS retains structured testthat
results; the CSV summarizes counts without requiring R object inspection.

## Before publication

Complete the real author name and maintainer email in DESCRIPTION and configure
the repository. maintainer@example.org is a declared placeholder. No GitHub or
CRAN publication was performed. Survey-specific statistical and performance
validation remains to be done with the intended real workload.
