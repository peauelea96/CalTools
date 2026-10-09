# CalTools 0.2.1

* Added `closest_consistent_totals()` to project incompatible total vectors onto
  the algebraically consistent subspace using Euclidean or user-scaled distances.
  Separate from bounded-weight calibration and from the existing diagnostic
  projection, which uses column normalization.
* Added `generate_calibration_problem()` for reproducible categorical, binary,
  continuous, zero-column and dependent-column examples, with guaranteed
  bounded-feasibility regimes and optional algebraic inconsistencies.
* Added tests, help pages, README examples and a repeatable stress-test script.

# CalTools 0.2.0

* Changed public identifiers, arguments, comments, messages and all documentation
  to English; added a complete migration table for this breaking API change.
* Added explicit null-space relations and total discrepancies to
  analyze_inconsistency(details=TRUE).
* Added calibration_preflight and initial screening in both search workflows.
  Inconsistent equalities stop by default; merging requires explicit opt-in.
* Added alternating affine/box projections as an optional feasibility method.
  Iteration exhaustion remains inconclusive or invokes the LP fallback.
* Added regression tests specifically for negative totals with positive starting
  weights, vector bounds, zero totals and sign invariance.
* Expanded roxygen help with runnable examples and rebuilt the English README.Rmd.
* Added mathematical derivations and a sourced assessment of projection methods.

# CalTools 0.1.0

* Rebuilt vector-bounded calibration using a convex final-weight formulation.
* Signed and zero initial weights supported with an explicit positive metric.
* Replaced permanent clipping by a reversible active-set dual formulation.
* Added LP feasibility, witness checks, and directional minimal expansion.
* Corrected raking/logit derivatives, objective evaluation and multiplier returns.
* Rebuilt margin aggregation, grid search, BFS and interactive selection.
* Removed external sourceCpp paths and session-wide library()/plot side effects.
* See inst/doc/AUDIT.md for the complete mathematical and migration review.
