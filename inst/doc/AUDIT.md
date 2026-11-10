> Historical 0.2.x notes. For the current 0.3.0 algorithms and solver dispatch,
> see [REVIEW_0.3.0.md](REVIEW_0.3.0.md) and the expanded project README.

# Prototype audit and release review

Review completed 2026-10-09 for CalTools 0.2.0. Sources supplied by the user:
calibra_tools.R, Herramientas_calibrado.rmd, its PDF reading version, pruebas.R,
and fastMatProd.cpp. Mathematical review uses the formulas and claims in the
R Markdown source. No administrative microdata were supplied or used.

## Findings and treatment

| Component | Finding | Current treatment |
|---|---|---|
| Quadratic objective | The exposition omitted d from the denominator but used formulas for scale d | Explicit positive metric_scale; distinguish scale d and scale 1 |
| Permanent clipping proof | Convexity does not justify fixing all exceeded coordinates; feasible problems can be rejected | Counterexample, reversible dual bounds, optional alternating projections |
| Signed starting weights | Factor bounds reverse under negative d; zero d cannot define factors; signed quadratic denominators break convexity | Final-weight formulation, explicit positive metric, NA factors at zero d |
| Negative/zero totals | Signed tolerance widths reverse inequalities; zero denominators are invalid | max(1,abs(totals)) scaling in residuals and LP constraints |
| Input validation | Recycling, invalid q, vector defaults and nonfinite inputs | Dimensions, names, finiteness, positivity and method validation |
| Iteration failure | Nonconverged output could survive as a usable result | Check original equations and bounds; distinguish numerical failure and infeasibility |
| Raking | Missing q in Jacobian; final multipliers overwritten | Complete derivative and correct returned multipliers |
| Logit | Incorrect derivative, unnecessary constant-factor-bound restriction, unstable exponential ratios | Stable vector-bound logistic formulas and exact derivative |
| Objective evaluation | Misspelled method branch, invalid conditional and zero-log-zero | Rewritten objectives with continuous endpoint values |
| Farkas LP | Dense matrices scaling quadratically with sample size, signed tolerances, incomplete status handling | Shifted primal LP, triplets, checked witnesses and explicit statuses |
| Farkas exposition | Missing normalization in written epsilon program | Clarify epsilon<=1, which the original code already included |
| Fast feasibility | Undefined state variable and same irreversible clipping flaw | Shared corrected reduced solver; no separate unvalidated fast path |
| Matrix-product C++ | Thin wrappers, mixed Eigen/Armadillo dependencies, by-value matrix inputs | BLAS-backed crossprod and one constraint-space reduction per solver call |
| Script loading | Local sourceCpp path and attaching packages inside functions | Declared package dependencies, no machine-specific paths |
| Bound grid | Absolute-weight lower branch multiplied by d; endpoints omitted | Correct absolute bounds, endpoint inclusion, nested grid and binary search |
| Rank inconsistency | Boolean rank comparison lacked explanation | Scaled SVD and optional checkable conflicting relations |
| Collinearity | Gram matrix squares conditioning; spectral indices are not redundant columns | Scaled QR, explicit relations and redundant-column indices |
| Margin changes | Dimension dropping, all-column edge cases and name collisions | Preserved matrix shape, empty-column support and name validation |
| Modified feasibility | Starting vector bounds ignored | Check actual box before common-bound grid |
| Sequential changes | Empty loops, unchecked action lengths and ambiguous final state | Validated lists, defined empty case, last attempted problem and history |
| Zero sample | Nonzero totals could be deleted silently | Diagnose, never automatically delete an equation |
| Adjacent elements | Singleton produced NA | Valid empty neighborhood |
| Merge recommendation | Ties could duplicate a name; multipliers depend on scale/representation | Distinct pairs, families, scale adjustment and documented heuristic scope |
| Automatic search | Premature termination, incomplete family updates and inconsistent outputs | BFS, deduplication, homogeneous candidates and explicit limits |
| Initial screening | Workflows attempted optimization before explaining contradictory totals | Preflight report, default stop, explicit opt-in to weakened merged equations |
| Interactive workflow | Unbounded menus and graphical side effects | Finite candidate selection, scripted choose callback and preflight report |
| Random probabilities | Endpoint loops and missing validation | Normalized random sizes with certainty units and defined endpoints |
| Exploratory tests | Local paths and no mathematical assertions | Reproducible tests, analytic/KKT constructions, QP and LP comparisons |
| Language and help | Mixed identifiers and incomplete usage documentation | English API, messages, comments, README.Rmd, roxygen and generated help |

## Scope

The suite checks analytic solutions, independently constructed KKT optima,
reference QP/LP results, regression examples and edge cases. This is synthetic
validation, not a proof that no numerical issue can occur. In particular, extreme
scaling, near dependence and boundary-only nonlinear targets need care. Statistical
validity of a chosen aggregation depends on the survey design and auxiliary data.
The package does not implement design-based variance estimation.

No novelty claim is made for standard numerical methods. Integration, correction
and diagnostics are the contribution of this release. A sparse numerical backend,
factorization caching across search nodes, warm starts and a compiled whole-solver
backend are possible future work, to be justified with actual workload profiles.

## Credits and proposed license

The supplied prototype states that its calibration function adapts sampling::calib.
The sampling package is by Yves Tille and Alina Matei and uses GPL (>= 2).
CalTools proposes the compatible GPL (>= 2) license and retains these credits.
The numerical implementation has been rewritten; there is no hidden fallback to
the original clipping algorithm. Complete author information and maintainer email
before publishing; maintainer@example.org is an explicit placeholder.

Official interface references:

- https://cran.r-project.org/web/packages/sampling/sampling.pdf
- https://search.r-project.org/CRAN/refmans/quadprog/html/solve.QP.html
- https://search.r-project.org/CRAN/refmans/lpSolve/html/lp.html

See MIGRATION.md for all public-name changes, MATHEMATICS.md for derivations,
and PROJECTIONS_AND_LITERATURE.md for the projection-method assessment.
