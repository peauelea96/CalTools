# Mathematical and implementation review, 0.3.0

## Scope and findings

Reviewed the supplied advanced mathematical foundations PDF, the original
Herramientas notes (PDF/R Markdown), prototype R functions, trial scripts and
CalTools 0.2.1 project. The advanced draft's distinction between feasibility,
optimality, convergence and speed is essential. The code implements that
separation; it does not assume the original permanent-fixing proof is valid.

1. The weighted affine projection is correct for a positive metric and consistent
   equalities. A generalized inverse does not make incompatible totals compatible.
2. Permanently fixing every violated coordinate is not a correct general QP
   algorithm or an infeasibility decision procedure. The five-unit objective-gap
   example and the centre-initialized false-negative example are reproduced in
   the tests. A bound on the number of activations is not a correctness theorem.
3. The box KKT characterization and clipped-multiplier identity are necessary
   and sufficient for the strictly convex quadratic problem. Fixed coordinates
   impose no sign condition. Multipliers may be nonunique, including when there
   are no free coordinates; optional LP recovery is needed for a complete search
   for a certificate rather than checking only one least-squares multiplier.
4. The duality gap and the weighted KKT error bounds in the advanced draft hold
   for exactly primal-feasible candidates. Returning them for numerical inputs
   requires this qualification. They are diagnostic estimates unless exact primal
   feasibility and arithmetic have been established separately.
5. The residual-minimization feasibility proposal is correct. Its objective is
   smooth convex with Lipschitz gradient; accelerated projected iterations have
   the stated objective rate. The projection variational inequality bounds the
   distance from Aw_k to the nearest point in the zonotope. Applying the support
   score's Lipschitz bound yields eventual strict separation when the infeasibility
   distance is positive. This does not imply a uniform n-step decision bound.
6. Quadratic calibration with linear inequalities is a strictly convex QP.
   Sign-correct multipliers and complementary slackness give the corresponding
   global certificate. This release selects one active-set QP implementation;
   it does not additionally implement Dykstra, OSQP or HiGHS.

## Implementation choices

- `solver="projection"` performs irreversible fixing centred on the original d.
  Every feasible return carries an independent certificate. Failure on a face
  never establishes that the original problem is infeasible.
- `method="fixing"` in feasibility uses the box centre but a positive Euclidean
  metric. The prototype's implicit centre-as-metric can be zero or negative for
  arbitrary final-weight bounds, so it is not suitable as a general metric.
- `method="pfc"` / `"auto"` first tries the fixing witness, then starts accelerated
  projected residual minimization at the box centre, as assumed in the stated
  convergence bound. It does not initialize FISTA from the failed heuristic.
- Control scaling uses max(1,abs(total)), not signed totals. Spectral-norm step
  size is computed for this scaled matrix. Distance brackets are in scaled units.
- Separation directions are normalized for numerical stability; the returned
  original-coordinate direction permits independent recomputation of support.
- A rounding guard and the directional primal-tolerance margin are subtracted
  before reporting numerical infeasibility. There is no interval-arithmetic
  or rational verification claim.
- Existing alternating projections remain available for feasibility. No claim
  is made that their witness minimizes the calibration distance.
- Mixed controls use `quadprog` for fitting and `lpSolve` for independent
  feasibility. Only the equality block is rank-reduced. Duplicate inequalities
  and interval endpoints are not incorrectly treated as exact equalities.
- Mixed LP uses triplets rather than allocating a dense identity, including
  when unrestricted variables must be split. Mixed QP is dense and capped at
  2000 observations. HiGHS is a future scalability option, not a dependency.
- The equality multipliers returned by `quadprog` can reflect internal sign
  flips. They are reconstructed from stationarity in the original reduced basis.
- Automatic aggregation continues to modify exact controls only. Mixed fitting
  and mixed feasibility are integrated in the general functions; inequality-aware
  aggregation policies are not silently inferred.

## Verification

The 0.3.0 validation log and test-results CSV contain the complete results.
The suite includes 478 passing assertions after implementation, with no failed,
warning or skipped expectations. It covers the mathematical counterexamples,
nonunique multipliers, fixed weights, zero controls, signed weights, mixed
constraints, random comparison with QP and workflow propagation. Runnable README
examples are in `inst/examples/advanced-examples.R`.

The advanced benchmark records 75 scenario/method combinations, five measured
repetitions each, with three calls per timing batch. The raw data retain failures,
uncertified solutions and unknown feasibility states. Small timings must be
interpreted at the clock resolution; no universal speed claim follows.

The original 0.2.x validation files are historical; current results use the
0.3.0 suffix or the `advanced-` prefix.
