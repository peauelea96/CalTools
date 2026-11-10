> Historical 0.2.x notes. For the current 0.3.0 algorithms and solver dispatch,
> see [REVIEW_0.3.0.md](REVIEW_0.3.0.md) and the expanded project README.

# Projection methods: efficiency, correctness and literature

CalTools 0.2.0, 2026-10-09.

## Assessment of the supplied algorithm

The computational motivation is sound: when n is large and the constraint rank
r is small, reuse a reduced factorization and work with r-dimensional systems.
This can be much cheaper than constructing a dense general-purpose LP or QP in
sample-size dimensions. This is an assessment of the algebra and implementation,
not a claim of universally faster convergence.

The original permanent-clipping rule, however, is not a correct general
feasibility test. Its selected face may exclude all solutions even when the box
contains feasible weights. The explicit four-unit example in MATHEMATICS.md
uses positive midpoint weights and refutes the proposed general justification.
A fast false negative is not a valid improvement in a feasibility algorithm.

CalTools offers three distinct paths:

| Method | Work after reduction | Meaning of success | Main limitation |
|---|---|---|---|
| Alternating affine/box projections | O(nr) per iteration | Checked feasible weights | May converge slowly; exhaustion is inconclusive |
| Reduced quadratic dual | Roughly O(nr^2+r^3) per Newton step, plus line search | Checked feasible weights; for calibration, KKT optimum | Saturation can make the Hessian singular |
| Shifted primal LP | Solver- and sparsity-dependent | Feasibility witness or solver infeasibility status | Potentially higher cost; numerical tolerances still matter |

For n>=p, the dense SVD preparation is roughly O(np^2). The stored basis is
O(nr). Detailed inconsistency reporting also creates a p-by-p projector. These
costs matter when there are many margins or many repeated modified problems.
The package has not yet implemented a sparse-matrix backend or warm starts
between search nodes. The C++ wrappers in the prototype only wrap matrix
products; they do not remove repeated factorizations or improve the mathematical
algorithm. crossprod and matrix multiplication already use R's BLAS backend.

## What is established in the literature?

Bauschke and Borwein (1996) survey projection algorithms for convex feasibility,
including convergence and the role of regularity. The general idea of solving
intersecting affine and convex constraints by projections is established.

Boyle and Dykstra (1986) address finding the nearest point in an intersection
using projections onto simpler convex sets. Dykstra's correction terms distinguish
best approximation from ordinary alternating projections. CalTools implements
ordinary alternating projections only as a feasibility option; it does not label
that output the closest calibrated weight vector.

For survey calibration, Devaud and Tille (2019), Section 3.2.2, are particularly
relevant. They discuss reduced Newton systems, loss of invertibility when truncated
calibration saturates at bounds, and a general primal-dual interior-point
alternative. Their computational discussion explains the appeal of low-dimensional
Newton calibration, while also emphasizing possible convergence failures.
It does not validate permanent clipping of all exceeded coordinates.

These references establish related method families. They do not establish that
the exact supplied procedure is a correct new algorithm, nor that it is novel.
The counterexample rules out its claimed unrestricted correctness. Any restricted
correctness claim would need explicit structural assumptions and a separate proof.

## How to interpret a comparison with Farkas

Farkas' lemma is a theorem of alternatives, not one particular solver. The
prototype encoded its certificate search in a dense LP. Avoiding that dense
construction can substantially improve memory and time without contradicting
or improving the theorem. A fair practical comparison must identify:

- The exact LP formulation, solver, tolerances and sparse/dense representation.
- Whether the problem is feasible, infeasible, or very near a boundary.
- The number of margins, numerical rank, conditioning and box tightness.
- Whether factorization and fallback costs are included.
- Whether every positive decision has a checked witness and iteration limits
  remain inconclusive rather than becoming false negatives.

The included benchmark compares complete calls to projection, auto and LP on
identical synthetic feasible problems, includes decomposition time, and reports
solver/status as well as elapsed time. It is a reproducible illustration, not
a representative survey workload or universal speed result. See VALIDATION.md.

## References and primary sources

1. Bauschke, H. H., and Borwein, J. M. (1996). On projection algorithms for solving
   convex feasibility problems. *SIAM Review*, 38(3), 367-426.
   [Author's publication list, 1996 section](https://bauschke.ca/publications).
2. Boyle, J. P., and Dykstra, R. L. (1986). A Method for Finding Projections onto
   the Intersection of Convex Sets in Hilbert Spaces. In *Advances in Order
   Restricted Statistical Inference*, Lecture Notes in Statistics 37, 28-47.
   [Publisher and DOI](https://doi.org/10.1007/978-1-4613-9940-7_3).
3. Devaud, D., and Tille, Y. (2019). Deville and Sarndal's calibration: revisiting
   a 25-years-old successful optimization problem. *TEST*, 28, 1033-1065.
   [DOI](https://doi.org/10.1007/s11749-019-00681-3).
   [Author manuscript](https://libra.unine.ch/server/api/core/bitstreams/9e007fb7-bad7-45c1-83b2-eb9ab7a8d289/content).
4. Deville, J.-C., and Sarndal, C.-E. (1992). Calibration Estimators in Survey
   Sampling. *Journal of the American Statistical Association*, 87(418), 376-382.
   [DOI](https://doi.org/10.1080/01621459.1992.10475217).
