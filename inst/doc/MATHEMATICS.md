> Historical 0.2.x notes. For the current 0.3.0 algorithms and solver dispatch,
> see [REVIEW_0.3.0.md](REVIEW_0.3.0.md) and the expanded project README.

# Mathematical specification and review

CalTools 0.2.0. Reviewed 2026-10-09. These derivations describe the implemented
objectives, numerical tolerances and limits. Literature is listed separately in
[PROJECTIONS_AND_LITERATURE.md](PROJECTIONS_AND_LITERATURE.md).

## 1. Problem, scale and existence

Let X have n rows and p columns, d be initial weights, and t be known totals.
Calibration requires X' w = t. For quadratic methods define v_i = s_i q_i > 0
and minimize

    Q(w) = (1/2) sum_i (w_i-d_i)^2 / v_i,
    subject to X'w=t and, when supplied, L<=w<=U.

The default s=d gives the usual chi-square calibration distance and requires
positive d. With s=1, this is a different weighted Euclidean objective. The
prototype exposition omitted d from the denominator while using formulas that
correspond to s=d. Both choices are legitimate, but must not be mixed.

Q has positive definite Hessian diag(1/v), is strictly convex and coercive.
Consequently it has a unique minimizer on any nonempty closed feasible set.
This does not establish feasibility. Negative or zero starting weights can be
used with an explicit positive s; they cannot be used as positive denominators.
Raking and logit require positive d. Factor weights g=w/d are undefined at d=0;
CalTools returns NA there while retaining the valid final weights.

## 2. Linear calibration and constraint reduction

Using Lagrangian Q(w)-lambda'(X'w-t), stationarity gives

    w = d + V X lambda,
    X' V X lambda = t-X'd,  V=diag(v).

A generalized inverse may be used for dependent constraints, provided t is
consistent. Multipliers need not be unique; the minimizing weights are unique.

For numerical stability, scale each column by its maximum absolute entry
(replace zero scales by one): Z=X D^-1. A thin SVD Z=U Sigma V' supplies rank r
using a relative singular-value threshold rank_tol. With retained factors,

    A = U_r,    M = D^-1 V_r Sigma_r^-1,    b = M't.

The independent equations are A'w=b, with A'A=I. Original multipliers are
M times reduced multipliers. This avoids repeatedly decomposing the same
constraint matrix. Every final solution is checked against X and t, not merely
the reduced equations. Numerical rank is a tolerance-dependent decision.

## 3. Explicit incompatible-total relations

In exact arithmetic, X'w=t is solvable if and only if appending t' to X does
not increase its row-space rank. Equivalently Xc=0 implies c't=0.

In scaled coordinates let P=I-V_r V_r'. For each coordinate j, define

    c_j = D^-1 P e_j.

Then X c_j=0 and c_j't is the jth residual of the projection of D^-1 t onto
the retained row space. CalTools normalizes these coefficients, gives readable
margin names, and reports the observed discrepancy and corresponding tolerance.
Approximate dependencies have a small, rather than exactly zero, sample residual;
that residual is included so the relation can be inspected.

The consistency criterion in original units is

    max_j |t_j - [D V_r V_r' D^-1 t]_j| / max(1,|t_j|) <= tol.

rank_X uses rank_tol; rank_augmented increases by one when this criterion fails.
This is a numerical extension-of-row-space diagnostic, not two unscaled rank
calls with an arbitrary shared threshold. The coefficient thresholds are
transformed with their normalization. Duplicate coefficient rows are collapsed.
The reported set can be dependent and is not guaranteed to be a minimal conflict
set. A relation identifies conflicting totals, not which source supplied a wrong
number. projected_totals are diagnostic; the input totals are never overwritten.
No bound expansion can repair a genuinely inconsistent equality system.

## 4. Correct bounded quadratic optimality conditions

For fixed lambda, minimize Q(w)-lambda'X'w over the box. The unique minimizer is

    w_i(lambda) = clip(d_i + v_i x_i'lambda, L_i, U_i).

If X'w(lambda)=t, these weights satisfy the KKT conditions and are the constrained
minimum. The negative dual has gradient X'w(lambda)-t and, away from breakpoints,
Hessian X' diag(v_i I{L_i<d_i+v_i x_i'lambda<U_i}) X.

The implementation works in reduced coordinates. It recomputes bound activity at
every evaluation, uses damped Newton steps with a generalized inverse and a BFGS
fallback. Saturation can make its Hessian singular. A finite optimization budget
does not guarantee convergence. A dense independent QP fallback is available up
to 2000 rows; its quadratic memory explains this explicit limit. Accepted weights
must meet original equations and bounds at the requested numerical tolerance.

## 5. Why permanent clipping is incorrect

The prototype projects, fixes every exceeded coordinate permanently at its bound,
and repeats on the remaining coordinates. Convexity does not justify selecting
that face permanently. Consider

    X = rbind(c(-2,3), c(0,0), c(3,0), c(0,1))
    t = c(1.9, 7.75)
    d = c(1,1,1,1),  L=0, U=2.

The unconstrained Euclidean projection is approximately

    (2.05797872, 1, 2.00531915, 1.57606383).

The prototype fixes both w_1 and w_3 at 2. The first equation then equals 2,
not 1.9, regardless of the two remaining coordinates. It declares infeasibility.
However,

    w = (39/20, 1, 29/15, 19/10)

lies in the box and satisfies both equations exactly. Thus a feasible problem
can be rejected even with positive equal starting weights at the box midpoint.
The regression suite contains this example. A claimed universal bound of n
iterations for this rule cannot establish a correct general feasibility solver.

## 6. Alternating projections and feasible witnesses

Let C={w:A'w=b}, B=[L,U]. Since A'A=I,

    P_C(w) = w + A(b-A'w),
    P_B(w) = clip(w,L,U).

The optional projection method iterates P_B(P_C(w)) from the box midpoint.
It never fixes a bound permanently. For nonempty intersections of these closed
convex sets in finite dimensions, the classical alternating-projection sequence
converges to an intersection point. This is an asymptotic statement: the count
needed for a given tolerance depends on geometry and conditioning.

This method supplies a feasible witness, not in general the closest point in
the intersection to its starting point. Dykstra's algorithm adds correction
vectors for the latter best-approximation problem; it is not implemented here.

After reduction each projection iteration costs O(nr), with O(nr) stored basis.
Exhausting max_iter returns feasible=NA when fallback=FALSE, or invokes the LP
when fallback=TRUE. Failure to find a witness is not a separating certificate.
The default auto method uses the quadratic dual and LP fallback instead.

## 7. LP feasibility, Farkas and signed totals

For finite bounds, write w=L+z. Feasibility becomes

    X'z=t-X'L,    0<=z<=U-L.

This preserves negative final weights while respecting lpSolve's nonnegative
variables. Constraint triplets avoid materializing a dense n-by-n identity.
The LP solves the primal feasibility problem and returns a witness when possible.
A box-support check first bounds each total using

    min_j = sum_i min(x_ij L_i, x_ij U_i),
    max_j = sum_i max(x_ij L_i, x_ij U_i).

A violated bound certifies infeasibility. Passing all separate checks is not
sufficient for joint feasibility. The rank test, dual witness or LP completes
the decision. LP infeasibility is a floating-point solver status, not an exact
rational certificate. Other solver failures or inaccurate witnesses return NA.

Farkas' lemma supplies an alternative: for Aw<=b with unrestricted w, either a
solution exists or there is y>=0 with A'y=0 and b'y<0. It is a theorem, not a
runtime algorithm. The prototype's dense LP implementing it has substantial
avoidable storage. Its written epsilon maximization needs normalization or
epsilon<=1; the original code did include this bound, although the exposition
did not. Comparing timings against that implementation is not a general theorem
that projections outperform linear programming.

All totals use the denominator max(1,abs(t_j)). A signed total must never be used
as a tolerance width: t_j*(1-tol) and t_j*(1+tol) reverse order for t_j<0, and
division by t_j fails at zero. Tests cover negative totals with positive starting
weights and observation-specific bounds, zero totals, all-negative targets, and
invariance to negating a column and its target together. An impossible negative
total on a nonnegative column under nonnegative bounds remains infeasible.

## 8. Raking and vector-bound logit

For d_i>0, raking minimizes

    sum_i d_i/q_i * (g_i log(g_i)-g_i+1),  g_i=w_i/d_i>=0,

with 0 log 0=0. Interior stationarity yields w_i=d_i exp(q_i eta_i), and the
dual Hessian is X'diag(q_i w_i)X. The prototype omitted q_i in the derivative.
The exponential parametrization seeks positive solutions. Boundary-only targets
can make its multipliers diverge; nonconvergence is reported explicitly.

For logit, let l_i=L_i/d_i<1<u_i=U_i/d_i and

    a_i=(u_i-l_i)/[(1-l_i)(u_i-1)].
    Q=sum_i d_i/(q_i a_i) * [(g_i-l_i) log((g_i-l_i)/(1-l_i))
                            +(u_i-g_i) log((u_i-g_i)/(u_i-1))].

Stationarity gives

    p_i=plogis(logit((1-l_i)/(u_i-l_i))+a_i q_i eta_i),
    g_i=l_i+(u_i-l_i)*p_i,
    dw_i/deta_i=d_i (u_i-l_i) a_i q_i p_i(1-p_i).

Thus there is no mathematical requirement that all factor bounds be identical.
Stable logistic and softplus formulas avoid direct exponential ratios. The
objective is extended continuously to finite endpoints. As with raking, the
finite-multiplier algorithm does not guarantee attainment of boundary-only targets.

## 9. Directional bound relaxation

Given nonnegative rates a,b, minimize rho>=0 under X'w=t and
L-rho*a<=w<=U+rho*b. The substitution z=w-L+rho*a>=0 gives

    X'z-rho*X'a=t-X'L,
    z-rho*(a+b)<=U-L.

This LP in n+1 nonnegative variables is equivalent to the stated problem.
It finds the minimal expansion along this direction, not under every distance
between boxes. Some directions may never permit feasibility. The common-bound
grid defines nested boxes; monotonicity justifies binary search for the first
feasible grid member. With method='both' the endpoints expand simultaneously,
not over every possible pair of lower and upper limits.

## 10. Aggregation and search

Replacing columns in J by their sum, with total sum_j t_j, weakens constraints:
every original solution satisfies the aggregate; the converse need not hold.
Deleting a margin also weakens constraints. Neither operation proves a reduction
in statistical bias or variance. An external survey-design rationale is needed.

Automatic and interactive calibration run an initial preflight. By default,
inconsistent totals stop the search with their explanatory report. Explicit
on_inconsistency='merge' allows changing equations through authorized merges.
Each node checks consistency before any bound search. Candidate outputs retain
residuals against all original totals, preventing a changed problem from being
mistaken for a solution to the original one.

Breadth-first search updates families, deduplicates partitions and respects
solution, depth and node limits. Ordered families allow adjacent merges only.
No nonzero total is discarded merely because its sample column is zero.
Multiplier-based merge scores are heuristic and representation-dependent under
collinearity, even after column-scale adjustment.
