# Migration from CalTools 0.1.0 and the supplied prototype

Version 0.2.0 deliberately changes public identifiers to a consistent English API.
Old aliases are not exported. Update named arguments and result-field access as
well as function calls. Mathematical notation such as X, d, q, n and p is retained.

| Previous function | Version 0.2.0 |
|---|---|
| calibra2 | calibrate (or the main calibrate_weights API) |
| squared_sum | calibration_distance_from_factors |
| feasibility_Farkas | is_feasible_lp |
| feasibility_proj | is_feasible |
| feasibility_proj_fast | is_feasible_fast |
| analyze_multicol | analyze_collinearity |
| suppress_marginal_variable | drop_margins |
| join_marginal_variable | merge_margins |
| get_feasibility | check_modified_feasibility |
| get_feasibility_list | check_sequential_feasibility |
| get_lambda1 | linear_multipliers |
| recommended_var_2_join | recommend_margin_merges |
| automatic_calib | automatic_calibration |
| interactive_calib | interactive_calibration |
| obtener_elementos_circundantes | adjacent_elements |
| def_Pi_i_distintas_v4 | random_inclusion_probabilities |

The existing names calibrate_weights, calibration_distance, calibration_feasibility,
relax_bounds, analyze_inconsistency, analyze_zero_sample and analyze_feasibility_bounds
remain. calibration_preflight is new.

| Previous argument/field | New argument/field |
|---|---|
| Xs | X |
| total | totals |
| L_bounds, U_bounds | lower, upper |
| EPS | rank_tol (collinearity); use tol for total consistency |
| EPS1 | tol |
| v_names, v_names_list | margin_names, margin_names_list |
| act, act_v | action, actions |
| join_name, join_name_l | merged_name, merged_names |
| mar_2_join | margin_groups |
| hierarchy | ordered |
| method_join | criterion |
| n_sol | max_solutions |
| vect, posicion | values, position |

Example:

```r
# Previous call:
# calibra2(Xs, d, total, method="truncated", L_bounds=L, U_bounds=U)
fit <- calibrate(X, d, totals, method="truncated", lower=L, upper=U)

report <- analyze_inconsistency(X, totals, details=TRUE)
report$incompatibilities
```

## Additional 0.2.0 behavior changes

- Automatic and interactive calibration run calibration_preflight. The default
  on_inconsistency="stop" returns an empty caltools_search with diagnostics if
  equalities conflict. Explicit "merge" permits altering those equations.
- Search objects have a preflight attribute. A selected interactive candidate
  also has a preflight field. No candidates returns the empty search;
  cancellation returns NULL.
- analyze_inconsistency remains logical by default; details=TRUE returns
  explanatory relations, coefficients, discrepancies and projected totals.
- calibration_feasibility(method="projection") performs actual alternating
  affine/box projections. With fallback=FALSE, an iteration limit returns NA.
- is_feasible and is_feasible_fast both use the validated auto engine.
  Their names do not promise distinct algorithms or a relative speed ranking.

## Important behavior inherited from 0.1.0

- Bounds apply to final weights. Signed/zero starting weights need an explicit
  positive metric_scale for quadratic methods; raking/logit require positive d.
- Fit fields are named; the first two remain g and lambda. Zero d produces NA
  factors, so use weights. linear_multipliers returns a named vector.
- Collinearity indices identify redundant columns, not eigenvalue positions.
- Margin results expose X and totals. Searches return homogeneous candidates
  with history and residuals against original totals.
- Sequential modification returns the last attempted problem after failure.
- The both-bound grid is a nested simultaneous expansion path, not a global
  search over all pairs. relax_bounds is minimal only along its stated direction.
- The probability generator preserves bounds and sum, not the old random law.
- C++ matrix-product helpers are not exported; use matrix multiplication or
  crossprod. No local sourceCpp path is required.
