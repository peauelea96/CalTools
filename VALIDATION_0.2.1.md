# Validation note for CalTools 0.2.1 additions

The supplied 0.2.0 package already includes prior test and R CMD check logs.
Those logs do **not** validate the newly added 0.2.1 functions.

New test file: `tests/testthat/test-generation-and-projection.R`.
It specifies the following regression cases:

- Analytic Euclidean projection and weighted projection with nontrivial scale.
- Closest point orthogonality, naming, rank-zero, zero-column and full-rank cases.
- Generator shape, column families, duplicated and empty margins, seed isolation.
- Guaranteed feasible tests below the difficulty threshold.
- Verified box-support contradictions for hard local/global stress cases.
- Deliberate inconsistent linear and zero-margin targets.
- Smoke test passing generated inputs to `automatic_calibration()`.

**Execution status:** The editing environment does not include an R executable or
installed R package dependencies. These new tests have therefore not been run.
Do not treat the previous 0.2.0 validation logs as evidence of a passing 0.2.1
R CMD check.

Execute in an environment with R >= 4.1, lpSolve, quadprog, testthat:

```r
install.packages(c("lpSolve", "quadprog", "testthat", "roxygen2"))
roxygen2::roxygenise("CalTools")
testthat::test_local("CalTools")
```

Then, from a shell, run:

```bash
R CMD build CalTools
R CMD check --no-manual CalTools_0.2.1.tar.gz
```
