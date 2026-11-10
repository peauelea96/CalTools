# CalTools 0.3.0 validation

Release checked on 20 October 2026 with R 4.3.3 on Ubuntu 24.04 (x86_64).

- Package installation from source: passed.
- Unit/regression tests: 478 passing assertions; 0 failures, warnings or skips.
- Advanced README examples: executed successfully. The five-unit projection
  example deliberately emits the documented uncertified-optimality warning.
- Advanced benchmarks: 75 scenario/method combinations, five timing repetitions
  per combination, each averaging three calls. Raw data and solver statuses retained.
- R CMD build: passed.
- R CMD check --no-manual: 0 errors, 0 warnings, 1 note.
- The note reports optional suggested packages unavailable in this environment:
  rmarkdown and sampling. The check used _R_CHECK_FORCE_SUGGESTS_=false. Neither
  is required by the new numerical methods. No PDF manual was built.

The runnable examples, documentation examples and tests check this version's
code. Earlier unsuffixed validation logs and 0.2.1 records are historical.

Files:

- validation/test-results-0.3.0.csv
- validation/tests-0.3.0.log
- validation/advanced-examples.log
- validation/advanced-benchmarks.csv
- validation/advanced-benchmark-summary.csv
- validation/benchmark-table.md
- validation/session-info-0.3.0.txt
- validation/R-CMD-check-0.3.0.log

Important scope: mixed calibration is quadratic only and its single dense QP
backend is limited to 2000 observations. The feasibility LP is independent.
Fixing alone can be inconclusive. Floating-point KKT and separation checks are
numerical certificates, not exact-arithmetic verification.
