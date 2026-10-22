# Exploratory release checks

## Executed locally

- Wheel built successfully with pip wheel, without build isolation, using the existing setuptools installation.
- Wheel installed in a separate Python 3.12 virtual environment; imports and the CLI ran from /tmp, outside the checkout. The package loaded from that environment's site-packages, rather than the source directory.
- That environment inherited third-party scientific dependencies using system-site-packages. Fully isolated dependency resolution was not tested locally.
- Subprocess tests cover zero permutations, stratified analysis, no eligible terms and unknown GO IDs with failed metadata preserved.
- Python 3.8 compileall succeeded for package, benchmark and tests. Numerical execution on Python 3.8 was not performed because its scientific dependencies are absent locally.

## CI configuration

.github/workflows/tests.yml builds source and wheel distributions, installs the wheel with dependencies, runs tests, checks imports and CLI outside the checkout and runs a small benchmark on Python 3.8 and 3.12. R/Bioconductor comparators are optional external integration checks and are not installed by this workflow. Remote CI results are pending; workflow creation is not evidence that the matrix passed.

## Publication gates

- The full 1,000-replicate / 9,999-permutation validation (63,000 analyses across 7 scenarios) was executed and consolidated in `benchmark_results/validation/` (see `docs/benchmark-protocol.md`).
- Gene-level comparisons, independent biological studies and paired uncertainty of bin differences remain open.
- Review the remote CI matrix and a fully isolated dependency installation before publishing a package.
- Preserve exploratory wording and the migration notes. Software test success does not establish statistical calibration.

No version bump, merge or external publication is performed by these checks.
