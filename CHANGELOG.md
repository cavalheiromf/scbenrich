# Changes

## Unreleased — exploratory analysis

- Independent marginal and stratified decisions; explicit support and balance states.
- BH-adjusted stratified p-values and separate discovery lists.
- Marginal OR and intervals computed using SciPy; explicit non-identifiable status.
- Shared input validation, string identifiers, configurable bins, term support and Monte Carlo diagnostics.
- GO replacement-chain resolution, ambiguity and cycle detection, separate summary maps.
- Empty-hypothesis outputs, protected output directories and success/failure provenance.
- Reproducible simulation runners and audited family-level external comparators.
- Precomputed deterministic-stratum contributions reduce permutation runtime; seeded streams may change (see migration notes).
- Python 3.8/3.12 CI configuration and subprocess CLI tests.

See docs/migration.md for renamed columns and files. The large-scale validation (63,000 analyses, 1,000 replicates across 7 scenarios) was executed and consolidated; biological case studies and external CI execution remain open. No package was published as part of these changes.
