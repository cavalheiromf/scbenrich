# Exploratory output migration

- `Evidence_Class` has been replaced by `marginal_significant`, `stratified_significant`, `analysis_status` and `balance_status`.
- `Odds_Ratio` is now `Marginal_OR`; its confidence intervals are `Marginal_CI_Low_95` and `Marginal_CI_High_95`.
- Numeric IDs are read as text, preserving leading zeros.
- The permutation API requires matching family universes in the matrix and covariates, binary matrix values and a nonempty foreground that is a proper subset of the universe.
- Missing or unknown foreground indicators and annotations outside the universe raise errors.
- `--bins` sets the requested number of bins per covariate.
- `run_analysis` returns the `terms`, `strata` and `balance` tables; `run_permutations` preserves its previous return format.

- Output directories must be empty. Failed runs retain `run_metadata.json` with error details.
- GO replacements now map to sets; ambiguous destinations and cycles are errors.
- GO summaries and removal maps have separate `_marginal` and `_stratified` filenames.
- Empty hypothesis sets produce tables with headers and `no_testable_terms` metadata.

- Sampling now precomputes contributions from fully selected strata. The test distribution is preserved, but the random sequence for a given seed may change. For reproduction, also preserve the code version and hashes; previous results and checkpoints must not be mixed with new ones.
