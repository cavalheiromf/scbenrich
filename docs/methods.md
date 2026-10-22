# Inference contracts

The unit of analysis is the family. Fisher's test and the odds ratio conditional on the margins assess marginal association; the permutation test assesses the annotated count under exchangeable selection within strata. BH is applied separately to each set of p-values.

`marginal_significant` and `stratified_significant` are independent decisions. `balance_status` is descriptive and does not select discoveries. `analysis_status` distinguishes `ok`, `constant_statistic`, `no_randomization` and `not_run`.

The support of each term's count is the sum of the hypergeometric bounds across strata. Many possible selections do not guarantee that a term's statistic varies. A constant statistic does not produce a discovery.

`tail_ci_low` and `tail_ci_high` are binomial intervals for the Monte Carlo exceedance probability. They are not intervals for the biological effect or the adjusted p-value. The quantiles in `balance_diagnostics.tsv` describe the distribution of covariate means across resampled sets, without guaranteeing the absence of confounding.
