# Interpretation of the historical pilot

The config.json, replicates.tsv and summary.tsv files have been preserved. They use the old metric schema.

In confounded scenarios, null_* terms have an expected marginal association with the foreground and no additional effect conditional on size. Therefore, the historical false_discoveries, false_discovery_proportion and empirical_fdr values for Fisher in these scenarios quantify disagreement with the conditional target. They do not demonstrate type I error or FDR inflation for Fisher's marginal hypothesis.

The 60% in confounded_null indicates that Fisher detected at least one of these marginal associations in 18 of 30 replicates. It does not justify claiming that the marginal test had 60% false positives.

rejection_rate is the frequency of at least one discovery defined as false under the old schema. Its binomial intervals are not FDR intervals in a mixture of null and alternative terms. For stratification, the target remains conditional on continuous size, whereas the implementation uses bins and may retain residual confounding.

New results use metric_schema=2, explicit marginal and conditional truths, and probability_any_false_discovery. Do not mix summaries from the two schemas.
