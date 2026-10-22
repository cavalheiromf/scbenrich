# Functional smoke test of expanded validation

These files were produced by the smoke-test command documented in docs/benchmark-protocol.md. The run used five replicates, 99 permutations, seven scenarios, OR 1/4 and bins 2/4. These results verify execution and metrics; they do not demonstrate nominal FDR control. The manifest preserves the configuration, code hashes and versions at execution time. Batches, term-level data and figures from this run are in /tmp/scbenrich-expanded-smoke.

FDR uses bootstrap resampling by replicate; intervals equal to zero in small samples are not reliable upper bounds. Coverage uses a binomial confidence interval and is calculated only for the known marginal OR in the unconfounded scenario. Missing values indicate metrics with unknown truth or no alternative effect.
