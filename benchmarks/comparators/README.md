# Comparing implementations with harmonized inputs

Run from the repository root:

```bash
python -m benchmarks.comparators.run \
  --master-table examples/mock_master_table.tsv \
  --annotations examples/mock_annotations.tsv \
  --foreground is_expanded --obo examples/mock_go.obo \
  --out-dir benchmark_results/comparators
```

GOATOOLS is optional; the R comparators require Rscript, topGO, GO.db and clusterProfiler to be installed. Missing tools are recorded. The program does not install packages automatically. An R failure is preserved in R.log and manifest status fields; incomplete results do not represent successful comparisons.

The unit of analysis is the family. The exporter writes the universe, foreground, binary matrix, family-to-term associations, reference results and hashes. SCB-Enrich performs is_a/part_of propagation before export, followed by the same size filters for all adapters. Each adapter records effective counts and denominators. topGO uses GO.db and may propagate annotations using a different ontology version: versions are recorded, and count differences invalidate interpretation as equivalent tests. The current topGO adapter is restricted to BP.

## Test alternatives and adjustments

- GOATOOLS: native two-sided Fisher test, without additional propagation; compare with the two-sided reference.
- topGO classic: one-sided Fisher overrepresentation test; compare with the one-sided reference.
- topGO elim: Fisher test with graph elimination. Exported counts are the original descriptive counts; p-values reflect the algorithm. Equality with classic Fisher is not expected.
- clusterProfiler: enricher with TERM2GENE, minGSSize=1, maxGSSize=N and cutoffs=1. It does not retrieve annotations from a remote organism database. Terms with no foreground overlap are omitted by the program; the adapter adds them with p=1, marked as source=reconstructed_zero_overlap, for the one-sided test. Other omissions remain explicit.

The audit_* files contain count comparisons, p-values and BH recalculated across the same retained terms. Native adjustments across namespaces or after filtering are not presented as equivalent. Matching p-values are not accepted if counts or denominators differ. The q_harmonized column is not calculated when terms are missing; for elim, the adjustment is only a common summary and does not establish FDR control.

The clusterProfiler universe and topGO feasible set may exclude unannotated families. The *_effective_universe.txt files make this exclusion auditable; denominators are compared with the reference. Do not force agreement by silently removing these families.

## Executed comparison

The supplied examples contain 20 families and three retained BP terms. In the run preserved in example_results, GOATOOLS 1.6.4, topGO classic 2.62.0 and clusterProfiler 4.18.4 agreed with their references on counts and p-values. topGO elim showed the p-value differences expected from its algorithm. GO.db 3.22.0 and R 4.5.1 were recorded.

This verifies integration and equivalence of the baseline tests for this example. It does not demonstrate superiority of stratification, calibration on real data or equivalence under other ontologies. Gene-level comparisons and real biological case studies still require their own data and protocols; they are not mixed with this family-level analysis.
