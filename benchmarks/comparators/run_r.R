args <- commandArgs(trailingOnly=TRUE)
out <- args[[1]]
universe <- readLines(file.path(out, 'universe.txt'))
foreground <- readLines(file.path(out, 'foreground.txt'))
pairs <- read.delim(file.path(out, 'term2family.tsv'), stringsAsFactors=FALSE)
reference <- read.delim(file.path(out, 'reference.tsv'), stringsAsFactors=FALSE)
versions <- data.frame(package='R', version=as.character(getRversion()))
for (pkg in c('topGO', 'GO.db', 'clusterProfiler', 'DOSE')) {
  versions <- rbind(versions, data.frame(package=pkg, version=if (requireNamespace(pkg, quietly=TRUE)) as.character(packageVersion(pkg)) else 'unavailable'))
}
write.table(versions, file.path(out, 'r_versions.tsv'), sep='\t', row.names=FALSE, quote=FALSE)
if (requireNamespace('topGO', quietly=TRUE)) {
  suppressPackageStartupMessages(library(topGO))
  mapping <- split(pairs$Term_ID, pairs$Family_ID)
  selected <- factor(as.integer(universe %in% foreground), levels=c(0,1))
  names(selected) <- universe
  go <- methods::new('topGOdata', ontology='BP', allGenes=selected,
                     annot=topGO::annFUN.gene2GO, gene2GO=mapping, nodeSize=1L)
  effective <- go@allGenes[go@feasible]
  writeLines(effective, file.path(out, 'topgo_effective_universe.txt'))
  counts <- topGO::genesInTerm(go, reference$Term_ID)
  for (algorithm in c('classic', 'elim')) {
    result <- topGO::runTest(go, algorithm=algorithm, statistic='fisher')
    scores <- topGO::score(result)
    available <- intersect(reference$Term_ID, names(scores))
    frame <- data.frame(Term_ID=available,
      foreground_count=vapply(available, function(t) length(intersect(counts[[t]], foreground)), integer(1)),
      universe_count=vapply(available, function(t) length(counts[[t]]), integer(1)),
      foreground_total=length(intersect(effective, foreground)), universe_total=length(effective),
      p_value=unname(scores[available]))
    write.table(frame, file.path(out, paste0('topgo_',algorithm,'.tsv')), sep='\t', row.names=FALSE, quote=FALSE)
  }
}
if (requireNamespace('clusterProfiler', quietly=TRUE)) {
  effective <- intersect(universe, unique(pairs$Family_ID))
  writeLines(effective, file.path(out, 'clusterprofiler_effective_universe.txt'))
  # enricher uses TERM2GENE rather than a live organism database.
  result <- clusterProfiler::enricher(foreground, universe=universe,
    TERM2GENE=pairs[, c('Term_ID','Family_ID')], pvalueCutoff=1, qvalueCutoff=1,
    pAdjustMethod='BH', minGSSize=1, maxGSSize=length(universe))
  frame <- as.data.frame(result)
  available <- intersect(reference$Term_ID, frame$ID)
  frame <- frame[match(available, frame$ID), , drop=FALSE]
  sizes <- table(pairs$Term_ID)
  exported <- data.frame(Term_ID=available, foreground_count=frame$Count,
    universe_count=as.integer(sizes[available]), foreground_total=length(intersect(effective, foreground)),
    universe_total=length(effective), p_value=frame$pvalue, source='native')
  # enricher omits terms with zero study overlap. Their upper-tail p is exactly 1.
  omitted <- setdiff(reference$Term_ID, available)
  for (term in omitted) {
    members <- pairs$Family_ID[pairs$Term_ID == term]
    if (length(intersect(members, foreground)) == 0) {
      exported <- rbind(exported, data.frame(Term_ID=term, foreground_count=0L,
        universe_count=length(unique(members)), foreground_total=length(intersect(effective, foreground)),
        universe_total=length(effective), p_value=1, source='reconstructed_zero_overlap'))
    }
  }
  write.table(exported, file.path(out, 'clusterprofiler.tsv'), sep='\t', row.names=FALSE, quote=FALSE)
}
