# Comparação de implementações com entradas harmonizadas

Executar na raiz do repositório:

```bash
python -m benchmarks.comparators.run \
  --master-table examples/mock_master_table.tsv \
  --annotations examples/mock_annotations.tsv \
  --foreground is_expanded --obo examples/mock_go.obo \
  --out-dir benchmark_results/comparators
```

GOATOOLS é opcional; os comparadores R exigem Rscript, topGO, GO.db e clusterProfiler instalados. A ausência de ferramenta é registrada. O programa não instala pacotes automaticamente. Uma falha R é preservada em R.log e nos estados do manifesto; resultados incompletos não representam comparações bem-sucedidas.

A unidade é a família. O exportador grava universo, foreground, matriz binária, associações família–termo, referência e hashes. A propagação is_a/part_of é feita pela SCB-Enrich antes da exportação, seguida dos mesmos filtros de tamanho para todos os adaptadores. Cada adaptador grava contagens e denominadores efetivos. topGO utiliza GO.db e pode propagar segundo uma versão diferente da ontologia: versões são registradas, e diferenças de contagem invalidam a interpretação como teste equivalente. O adaptador topGO atual é restrito a BP.

## Alternativas e ajustes

- GOATOOLS: Fisher bilateral nativo, sem propagação adicional; comparar com a referência bilateral.
- topGO classic: Fisher unilateral de sobre-representação; comparar com a referência unilateral.
- topGO elim: Fisher com eliminação no grafo. As contagens exportadas são as contagens descritivas originais; os p-valores refletem o algoritmo. Igualdade com Fisher clássico não é esperada.
- clusterProfiler: enricher com TERM2GENE, minGSSize=1, maxGSSize=N e cutoffs=1. Não usa busca de anotações num organismo remoto. Termos sem sobreposição com o foreground são omitidos pelo programa; o adaptador os adiciona com p=1, marcado como source=reconstructed_zero_overlap, para o teste unilateral. Outras omissões permanecem explícitas.

Os arquivos audit_* contêm correspondência de contagens, p-valores e BH recalculado sobre os mesmos termos retidos. Ajustes nativos entre namespaces ou após filtragens não são apresentados como equivalentes. Igualdade de p não é aceita se contagens/denominadores diferirem. A coluna q_harmonized não é calculada se há termos ausentes; para elim, o ajuste é apenas um resumo comum e não estabelece controle de FDR.

O universo de clusterProfiler e o conjunto viável de topGO podem excluir famílias sem anotação. Os arquivos *_effective_universe.txt tornam essa exclusão auditável; os denominadores são comparados com a referência. Não forçar concordância removendo essas famílias silenciosamente.

## Resultado executado

Os exemplos fornecidos têm 20 famílias e três termos BP retidos. Na execução preservada em example_results, GOATOOLS 1.6.4, topGO classic 2.62.0 e clusterProfiler 4.18.4 concordaram com suas referências em contagens e p-valores. topGO elim apresentou diferenças de p esperadas pelo algoritmo. GO.db 3.22.0 e R 4.5.1 foram registrados.

Isso verifica integração e equivalência dos testes de base neste exemplo. Não demonstra superioridade da estratificação, calibração em dados reais ou equivalência sob outras ontologias. Comparações por gene e estudos de caso reais ainda precisam de dados e protocolos próprios; não são misturados a esta análise familiar.
