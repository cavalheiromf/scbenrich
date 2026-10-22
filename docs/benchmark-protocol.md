# Benchmark: hipóteses, métricas e interpretação

## Verdade de referência — esquema 2

O gerador define a verdade pela distribuição populacional, sem consultar os p-valores observados. `term_truth.tsv` exporta `marginal_null` e `conditional_null` para cada termo e cenário. O alvo condicionado é o tamanho contínuo; o teste usa bins e pode apresentar confundimento residual dentro deles.

| Cenário | Termos `null_*`: nulo marginal | Termos `null_*`: nulo condicionado | Termo `effect` |
|---|---|---|---|
| simple_null | verdadeiro | verdadeiro | ausente |
| confounded_null | falso | verdadeiro | ausente |
| simple_positive | verdadeiro | verdadeiro | alternativa marginal e condicionada |
| confounded_positive | falso | verdadeiro | alternativa marginal e condicionada |

O nome histórico `confounded_null` significa nulo condicionado. Seleção e anotação aumentam com o tamanho, gerando covariância marginal positiva. Uma descoberta de Fisher nesses termos pode representar associação marginal real.

## Métricas

- `target`: hipótese marginal para Fisher; condicionada ao tamanho para estratificação.
- `false_discoveries`: descobertas de termos nulos para o alvo do método.
- `n_null_terms`: denominador conceitual que identifica cenários sem hipóteses nulas para aquele alvo. FDR zero quando não há nulos é trivial e não demonstra calibração.
- `conditional_target_mismatch`: quantidade de termos condicionadamente nulos descobertos por Fisher nos cenários confundidos. É uma divergência de objetivo, não erro tipo I marginal.
- `probability_any_false_discovery`: frequência de réplicas com pelo menos uma falsa descoberta contra o alvo correspondente. Coincide com FDR apenas sob nulo global.
- `empirical_fdr`: média da proporção de falsas descobertas por réplica, considerando zero quando não há descobertas.
- `power`: frequência de detecção do termo `effect`, quando presente. Comparação direta de poder entre métodos exige alvos coincidentes; nos cenários confundidos, ler como detecção do termo predefinido para hipóteses diferentes.
- `any_false_discovery_ci_low/high`: IC binomial para o evento de ao menos uma falsa descoberta; não é IC do FDR nos cenários mistos.
- `ci_coverage`: cobertura da OR marginal conhecida apenas no cenário simple_positive, sem seleção prévia por significância.

## Protocolo do benchmark atual

Configuração padrão: seed 42, 100 réplicas, 999 permutações, bins 2/4/8, 120 famílias, oito termos `null_*`, log-OR do efeito 1.5. Os cenários são avaliados sobre a mesma matriz de termos para os dois métodos. BH e limiar estrito de 0.05 são aplicados separadamente. O binning não é escolhido retrospectivamente pelo melhor resultado.

Esse experimento não modela ainda três covariáveis, anotações GO incompletas ou estrutura filogenética. Também não fornece IC de FDR em cenário misto. Esses componentes pertencem à validação ampliada. A resolução e os intervalos de incerteza devem acompanhar qualquer conclusão.

## Piloto histórico

`benchmarks/pilot` é preservado sem alteração dos arquivos numéricos. Seus rótulos anteriores mediam descobertas contra o nulo condicionado para ambos os métodos. Nos cenários confundidos, as colunas de falsas descobertas de Fisher não podem ser interpretadas como erros de seu teste marginal. Ver `benchmarks/pilot/INTERPRETATION.md`.

## Validação ampliada

`benchmarks.expanded` acrescenta sete cenários e três covariáveis. Nos cenários discretos, valores são balanceados para que bins de quantis separem exatamente cada covariável binária. No cenário contínuo, os bins aproximam o condicionamento. `missing_size` aplica perda de anotação dependente do tamanho; `dependent` cria dois filhos e um ancestral com união dos filhos; `sparse` usa seleção quase determinística; `derived_degree` calcula o grau da própria matriz. `unconfounded` oferece OR marginal conhecida para medir cobertura.

A verdade condicionada de `derived_degree` é desconhecida: condicionar no grau calculado do desfecho pode mudar a hipótese. A verdade marginal de `missing_size` também é desconhecida, pois perda diferencial pode reverter ou cancelar a associação. Essas células não recebem FDR. OR=1 representa nulo; OR>1 representa efeito positivo. A taxa de detecção de efeito não deve ser comparada entre métodos com alvos diferentes. FDR zero em cenários sem termos nulos marginais é trivial.

Execução curta para verificar funcionamento:

```bash
python -m benchmarks.expanded --config benchmarks/configs/smoke.json --out-dir benchmark_results/smoke
python -m benchmarks.report benchmark_results/smoke --plot
```

Para a configuração ampliada:

```bash
python -m benchmarks.expanded --config benchmarks/configs/validation.json --out-dir benchmark_results/validation
```

A configuração ampliada define 1.000 réplicas, 9.999 permutações, 240 famílias, bins 2/4/8 e OR 1/2/4. Essa configuração é computacionalmente maior e não foi executada nesta entrega. O ensaio curto executado usa cinco réplicas, 99 permutações, 120 famílias, bins 2/4 e OR 1/4; demonstra execução, não calibração.

Cada lote é gravado atomicamente. Repetir o comando reutiliza lotes concluídos se configuração, versões e hashes do código coincidirem. Mudanças exigem outra pasta. `SeedSequence` recebe cenário, efeito, réplica, propósito e bins: reordenar cenários não muda os dados. Em retomada, o lote incompleto é refeito de forma determinística. IDs de cenário são estáveis e não devem ser renumerados.

`terms.tsv` preserva p-valores, ajustes, verdade de referência, suporte e excedências. `replicates.tsv` preserva descobertas, proporção de falsas descobertas, detecção do efeito, cobertura e suporte. `summary.tsv` calcula FDR com bootstrap de réplicas e intervalos binomiais para detecção/cobertura. Células com menos de duas réplicas não recebem intervalos. Intervalos bootstrap degenerados em zero com poucas réplicas não demonstram ausência de erro; o ensaio curto não serve como validação nominal.

`calibration.tsv` registra, para cada limiar, a média entre réplicas das frações de rejeição entre termos nulos conhecidos. As figuras SVG são curvas de calibração descritivas, sem intervalos simultâneos. Discreção dos p-valores pode produzir conservadorismo. Não há ainda bootstrap pareado das diferenças entre bins, incerteza simultânea das curvas, simulação filogenética nem avaliação de políticas alternativas de evidência IEA.

## Comparadores externos executados

O runner `benchmarks.comparators.run` harmoniza entradas familiares e audita contagens antes de comparar p-valores. GOATOOLS bilateral é comparado com Fisher bilateral; topGO classic e clusterProfiler são comparados com Fisher unilateral. topGO elim é identificado como algoritmo diferente. As versões, omissões e ajustes BH comuns estão documentados em `benchmarks/comparators/README.md`. A execução com exemplos verifica integração; ainda não constitui benchmark de calibração dos comparadores ou estudo biológico real.

## Custo medido da execução completa

Uma medição de 21 análises com 9.999 permutações estimou aproximadamente 35 horas de execução sequencial após a otimização dos estratos determinísticos (73 horas antes) para as 63.000 análises da configuração completa. A extrapolação usa apenas OR=1 e uma réplica; não inclui todo o custo de geração, escrita e agregação. Consulte [tempos, suporte e reprodução](../benchmarks/timing/README.md). Os resultados da execução completa subsequente constam na seção abaixo.

## Execução completa concluída em 2026-10-22

A validação completa foi executada com a configuração estatística de `benchmarks/configs/validation_checkpointed.json` (1.000 réplicas, 9.999 permutações, 240 famílias, bins 2/4/8, Odds Ratios 1/2/4 e sete cenários).

Para viabilizar a execução dentro dos recursos do cluster Slurm:
- As análises foram distribuídas em um job array de 84 tarefas simultâneas (7 cenários × 3 Odds Ratios × 4 blocos de réplicas de 250 réplicas / 50 lotes).
- O tempo total de relógio (*wall-clock*) foi de aproximadamente **39 minutos** (reduzindo das ~35 horas da estimativa sequencial).
- Foram gerados e agregados com sucesso os **4.200 lotes** (`batch_*.json`), totalizando **63.000 análises**.

Os arquivos finais estão consolidados em `benchmark_results/validation/` (ignorados pelo Git):
- `summary.tsv`: resumo das 126 combinações com estimativas de FDR empírico, poder, cobertura e intervalos bootstrap/binomiais.
- `calibration.tsv`: taxas de rejeição observadas sob o nulo em limiares nominais para avaliação de calibração.
- `calibration_*.svg`: figuras vetoriais descritivas das curvas de calibração por cenário.
- `replicates.tsv` e `terms.tsv`: dados brutos por réplica e por termo.
