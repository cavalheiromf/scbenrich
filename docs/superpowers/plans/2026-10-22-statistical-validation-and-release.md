# SCB-Enrich Statistical Validation and Release Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking. Execute inline unless delegation is explicitly authorized.

**Goal:** Corrigir as pendências inferenciais e operacionais identificadas na revisão, produzir evidência de validação proporcional às alegações e preparar um release exploratório reproduzível.

**Architecture:** Preservar os testes por família e a interface básica de entrada. Separar decisão estatística, suporte da randomização e diagnóstico de covariáveis; centralizar validação e proveniência em módulos pequenos. A validação deve distinguir hipóteses marginais e condicionadas antes de comparar métodos.

**Tech Stack:** Python >=3.8, NumPy, pandas, SciPy >=1.10, pytest; ambiente `.venv` existente. Matplotlib como dependência opcional de benchmark para figuras.

## Global Constraints

- Hipótese principal: sobre-representação condicionada aos estratos especificados, sob permutabilidade dentro deles.
- Fisher e OR condicional às margens continuam sendo análises marginais complementares.
- Não interpretar perda de significância como artefato, causalidade ou ausência de efeito.
- Preservar as alterações locais existentes; `.ll` não pertence a este trabalho.
- Não executar `git add .`; selecionar os arquivos de cada tarefa.
- Não usar o resultado do benchmark para escolher retrospectivamente o default mais favorável.
- Não considerar os 16 testes atuais evidência de calibração geral.
- Não prometer equivalência de resultados entre versões de dependências apenas por usar a mesma semente.
- Mudanças de colunas e nomes de arquivos precisam de nota de migração.
- Commit, publicação e instalação global não fazem parte da criação deste plano.

## Situação inicial e escopo

Já existem BH para p-valores estratificados, tratamento parcial de SMD indefinido, relatório por estrato, OR via SciPy, normalização GO simples, mapa de redução, manifesto parcial e benchmark piloto. Esses componentes devem ser corrigidos e ampliados, sem reimplementação desnecessária.

O plano é dividido em entregas independentemente revisáveis:

1. Contratos inferenciais, diagnósticos e segurança das entradas.
2. GO, saídas e reprodução.
3. Validação metodológica e preparação do release.

## Mapa de arquivos

| Arquivo | Responsabilidade |
|---|---|
| `scbenrich/core/inference.py` (novo) | Campos independentes de decisão e estado. |
| `scbenrich/utils/validation.py` (novo) | Validação compartilhada pela CLI e API. |
| `scbenrich/core/stratified.py` | Amostragem, suporte por termo e diagnósticos das amostras nulas. |
| `scbenrich/core/stats.py` | Ajustes múltiplos e interface numérica da OR. |
| `scbenrich/core/obo.py` | Parsing, normalização de IDs e detecção de ciclos. |
| `scbenrich/core/dag.py` | Resumo de apresentação e mapa de remoções. |
| `scbenrich/utils/provenance.py` (novo) | Manifesto, hashes e identificação do código. |
| `scbenrich/cli.py` | Orquestração e exportação, sem duplicar validações. |
| `benchmarks/simulate.py` | Geradores com verdades marginais e condicionadas explícitas. |
| `benchmarks/report.py` (novo) | Resumos, incerteza e figuras a partir dos resultados brutos. |
| `tests/test_inference.py`, `tests/test_validation.py`, `tests/test_ontology.py`, `tests/test_provenance.py` (novos) | Contratos e regressões dos novos componentes. |
| `tests/test_benchmark.py` | Verdade simulada, métricas e reprodução. |
| `README.md`, `docs/methods.md`, `docs/migration.md` | Método, pressupostos, limitações e migração. |

## Tarefa 1 — Fixar a interpretação e o contrato das decisões

**Dependências:** nenhuma.

**Arquivos:** criar `scbenrich/core/inference.py`, `tests/test_inference.py`, `docs/methods.md`, `docs/migration.md`; modificar `scbenrich/cli.py`, `scbenrich/__init__.py`, `README.md`.

**Interface:** `decision_fields(q_marginal, q_stratified, analysis_status, balance_status, alpha=0.05) -> dict`.

Retornar `marginal_significant`, `stratified_significant`, `analysis_status`, `balance_status`. Estados da análise: `ok`, `no_randomization`, `constant_statistic`, `not_run`. Estados do balanceamento: `within_threshold`, `above_threshold`, `undefined`; esses campos não filtram a lista de descobertas estratificadas.

- [ ] Escrever o teste de independência entre significância e diagnóstico:

```python
result = decision_fields(0.8, 0.01, "ok", "undefined")
assert result["stratified_significant"] is True
assert result["marginal_significant"] is False
assert result["balance_status"] == "undefined"
assert not decision_fields(0.01, 0.01, "no_randomization", "undefined")["stratified_significant"]
```

- [ ] Executar `.venv/bin/python -m pytest tests/test_inference.py -q`; confirmar falha por interface ausente.
- [ ] Implementar a decisão com comparações em precisão completa: `q < alpha`; p-valor ausente nunca produz descoberta. Validar `0 < alpha < 1`.
- [ ] Substituir o uso decisório de `Evidence_Class`; documentar a remoção dessa coluna em `docs/migration.md`.
- [ ] Renomear a descrição da OR para `Marginal_OR` e manter os ICs explicitamente marginais; registrar a migração das colunas antigas.
- [ ] Remover das docstrings as promessas de controle explícito de todas as covariáveis e de nulos rigorosos sem pressupostos.
- [ ] Executar os testes da tarefa e os testes existentes da CLI; atualizar expectativas apenas para as mudanças intencionais.

**Aceite:** uma descoberta estratificada pode coexistir com diagnóstico indefinido sem contradição entre tabela completa e lista exportada; não há rótulo de confirmação ou artefato.

## Tarefa 2 — Validação compartilhada e contratos da API

**Dependência:** tarefa 1 para os estados exportados.

**Arquivos:** criar `scbenrich/utils/validation.py`, `tests/test_validation.py`; modificar `scbenrich/utils/io.py`, `scbenrich/core/stratified.py`, `scbenrich/cli.py`.

**Interfaces:** `validate_family_inputs(covariates, foreground_ids, term_matrix) -> None`; `validate_parameters(permutations, bins, alpha) -> None`.

- [ ] Criar testes parametrizados para IDs duplicados, nulos ou somente espaços, matriz não binária, IDs extras/ausentes na matriz, foreground vazio/completo, covariável não numérica/não finita, bins <2 e contagem negativa ou fracionária de permutações.

```python
with pytest.raises(ValueError, match="binary"):
    validate_family_inputs(cov, {"A"}, pd.DataFrame({"T": [1, 2]}, index=["A", "B"]))
```

- [ ] Rodar `.venv/bin/python -m pytest tests/test_validation.py -q` e verificar as falhas esperadas.
- [ ] Ler IDs como strings, preservando zeros iniciais; rejeitar valores vazios após inspeção de espaços. Não transformar silenciosamente IDs distintos em um único ID.
- [ ] Na API, alinhar a ordem da matriz por índice somente depois de validar igualdade dos conjuntos de IDs.
- [ ] No foreground por coluna, rejeitar ausências e tokens desconhecidos; definir explicitamente os valores falsos aceitos (`0`, `false`, `no`, `not_expanded`).
- [ ] Anotações fora do universo: rejeitar por padrão com contagem e exemplos; não descartá-las silenciosamente.
- [ ] Sem termos elegíveis: exportar tabelas vazias com cabeçalhos e estado `no_testable_terms`, sem erro de indexação.
- [ ] Executar testes novos e suíte completa.

**Aceite:** CLI e API rejeitam as mesmas entradas inválidas, preservam IDs e produzem resultado legível quando não há hipóteses testáveis.

## Tarefa 3 — Suporte por termo, bins e incerteza Monte Carlo

**Dependências:** tarefas 1–2.

**Arquivos:** modificar `scbenrich/core/stratified.py`, `scbenrich/cli.py`, `tests/test_delivery_one.py`; criar `tests/test_randomization.py`.

**Interface nova:** `run_analysis(foreground_ids, term_matrix, n_permutations=10000) -> dict`, contendo `terms` (DataFrame), `strata` (DataFrame), `balance` (DataFrame). Manter `run_permutations` como adaptador que devolve o dicionário antigo.

- [ ] Testar um universo de quatro famílias e dois selecionados contra as seis seleções enumeráveis; conferir estatística observada, limites possíveis e amostras sem duplicação.
- [ ] Testar um termo constante em todos os subconjuntos apesar de haver randomização global.
- [ ] Por estrato, com `N` famílias, `K` anotadas e `k` selecionadas, somar os limites `max(0, k-(N-K))` e `min(k,K)` para obter suporte da contagem por termo.
- [ ] Marcar `constant_statistic` quando limites coincidirem; retornar p=1 e nenhuma descoberta. Registrar `no_randomization` quando todas as seleções forem fixas.
- [ ] Exportar excedências, B realizado, `(b+1)/(B+1)` e resolução `1/(B+1)`.
- [ ] Calcular IC binomial para a probabilidade de excedência Monte Carlo usando `scipy.stats.binomtest`; nomear como incerteza da estimativa da cauda, não IC do efeito biológico.
- [ ] Expor `--bins` na CLI e preservar o default atual quando omitido. Registrar limites efetivos dos bins, incluindo colapso de limites repetidos.
- [ ] Exportar, por covariável, média e quantis do foreground observado e da distribuição das médias das réplicas; manter SMD por termo apenas como descrição separada.
- [ ] Testar determinismo com a mesma semente e avanço de estado ao reutilizar a instância.
- [ ] Rodar `.venv/bin/python -m pytest tests/test_randomization.py tests/test_delivery_one.py -q`.

**Aceite:** o usuário distingue falta de suporte, pouca precisão Monte Carlo e ausência de significância; a definição do nulo e a regra de parada permanecem fixas.

## Tarefa 4 — GO sem substituições ambíguas ou recursão cíclica

**Dependência:** tarefa 2.

**Arquivos:** modificar `scbenrich/core/obo.py`, `scbenrich/cli.py`; criar `tests/test_ontology.py`.

**Interface:** `normalize_go_id(term_id, terms, alt_map, replacements) -> str`; `replacements` passa a ser `dict[str, set[str]]` conceitualmente, usando anotações compatíveis com Python 3.8.

- [ ] Criar pequenos OBOs sintéticos com dois `replaced_by`, cadeia de substituição, destino inexistente, ciclo e múltiplos pais válidos.

```python
with pytest.raises(ValueError, match="ambiguous"):
    normalize_go_id("GO:0000001", terms, {}, {"GO:0000001": {"GO:0000002", "GO:0000003"}})
```

- [ ] Rodar `.venv/bin/python -m pytest tests/test_ontology.py -q` e confirmar falhas.
- [ ] Preservar todas as substituições no parser. Resolver cadeias iterativamente com conjunto de visitados; rejeitar ciclo, múltiplos destinos e destino ausente.
- [ ] Para ancestrais, usar estados não visitado/em processamento/concluído; uma aresta para nó em processamento deve produzir erro que identifique o ciclo.
- [ ] Distinguir GO desconhecido de termo não GO; com OBO, rejeitar GO desconhecido e preservar termos não GO sem propagação.
- [ ] Exportar caminho completo de normalização e motivo de cada transformação. Atualizar nomes e categorias com o ID canônico.
- [ ] Manter somente `is_a` e `part_of`; testar explicitamente que `has_part` e `regulates` não propagam.
- [ ] Rodar testes de ontologia e integração da CLI.

**Aceite:** nenhuma ambiguidade é resolvida pela ordem das linhas do OBO; grafos cíclicos falham com mensagem identificável.

## Tarefa 5 — Resumos e proveniência sem mistura de execuções

**Dependências:** tarefas 1–4.

**Arquivos:** criar `scbenrich/utils/provenance.py`, `tests/test_provenance.py`; modificar `scbenrich/cli.py`, `scbenrich/core/dag.py`, `README.md`.

**Interfaces:** `build_manifest(args, inputs, outputs, status) -> dict`; `write_manifest(path, manifest) -> None`.

- [ ] Testar que duas colunas diferentes de foreground geram manifestos diferentes mesmo com os mesmos hashes de entrada.
- [ ] Testar que uma pasta com artefatos de outra execução é rejeitada, preservando todos os arquivos existentes.
- [ ] Criar diretório exclusivo por execução ou exigir pasta vazia. Nunca limpar silenciosamente uma pasta fornecida pelo usuário.
- [ ] Registrar argumentos normalizados, foreground (arquivo ou coluna), hashes, versões, gerador RNG, seed, B, bins, método de ajuste, commit, indicador de mudanças locais e hashes dos módulos utilizados.
- [ ] Tratar ausência de Git como estado explícito; versão/hash do código continua obrigatório.
- [ ] Gravar manifesto inicialmente com `running`, depois `completed` ou `failed`; em falha registrar tipo e mensagem, sem substituir o erro original.
- [ ] Exportar resumos GO marginais e estratificados em arquivos diferentes. A redução deve usar exatamente a lista declarada no nome do arquivo.
- [ ] Manter mapas separados para cada resumo. Informar regra, limiar e representante de cada remoção.
- [ ] Testar termo significativo apenas na análise estratificada: deve estar no resumo estratificado e ausente no marginal.
- [ ] Rodar `.venv/bin/python -m pytest tests/test_provenance.py tests/test_delivery_two.py -q`.

**Aceite:** cada arquivo pertence a uma execução identificável; a seleção usada pode ser reconstruída; resumos não omitem silenciosamente descobertas da análise que representam.

## Tarefa 6 — Corrigir a verdade de referência do benchmark

**Dependências:** tarefa 1; executável antes das demais tarefas operacionais.

**Arquivos:** modificar `benchmarks/simulate.py`, `tests/test_benchmark.py`; criar `docs/benchmark-protocol.md`; atualizar `README.md`.

**Contrato:** cada termo simulado possui `marginal_null` e `conditional_null`. As duas verdades não são inferidas do p-valor obtido.

- [ ] Criar teste determinístico para o cenário em que tamanho aumenta seleção e anotação: a covariância marginal esperada é positiva, embora a associação condicional adicional seja zero.
- [ ] Rotular os resultados de Fisher contra o nulo condicionado como `conditional_target_mismatch`, e não como erro tipo I de Fisher.
- [ ] Manter métricas de erro de cada método contra seu próprio alvo; comparar poder diretamente apenas nos cenários de alvos coincidentes.
- [ ] Preservar o piloto original como artefato histórico, acompanhado de nota de interpretação; não sobrescrever resultados antigos com novas definições.
- [ ] Renomear `rejection_rate` para `probability_any_false_discovery`. Em nulo global ela coincide com FDR; em cenário misto não coincide.
- [ ] Especificar em `docs/benchmark-protocol.md` cenários, sementes, bins, métricas e regras de análise antes da execução ampliada.
- [ ] Rodar `.venv/bin/python -m pytest tests/test_benchmark.py -q`.

**Aceite:** nenhuma conclusão de superioridade depende de chamar associação marginal verdadeira de falso positivo marginal.

## Tarefa 7 — Ampliar simulações e incerteza

**Dependências:** tarefas 3 e 6.

**Arquivos:** modificar `benchmarks/simulate.py`, `tests/test_benchmark.py`; criar `benchmarks/report.py`, `benchmarks/configs/validation.json`, `tests/test_benchmark_report.py`.

**Saídas:** resultados por termo, por réplica, resumo por cenário, configuração, hashes do código e versões. Todas as unidades de erro devem estar identificadas.

- [ ] Acrescentar tamanho, homólogos de referência e grau anotacional. Incluir variantes em que grau é externo e em que é calculado da matriz testada.
- [ ] Incluir cinco famílias de cenários: permutabilidade exata; confundimento contínuo residual; anotação incompleta dependente do tamanho; termos dependentes por estrutura ancestral; estratos esparsos e seleções quase determinísticas.
- [ ] Fixar grids iniciais de bins `[2, 4, 8]`, efeitos OR `[1, 2, 4]`, taxas de ausência `[0, 0.3, 0.6]`; executar primeiro combinações predefinidas que isolem um fator por vez, evitando expansão combinatória sem objetivo.
- [ ] Executar um smoke test com 5 réplicas e 99 permutações para cada cenário, verificando arquivos, estados e reprodução.
- [ ] Executar validação com 1.000 réplicas por cenário selecionado e 9.999 permutações, em lotes retomáveis de 50 réplicas; registrar tempo e falhas de cada lote.
- [ ] Registrar sementes independentes por cenário, réplica e método via `SeedSequence`; adicionar um cenário não deve alterar as simulações anteriores.
- [ ] Reportar FDR como média da proporção de falsas descobertas por réplica; usar bootstrap por réplica para IC de FDR e diferenças pareadas. Para poder binário e cobertura, usar IC binomial quando a unidade for uma réplica independente.
- [ ] Gerar curvas de calibração por termo nulo, gráficos de poder/FDR por bins e cobertura de ICs. Não tratar milhares de termos da mesma réplica como observações independentes.
- [ ] Calcular cobertura da OR somente quando seu alvo marginal é conhecido; não comparar IC marginal com OR condicional ao tamanho.
- [ ] Não trocar defaults após observar qual configuração teve melhor resultado; resultados fora do nominal motivam investigação e nova validação com protocolo separado.

**Aceite:** relatório reproduzível com incerteza, pressupostos e situações de falha; ausência de diferenças significativas não é prova de equivalência.

## Tarefa 8 — Comparadores externos e estudos de caso

**Dependência:** tarefa 6; análises ampliadas dependem da tarefa 7.

**Arquivos:** criar `benchmarks/comparators/README.md` e adaptadores específicos apenas para ferramentas realmente instaladas e fixadas por versão.

- [ ] Definir duas comparações: unidade gene e unidade família, sem misturá-las na mesma estimativa de erro.
- [ ] Exportar matrizes e universos idênticos para Fisher, GOATOOLS, clusterProfiler e algoritmos topGO explicitamente identificados.
- [ ] Registrar propagação, filtros, relações GO, alternativas dos testes e ajuste múltiplo de cada comparador.
- [ ] Confirmar que dois implementadores de Fisher recuperam os mesmos resultados no mesmo problema; diferenças devem ser investigadas antes de comparar métodos distintos.
- [ ] Acrescentar estudos de caso com quinases, P450s e NLRs somente com dados e versões recuperáveis; tratá-los como plausibilidade biológica, não verdade automática de calibração.
- [ ] Caso os dados reais não estejam disponíveis, marcar a seção como não executada e concluir separadamente a validação simulada; não inventar controles.

**Aceite:** comparação deixa explícita a hipótese de cada ferramenta e permite reproduzir seu universo e suas anotações.

## Tarefa 9 — Revisão integrada e release exploratório

**Dependências:** tarefas 1–6 para release exploratório; tarefas 7–8 para alegações metodológicas mais amplas.

**Arquivos:** `README.md`, `docs/methods.md`, `docs/migration.md`, `pyproject.toml`, `.github/workflows/tests.yml` (novo).

- [ ] Executar `.venv/bin/python -m pytest -q` e `git diff --check`; exigir sucesso.
- [ ] Executar CLI com os exemplos, sem permutações, com permutações, sem termos elegíveis e com erro GO proposital; conferir estados e arquivos de cada execução.
- [ ] Testar instalação limpa e Python mínimo anunciado, além de uma versão moderna; atualizar o suporte anunciado se a matriz não puder ser atendida.
- [ ] Revisar todo texto que contenha `rigorous`, `robust`, `confirmatory`, `artifact` e equivalentes; cada afirmação deve apontar para evidência ou ser reformulada.
- [ ] Revisar resultados do benchmark separadamente da aprovação dos testes de software.
- [ ] Preparar notas de migração para colunas, arquivos, regras de entrada e substituições GO.
- [ ] Criar commits por entrega revisada, selecionando arquivos explicitamente. Publicação de pacote é uma ação posterior e separada.

**Aceite:** release descrito como exploratório, sem alegação de controle geral de FDR; limitações conhecidas não ficam ocultas por testes aprovados.

## Ordem de execução e pontos de revisão

1. Tarefas 1–3: fechar contrato inferencial e diagnósticos.
2. Tarefas 4–5: fechar GO e rastreabilidade.
3. Tarefa 6: corrigir a interpretação do piloto antes de ampliar suas conclusões.
4. Tarefa 7: executar validação em lotes com atualizações de progresso.
5. Tarefa 8: executar comparações com dados e ambientes disponíveis.
6. Tarefa 9: revisar o conjunto e preparar o release.

A criação deste documento não implementa essas tarefas. Cada checkbox só deve ser marcado após código, testes e documentação correspondentes estarem verificados.

## Execution checkpoint

First implementation block verified: 23 tests passed. Decision fields, marginal OR names, shared API validation, ID string preservation, foreground token validation, term support, Monte Carlo intervals, null-mean diagnostics and CLI bins are implemented. Tasks 1–3 remain partially open: empty-term output, exact bin-boundary export, broader parameterized validation and distribution diagnostics require completion. Tasks 4–9 have not been executed under this plan. No release or calibration claim is complete.

## Second execution checkpoint

Implemented: empty-term output, actual bin-edge export, GO replacement sets/chains/ambiguity rejection, cycle detection, canonical GO names, complete mapping paths, empty output-directory enforcement, running/failed manifests, normalized arguments, Git state and source hashes, and separate marginal/stratified summaries. Tests cover these contracts; remaining checklist items are not automatically marked complete. Tasks 6–9 (benchmark interpretation, expanded validation, external comparators and release) remain open. Distribution diagnostics currently summarize means; quantiles of foreground values and wider validation tests remain to expand.

## Benchmark interpretation checkpoint

Task 6 implemented: population truth exported per term for marginal and continuous-size-conditional targets; method-specific false-discovery metrics, explicit conditional_target_mismatch, renamed probability_any_false_discovery, historical pilot interpretation preserved, and benchmark protocol documented. Tests verify positive population covariance under conditional null and method-specific error counting. Expanded simulations and external comparator work remain open.

## Expanded validation checkpoint

Task 7 infrastructure implemented in benchmarks/expanded.py and benchmarks/report.py: seven scenarios, three covariates, unknown truth preserved, independent SeedSequence streams, atomic restartable batches, source/version validation, term-level outputs, replicate bootstrap FDR intervals, binomial power/coverage intervals and calibration figures. Executed smoke configuration: five replicates, 99 permutations, OR 1/4, bins 2/4. Summary artifacts are preserved in benchmarks/expanded_smoke. Full 1,000-replicate/9,999-permutation configuration is supplied but not executed. Coverage assertions use the unconfounded marginal target only. Paired difference intervals and fuller distribution diagnostics remain open; these tasks are not marked as fully complete. Suite: 35 tests passed.

## External comparator checkpoint

Task 8 family-level adapters implemented and executed with installed GOATOOLS 1.6.4, topGO 2.62.0 classic/elim and clusterProfiler 4.18.4. Exported shared matrix/universe and audited counts, denominators, alternatives, p-values and harmonized BH. GOATOOLS two-sided, topGO classic greater and clusterProfiler greater agreed with matched Fisher references on mock data; elim differed as expected. clusterProfiler zero-overlap omissions are reconstructed as upper-tail p=1 with explicit source marking. Results and versions are preserved in benchmarks/comparators/example_results. Gene-level comparators and biological case studies remain open. Suite: 37 tests passed. No release or completed large-scale validation is claimed.

## Exploratory release checkpoint

Task 9 preparation implemented: Python 3.8/3.12 CI workflow, wheel build and installation outside the checkout on Python 3.12, four subprocess CLI cases, corrected output names and exploratory wording, CHANGELOG and release-check limits. The local wheel environment inherited scientific dependencies; fully isolated resolution and remote CI remain unverified. Python 3.8 syntax compilation passed; numerical execution is pending dependencies/CI. Suite: 41 tests passed; git diff --check passed. No release publication or version bump occurred. Full simulation configuration, gene-level comparisons, biological studies and other previously listed validation gaps remain open.
