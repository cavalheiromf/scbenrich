# Migração da saída exploratória

- `Evidence_Class` foi substituído por `marginal_significant`, `stratified_significant`, `analysis_status` e `balance_status`.
- `Odds_Ratio` passa a ser `Marginal_OR`; os ICs são `Marginal_CI_Low_95` e `Marginal_CI_High_95`.
- IDs numéricos são lidos como texto, preservando zeros iniciais.
- API de permutação exige igualdade dos universos da matriz e das covariáveis, valores binários e foreground próprio não vazio.
- Indicadores de foreground ausentes ou desconhecidos e anotações fora do universo geram erro.
- `--bins` configura a quantidade solicitada de bins por covariável.
- `run_analysis` retorna tabelas `terms`, `strata` e `balance`; `run_permutations` preserva o retorno antigo.

- Output directories must be empty. Failed runs retain `run_metadata.json` with error details.
- GO replacements now map to sets; ambiguous destinations and cycles are errors.
- GO summaries and removal maps have separate `_marginal` and `_stratified` filenames.
- Empty hypothesis sets produce tables with headers and `no_testable_terms` metadata.

- A amostragem agora pré-calcula estratos totalmente selecionados. A distribuição do teste é preservada, mas a sequência aleatória para uma mesma semente pode mudar. Para reprodução, preserve também a versão e os hashes do código; resultados e checkpoints anteriores não devem ser misturados com os novos.
