# Contratos inferenciais

A unidade é a família. Fisher e a OR condicional às margens avaliam associação marginal; a permutação avalia a contagem anotada sob seleção permutável dentro dos estratos. BH é aplicado separadamente a cada conjunto de p-valores.

`marginal_significant` e `stratified_significant` são decisões independentes. `balance_status` é descritivo e não seleciona descobertas. `analysis_status` distingue `ok`, `constant_statistic`, `no_randomization` e `not_run`.

O suporte por termo é a soma dos limites hipergeométricos de cada estrato. Muitas seleções possíveis não garantem que a estatística de um termo varie. Uma estatística constante não produz descoberta.

`tail_ci_low` e `tail_ci_high` são intervalos binomiais para a probabilidade de excedência Monte Carlo. Não são intervalos do efeito biológico ou do p-valor ajustado. Os quantis em `balance_diagnostics.tsv` descrevem a distribuição das médias das covariáveis nas réplicas, sem garantir ausência de confundimento.
