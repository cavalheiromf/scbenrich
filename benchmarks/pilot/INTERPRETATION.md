# Interpretação do piloto histórico

Os arquivos config.json, replicates.tsv e summary.tsv foram preservados. Usam o esquema antigo de métricas.

Nos cenários confundidos, os termos null_* têm associação marginal esperada com o foreground e nenhum efeito adicional condicionado ao tamanho. Portanto, os valores históricos false_discoveries, false_discovery_proportion e empirical_fdr de Fisher nesses cenários quantificam divergência em relação ao alvo condicionado. Não demonstram erro tipo I ou inflação de FDR para a hipótese marginal do Fisher.

Os 60% em confounded_null indicam que Fisher detectou pelo menos uma dessas associações marginais em 18 de 30 réplicas. Não autorizam afirmar que o teste marginal teve 60% de falsos positivos.

rejection_rate é a frequência de ao menos uma descoberta definida como falsa pelo esquema antigo. Seus intervalos binomiais não são intervalos para FDR no cenário misto. Para estratificação, o alvo continua condicionado ao tamanho contínuo, enquanto a implementação usa bins e pode reter confundimento residual.

Resultados novos usam metric_schema=2, verdades marginais e condicionadas explícitas e probability_any_false_discovery. Não misturar resumos dos dois esquemas.
