# Custo da validação ampliada após otimização

Medição em 2026-10-22: Python 3.12.12, 240 famílias, 9.999 permutações, sete cenários e bins 2/4/8. Foram executadas 21 análises por versão, com OR=1 e réplica 0.

A estimativa sequencial para as 63.000 análises da configuração completa caiu de **73.0 para 34.7 horas**, um ganho de **2.10×**. O tempo acumulado das análises caiu de 87.6 para 41.7 segundos.

A contribuição dos estratos totalmente selecionados agora é pré-calculada; somente estratos com seleção parcial usam o gerador aleatório. A distribuição de amostragem sem reposição e o compartilhamento dos sorteios entre termos são preservados. A sequência para uma mesma semente pode mudar em relação à implementação anterior. Hashes do código diferenciam as execuções; checkpoints antigos exigem outra pasta.

| Cenário | Bins | Antes (s) | Depois (s) | Fração de estratos informativos |
|---|---:|---:|---:|---:|
| exact | 2 | 1.74 | 1.42 | 1.000 |
| exact | 4 | 1.47 | 1.25 | 1.000 |
| exact | 8 | 1.46 | 1.25 | 1.000 |
| continuous | 2 | 1.46 | 1.22 | 1.000 |
| continuous | 4 | 5.82 | 4.17 | 0.603 |
| continuous | 8 | 10.82 | 2.16 | 0.093 |
| missing_size | 2 | 1.43 | 1.20 | 1.000 |
| missing_size | 4 | 5.83 | 4.26 | 0.619 |
| missing_size | 8 | 10.24 | 3.28 | 0.102 |
| dependent | 2 | 1.51 | 1.30 | 1.000 |
| dependent | 4 | 1.55 | 1.30 | 1.000 |
| dependent | 8 | 1.60 | 1.30 | 1.000 |
| sparse | 2 | 1.49 | 1.08 | 0.750 |
| sparse | 4 | 6.20 | 2.53 | 0.333 |
| sparse | 8 | 13.53 | 1.09 | 0.037 |
| derived_degree | 2 | 1.44 | 1.24 | 1.000 |
| derived_degree | 4 | 5.81 | 4.25 | 0.655 |
| derived_degree | 8 | 10.00 | 3.77 | 0.219 |
| unconfounded | 2 | 1.41 | 1.21 | 1.000 |
| unconfounded | 4 | 1.41 | 1.21 | 1.000 |
| unconfounded | 8 | 1.41 | 1.21 | 1.000 |

Estimativas pontuais, sem intervalo de incerteza. OR=2/4, outras réplicas e carga da máquina podem alterar os tempos. Geração dos dados, checkpoints, agregação e relatórios ficam fora da extrapolação. Os testes ocorreram parcialmente em paralelo com a medição. Esta comparação mede o tempo computacional; a validação estatística completa (63.000 análises) foi subsequentemente executada e consolidada em `benchmark_results/validation/` (ver [protocolo de benchmark](../../docs/benchmark-protocol.md)).

Com oito bins, alguns cenários apresentam poucos estratos informativos, aspecto a examinar na análise de poder e suporte.

Verificação: 44 testes passaram, incluindo distribuição exata de um universo com seis pares possíveis, sorteios compartilhados entre termos duplicados e estratos com seleção fixa e variável.

Reprodução na raiz do repositório:

```bash
PYTHONPATH=. .venv/bin/python benchmarks/timing/measure.py
```

`runtime.json` contém os novos tempos, configuração, versões e hashes; `baseline/runtime.json` preserva a medição anterior. `measurements.json` é atualizado após cada análise. Reexecutar substitui a medição atual.
