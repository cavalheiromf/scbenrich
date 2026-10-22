# Ensaio funcional da validação ampliada

Arquivos produzidos pelo comando de smoke test documentado em docs/benchmark-protocol.md. Cinco réplicas, 99 permutações, sete cenários, OR 1/4 e bins 2/4. Estes resultados verificam a execução e as métricas; não demonstram controle nominal de FDR. O manifesto preserva a configuração, hashes do código e versões do momento da execução. Os lotes, dados por termo e figuras desta execução estão em /tmp/scbenrich-expanded-smoke.

FDR usa bootstrap por réplica; intervalos iguais a zero em amostras pequenas não são limites superiores confiáveis. Cobertura usa IC binomial e é calculada apenas para a OR marginal conhecida no cenário unconfounded. Valores ausentes indicam métricas sem verdade conhecida ou sem efeito alternativo.
