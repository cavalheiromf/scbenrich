#!/bin/bash
#SBATCH --job-name=scbenrich-report
#SBATCH --partition=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=benchmark_results/slurm/report-%j.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
.venv/bin/python -u -m benchmarks.slurm collect --config benchmarks/configs/validation_checkpointed.json --out-dir benchmark_results/validation
.venv/bin/python -m benchmarks.report benchmark_results/validation
