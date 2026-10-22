#!/bin/bash
#SBATCH --job-name=scbenrich-validation
#SBATCH --partition=long
#SBATCH --array=0-83%84
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=04:00:00
#SBATCH --output=benchmark_results/slurm/validation-%A_%a.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
.venv/bin/python -u -m benchmarks.slurm task --config benchmarks/configs/validation_checkpointed.json --out-dir benchmark_results/validation --index "${SLURM_ARRAY_TASK_ID}"
