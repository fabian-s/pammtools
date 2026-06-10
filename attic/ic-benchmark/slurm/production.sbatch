#!/bin/bash
#SBATCH -J ic-bench
#SBATCH -o attic/ic-benchmark/logs/%x.%A_%a.out
#SBATCH -e attic/ic-benchmark/logs/%x.%A_%a.err
#SBATCH -D /dss/dsshome1/lxc04/ri68suz2/pammtools-icbench
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --array=1-39
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=24000mb
#SBATCH --time=04:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fabian.scheipl@stat.uni-muenchen.de
#SBATCH --export=NONE
# NOTE: --time is a placeholder cap; after the P3 pilot, set per-task walltime
# (pilot s/rep x 500 / 16 x 2) by resubmitting slow cells with --time overrides:
#   sbatch --time=HH:MM:SS --array=<id> attic/ic-benchmark/slurm/production.sbatch

module load slurm_setup
module load r/4.3.3-gcc13-mkl

export R_LIBS_USER="${HOME}/R/ic-bench-lib:${HOME}/R/x86_64-pc-linux-gnu-library/4.3"
export IC_BENCH_LIB="${HOME}/R/ic-bench-lib"
# MKL x mclapply oversubscription guard
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MKL_DYNAMIC=FALSE

Rscript attic/ic-benchmark/run-task.R "$SLURM_ARRAY_TASK_ID" "$SLURM_CPUS_PER_TASK"
