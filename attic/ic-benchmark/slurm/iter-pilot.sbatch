#!/bin/bash
#SBATCH -J ic-iter
#SBATCH -o attic/ic-benchmark/logs/%x.%A_%a.out
#SBATCH -e attic/ic-benchmark/logs/%x.%A_%a.err
#SBATCH -D /dss/dsshome1/lxc04/ri68suz2/pammtools-icbench
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --array=1-10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=24000mb
#SBATCH --time=02:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fabian.scheipl@stat.uni-muenchen.de
#SBATCH --export=NONE
# Iterated-MI extension pilot: mi_iter (pamm_ic iter = 2) only, 50 reps on the
# 10 cells listed in run-iter-pilot.R, same per-rep seeds as production.

module load slurm_setup
module load r/4.3.3-gcc13-mkl

export R_LIBS_USER="${HOME}/R/ic-bench-lib:${HOME}/R/x86_64-pc-linux-gnu-library/4.3"
export IC_BENCH_LIB="${HOME}/R/ic-bench-lib"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MKL_DYNAMIC=FALSE
export ITER_PILOT_REPS=50

Rscript attic/ic-benchmark/run-iter-pilot.R "$SLURM_ARRAY_TASK_ID" "$SLURM_CPUS_PER_TASK"
