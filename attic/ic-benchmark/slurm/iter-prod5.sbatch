#!/bin/bash
#SBATCH -J ic-it5
#SBATCH -o attic/ic-benchmark/logs/%x.%A_%a.out
#SBATCH -e attic/ic-benchmark/logs/%x.%A_%a.err
#SBATCH -D /dss/dsshome1/lxc04/ri68suz2/pammtools-icbench
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --array=1-8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=24000mb
#SBATCH --time=06:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fabian.scheipl@stat.uni-muenchen.de
#SBATCH --export=NONE
# Iterated-MI convergence exhibit: mi_iter5 (pamm_ic iter = 5), 500 reps on
# the 8 rate-0.3 core cells (where iteration matters). Production data seeds.

module load slurm_setup
module load r/4.3.3-gcc13-mkl

export R_LIBS_USER="${HOME}/R/ic-bench-lib:${HOME}/R/x86_64-pc-linux-gnu-library/4.3"
export IC_BENCH_LIB="${HOME}/R/ic-bench-lib"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MKL_DYNAMIC=FALSE
export ITER_PILOT_REPS=500
export ITER_PILOT_METHODS="mi_iter5"
export ITER_PILOT_OUTDIR="iter-prod5"
export ITER_PILOT_CELLSET="r03"

Rscript attic/ic-benchmark/run-iter-pilot.R "$SLURM_ARRAY_TASK_ID" "$SLURM_CPUS_PER_TASK"
