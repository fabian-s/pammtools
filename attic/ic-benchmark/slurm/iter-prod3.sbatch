#!/bin/bash
#SBATCH -J ic-it3
#SBATCH -o attic/ic-benchmark/logs/%x.%A_%a.out
#SBATCH -e attic/ic-benchmark/logs/%x.%A_%a.err
#SBATCH -D /dss/dsshome1/lxc04/ri68suz2/pammtools-icbench
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --array=1-24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=24000mb
#SBATCH --time=04:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fabian.scheipl@stat.uni-muenchen.de
#SBATCH --export=NONE
# Iterated-MI production arm: mi_iter3 (pamm_ic iter = 3, the recommended
# sparse-inspection setting), 500 reps on ALL 24 core cells -- sparse cells for
# the headline, dense cells to verify no harm. Production data seeds (paired).

module load slurm_setup
module load r/4.3.3-gcc13-mkl

export R_LIBS_USER="${HOME}/R/ic-bench-lib:${HOME}/R/x86_64-pc-linux-gnu-library/4.3"
export IC_BENCH_LIB="${HOME}/R/ic-bench-lib"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MKL_DYNAMIC=FALSE
export ITER_PILOT_REPS=500
export ITER_PILOT_METHODS="mi_iter3"
export ITER_PILOT_OUTDIR="iter-prod3"
export ITER_PILOT_CELLSET="core24"

Rscript attic/ic-benchmark/run-iter-pilot.R "$SLURM_ARRAY_TASK_ID" "$SLURM_CPUS_PER_TASK"
