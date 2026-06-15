#!/bin/bash
#SBATCH -J ic-itsens
#SBATCH -o attic/ic-benchmark/logs/%x.%A_%a.out
#SBATCH -e attic/ic-benchmark/logs/%x.%A_%a.err
#SBATCH -D /dss/dsshome1/lxc04/ri68suz2/pammtools-icbench
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --array=1,2,4,6
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=24000mb
#SBATCH --time=02:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fabian.scheipl@stat.uni-muenchen.de
#SBATCH --export=NONE
# iter-sensitivity: mi_iter3 + mi_iter5 on the 4 worst pilot cells
# (array ids index ITER_PILOT_CELLS: 1 = const-ph n300, 2 = const-ph n1000,
#  4 = peaked-ph n1000, 6 = const-tv n1000; all rate 0.3), 50 reps,
# production data seeds. Answers converge-vs-bias-floor.

module load slurm_setup
module load r/4.3.3-gcc13-mkl

export R_LIBS_USER="${HOME}/R/ic-bench-lib:${HOME}/R/x86_64-pc-linux-gnu-library/4.3"
export IC_BENCH_LIB="${HOME}/R/ic-bench-lib"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MKL_DYNAMIC=FALSE
export ITER_PILOT_REPS=50
export ITER_PILOT_METHODS="mi_iter3,mi_iter5"
export ITER_PILOT_OUTDIR="iter-sens"

Rscript attic/ic-benchmark/run-iter-pilot.R "$SLURM_ARRAY_TASK_ID" "$SLURM_CPUS_PER_TASK"
