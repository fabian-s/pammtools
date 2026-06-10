#!/bin/bash
#SBATCH -J ic-nsim
#SBATCH -o attic/ic-benchmark/logs/%x.%j.out
#SBATCH -e attic/ic-benchmark/logs/%x.%j.err
#SBATCH -D /dss/dsshome1/lxc04/ri68suz2/pammtools-icbench
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=24000mb
#SBATCH --time=00:30:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fabian.scheipl@stat.uni-muenchen.de
#SBATCH --export=NONE

module load slurm_setup
module load r/4.3.3-gcc13-mkl

export R_LIBS_USER="${HOME}/R/ic-bench-lib:${HOME}/R/x86_64-pc-linux-gnu-library/4.3"
export IC_BENCH_LIB="${HOME}/R/ic-bench-lib"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MKL_DYNAMIC=FALSE

Rscript attic/ic-benchmark/nsim-sensitivity.R 20 16
