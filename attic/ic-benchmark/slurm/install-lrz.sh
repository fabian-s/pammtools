#!/bin/bash
# ===========================================================================
# IC benchmark: one-time / per-code-change setup on the LRZ login node.
# Run from anywhere on LRZ:  bash ~/pammtools/attic/ic-benchmark/slurm/install-lrz.sh
# Light work only (package install) -- no computation here.
# ===========================================================================
set -euo pipefail

module load slurm_setup
module load r/4.3.3-gcc13-mkl

REPO="$HOME/pammtools"
RLIB="$HOME/R/x86_64-pc-linux-gnu-library/4.3"
BRANCH="claude/keen-dirac-ZfYvO"

mkdir -p "$RLIB"
cd "$REPO"
git fetch origin
git checkout "$BRANCH"
git pull origin "$BRANCH"

# CRAN dependencies (pak resolves conflicts correctly)
Rscript -e "if (!requireNamespace('pak', quietly=TRUE)) install.packages('pak', repos='https://cloud.r-project.org', lib='$RLIB')"
Rscript -e "pak::pkg_install(c('icenReg','mvtnorm','dplyr','tidyr','purrr','tibble','mgcv','survival'), lib='$RLIB', ask=FALSE)"

# the branch itself
R CMD INSTALL --library="$RLIB" .

# sanity: branch API present, icenReg loads (compiled C++ against current toolchain)
Rscript -e "library(pammtools); stopifnot(exists('pamm_ic'), exists('add_inspections')); cat('pammtools', as.character(packageVersion('pammtools')), 'ok\n')"
Rscript -e "library(icenReg); cat('icenReg', as.character(packageVersion('icenReg')), 'ok\n')"

# partition check: serial_std must grant 16 CPUs on one node
scontrol -M serial show partition serial_std | grep -E "MaxNode|MaxCPUsPerNode|MaxTime" || true

mkdir -p "$REPO/attic/ic-benchmark/logs" "$REPO/attic/ic-benchmark/results/raw"
echo "install-lrz.sh: done"
