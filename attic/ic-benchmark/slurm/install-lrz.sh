#!/bin/bash
# ===========================================================================
# IC benchmark: one-time / per-code-change setup on the LRZ login node.
# Run from anywhere on LRZ:
#   bash ~/pammtools-icbench/attic/ic-benchmark/slurm/install-lrz.sh
# Light work only (package install) -- no computation here.
#
# Layout (deliberately isolated from other pammtools work on this account --
# ~/pammtools is a DIFFERENT project and must not be touched):
#   repo:      ~/pammtools-icbench   (clone of fabian-s/pammtools, this branch)
#   pammtools: ~/R/ic-bench-lib      (prepended via R_LIBS_USER in the sbatch)
#   CRAN deps: ~/R/x86_64-pc-linux-gnu-library/4.3  (shared user lib)
# ===========================================================================
set -euo pipefail

module load slurm_setup
module load r/4.3.3-gcc13-mkl

REPO="$HOME/pammtools-icbench"
RLIB="$HOME/R/ic-bench-lib"
RLIB_BASE="$HOME/R/x86_64-pc-linux-gnu-library/4.3"
BRANCH="claude/keen-dirac-ZfYvO"

mkdir -p "$RLIB" "$RLIB_BASE"
if [ ! -d "$REPO/.git" ]; then
  git clone --branch "$BRANCH" https://github.com/fabian-s/pammtools.git "$REPO"
fi
cd "$REPO"
git fetch origin
git checkout "$BRANCH"
git pull origin "$BRANCH"

# CRAN dependencies into the SHARED lib. NOTE: pak's solver fails on this R
# 4.3.3 setup (current rms needs R >= 4.4 and pak refuses partial solutions),
# so: pin rms 6.7-1 via pak once, then plain install.packages for the rest,
# installing only what is missing.
Rscript -e "if (!requireNamespace('pak', quietly=TRUE)) install.packages('pak', repos='https://cloud.r-project.org', lib='$RLIB_BASE')"
Rscript -e "if (!requireNamespace('rms', quietly=TRUE)) pak::pkg_install('rms@6.7-1', lib='$RLIB_BASE', ask=FALSE)"
Rscript -e "pkgs <- c('icenReg','mvtnorm','dplyr','tidyr','purrr','tibble','mgcv','survival','Formula','pec','scam','checkmate','magrittr','rlang','ggplot2','lazyeval','vctrs'); miss <- pkgs[!vapply(pkgs, requireNamespace, TRUE, quietly=TRUE)]; if (length(miss)) install.packages(miss, repos='https://cloud.r-project.org', lib='$RLIB_BASE') else cat('all deps present\n')"

# the branch itself, into the ISOLATED lib
export R_LIBS_USER="$RLIB:$RLIB_BASE"
R CMD INSTALL --library="$RLIB" .

# sanity: branch API present from the isolated lib; icenReg loads
Rscript -e "library(pammtools, lib.loc='$RLIB'); stopifnot(exists('pamm_ic'), exists('add_inspections')); cat('pammtools', as.character(packageVersion('pammtools')), 'ok from', '$RLIB', '\n')"
Rscript -e "library(icenReg); cat('icenReg', as.character(packageVersion('icenReg')), 'ok\n')"

# partition check: serial_std must grant 16 CPUs on one node
scontrol -M serial show partition serial_std | grep -E "MaxNode|MaxCPUsPerNode|MaxTime" || true

mkdir -p "$REPO/attic/ic-benchmark/logs" "$REPO/attic/ic-benchmark/results/raw"
echo "install-lrz.sh: done"
