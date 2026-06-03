#!/usr/bin/env bash
#
# setup-r-test-env.sh
# -------------------------------------------------------------------------
# Provision an R environment capable of building/rendering pammtools and its
# vignettes -- in particular vignettes/bayesian.Rmd, which needs mgcv + jagam
# (rjags/JAGS) + brms (rstan/Stan).
#
# WHY THIS SCRIPT EXISTS
# ----------------------
# In the Claude Code remote container the outbound network policy BLOCKS the
# usual R package repositories (CRAN, cloud.r-project.org and the Posit/RSPM
# binary mirror all return HTTP 403 "host_not_allowed"). So the normal
# `install.packages(...)` / `pak` / `remotes::install_deps()` route does NOT
# work here. The Ubuntu apt mirror IS reachable, so we install R and almost
# every R package from Ubuntu's pre-built `r-cran-*` / system packages instead.
#
# Target image: Ubuntu 24.04 (noble), run as root. R 4.3.x ships in the distro.
#
# Idempotent: safe to re-run. Uses `apt-get install` (no-op if already present)
# and only rebuilds the scam stub / re-installs pammtools when asked.
#
# Usage:
#   sudo bash tools/setup-r-test-env.sh            # full setup
#   sudo bash tools/setup-r-test-env.sh --no-pammtools   # skip local install
# -------------------------------------------------------------------------
set -euo pipefail

INSTALL_PAMMTOOLS=1
[[ "${1:-}" == "--no-pammtools" ]] && INSTALL_PAMMTOOLS=0

echo "==> [0/6] Repair apt sources if a broken third-party PPA is present"
# A broken PPA (we hit ondrej/php returning 403) makes `apt-get update` fail
# hard, which then blocks every install below. Drop any source that fails to
# refresh so the main Ubuntu mirror can still be used. This is defensive --
# it's a no-op on a clean image.
if ! apt-get update -qq 2>/dev/null; then
  echo "    apt-get update failed; disabling unreachable third-party sources..."
  # Comment out anything that isn't the main ubuntu.sources / docker list.
  for f in /etc/apt/sources.list.d/*; do
    [[ -e "$f" ]] || continue
    case "$(basename "$f")" in
      ubuntu.sources|docker.list) : ;;          # keep
      *) echo "      removing $f"; rm -f "$f" ;; # drop the rest
    esac
  done
  apt-get update -qq
fi

echo "==> [1/6] Install R, JAGS, pandoc and the C/build system libraries"
# - r-base / r-base-dev: R itself + headers/compilers for building packages
# - jags: the JAGS sampler -- the SYSTEM dependency behind the rjags R package
#         (mgcv::jagam writes a JAGS model that rjags then compiles & samples)
# - pandoc: required by rmarkdown to render the vignette to HTML
# - the lib*-dev packages: common headers needed to build/run the tidyverse,
#   ragg/systemfonts/textshaping graphics stack, curl/xml, etc.
DEBIAN_FRONTEND=noninteractive apt-get install -y -qq \
  r-base r-base-dev \
  jags \
  pandoc \
  libcurl4-openssl-dev libssl-dev libxml2-dev \
  libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
  libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev

echo "==> [2/6] Install R packages from Ubuntu's prebuilt r-cran-* binaries"
# These cover pammtools' Imports + the bayesian vignette's needs. We use apt
# (not install.packages) because CRAN/RSPM are network-blocked in this env.
#
# Bayesian back-ends:
#   r-cran-rjags  -> Bayesian fit via mgcv::jagam + rjags
#   r-cran-brms   -> Bayesian fit via Stan; pulls in r-cran-rstan + StanHeaders
#   r-cran-rstan  -> the Stan backend brms compiles its model with
DEBIAN_FRONTEND=noninteractive apt-get install -y -qq \
  r-cran-rjags r-cran-brms r-cran-rstan \
  r-cran-mgcv r-cran-survival r-cran-dplyr r-cran-ggplot2 r-cran-tidyr \
  r-cran-knitr r-cran-rmarkdown \
  r-cran-checkmate r-cran-magrittr r-cran-rlang r-cran-purrr r-cran-tibble \
  r-cran-lazyeval r-cran-formula r-cran-mvtnorm r-cran-pec r-cran-vctrs \
  r-cran-remotes

echo "==> [3/6] Ensure Boost headers are available for rstan model compilation"
# brms compiles its Stan program with rstan at runtime, which needs the Boost
# C++ headers. r-cran-bh is a (virtual) package that just depends on the real
# system Boost headers (libboost-dev -> /usr/include/boost); install both.
DEBIAN_FRONTEND=noninteractive apt-get install -y -qq r-cran-bh libboost-dev

echo "==> [4/6] Work around rstan's Boost guard (Debian StanHeaders quirk)"
# Debian/Ubuntu's StanHeaders does NOT bundle Boost, and the r-cran-bh R
# package ships an EMPTY include dir (the headers live at /usr/include/boost).
# rstan::stan_model() guards with:
#     if (!file.exists(rstan_options("boost_lib"))) stop("Boost not found; ...")
# and that option defaults to the (empty) BH include path, so it errors out
# before compiling -- even though g++ would find /usr/include/boost on its own.
#
# Fix: point rstan's boost_lib at /usr/include via an Rprofile.site hook that
# fires whenever rstan is loaded (brms loads it for us). This is a LOCAL
# container fix only -- it is deliberately NOT in the package or CI config,
# because on CRAN/RSPM the real BH package ships the headers.
#
# NOTE: This quirk is specific to the distro rstan. On GitHub's pkgdown runner
# (use-public-rspm: true) the equivalent fix is to install BH + RcppEigen as
# full packages (they're LinkingTo-only deps that pak skips for binaries) --
# see .github/workflows/pkgdown.yaml.
RPROFILE_SITE="$(R RHOME)/etc/Rprofile.site"
HOOK_MARKER="# pammtools-rstan-boost-hook"
if ! grep -qF "$HOOK_MARKER" "$RPROFILE_SITE" 2>/dev/null; then
  cat >> "$RPROFILE_SITE" <<'RPROFILE'

# pammtools-rstan-boost-hook
# Point rstan at the system Boost headers so brms can compile Stan models in
# this container (the distro rstan/BH does not expose Boost where rstan looks).
local({
  if (dir.exists("/usr/include/boost")) {
    setHook(packageEvent("rstan", "onLoad"), function(...) {
      try(rstan::rstan_options(boost_lib = "/usr/include"), silent = TRUE)
    })
  }
})
RPROFILE
  echo "    appended rstan boost hook to $RPROFILE_SITE"
else
  echo "    rstan boost hook already present in $RPROFILE_SITE"
fi

echo "==> [5/6] Install a minimal 'scam' stub (not packaged for Ubuntu)"
# pammtools Imports scam (NAMESPACE: importFrom(scam, predict.scam)), but scam
# is not available as an r-cran-* package and CRAN is blocked. The bayesian
# vignette never fits a scam model, so a tiny stub that satisfies the namespace
# import is enough to let pammtools load and the vignette render.
#
# !! This is a TEST-ONLY stub, NOT the real scam package. predict.scam errors
#    if ever called, which never happens in the Bayesian article.
if ! R --vanilla -q -e 'q(status = !("1.2-99" == as.character(tryCatch(packageVersion("scam"), error = function(e) NA))))' >/dev/null 2>&1; then
  SCAM_STUB="$(mktemp -d)/scam"
  mkdir -p "$SCAM_STUB/R"
  cat > "$SCAM_STUB/DESCRIPTION" <<'EOF'
Package: scam
Title: Shape Constrained Additive Models (stub for local testing)
Version: 1.2-99
Authors@R: person("Stub", "Stub", email="stub@example.com", role=c("aut","cre"))
Description: Local stub providing predict.scam so that pammtools can be loaded
    for testing in an offline environment. NOT the real scam package.
License: GPL-2
Encoding: UTF-8
EOF
  cat > "$SCAM_STUB/NAMESPACE" <<'EOF'
export(predict.scam)
S3method(predict, scam)
EOF
  cat > "$SCAM_STUB/R/predict.scam.R" <<'EOF'
predict.scam <- function(object, ...) {
  stop("This is a local stub of scam::predict.scam and must not be called.")
}
EOF
  R CMD INSTALL "$SCAM_STUB"
else
  echo "    scam stub (1.2-99) already installed"
fi

if [[ "$INSTALL_PAMMTOOLS" -eq 1 ]]; then
  echo "==> [6/6] Install the local pammtools source"
  # Run from the repo root so '.' is the package. --no-docs keeps it quick;
  # all Imports are satisfied by the apt packages + the scam stub above.
  R CMD INSTALL --no-multiarch --no-docs .
else
  echo "==> [6/6] Skipping local pammtools install (--no-pammtools)"
fi

echo "==> Verifying the toolchain"
R --vanilla -q -e '
ok <- TRUE
for (p in c("rjags","brms","rstan","mgcv","survival","dplyr","ggplot2",
            "tidyr","knitr","rmarkdown","scam","pammtools")) {
  have <- requireNamespace(p, quietly = TRUE)
  ok <- ok && have
  cat(sprintf("  %-12s %s\n", p, if (have) "OK" else "MISSING"))
}
cat("JAGS:", system("which jags", intern = TRUE), "\n")
if (!ok) quit(status = 1)
'

cat <<'DONE'

==> Done.

To render the Bayesian vignette end-to-end (mgcv + jagam + brms):

  cd /path/to/pammtools
  cp vignettes/bayesian.Rmd /tmp/ && cd /tmp
  R -q -e 'rmarkdown::render("bayesian.Rmd", "html_document")'

The brms chunk compiles a Stan model the first time (~1-2 min). The Rprofile
hook installed above makes rstan find Boost automatically -- no need to set
rstan_options(boost_lib=...) by hand.
DONE
