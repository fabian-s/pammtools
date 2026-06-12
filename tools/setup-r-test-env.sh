#!/usr/bin/env bash
set -euo pipefail
INSTALL_PAMMTOOLS=1
[[ "${1:-}" == "--no-pammtools" ]] && INSTALL_PAMMTOOLS=0
echo "==> [0/6] Repair apt sources if a broken third-party PPA is present"
if ! apt-get update -qq 2>/dev/null; then
  echo "    apt-get update failed; disabling unreachable third-party sources..."
  for f in /etc/apt/sources.list.d/*; do
    [[ -e "$f" ]] || continue
    case "$(basename "$f")" in
      ubuntu.sources|docker.list) : ;;
      *) echo "      removing $f"; rm -f "$f" ;;
    esac
  done
  apt-get update -qq
fi
echo "==> [1/6] Install R, JAGS, pandoc and the C/build system libraries"
DEBIAN_FRONTEND=noninteractive apt-get install -y -qq \
  r-base r-base-dev jags pandoc \
  libcurl4-openssl-dev libssl-dev libxml2-dev \
  libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
  libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
echo "==> [2/6] Install R packages from Ubuntu's prebuilt r-cran-* binaries"
DEBIAN_FRONTEND=noninteractive apt-get install -y -qq \
  r-cran-rjags r-cran-brms r-cran-rstan \
  r-cran-mgcv r-cran-survival r-cran-dplyr r-cran-ggplot2 r-cran-tidyr \
  r-cran-knitr r-cran-rmarkdown \
  r-cran-checkmate r-cran-magrittr r-cran-rlang r-cran-purrr r-cran-tibble \
  r-cran-lazyeval r-cran-formula r-cran-mvtnorm r-cran-pec r-cran-vctrs \
  r-cran-remotes
echo "==> [3/6] Ensure Boost headers are available for rstan model compilation"
DEBIAN_FRONTEND=noninteractive apt-get install -y -qq r-cran-bh libboost-dev
echo "==> [4/6] Work around rstan's Boost guard (Debian StanHeaders quirk)"
RPROFILE_SITE="$(R RHOME)/etc/Rprofile.site"
HOOK_MARKER="# pammtools-rstan-boost-hook"
if ! grep -qF "$HOOK_MARKER" "$RPROFILE_SITE" 2>/dev/null; then
  cat >> "$RPROFILE_SITE" <<'RPROFILE'

# pammtools-rstan-boost-hook
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
if ! R --vanilla -q -e 'q(status = !("1.2-99" == as.character(tryCatch(packageVersion("scam"), error = function(e) NA))))' >/dev/null 2>&1; then
  SCAM_STUB="$(mktemp -d)/scam"
  mkdir -p "$SCAM_STUB/R"
  cat > "$SCAM_STUB/DESCRIPTION" <<'EOF'
Package: scam
Title: Shape Constrained Additive Models (stub for local testing)
Version: 1.2-99
Authors@R: person("Stub", "Stub", email="stub@example.com", role=c("aut","cre"))
Description: Local stub providing predict.scam so that pammtools can be loaded.
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
  R CMD INSTALL --no-multiarch --no-docs .
else
  echo "==> [6/6] Skipping local pammtools install (--no-pammtools)"
fi
echo "==> Verifying the toolchain"
R --vanilla -q -e '
ok <- TRUE
for (p in c("mgcv","survival","dplyr","ggplot2","tidyr","knitr","rmarkdown","scam","pammtools")) {
  have <- requireNamespace(p, quietly = TRUE)
  ok <- ok && have
  cat(sprintf("  %-12s %s\n", p, if (have) "OK" else "MISSING"))
}
if (!ok) quit(status = 1)
'
