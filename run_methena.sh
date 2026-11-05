#!/usr/bin/env bash
# =============================================================================
# methena - one-command wrapper
#   ./run_methena.sh setup                     install R packages + methena (once, ~30-60 min)
#   ./run_methena.sh demo                      run a 2-minute simulated example to check the install
#   ./run_methena.sh check     my_config.yaml  validate the config and sample sheet only
#   ./run_methena.sh preprocess my_config.yaml QC + normalisation  -> <output_dir>/preprocess
#   ./run_methena.sh analyze   my_config.yaml  statistics + report -> <output_dir>/analysis
#   ./run_methena.sh all       my_config.yaml  preprocess + analyze
#   ./run_methena.sh shiny     <output_dir>/analysis   open the interactive results browser
#   ./run_methena.sh test                      run the automated test suite
# Environment variables: METHENA_LIB (R library folder, default ./Rlib), RSCRIPT (path to Rscript),
#                        METHENA_DEBUG=1 (print call stack on errors)
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export METHENA_LIB="${METHENA_LIB:-$HERE/Rlib}"
export R_LIBS_USER="$METHENA_LIB"
RSCRIPT="${RSCRIPT:-$(command -v Rscript || true)}"
if [ -z "$RSCRIPT" ] && [ -x /Library/Frameworks/R.framework/Resources/bin/Rscript ]; then RSCRIPT=/Library/Frameworks/R.framework/Resources/bin/Rscript; fi
[ -n "$RSCRIPT" ] || { echo "ERROR: R is not installed. Install R (>= 4.4) from https://cran.r-project.org and re-run."; exit 1; }
CLI="$HERE/inst/scripts/methena.R"
cmd="${1:-help}"; shift || true
case "$cmd" in
  setup)
    mkdir -p "$METHENA_LIB"
    "$RSCRIPT" "$HERE/inst/scripts/install_dependencies.R" --lib "$METHENA_LIB" "$@"
    "$(dirname "$RSCRIPT")/R" CMD INSTALL --library="$METHENA_LIB" "$HERE"
    echo; echo "Setup complete. Try:  ./run_methena.sh demo" ;;
  demo)       "$RSCRIPT" "$CLI" demo --out "${1:-$HERE/demo_results}" ;;
  check|preprocess|analyze|all)
    [ $# -ge 1 ] || { echo "Usage: $0 $cmd <config.yaml>"; exit 2; }
    "$RSCRIPT" "$CLI" "$cmd" --config "$1" ;;
  shiny)
    [ $# -ge 1 ] || { echo "Usage: $0 shiny <output_dir>/analysis [port]"; exit 2; }
    "$RSCRIPT" "$CLI" shiny --results "$1" --port "${2:-3838}" ;;
  test)
    "$RSCRIPT" -e "library(testthat); library(methena); res <- test_dir('$HERE/tests/testthat', reporter = 'summary', stop_on_failure = TRUE)" ;;
  help|-h|--help) sed -n '3,13p' "$0" ;;
  *) echo "Unknown command '$cmd'"; sed -n '3,13p' "$0"; exit 2 ;;
esac
