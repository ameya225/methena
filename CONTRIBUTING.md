# Contributing

* Report bugs with the full error message; re-run with `METHYLFLOW_DEBUG=1` to include the call stack.
* Before opening a pull request run `./run_methylflow.sh test` (or `R CMD check`); CI runs the same checks
  on Linux and macOS.
* New functions go in `R/` with roxygen comments (`roxygen2::roxygenise()` regenerates `man/`).
  Exported functions use the `mf_` prefix.
* Never commit patient-level data; `data/raw`, `data/idat`, `data/rrbs` and `results/` are git-ignored.
