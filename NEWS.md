# methena 1.1.0

* Renamed from **methylflow** to **methena** (package, repository, command-line script `inst/scripts/methena.R`,
  wrapper `run_methena.sh`, environment variables `METHENA_LIB` / `METHENA_DEBUG`, dataset file
  `methena_dataset.rds`). Exported functions keep the `mf_` prefix, so analysis code only needs
  `library(methylflow)` changed to `library(methena)`.
* Results produced by methylflow 1.0.0 can still be analysed: pass the old file explicitly with
  `./run_methena.sh analyze` after renaming `methylflow_dataset.rds` to `methena_dataset.rds`, or
  `mf_analyze(cfg, dataset = "path/to/methylflow_dataset.rds")`.
* CI runs R CMD check inside the Bioconductor container.

# methylflow 1.0.0

* First release: IDAT / intensity / beta / Bismark input, QC and normalisation, DMP/DMR/VMP analysis
  with bacon inflation correction, 17 epigenetic clocks, gene/pathway/regulatory/TF-motif mapping,
  HTML report and Shiny app.
