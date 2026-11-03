#!/usr/bin/env Rscript
# methylflow command-line interface
#   Rscript methylflow.R preprocess --config config.yaml
#   Rscript methylflow.R analyze    --config config.yaml
#   Rscript methylflow.R all        --config config.yaml
#   Rscript methylflow.R shiny      --results results/analysis [--port 3838]
#   Rscript methylflow.R check      --config config.yaml      (validate config + inputs only)
#   Rscript methylflow.R geo-sheet  --series GSE..._series_matrix.txt.gz --out samples.tsv
suppressPackageStartupMessages({
  if (!requireNamespace("optparse", quietly = TRUE)) stop("Package 'optparse' missing: run install_dependencies.R")
  library(optparse)
})
args <- commandArgs(TRUE)
cmd <- if (length(args)) args[1] else "help"
opts <- parse_args(OptionParser(option_list = list(
  make_option("--config", type = "character"), make_option("--results", type = "character"),
  make_option("--dataset", type = "character", help = "override dataset RDS for 'analyze'"),
  make_option("--port", type = "integer", default = 3838L), make_option("--series", type = "character"),
  make_option("--out", type = "character"), make_option("--no-browser", action = "store_true", default = FALSE))),
  args = args[-1])
if (cmd %in% c("help", "-h", "--help")) { cat(readLines(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))[2:8], sep = "\n"); quit(status = 0) }
suppressPackageStartupMessages(library(methylflow))
`%||%` <- function(a, b) if (is.null(a)) b else a
# set METHYLFLOW_DEBUG=1 to print the call stack on errors
run <- function(expr) tryCatch(withCallingHandlers(expr,
  error = function(e) if (identical(Sys.getenv("METHYLFLOW_DEBUG"), "1")) { message("Call stack:"); print(utils::tail(sys.calls(), 12)) }),
  error = function(e) { message("\nERROR: ", conditionMessage(e)); quit(status = 1) })
need_cfg <- function() { if (is.null(opts$config)) { message("ERROR: --config is required"); quit(status = 2) }; run(mf_read_config(opts$config)) }
switch(cmd,
  check = { cfg <- need_cfg(); ss <- run(mf_read_samples(cfg)); message("Config OK. ", nrow(ss), " samples; columns: ", paste(names(ss), collapse = ", "))
            if (!is.null(ss$group)) print(table(group = ss$group, useNA = "ifany")) },
  preprocess = { cfg <- need_cfg(); run(mf_preprocess(cfg)) },
  analyze = { cfg <- need_cfg(); run(mf_analyze(cfg, opts$dataset)) },
  all = { cfg <- need_cfg(); run(mf_preprocess(cfg)); run(mf_analyze(cfg)) },
  shiny = { if (is.null(opts$results)) { message("ERROR: --results <analysis folder> is required"); quit(status = 2) }
            run(mf_shiny(opts$results, port = opts$port, launch.browser = !isTRUE(opts$`no-browser`))) },
  demo = { out <- opts$out %||% "methylflow_demo"; dir.create(out, showWarnings = FALSE, recursive = TRUE); out <- normalizePath(out)
           message("Simulating a small case-control 450K dataset in ", out)
           sim <- run(mf_simulate(file.path(out, "data"), "intensity", n_per_group = 15, n_cpg = 20000, n_dmp = 200))
           cfg <- list(project = "methylflow_demo", output_dir = file.path(out, "results"), seed = 1L,
                       input = list(type = "intensity", sample_sheet = sim$sample_sheet, sample_id_col = "sample_id", data_id_col = "data_id",
                                    meth_file = sim$meth_file, unmeth_file = sim$unmeth_file, detp_file = sim$detp_file, tissue = "other"),
                       phenotype = list(group_col = "group", control_level = "control", case_level = "case", age_col = "age", sex_col = "sex"),
                       normalization = list(method = "quantile"), covariates = list(cell_type = "none"),
                       analysis = list(covariates = c("age", "sex"), clocks = FALSE, enrichment = TRUE))
           yaml::write_yaml(cfg, file.path(out, "demo_config.yaml"))
           cfg <- run(mf_read_config(file.path(out, "demo_config.yaml")))
           run(mf_preprocess(cfg)); r <- run(mf_analyze(cfg))
           t <- r$dmp$table; found <- mean(intersect(sim$planted, t$cpg) %in% t$cpg[t$FDR < 0.05])
           message(sprintf("\nDEMO OK: recovered %.0f%% of the 200 planted differences. Open %s", 100 * found, file.path(out, "results", "analysis", "analysis_report.html")))
           if (found < 0.8) quit(status = 1) },
  `geo-sheet` = { if (is.null(opts$series) || is.null(opts$out)) { message("ERROR: --series and --out are required"); quit(status = 2) }
                  ph <- run(mf_read_geo_pheno(opts$series)); mf_write_tsv(ph, opts$out); message("Wrote ", nrow(ph), " samples x ", ncol(ph), " columns to ", opts$out) },
  { message("Unknown command '", cmd, "'. Use: preprocess | analyze | all | shiny | check | demo | geo-sheet"); quit(status = 2) })
