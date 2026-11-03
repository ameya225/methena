#!/usr/bin/env Rscript
# Installs every dependency of methylflow. Safe to re-run (skips what is installed).
# Usage: Rscript install_dependencies.R [--lib /path/to/library] [--minimal]
args <- commandArgs(TRUE)
lib <- if ("--lib" %in% args) args[which(args == "--lib") + 1] else NULL
minimal <- "--minimal" %in% args
if (!is.null(lib)) { dir.create(lib, showWarnings = FALSE, recursive = TRUE); .libPaths(c(lib, .libPaths())) }
ncores <- suppressWarnings(parallel::detectCores()); if (is.na(ncores)) ncores <- 2L
if (is.null(getOption("repos")) || getOption("repos")[["CRAN"]] %in% c("@CRAN@", "")) options(repos = c(CRAN = "https://cloud.r-project.org"))
options(Ncpus = max(1L, ncores - 1L),
        timeout = 3600)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", lib = .libPaths()[1])
cran <- c("data.table", "yaml", "optparse", "ggplot2", "ggrepel", "patchwork", "rmarkdown", "knitr",
          "testthat", "shiny", "DT", "plotly", "msigdbr", "glmnet", "matrixStats", "R.utils", "dplyr",
          "tidyr", "remotes", "jsonlite", "digest", "htmltools", "scales", "BiasedUrn", "RPMM",
          "bslib", "lmtest", "roxygen2", "base64enc", "withr")
bioc <- c("minfi", "limma", "DMRcate", "missMethyl", "wateRmelon", "EpiDISH", "methylKit", "bsseq",
          "DSS", "GEOquery", "bacon", "HiBED", "preprocessCore",
          "IlluminaHumanMethylation450kanno.ilmn12.hg19", "IlluminaHumanMethylation450kmanifest",
          "IlluminaHumanMethylationEPICanno.ilm10b4.hg19", "IlluminaHumanMethylationEPICmanifest",
          "IlluminaHumanMethylationEPICv2anno.20a1.hg38", "IlluminaHumanMethylationEPICv2manifest",
          "UCSC.utils", "GenomeInfoDb", "GO.db", "TxDb.Hsapiens.UCSC.hg19.knownGene", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db", "GenomicRanges", "rtracklayer")
motif <- c("monaLisa", "JASPAR2024", "JASPAR2020", "TFBSTools", "BSgenome.Hsapiens.UCSC.hg19", "RSQLite")
github <- c(dnaMethyAge = "yiluyucheng/dnaMethyAge")
want <- c(cran, bioc, if (!minimal) motif)
missing <- setdiff(want, rownames(installed.packages()))
message("Installing ", length(missing), " packages: ", paste(missing, collapse = ", "))
if (length(missing)) BiocManager::install(missing, update = FALSE, ask = FALSE, lib = .libPaths()[1])
for (p in names(github)) if (!requireNamespace(p, quietly = TRUE))
  remotes::install_github(github[[p]], upgrade = "never", dependencies = FALSE, lib = .libPaths()[1])
loads <- function(p) tryCatch({ suppressMessages(loadNamespace(p)); TRUE }, error = function(e) FALSE)
# second pass: anything that is present but cannot be loaded (e.g. a dependency failed in parallel install)
broken <- Filter(function(p) !loads(p), c(want, names(github)))
if (length(broken)) {
  message("Re-installing packages that do not load: ", paste(broken, collapse = ", "))
  BiocManager::install(intersect(broken, c(cran, bioc, motif)), update = FALSE, ask = FALSE, force = TRUE, lib = .libPaths()[1], Ncpus = 1L)
}
still <- Filter(function(p) !loads(p), c(want, names(github)))
if (length(still)) { msg <- paste0("FAILED to install: ", paste(still, collapse = ", "))
  if (interactive()) stop(msg) else { message(msg); quit(status = 1) } }
message("All methylflow dependencies are installed.")
