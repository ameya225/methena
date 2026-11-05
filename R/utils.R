# ---- logging ---------------------------------------------------------------
.mf_env <- new.env(parent = emptyenv())
.mf_env$log_file <- NULL
.mf_env$timings <- list()

#' Write a timestamped log message to the console and (optionally) a log file
#' @param ... message parts
#' @param level one of INFO, WARN, ERROR
#' @export
mf_log <- function(..., level = "INFO") {
  msg <- sprintf("[%s] %-5s %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, paste0(...))
  message(msg)
  if (!is.null(.mf_env$log_file))
    tryCatch(suppressWarnings(cat(msg, "\n", file = .mf_env$log_file, append = TRUE, sep = "")), error = function(e) .mf_env$log_file <- NULL)
  invisible(msg)
}

#' Start writing log messages to a file
#' @param path log file path
#' @export
mf_log_to <- function(path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  .mf_env$log_file <- path
  invisible(path)
}

#' Time an expression and record it under a step name
#' @param step step name
#' @param expr expression
#' @export
mf_step <- function(step, expr) {
  mf_log(">> ", step)
  t0 <- proc.time()[["elapsed"]]
  res <- force(expr)
  dt <- proc.time()[["elapsed"]] - t0
  .mf_env$timings[[step]] <- round(dt, 2)
  mf_log("<< ", step, " (", round(dt, 1), " s)")
  res
}

#' Run an optional analysis module; on error log a warning and return NULL instead of aborting
#' @param step step name
#' @param expr expression
#' @export
mf_try_step <- function(step, expr) {
  tryCatch(mf_step(step, expr), error = function(e) {
    mf_log("Optional step '", step, "' failed and was skipped: ", conditionMessage(e), level = "WARN"); NULL })
}

#' Log an error and stop with a readable message
#' @param ... message parts
#' @keywords internal
#' @export
mf_stop <- function(...) {
  msg <- paste0(...)
  mf_log(msg, level = "ERROR")
  stop(msg, call. = FALSE)
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

need_pkg <- function(pkg, why = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE))
    mf_stop("Package '", pkg, "' is required", if (!is.null(why)) paste0(" for ", why),
            ". Run inst/scripts/install_dependencies.R")
  invisible(TRUE)
}

# ---- configuration ---------------------------------------------------------

#' Default configuration (every option documented in inst/config/config_template.yaml)
#' @export
mf_default_config <- function() {
  list(
    project = "methylation_study",
    output_dir = "methena_results",
    seed = 20240913L,
    threads = 1L,
    input = list(type = NULL, array = "auto", genome = "hg19",
                 sample_sheet = NULL, geo_series_matrix = NULL, sample_id_col = NULL,
                 data_id_col = NULL, idat_dir = NULL, basename_col = NULL,
                 meth_file = NULL, unmeth_file = NULL, detp_file = NULL, beta_file = NULL,
                 bismark_dir = NULL, bismark_pattern = "\\.cov(\\.gz)?$", bismark_file_col = NULL,
                 tissue = "blood"),
    phenotype = list(group_col = NULL, case_level = NULL, control_level = NULL,
                     age_col = NULL, sex_col = NULL, batch_cols = NULL, extra_covariates = NULL,
                     recode = NULL, keep_samples_query = NULL),
    qc = list(detp_threshold = 0.01, sample_max_failed_frac = 0.05, probe_max_failed_frac = 0.05,
              min_beadcount = 3L, sex_check = TRUE, drop_sex_mismatch = TRUE,
              outlier_sd = 4, drop_outliers = TRUE, remove_snp_probes = TRUE, snp_maf = 0.05,
              remove_cross_reactive = TRUE, remove_xy = TRUE,
              bs_min_coverage = 10L, bs_min_sample_frac = 0.8, bs_max_coverage_quantile = 0.999,
              bs_min_cpgs_per_sample = 1000L),
    normalization = list(method = "auto"),
    covariates = list(cell_type = "auto", smoking_score = TRUE, predict_sex = TRUE, dnam_age = TRUE),
    analysis = list(covariates = c("age", "sex", "cells", "smoking"), batch_cols = NULL,
                    n_pcs = 0L, inflation_correction = "bacon", fdr = 0.05, min_delta_beta = 0,
                    dmr = TRUE, dmr_lambda = 1000, dmr_C = 2, dmr_min_cpgs = 3L,
                    vmp = TRUE, clocks = TRUE,
                    clock_set = c("HorvathS2013", "HannumG2013", "LevineM2018", "HorvathS2018", "ZhangQ2019",
                                  "ShirebyG2020", "McEwenL2019", "PCHorvathS2013", "PCHannumG2013",
                                  "PCHorvathS2018", "PCPhenoAge", "PCGrimAge", "DunedinPACE", "LuA2019", "YangZ2016",
                                  "epiTOC2", "ZhangY2017"),
                    enrichment = TRUE, enrichment_top_n = 1000L, msigdb_collections = c("H"), gene_set_files = NULL,
                    motif = FALSE, motif_top_n = 1000L, motif_flank = 100L,
                    compare_with = NULL, report = TRUE)
  )
}

merge_list <- function(base, over) {
  for (n in names(over)) {
    if (is.list(base[[n]]) && is.list(over[[n]]) && !is.null(names(over[[n]])))
      base[[n]] <- merge_list(base[[n]], over[[n]])
    else base[n] <- list(over[[n]])
  }
  base
}

#' Read a YAML config, fill defaults, resolve relative paths, validate
#' @param path YAML file
#' @param overrides optional named list merged on top
#' @export
mf_read_config <- function(path, overrides = NULL) {
  if (!file.exists(path)) mf_stop("Config file not found: ", path)
  user <- tryCatch(yaml::read_yaml(path), error = function(e) mf_stop("Config is not valid YAML: ", conditionMessage(e)))
  cfg <- merge_list(mf_default_config(), user %||% list())
  if (!is.null(overrides)) cfg <- merge_list(cfg, overrides)
  base <- dirname(normalizePath(path))
  cfg$config_dir <- base
  resolve <- function(p) if (is.null(p) || grepl("^(/|~|[A-Za-z]:)", p)) p else file.path(base, p)
  for (k in c("sample_sheet", "geo_series_matrix", "idat_dir", "meth_file", "unmeth_file",
              "detp_file", "beta_file", "bismark_dir"))
    cfg$input[[k]] <- resolve(cfg$input[[k]])
  cfg$output_dir <- resolve(cfg$output_dir)
  if (!is.null(cfg$analysis$compare_with))
    cfg$analysis$compare_with <- lapply(cfg$analysis$compare_with, function(x) { x$path <- resolve(x$path); x })
  mf_validate_config(cfg)
}

#' Validate a configuration list, stopping with a readable message on problems
#' @param cfg config list
#' @export
mf_validate_config <- function(cfg) {
  inp <- cfg$input
  types <- c("idat", "intensity", "beta", "bismark")
  if (is.null(inp$type) || !inp$type %in% types)
    mf_stop("input$type must be one of: ", paste(types, collapse = ", "))
  chk_file <- function(f, what) {
    if (is.null(f)) mf_stop("input$", what, " is required for input type '", inp$type, "'")
    if (!file.exists(f)) mf_stop("File/folder for input$", what, " does not exist: ", f)
  }
  if (is.null(inp$sample_sheet) && is.null(inp$geo_series_matrix))
    mf_stop("Provide input$sample_sheet (CSV/TSV) or input$geo_series_matrix (GEO *_series_matrix.txt.gz)")
  if (!is.null(inp$sample_sheet)) chk_file(inp$sample_sheet, "sample_sheet")
  if (!is.null(inp$geo_series_matrix)) chk_file(inp$geo_series_matrix, "geo_series_matrix")
  switch(inp$type,
         idat = chk_file(inp$idat_dir, "idat_dir"),
         intensity = { chk_file(inp$meth_file, "meth_file"); chk_file(inp$unmeth_file, "unmeth_file")
                       if (!is.null(inp$detp_file)) chk_file(inp$detp_file, "detp_file") },
         beta = chk_file(inp$beta_file, "beta_file"),
         bismark = chk_file(inp$bismark_dir, "bismark_dir"))
  if (!inp$array %in% c("auto", "450k", "epic", "epicv2")) mf_stop("input$array must be auto, 450k, epic or epicv2")
  if (!inp$genome %in% c("hg19", "hg38")) mf_stop("input$genome must be hg19 or hg38")
  nm <- cfg$normalization$method
  ok <- c("auto", "noob", "noob_bmiq", "funnorm", "quantile", "swan", "dasen", "bmiq", "none")
  if (!nm %in% ok) mf_stop("normalization$method must be one of: ", paste(ok, collapse = ", "))
  q <- cfg$qc
  for (k in c("detp_threshold", "sample_max_failed_frac", "probe_max_failed_frac", "bs_min_sample_frac"))
    if (!is.numeric(q[[k]]) || q[[k]] <= 0 || q[[k]] > 1) mf_stop("qc$", k, " must be a number in (0, 1]")
  if (!is.numeric(cfg$analysis$fdr) || cfg$analysis$fdr <= 0 || cfg$analysis$fdr >= 1)
    mf_stop("analysis$fdr must be in (0, 1)")
  if (!(cfg$analysis$inflation_correction %||% "none") %in% c("bacon", "none")) mf_stop("analysis$inflation_correction must be bacon or none")
  if (!cfg$covariates$cell_type %in% c("auto", "blood", "brain", "none"))
    mf_stop("covariates$cell_type must be auto, blood, brain or none")
  cfg
}

# ---- reproducibility -------------------------------------------------------

file_sha256 <- function(f) {
  if (!file.exists(f) || dir.exists(f)) return(NA_character_)
  digest::digest(file = f, algo = "sha256")
}

#' Write a JSON run manifest (inputs, outputs, checksums, parameters, session)
#' @param path output JSON
#' @param cfg config list
#' @param inputs character vector of input files
#' @param outputs character vector of output files
#' @param extra named list of extra fields
#' @export
mf_write_manifest <- function(path, cfg, inputs = character(), outputs = character(), extra = list()) {
  si <- utils::sessionInfo()
  pk <- c(names(si$otherPkgs), names(si$loadedOnly))
  vers <- vapply(sort(unique(pk)), function(p) as.character(utils::packageVersion(p)), "")
  cfg_clean <- cfg; cfg_clean$config_dir <- NULL
  man <- list(
    tool = "methena", version = as.character(utils::packageVersion("methena")),
    started = .mf_env$start_time %||% NA, finished = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    R = R.version.string, platform = R.version$platform,
    parameters = cfg_clean,
    parameter_digest = digest::digest(cfg_clean, algo = "sha256"),
    inputs = lapply(inputs[!is.na(inputs)], function(f) list(path = f, sha256 = file_sha256(f),
                                                             bytes = file.size(f))),
    outputs = lapply(outputs[file.exists(outputs)], function(f) list(path = basename(f), sha256 = file_sha256(f),
                                                                    bytes = file.size(f))),
    step_seconds = .mf_env$timings, packages = as.list(vers))
  man <- c(man, extra)
  jsonlite::write_json(man, path, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
  invisible(path)
}

#' Write a data.frame as a tab-separated file (optionally gzipped)
#' @param df data.frame
#' @param path output path
#' @export
mf_write_tsv <- function(df, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  data.table::fwrite(df, path, sep = "\t", quote = FALSE, na = "NA")
  invisible(path)
}

#' Reset per-run timing state
#' @keywords internal
#' @export
mf_reset <- function() {
  .mf_env$timings <- list()
  .mf_env$start_time <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
}
