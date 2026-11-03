# Clock catalogue (dnaMethyAge engine, coefficients bundled offline). 'age' clocks are in years
# and are benchmarked against chronological age; the others are rates, scores or telomere length.
.clock_info <- data.frame(
  clock = c("HorvathS2013", "HannumG2013", "LevineM2018", "HorvathS2018", "ZhangQ2019", "ShirebyG2020",
            "BernabeuE2023c", "McEwenL2019", "PCHorvathS2013", "PCHannumG2013", "PCHorvathS2018", "PCPhenoAge",
            "PCGrimAge", "LuA2023p1", "LuA2023p2", "LuA2023p3", "CBL_common", "Cortex_common",
            "DunedinPACE", "LuA2019", "YangZ2016", "epiTOC2", "ZhangY2017", "CBL_specific", "FuentealbaM2025"),
  label = c("Horvath pan-tissue (2013)", "Hannum blood (2013)", "PhenoAge (Levine 2018)", "Skin & blood (Horvath 2018)",
            "Zhang EN (2019)", "Cortical clock (Shireby 2020)", "cAge (Bernabeu 2023)", "PedBE paediatric (McEwen 2019)",
            "PC-Horvath", "PC-Hannum", "PC-SkinBlood", "PC-PhenoAge", "PC-GrimAge", "Universal clock 1 (Lu 2023)",
            "Universal clock 2 (Lu 2023)", "Universal clock 3 (Lu 2023)", "Cerebellum common", "Cortex common",
            "DunedinPACE (pace of ageing)", "DNAm telomere length (Lu 2019)", "epiTOC mitotic clock (Yang 2016)",
            "epiTOC2 mitotic clock", "Mortality risk score (Zhang 2017)", "Cerebellum specific", "Fuentealba 2025"),
  type = c(rep("age", 18), "pace", "telomere", "mitotic", "mitotic", "risk score", "age", "other"),
  stringsAsFactors = FALSE)
.age_clocks <- .clock_info$clock[.clock_info$type == "age"]

#' Available epigenetic clocks
#' @export
mf_available_clocks <- function() .clock_info

clock_beta_matrix <- function(se) {
  B <- SummarizedExperiment::assay(se, "beta")
  ann <- as.data.frame(SummarizedExperiment::rowData(se))
  if (!is.null(ann$array_probe)) {           # sequencing data mapped to array IDs
    keep <- !is.na(ann$array_probe)
    B <- B[keep, , drop = FALSE]; rownames(B) <- ann$array_probe[keep]
    B <- B[!duplicated(rownames(B)), , drop = FALSE]
  }
  extra <- S4Vectors::metadata(se)$beta_filtered_probes
  if (!is.null(extra) && nrow(extra)) B <- rbind(B, extra[setdiff(rownames(extra), rownames(B)), colnames(B), drop = FALSE])
  rownames(B) <- sub("_.*$", "", rownames(B))           # EPICv2 suffixes
  B[!duplicated(rownames(B)), , drop = FALSE]
}

#' Compute epigenetic clocks
#'
#' Uses the dnaMethyAge engine (coefficients bundled with the package, no internet needed):
#' first-generation (Horvath, Hannum, skin & blood, Zhang), second-generation (PhenoAge,
#' PC-GrimAge, DunedinPACE), principal-component versions with improved reliability,
#' tissue-specific (cortical, cerebellar), mitotic (epiTOC) and DNAm telomere length. Missing
#' clock CpGs are imputed from the package's reference medians; the fraction of clock CpGs
#' actually measured is reported so that low-coverage clocks can be discounted. SNP /
#' cross-reactive probes removed during QC are re-added because the clocks were trained with them.
#' @param se dataset
#' @param cfg config
#' @return list(values = data.frame per sample, coverage = data.frame per clock)
#' @export
mf_clocks <- function(se, cfg) {
  B <- clock_beta_matrix(se)
  B <- row_impute(B)
  cd0 <- as.data.frame(SummarizedExperiment::colData(se))
  sexv <- if (!is.null(cd0$sex) && any(!is.na(cd0$sex))) as.character(cd0$sex) else as.character(cd0$predicted_sex %||% rep(NA, nrow(cd0)))
  age_info <- if (!is.null(cd0$age)) data.frame(Sample = rownames(cd0), Age = cd0$age,
                                                Sex = ifelse(sexv == "F", "Female", ifelse(sexv == "M", "Male", NA))) else NA
  set <- cfg$analysis$clock_set
  out <- data.frame(sample_id = colnames(B), row.names = colnames(B))
  cov <- list()
  if (!requireNamespace("dnaMethyAge", quietly = TRUE)) { mf_log("dnaMethyAge not installed: clocks skipped", level = "WARN")
    return(list(values = out, coverage = data.frame())) }
  # dnaMethyAge resolves its coefficient tables from the attached package environment
  if (!"package:dnaMethyAge" %in% search()) suppressPackageStartupMessages(do.call("require", list("dnaMethyAge", character.only = TRUE)))
  unknown <- setdiff(set, .clock_info$clock)
  if (length(unknown)) mf_log("Unknown clocks ignored: ", paste(unknown, collapse = ", "), " (see mf_available_clocks())", level = "WARN")
  for (k in intersect(set, .clock_info$clock)) {
    e <- new.env()
    suppressWarnings(tryCatch(utils::data(list = k, package = "dnaMethyAge", envir = e), error = function(err) NULL))
    cf <- e$coefs
    pr <- if (is.data.frame(cf) && "Probe" %in% names(cf)) as.character(cf$Probe) else if (is.matrix(cf)) rownames(cf) else if (!is.null(e$epiTOCcpgs)) e$epiTOCcpgs else NULL
    if (length(pr)) {
      pr <- setdiff(pr, c("(Intercept)", "Intercept"))
      cov[[k]] <- data.frame(clock = k, n_cpgs = length(pr), n_missing = sum(!pr %in% rownames(B)),
                             frac_measured = mean(pr %in% rownames(B)))
    }
    r <- tryCatch({
      junk <- utils::capture.output(rr <- suppressWarnings(suppressMessages(
        dnaMethyAge::methyAge(B, clock = k, age_info = if (k %in% c("PCGrimAge", "BernabeuE2023c") && is.data.frame(age_info)) age_info[stats::complete.cases(age_info), ] else NA,
                              do_plot = FALSE, imputation = TRUE, use_cores = max(1L, as.integer(cfg$threads %||% 1L))))))
      rr }, error = function(err) { mf_log("Clock ", k, " failed: ", conditionMessage(err), level = "WARN"); NULL })
    if (is.data.frame(r)) out[[k]] <- r$mAge[match(out$sample_id, r$Sample)]
  }
  covdf <- if (length(cov)) do.call(rbind, cov) else data.frame(clock = character(), n_cpgs = integer(), n_missing = integer(), frac_measured = numeric())
  mf_log("Computed ", ncol(out) - 1, " clocks")
  list(values = out, coverage = covdf)
}

#' Benchmark clocks against chronological age and test age acceleration between groups
#'
#' Age acceleration = residual of clock ~ chronological age (fit in all samples). Group
#' differences are tested with lm(acceleration ~ group + sex + cells + smoking) using the
#' covariates available in the analysis design.
#' @param clocks result of [mf_clocks]
#' @param des design list from [mf_design]
#' @export
mf_clock_benchmark <- function(clocks, des) {
  v <- clocks$values; cd <- des$cd
  common <- intersect(rownames(cd), v$sample_id)
  v <- v[common, , drop = FALSE]; cd <- cd[common, , drop = FALSE]
  cl <- setdiff(names(v), "sample_id")
  cl <- cl[vapply(cl, function(c) sum(!is.na(v[[c]])) >= 5, TRUE)]
  age <- cd$age
  bench <- do.call(rbind, lapply(cl, function(c) {
    x <- v[[c]]; ok <- !is.na(x) & !is.na(age)
    is_age <- c %in% .age_clocks
    r <- if (sum(ok) > 3) stats::cor(x[ok], age[ok]) else NA
    ty <- .clock_info$type[match(c, .clock_info$clock)]
    data.frame(clock = c, label = .clock_info$label[match(c, .clock_info$clock)],
               type = switch(ty, age = "age estimate (years)", pace = "pace of ageing (years/yr)", telomere = "telomere length (kb)",
                             mitotic = "mitotic score", `risk score` = "risk score", "other"),
               n = sum(ok), pearson_r = r,
               MAE_years = if (is_age) mean(abs(x[ok] - age[ok])) else NA,
               median_AE_years = if (is_age) stats::median(abs(x[ok] - age[ok])) else NA,
               mean_bias_years = if (is_age) mean(x[ok] - age[ok]) else NA)
  }))
  covs <- intersect(des$covariates, names(cd)); covs <- setdiff(covs, c("age", "dnam_age", "age_or_dnam"))
  covs <- covs[!grepl("^resPC", covs)]
  accel <- data.frame(sample_id = common, group = cd$group)
  tests <- do.call(rbind, lapply(cl, function(c) {
    x <- v[[c]]; ok <- !is.na(x) & !is.na(age)
    if (sum(ok) < 10) return(NULL)
    acc <- rep(NA_real_, length(x)); acc[ok] <- stats::residuals(stats::lm(x[ok] ~ age[ok]))
    accel[[paste0(c, "_accel")]] <<- acc
    d <- data.frame(acc = acc, group = cd$group, cd[, covs, drop = FALSE])[ok, ]
    f <- stats::as.formula(paste("acc ~ group", if (length(covs)) paste("+", paste(sprintf("`%s`", covs), collapse = "+")) else ""))
    m <- tryCatch(summary(stats::lm(f, data = d))$coefficients, error = function(e) NULL)
    if (is.null(m) || nrow(m) < 2) return(NULL)
    m0 <- summary(stats::lm(acc ~ group, data = d))$coefficients
    d$acc <- d$acc / stats::sd(d$acc)
    ms <- summary(stats::lm(f, data = d))$coefficients
    data.frame(clock = c, effect_case_vs_control = m[2, 1], se = m[2, 2], effect_sd_units = ms[2, 1], se_sd_units = ms[2, 2],
               p_adjusted_model = m[2, 4],
               effect_unadjusted = m0[2, 1], p_unadjusted = m0[2, 4],
               mean_accel_control = mean(acc[ok][d$group == levels(d$group)[1]]),
               mean_accel_case = mean(acc[ok][d$group == levels(d$group)[2]]))
  }))
  if (!is.null(tests) && nrow(tests)) tests$fdr <- stats::p.adjust(tests$p_adjusted_model, "BH")
  bench <- merge(bench, clocks$coverage, by = "clock", all.x = TRUE)
  list(benchmark = bench[order(-abs(bench$pearson_r)), ], group_tests = tests, acceleration = accel,
       values = cbind(v, age = age, group = cd$group))
}
