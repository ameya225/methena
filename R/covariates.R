# Robust smoking-associated CpGs (all hypomethylated in current smokers; Zeilinger 2013,
# Joehanes 2016). Score = -mean of per-CpG z-scores (higher = more smoking-like).
.smoking_cpgs <- c(AHRR = "cg05575921", ALPP2_a = "cg21566642", ALPP2_b = "cg01940273",
                   F2RL3 = "cg03636183", GPR15 = "cg19859270", GFI1 = "cg09935388",
                   chr6p21 = "cg06126421", AHRR_b = "cg21161138", chr2q37 = "cg05951221",
                   RARA = "cg17739917")

#' Methylation-derived smoking score
#' @param beta beta matrix with array probe IDs (EPICv2 suffixes are stripped)
#' @export
mf_smoking_score <- function(beta) {
  ids <- sub("_.*$", "", rownames(beta))
  hit <- match(.smoking_cpgs, ids)
  if (sum(!is.na(hit)) < 3) return(NULL)
  b <- beta[stats::na.omit(hit), , drop = FALSE]
  z <- t(scale(t(b)))
  data.frame(smoking_score = -colMeans(z, na.rm = TRUE),
             AHRR_cg05575921 = if (!is.na(hit[1])) beta[hit[1], ] else NA_real_,
             smoking_n_cpgs = sum(!is.na(hit)))
}

#' Estimate cell-type proportions
#'
#' Blood: EpiDISH robust partial correlations with the 7-cell centDHSbloodDMC.m reference.
#' Brain: HiBED (if installed), else neuron/non-neuron with EpiDISH-style RPC on CETS markers is
#' not available and the step is skipped with a warning.
#' @param beta beta matrix (array probe IDs)
#' @param tissue 'blood' or 'brain'
#' @export
mf_cell_props <- function(beta, tissue) {
  rownames(beta) <- sub("_.*$", "", rownames(beta))
  beta <- beta[!duplicated(rownames(beta)), , drop = FALSE]
  beta <- row_impute(beta)
  if (tissue == "blood") {
    need_pkg("EpiDISH", "blood cell-type deconvolution")
    env <- new.env(); utils::data("centDHSbloodDMC.m", package = "EpiDISH", envir = env)
    ref <- env$centDHSbloodDMC.m
    if (sum(rownames(ref) %in% rownames(beta)) < 50) { mf_log("Too few reference CpGs for blood deconvolution", level = "WARN"); return(NULL) }
    est <- EpiDISH::epidish(beta.m = beta, ref.m = ref, method = "RPC")$estF
    colnames(est) <- paste0("cell_", colnames(est))
    return(as.data.frame(est))
  }
  if (tissue == "brain") {
    if (!requireNamespace("HiBED", quietly = TRUE)) { mf_log("HiBED not installed: skipping brain deconvolution", level = "WARN"); return(NULL) }
    est <- tryCatch(HiBED::HiBED_deconvolution(beta, h = 1), error = function(e) { mf_log("HiBED failed: ", conditionMessage(e), level = "WARN"); NULL })
    if (is.null(est)) return(NULL)
    est <- as.data.frame(est)
    if (max(colSums(est), na.rm = TRUE) > 2) est <- est / 100
    colnames(est) <- paste0("cell_", clean_names(colnames(est)))
    return(est)
  }
  NULL
}

#' Add estimated covariates (cell composition, smoking score) to sample data
#' @param beta normalised beta matrix
#' @param cd sample data.frame
#' @param cfg config
#' @param data_type input type
#' @export
mf_add_covariates <- function(beta, cd, cfg, data_type) {
  if (data_type == "bismark") return(cd)
  tissue <- cfg$covariates$cell_type
  if (tissue == "auto") tissue <- if (tolower(cfg$input$tissue %||% "") %in% c("blood", "brain")) tolower(cfg$input$tissue) else "none"
  if (tissue != "none") {
    cp <- mf_cell_props(beta, tissue)
    if (!is.null(cp)) { cd <- cbind(cd[, setdiff(names(cd), names(cp)), drop = FALSE], cp[rownames(cd), , drop = FALSE])
      mf_log("Estimated ", ncol(cp), " cell-type proportions (", tissue, ")") }
  }
  if (isTRUE(cfg$covariates$dnam_age %||% TRUE)) {
    clk <- switch(tissue, blood = "ZhangQ2019", brain = "ShirebyG2020", "HorvathS2013")
    da <- tryCatch(mf_dnam_age(beta, clk), error = function(e) { mf_log("DNAm age covariate failed: ", conditionMessage(e), level = "WARN"); NULL })
    if (!is.null(da)) { cd$dnam_age <- da[rownames(cd)]; cd$dnam_age_clock <- clk
      mf_log("Added DNAm-predicted age (", clk, ") as covariate 'dnam_age'") }
  }
  if (isTRUE(cfg$covariates$smoking_score)) {
    sm <- mf_smoking_score(beta)
    if (!is.null(sm)) cd <- cbind(cd[, setdiff(names(cd), names(sm)), drop = FALSE], sm[rownames(cd), , drop = FALSE])
  }
  cd
}

#' DNAm-predicted age from a single clock (used as a covariate when reported age is missing)
#' @param beta beta matrix with array probe IDs
#' @param clock dnaMethyAge clock name
#' @export
mf_dnam_age <- function(beta, clock = "ZhangQ2019") {
  need_pkg("dnaMethyAge", "DNAm age")
  if (!"package:dnaMethyAge" %in% search()) suppressPackageStartupMessages(do.call("require", list("dnaMethyAge", character.only = TRUE)))
  rownames(beta) <- sub("_.*$", "", rownames(beta)); beta <- row_impute(beta[!duplicated(rownames(beta)), , drop = FALSE])
  junk <- utils::capture.output(r <- suppressWarnings(suppressMessages(dnaMethyAge::methyAge(beta, clock = clock, do_plot = FALSE, use_cores = 1))))
  stats::setNames(r$mAge, r$Sample)
}
