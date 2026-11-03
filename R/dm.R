#' Load a methylflow dataset (RDS written by [mf_preprocess])
#' @param path path to methylflow_dataset.rds
#' @export
mf_load_dataset <- function(path) {
  if (!file.exists(path)) mf_stop("Dataset not found: ", path, ". Run the preprocessing step first.")
  se <- readRDS(path)
  if (!methods::is(se, "SummarizedExperiment") || !"beta" %in% SummarizedExperiment::assayNames(se))
    mf_stop(path, " is not a methylflow dataset")
  se
}

#' Resolve analysis covariates and build the design matrix
#'
#' Covariate keywords: 'age', 'sex', 'cells' (all estimated cell proportions except the
#' most abundant, to avoid collinearity), 'smoking', plus any sample-sheet column names.
#' Covariates that are absent, constant or entirely missing are skipped with a warning.
#' Samples with missing values in the group or any used covariate are dropped.
#' @param cd sample data.frame (must contain factor 'group' with 2 levels)
#' @param cfg config
#' @return list(design, cd, covariates, formula)
#' @export
mf_design <- function(cd, cfg) {
  a <- cfg$analysis
  if (is.null(cd$group)) mf_stop("No 'group' column: set phenotype$group_col in the config")
  cd$group <- droplevels(factor(cd$group))
  if (nlevels(cd$group) != 2) mf_stop("Case-control analysis needs exactly 2 groups; found: ", paste(levels(cd$group), collapse = ", "))
  used <- character()
  for (cv in unique(c(a$covariates, a$batch_cols, cfg$phenotype$extra_covariates))) {
    cols <- switch(cv,
      age = "age", sex = if (!is.null(cd$sex) && length(unique(stats::na.omit(cd$sex))) > 1) "sex" else if (!is.null(cd$predicted_sex)) "predicted_sex" else "sex",
      cells = { cc <- grep("^cell_", names(cd), value = TRUE)
                if (length(cc)) setdiff(cc, cc[which.max(colMeans(cd[, cc, drop = FALSE], na.rm = TRUE))]) else character() },
      smoking = "smoking_score",
      age_or_dnam = { if (!is.null(cd$age) && !is.null(cd$dnam_age)) { cd$age_or_dnam <- ifelse(is.na(cd$age), cd$dnam_age, cd$age) }; "age_or_dnam" },
      cv)
    for (cc in cols) {
      if (!cc %in% names(cd)) { mf_log("Covariate '", cc, "' not available - skipped", level = "WARN"); next }
      v <- cd[[cc]]
      if (all(is.na(v)) || length(unique(stats::na.omit(v))) < 2) { mf_log("Covariate '", cc, "' is constant/missing - skipped", level = "WARN"); next }
      if (cv %in% c(a$batch_cols) || is.character(v)) cd[[cc]] <- factor(v)
      used <- c(used, cc)
    }
  }
  used <- unique(used)
  keep <- stats::complete.cases(cd[, c("group", used), drop = FALSE])
  if (any(!keep)) mf_log("Dropping ", sum(!keep), " samples with missing group/covariates")
  cd <- cd[keep, , drop = FALSE]
  for (cc in used) if (is.factor(cd[[cc]])) cd[[cc]] <- droplevels(cd[[cc]])
  used <- used[vapply(used, function(cc) length(unique(cd[[cc]])) > 1, TRUE)]
  fml <- stats::as.formula(paste("~ group", if (length(used)) paste("+", paste(sprintf("`%s`", used), collapse = " + ")) else ""))
  design <- stats::model.matrix(fml, data = cd)
  if (qr(design)$rank < ncol(design)) {
    mf_log("Design is rank-deficient; removing aliased columns", level = "WARN")
    q <- qr(design); design <- design[, q$pivot[seq_len(q$rank)], drop = FALSE]
  }
  colnames(design) <- make.names(colnames(design))
  if (ncol(design) >= nrow(design) - 1) mf_stop("Too many covariates for ", nrow(design), " samples")
  mf_log("Model: ", paste(deparse(fml, width.cutoff = 500L), collapse = " "), " | n = ", nrow(cd), " (", paste(names(table(cd$group)), table(cd$group), sep = "=", collapse = ", "), ")")
  list(design = design, cd = cd, covariates = used, formula = paste(deparse(fml, width.cutoff = 500L), collapse = " "), coef = 2)
}

add_residual_pcs <- function(des, M, n) {
  if (n <= 0) return(des)
  fit <- limma::lmFit(M, des$design)
  res <- M - fit$coefficients %*% t(des$design)
  v <- matrixStats::rowVars(res); sel <- order(v, decreasing = TRUE)[seq_len(min(20000, nrow(res)))]
  pc <- stats::prcomp(t(res[sel, ]), center = TRUE, rank. = n)$x
  colnames(pc) <- paste0("resPC", seq_len(ncol(pc)))
  des$design <- cbind(des$design, pc); des$cd <- cbind(des$cd, pc); des$covariates <- c(des$covariates, colnames(pc))
  mf_log("Added ", n, " residual principal components")
  des
}

genomic_lambda <- function(p) { p <- p[!is.na(p)]; stats::median(stats::qchisq(p, 1, lower.tail = FALSE)) / stats::qchisq(0.5, 1) }

#' Differentially methylated positions with limma
#' @param se dataset
#' @param cfg config
#' @return list(table, design, lambda, M (imputed M-values), fit)
#' @export
mf_dmp <- function(se, cfg) {
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  des <- mf_design(cd, cfg)
  B <- SummarizedExperiment::assay(se, "beta")[, rownames(des$cd), drop = FALSE]
  nmiss <- sum(is.na(B))
  B <- row_impute(B)
  if (nmiss) mf_log("Imputed ", nmiss, " missing beta values with CpG medians for model fitting")
  w <- NULL
  if ("coverage" %in% SummarizedExperiment::assayNames(se)) {
    # sequencing: M-values from read counts with a 1-read offset, coverage as precision weight
    cv <- SummarizedExperiment::assay(se, "coverage")[, rownames(des$cd), drop = FALSE]
    cv[is.na(cv)] <- 0
    mc <- B * cv
    M <- log2((mc + 1) / (cv - mc + 1))
    w <- sqrt(pmax(cv, 1)); w <- w / mean(w)
  } else M <- mf_beta2m(B)
  des <- add_residual_pcs(des, M, cfg$analysis$n_pcs)
  fit <- limma::eBayes(limma::lmFit(M, des$design, weights = w), robust = TRUE)
  tt <- limma::topTable(fit, coef = des$coef, number = Inf, sort.by = "none")
  fitb <- limma::lmFit(B, des$design)
  g <- des$cd$group
  ann <- as.data.frame(SummarizedExperiment::rowData(se))[rownames(M), , drop = FALSE]
  res <- data.frame(cpg = rownames(M), chr = ann$chr, pos = ann$pos, gene = ann$gene, gene_group = ann$gene_group,
                    island = ann$island, enhancer = ann$enhancer,
                    logFC_M = tt$logFC, delta_beta_adj = fitb$coefficients[, des$coef],
                    mean_beta_control = rowMeans(B[, g == levels(g)[1], drop = FALSE]),
                    mean_beta_case = rowMeans(B[, g == levels(g)[2], drop = FALSE]),
                    t = tt$t, P.Value = tt$P.Value, FDR = tt$adj.P.Val,
                    bonferroni = pmin(1, tt$P.Value * nrow(tt)), stringsAsFactors = FALSE)
  res$delta_beta <- res$mean_beta_case - res$mean_beta_control
  res$direction <- ifelse(res$t > 0, "hyper", "hypo")
  lam <- genomic_lambda(res$P.Value)
  bac <- NULL
  if (identical(cfg$analysis$inflation_correction %||% "none", "bacon")) {
    if (requireNamespace("bacon", quietly = TRUE)) {
      # convert moderated t to z, then estimate the empirical null (van Iterson et al. 2017)
      df <- fit$df.total   # same row order as res (topTable sort.by = "none")
      # exact t -> z conversion on the log scale (avoids pt() rounding to 1 for large |t|)
      z <- -sign(res$t) * stats::qnorm(stats::pt(-abs(res$t), df, log.p = TRUE), log.p = TRUE)
      set.seed(cfg$seed)
      b <- bacon::bacon(z)
      bac <- c(inflation = unname(bacon::inflation(b)), bias = unname(bacon::bias(b)))
      if (bac[1] > 1) {
        zc <- as.vector(bacon::tstat(b))
        res$P.uncorrected <- res$P.Value; res$FDR.uncorrected <- res$FDR
        res$P.Value <- 2 * stats::pnorm(-abs(zc)); res$FDR <- stats::p.adjust(res$P.Value, "BH")
        res$bonferroni <- pmin(1, res$P.Value * nrow(res))
        mf_log(sprintf("bacon: inflation = %.3f, bias = %.3f; P-values corrected for test-statistic inflation", bac[1], bac[2]))
      } else {
        # never let the correction make results *more* significant (deflation estimates are unreliable in small studies)
        mf_log(sprintf("bacon: inflation = %.3f (<= 1), bias = %.3f; no correction applied (uncorrected P-values kept)", bac[1], bac[2]))
        bac <- NULL
      }
    } else mf_log("bacon not installed: no inflation correction applied", level = "WARN")
  }
  res$significant <- res$FDR < cfg$analysis$fdr & abs(res$delta_beta_adj) >= cfg$analysis$min_delta_beta
  if (!is.null(bac)) mf_log("DMPs before bacon correction: ", sum(res$FDR.uncorrected < cfg$analysis$fdr))
  mf_log("DMPs at FDR < ", cfg$analysis$fdr, ": ", sum(res$significant), " (", sum(res$significant & res$direction == "hyper"), " hyper / ",
         sum(res$significant & res$direction == "hypo"), " hypo); lambda = ", round(lam, 3))
  res <- res[order(res$P.Value), ]
  list(table = res, design = des, lambda = lam, bacon = bac, lambda_corrected = genomic_lambda(res$P.Value), M = M, beta = B, fit = fit, weights = w)
}

#' Variably methylated positions (DiffVar, missMethyl::varFit)
#' @param dmp result of [mf_dmp]
#' @param cfg config
#' @export
mf_vmp <- function(dmp, cfg) {
  need_pkg("missMethyl", "variance testing")
  des <- dmp$design$design
  fv <- missMethyl::varFit(dmp$M, design = des, coef = c(1, 2))
  tv <- missMethyl::topVar(fv, coef = 2, number = nrow(dmp$M), sort = FALSE)
  res <- data.frame(cpg = rownames(tv), log_var_ratio = tv$LogVarRatio, t = tv$t, P.Value = tv$P.Value, FDR = tv$Adj.P.Value,
                    sd_control = matrixStats::rowSds(dmp$beta[rownames(tv), dmp$design$cd$group == levels(dmp$design$cd$group)[1]]),
                    sd_case = matrixStats::rowSds(dmp$beta[rownames(tv), dmp$design$cd$group == levels(dmp$design$cd$group)[2]]))
  res$direction <- ifelse(res$log_var_ratio > 0, "more_variable_in_cases", "less_variable_in_cases")
  res$significant <- res$FDR < cfg$analysis$fdr
  a <- dmp$table[match(res$cpg, dmp$table$cpg), c("chr", "pos", "gene", "gene_group", "island")]
  res <- cbind(res[, 1, drop = FALSE], a, res[, -1])
  mf_log("VMPs at FDR < ", cfg$analysis$fdr, ": ", sum(res$significant), "; lambda = ", round(genomic_lambda(res$P.Value), 3))
  res[order(res$P.Value), ]
}

#' Differentially methylated regions (DMRcate kernel smoothing on limma statistics)
#'
#' Works for arrays and sequencing data: the per-CpG limma t-statistics from [mf_dmp]
#' are passed to DMRcate as a CpGannotated object.
#' @param dmp result of [mf_dmp]
#' @param cfg config
#' @param genome hg19 or hg38
#' @export
mf_dmr <- function(dmp, cfg, genome = "hg19") {
  need_pkg("DMRcate", "DMR calling")
  tt <- dmp$table[!is.na(dmp$table$chr), ]
  tt <- tt[order(tt$chr, tt$pos), ]
  nsig <- sum(tt$FDR < cfg$analysis$fdr)
  if (nsig == 0) { mf_log("No individually significant CpGs: DMRcate cannot seed regions; skipping DMRs", level = "WARN"); return(NULL) }
  gr <- GenomicRanges::GRanges(tt$chr, IRanges::IRanges(tt$pos, width = 1), stat = tt$t, rawpval = tt$P.Value,
                               diff = tt$delta_beta_adj, ind.fdr = tt$FDR, is.sig = tt$FDR < cfg$analysis$fdr)
  names(gr) <- tt$cpg
  ann <- methods::new("CpGannotated", ranges = gr)
  d <- tryCatch(DMRcate::dmrcate(ann, lambda = cfg$analysis$dmr_lambda, C = cfg$analysis$dmr_C, min.cpgs = cfg$analysis$dmr_min_cpgs),
                error = function(e) { mf_log("dmrcate: ", conditionMessage(e), level = "WARN"); NULL })
  if (is.null(d)) return(NULL)
  # build the region table directly from the DMResults slots (no ExperimentHub download needed)
  co <- d@coord
  df <- data.frame(chr = sub(":.*", "", co), start = as.integer(sub(".*:(\\d+)-.*", "\\1", co)),
                   end = as.integer(sub(".*-", "", co)), no.cpgs = d@no.cpgs, min_smoothed_fdr = d@min_smoothed_fdr,
                   Stouffer = d@Stouffer, HMFDR = d@HMFDR, Fisher = d@Fisher, maxdiff = d@maxdiff, meandiff = d@meandiff)
  df$width <- df$end - df$start + 1
  df$overlapping.genes <- region_genes(df, genome)
  df$region <- paste0(df$chr, ":", df$start, "-", df$end)
  df$direction <- ifelse(df$meandiff > 0, "hyper", "hypo")
  df <- df[order(df$min_smoothed_fdr), ]
  mf_log("DMRs: ", nrow(df))
  df
}

#' Gene-level summary with predicted transcriptional consequence
#'
#' Promoter CpGs (TSS1500, TSS200, 5'UTR, 1stExon) gaining methylation are scored as
#' 'promoter hypermethylation (likely reduced transcription)'; losing methylation as
#' 'promoter hypomethylation (likely increased transcription)'. Gene-body changes are reported
#' separately because body methylation correlates positively with expression.
#' @param dmp_table DMP table
#' @param fdr significance threshold
#' @param top_n if fewer than 20 CpGs are significant, use the top_n CpGs by P
#' @export
mf_gene_level <- function(dmp_table, fdr = 0.05, top_n = 1000) {
  sig <- dmp_table[dmp_table$FDR < fdr, ]
  basis <- "FDR"
  if (nrow(sig) < 20) { sig <- utils::head(dmp_table[order(dmp_table$P.Value), ], top_n); basis <- paste0("top", top_n) }
  sig <- sig[!is.na(sig$gene) & sig$gene != "", ]
  if (!nrow(sig)) return(data.frame())
  genes <- strsplit(sig$gene, ";"); groups <- strsplit(ifelse(is.na(sig$gene_group), "", sig$gene_group), ";")
  long <- do.call(rbind, lapply(seq_along(genes), function(i) {
    g <- genes[[i]]; gr <- groups[[i]]; if (length(gr) < length(g)) gr <- rep_len(c(gr, ""), length(g))
    data.frame(gene = g, region = gr[seq_along(g)], cpg = sig$cpg[i], db = sig$delta_beta_adj[i], t = sig$t[i], P = sig$P.Value[i])
  }))
  long <- unique(long)
  long$promoter <- long$region %in% c("TSS1500", "TSS200", "5'UTR", "1stExon", "5URT", "ExonBnd")
  agg <- do.call(rbind, lapply(split(long, long$gene), function(d) {
    pr <- d[d$promoter, ]; bo <- d[!d$promoter, ]
    data.frame(gene = d$gene[1], n_cpgs = length(unique(d$cpg)), n_promoter = nrow(pr), n_body_other = nrow(bo),
               promoter_mean_delta_beta = if (nrow(pr)) mean(pr$db) else NA_real_,
               body_mean_delta_beta = if (nrow(bo)) mean(bo$db) else NA_real_,
               min_p = min(d$P), cpgs = paste(utils::head(unique(d$cpg[order(d$P)]), 5), collapse = ";"))
  }))
  agg$predicted_effect <- ifelse(is.na(agg$promoter_mean_delta_beta), "gene body / other only (weak or positive link to expression)",
                                 ifelse(agg$promoter_mean_delta_beta > 0, "promoter hypermethylation: likely reduced transcription",
                                        "promoter hypomethylation: likely increased transcription"))
  agg$basis <- basis
  agg[order(agg$min_p), ]
}

region_genes <- function(df, genome) {
  txpkg <- paste0("TxDb.Hsapiens.UCSC.", genome, ".knownGene")
  if (!nrow(df) || !requireNamespace(txpkg, quietly = TRUE) || !requireNamespace("org.Hs.eg.db", quietly = TRUE)) return(rep(NA_character_, nrow(df)))
  genes <- suppressMessages(GenomicFeatures::genes(get(txpkg, envir = asNamespace(txpkg)), single.strand.genes.only = TRUE))
  prom <- GenomicRanges::promoters(genes, upstream = 1500, downstream = 0)
  span <- GenomicRanges::punion(genes, prom, fill.gap = TRUE, ignore.strand = TRUE)
  names(span) <- names(genes)
  sym <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, names(genes), "SYMBOL", "ENTREZID"))
  gr <- GenomicRanges::GRanges(df$chr, IRanges::IRanges(df$start, df$end))
  ov <- GenomicRanges::findOverlaps(gr, span, ignore.strand = TRUE)
  s <- split(unname(sym[names(span)[S4Vectors::subjectHits(ov)]]), S4Vectors::queryHits(ov))
  out <- rep(NA_character_, nrow(df))
  out[as.integer(names(s))] <- vapply(s, function(x) paste(unique(stats::na.omit(x)), collapse = ", "), "")
  out
}

#' Association of top principal components with group and covariates (confounding diagnostic)
#' @param beta beta matrix
#' @param cd sample data (design$cd)
#' @param vars columns to test
#' @param k number of PCs
#' @export
mf_pc_association <- function(beta, cd, vars, k = 10) {
  p <- top_var_pca(beta[, rownames(cd), drop = FALSE], n = 20000, k = k)
  vars <- intersect(unique(c("group", vars)), names(cd))
  out <- do.call(rbind, lapply(vars, function(v) do.call(rbind, lapply(seq_len(ncol(p$pcs)), function(i) {
    x <- cd[[v]]; y <- p$pcs[, i]
    pv <- tryCatch(if (is.numeric(x)) stats::cor.test(x, y)$p.value else stats::kruskal.test(y, factor(x))$p.value, error = function(e) NA)
    data.frame(variable = v, pc = colnames(p$pcs)[i], var_explained = p$var_explained[i], p = pv) }))))
  out
}
