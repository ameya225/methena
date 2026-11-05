chr_order <- function(chr) { x <- sub("^chr", "", chr); x[x == "X"] <- "23"; x[x == "Y"] <- "24"; suppressWarnings(as.integer(x)) }

dm_figures <- function(res, fd, fdr, title = "DMP") {
  f <- list(); tab <- res
  n <- nrow(tab); exp <- -log10(stats::ppoints(n)); obs <- -log10(sort(tab$P.Value))
  lam <- genomic_lambda(tab$P.Value)
  thin <- c(seq_len(min(n, 5000)), if (n > 5000) sort(sample(seq(5001, n), min(n - 5000, 20000))))
  f$qq <- mf_save_plot(ggplot2::ggplot(data.frame(e = exp[thin], o = obs[thin]), ggplot2::aes(e, o)) + ggplot2::geom_point(size = 0.6) +
    ggplot2::geom_abline(colour = "red") + ggplot2::theme_bw() +
    ggplot2::labs(title = sprintf("%s Q-Q plot (lambda = %.3f)", title, lam), x = "Expected -log10 P", y = "Observed -log10 P"),
    file.path(fd, paste0(tolower(title), "_qq.png")), 5, 5)
  if (!is.null(tab$delta_beta_adj)) {
    set.seed(1); keep <- tab$P.Value < 0.1 | stats::runif(n) < 0.25
    v <- tab[keep, ]; v$sig <- ifelse(v$FDR < fdr, ifelse(v$t > 0, "Hyper (FDR)", "Hypo (FDR)"), "n.s.")
    lab <- utils::head(v[v$FDR < fdr & !is.na(v$gene) & v$gene != "", ], 12); lab$g <- sub(";.*", "", lab$gene)
    f$volcano <- mf_save_plot(ggplot2::ggplot(v, ggplot2::aes(delta_beta_adj, -log10(P.Value), colour = sig)) +
      ggplot2::geom_point(size = 0.5, alpha = 0.6) + ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = g), size = 3, colour = "black", max.overlaps = 30) +
      ggplot2::scale_colour_manual(values = c("Hyper (FDR)" = "#d7301f", "Hypo (FDR)" = "#2b8cbe", "n.s." = "grey70")) + ggplot2::theme_bw() +
      ggplot2::labs(title = "Volcano plot", x = "Adjusted difference in beta (case - control)", y = "-log10 P", colour = NULL),
      file.path(fd, "dmp_volcano.png"), 7, 5)
  }
  m <- tab[!is.na(tab$chr), ]; m$cn <- chr_order(m$chr); m <- m[!is.na(m$cn), ]
  set.seed(2); m <- m[m$P.Value < 0.1 | stats::runif(nrow(m)) < 0.25, ]
  m <- m[order(m$cn, m$pos), ]
  off <- tapply(m$pos, m$cn, max); cum <- c(0, cumsum(as.numeric(off)))[seq_along(off)]; names(cum) <- names(off)
  m$x <- m$pos + cum[as.character(m$cn)]
  mid <- tapply(m$x, m$cn, function(z) mean(range(z)))
  thr <- if (any(tab$FDR < fdr)) max(tab$P.Value[tab$FDR < fdr]) else NA
  p <- ggplot2::ggplot(m, ggplot2::aes(x, -log10(P.Value), colour = factor(cn %% 2))) + ggplot2::geom_point(size = 0.4) +
    ggplot2::scale_colour_manual(values = c("#08519c", "#6baed6"), guide = "none") +
    ggplot2::scale_x_continuous(breaks = mid, labels = names(mid)) + ggplot2::theme_bw() +
    ggplot2::labs(title = paste(title, "Manhattan plot"), x = "Chromosome", y = "-log10 P")
  if (!is.na(thr)) p <- p + ggplot2::geom_hline(yintercept = -log10(thr), linetype = 2, colour = "red")
  f$manhattan <- mf_save_plot(p, file.path(fd, paste0(tolower(title), "_manhattan.png")), 10, 4)
  f
}

pcassoc_figure <- function(pa, fd) {
  pa$pc <- factor(pa$pc, levels = unique(pa$pc)); pa$lp <- pmin(-log10(pa$p), 20)
  lab <- paste0(levels(pa$pc), "\n", round(100 * pa$var_explained[match(levels(pa$pc), pa$pc)], 1), "%")
  mf_save_plot(ggplot2::ggplot(pa, ggplot2::aes(pc, variable, fill = lp)) + ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(p < 0.05, formatC(p, format = "e", digits = 0), "")), size = 2.6) +
    ggplot2::scale_fill_gradient(low = "white", high = "#b30000", name = "-log10 P") + ggplot2::scale_x_discrete(labels = lab) +
    ggplot2::theme_minimal() + ggplot2::labs(title = "Which variables drive the main axes of methylation variation?",
    subtitle = "Association of the top principal components with group and covariates", x = NULL, y = NULL),
    file.path(fd, "pc_covariate_association.png"), 10, 0.5 + 0.35 * length(unique(pa$variable)) + 1.5)
}

top_boxplots <- function(dmp, fd, n = 6) {
  top <- utils::head(dmp$table, n)
  d <- do.call(rbind, lapply(seq_len(nrow(top)), function(i) data.frame(
    cpg = paste0(top$cpg[i], "\n", sub(";.*", "", ifelse(is.na(top$gene[i]) | top$gene[i] == "", "intergenic", top$gene[i])),
                 "  P=", signif(top$P.Value[i], 2)),
    beta = dmp$beta[top$cpg[i], ], group = dmp$design$cd$group)))
  mf_save_plot(ggplot2::ggplot(d, ggplot2::aes(group, beta, fill = group)) + ggplot2::geom_boxplot(outlier.size = 0.3) +
    ggplot2::facet_wrap(~cpg, scales = "free_y") + ggplot2::theme_bw() + ggplot2::theme(legend.position = "none") +
    ggplot2::labs(title = "Top differentially methylated CpGs", x = NULL, y = "Beta"), file.path(fd, "top_dmp_boxplots.png"), 9, 6)
}

clock_figures <- function(cb, fd) {
  f <- list(); v <- cb$values
  ac <- intersect(cb$benchmark$clock[cb$benchmark$type == "age estimate (years)"], names(v))
  if (length(ac) && any(!is.na(v$age))) {
    d <- do.call(rbind, lapply(ac, function(c) data.frame(clock = c, age = v$age, pred = v[[c]], group = v$group)))
    d <- d[!is.na(d$pred), ]
    f$scatter <- mf_save_plot(ggplot2::ggplot(d, ggplot2::aes(age, pred, colour = group)) + ggplot2::geom_point(size = 0.6, alpha = 0.6) +
      ggplot2::geom_abline(linetype = 2) + ggplot2::facet_wrap(~clock, scales = "free_y") + ggplot2::theme_bw() +
      ggplot2::labs(title = "Epigenetic age vs chronological age", x = "Chronological age (years)", y = "DNAm age (years)", colour = "Group"),
      file.path(fd, "clock_vs_age.png"), 10, 7)
  }
  if (!is.null(cb$group_tests) && nrow(cb$group_tests)) {
    g <- cb$group_tests; g$clock <- factor(g$clock, levels = g$clock[order(g$effect_sd_units)])
    f$forest <- mf_save_plot(ggplot2::ggplot(g, ggplot2::aes(effect_sd_units, clock)) +
      ggplot2::geom_vline(xintercept = 0, linetype = 2) + ggplot2::geom_point(ggplot2::aes(colour = p_adjusted_model < 0.05)) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = effect_sd_units - 1.96 * se_sd_units, xmax = effect_sd_units + 1.96 * se_sd_units), height = 0.2) +
      ggplot2::scale_colour_manual(values = c(`TRUE` = "#d7301f", `FALSE` = "black"), name = "P < 0.05") +
      ggplot2::theme_bw() + ggplot2::labs(title = "Age acceleration, case vs control", subtitle = "covariate-adjusted difference, 95% CI",
      x = "Difference in acceleration (SD units)", y = NULL), file.path(fd, "clock_group_forest.png"), 7, 4.5)
  }
  f
}

enrich_figure <- function(en, fd, label) {
  all <- do.call(rbind, lapply(setdiff(names(en), "basis"), function(k) {
    x <- en[[k]]; if (is.null(x)) return(NULL); x <- utils::head(x[order(x$P.DE), ], 8)
    data.frame(term = substr(paste0("[", k, "] ", x$term), 1, 70), p = x$P.DE, fdr = x$FDR) }))
  if (is.null(all) || !nrow(all)) return(NULL)
  all$term <- factor(all$term, levels = rev(unique(all$term)))
  mf_save_plot(ggplot2::ggplot(all, ggplot2::aes(-log10(p), term, fill = fdr < 0.05)) + ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#d95f0e", `FALSE` = "grey60")) + ggplot2::theme_bw() +
    ggplot2::labs(title = paste(label, "gene-set enrichment (probe-bias corrected)"), x = "-log10 P", y = NULL, fill = "FDR < 0.05"),
    file.path(fd, paste0(tolower(label), "_enrichment.png")), 9, 6)
}

reg_figure <- function(reg, fd) {
  if (is.null(reg) || !nrow(reg)) return(NULL)
  reg <- reg[is.finite(log2(reg$odds_ratio)), ]; reg$odds_ratio <- pmin(pmax(reg$odds_ratio, 2^-6), 2^6)
  reg$feature <- factor(reg$feature, levels = rev(unique(reg$feature)))
  mf_save_plot(ggplot2::ggplot(reg, ggplot2::aes(log2(odds_ratio), feature, colour = direction)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_pointrange(ggplot2::aes(xmin = pmax(log2(ci_low), -6), xmax = pmin(log2(ci_high), 6)), position = ggplot2::position_dodge(width = 0.5)) +
    ggplot2::scale_colour_manual(values = c(hyper = "#d7301f", hypo = "#2b8cbe")) + ggplot2::theme_bw() +
    ggplot2::labs(title = "Genomic / regulatory context of DMPs vs all tested CpGs", x = "log2 odds ratio (95% CI)", y = NULL, colour = NULL),
    file.path(fd, "regulatory_context.png"), 10, 5)
}

motif_figure <- function(mo, fd) {
  if (is.null(mo) || !nrow(mo)) return(NULL)
  top <- unique(utils::head(mo$tf[order(-mo$neg_log10_padj)], 20))
  d <- mo[mo$tf %in% top, ]; d <- d[!duplicated(paste(d$tf, d$bin)), ]
  mf_save_plot(ggplot2::ggplot(d, ggplot2::aes(bin, tf, colour = log2_enrichment, size = neg_log10_padj)) + ggplot2::geom_point() +
    ggplot2::scale_colour_gradient2(low = "#2b8cbe", high = "#d7301f") + ggplot2::theme_bw() +
    ggplot2::labs(title = "TF motif enrichment near DMPs", subtitle = "JASPAR CORE vertebrate motifs vs GC-matched background CpGs", x = NULL, y = NULL,
                  colour = "log2 enrichment", size = "-log10 adj. P"), file.path(fd, "tf_motifs.png"), 7, 7)
}

fmtp <- function(p) ifelse(is.na(p), "NA", ifelse(p < 1e-3, formatC(p, format = "e", digits = 1), formatC(p, format = "f", digits = 3)))

insight_text <- function(r, cfg) {
  b <- character(); add <- function(...) b <<- c(b, paste0(...))
  cd <- r$dmp$design$cd; lv <- levels(cd$group); tb <- table(cd$group)
  add(sprintf("<b>Cohort.</b> %d samples (%s: %d, %s: %d) analysed with the model <code>%s</code>.", nrow(cd), lv[1], tb[1], lv[2], tb[2], htmltools::htmlEscape(r$dmp$design$formula)))
  if (!is.null(cd$pca_outlier) && any(cd$pca_outlier %in% TRUE)) {
    ot <- table(cd$group[cd$pca_outlier %in% TRUE])
    add(sprintf("<b>Caution.</b> %d sample(s) flagged as PCA/profile outliers were kept (%s). Check whether results persist without them (qc: force_drop_outliers: true).",
                sum(cd$pca_outlier %in% TRUE), paste(names(ot), ot, sep = ": ", collapse = ", ")))
  }
  t <- r$dmp$table; ns <- sum(t$FDR < cfg$analysis$fdr)
  add(sprintf("<b>Differential methylation.</b> %s of %s CpGs differ at FDR &lt; %s (%d hyper-, %d hypomethylated in %s); genomic inflation &lambda; = %.2f%s.",
              format(ns, big.mark = ","), format(nrow(t), big.mark = ","), cfg$analysis$fdr, sum(t$FDR < cfg$analysis$fdr & t$t > 0),
              sum(t$FDR < cfg$analysis$fdr & t$t < 0), lv[2], r$dmp$lambda,
              if (r$dmp$lambda > 1.2) " (inflated: residual confounding or a broad polygenic signal; see the PC diagnostics and consider residual PCs)" else ""))
  if (!is.null(r$dmp$bacon)) add(sprintf("Counts above use bacon-corrected P-values (estimated inflation %.2f, bias %.2f; &lambda; after correction = %.2f). Without correction, %s CpGs would pass FDR &lt; %s.",
              r$dmp$bacon[1], r$dmp$bacon[2], r$dmp$lambda_corrected, format(sum(t$FDR.uncorrected < cfg$analysis$fdr), big.mark = ","), cfg$analysis$fdr))
  if (!is.null(r$pcassoc)) { g <- r$pcassoc[r$pcassoc$variable == "group" & r$pcassoc$p < 0.05, ]
    oth <- r$pcassoc[r$pcassoc$variable != "group" & r$pcassoc$p < 1e-3 & r$pcassoc$pc %in% c("PC1", "PC2", "PC3"), ]
    add("<b>Sources of variation.</b> ", if (nrow(g)) paste0("Group is associated with ", paste(g$pc, collapse = ", "), ". ") else "Group is not associated with any of the top 10 principal components. ",
        if (nrow(oth)) paste0("The leading components track ", paste(unique(oth$variable), collapse = ", "), ", which the model adjusts for.") else "") }
  if (ns > 0) {
    top <- utils::head(t[t$FDR < cfg$analysis$fdr, ], 10)
    add("Strongest positions: ", paste(sprintf("%s (%s, &Delta;&beta; = %+.3f)", top$cpg, ifelse(is.na(top$gene) | top$gene == "", "intergenic", sub(";.*", "", top$gene)), top$delta_beta_adj), collapse = "; "), ".")
    add(sprintf("Median absolute adjusted &Delta;&beta; among significant CpGs is %.3f, i.e. effects are %s.", stats::median(abs(t$delta_beta_adj[t$FDR < cfg$analysis$fdr])),
                if (stats::median(abs(t$delta_beta_adj[t$FDR < cfg$analysis$fdr])) < 0.02) "small (typical of blood EWAS of psychiatric traits)" else "moderate-to-large"))
  }
  if (!is.null(r$dmr) && nrow(r$dmr)) add(sprintf("<b>Regions.</b> %d DMRs (%d hyper / %d hypo); top: %s.", nrow(r$dmr), sum(r$dmr$direction == "hyper"), sum(r$dmr$direction == "hypo"),
                                                  paste(utils::head(sprintf("%s (%s)", r$dmr$region, ifelse(is.na(r$dmr$overlapping.genes), "-", r$dmr$overlapping.genes)), 5), collapse = "; ")))
  if (!is.null(r$vmp)) {
    nv <- sum(r$vmp$significant); ov <- length(intersect(r$vmp$cpg[r$vmp$significant], t$cpg[t$FDR < cfg$analysis$fdr]))
    add(sprintf("<b>Variability.</b> %d variably methylated positions (FDR &lt; %s; %d more variable and %d less variable in %s). %d of them are also DMPs, so variance and mean effects largely mark %s CpGs.",
                nv, cfg$analysis$fdr, sum(r$vmp$significant & r$vmp$log_var_ratio > 0), sum(r$vmp$significant & r$vmp$log_var_ratio < 0), lv[2], ov,
                if (nv > 0 && ov / max(nv, 1) < 0.2) "different" else "overlapping"))
  }
  if (!is.null(r$reg)) {
    s <- r$reg[r$reg$fdr < 0.05, ]
    if (nrow(s)) add("<b>Regulatory context.</b> ", paste(sprintf("%s DMPs are %s in %s (OR = %.2f, FDR = %s)", s$direction, ifelse(s$odds_ratio > 1, "enriched", "depleted"), s$feature, s$odds_ratio, fmtp(s$fdr)), collapse = "; "), ".")
    else add("<b>Regulatory context.</b> No genomic/regulatory feature is significantly over- or under-represented among DMPs (FDR &lt; 0.05).")
  }
  if (!is.null(r$genes) && nrow(r$genes)) {
    up <- r$genes[grepl("hypermethylation", r$genes$predicted_effect), ]; dn <- r$genes[grepl("hypomethylation", r$genes$predicted_effect), ]
    add(sprintf("<b>Transcriptional consequences.</b> %d genes carry promoter hypermethylation (predicted reduced transcription; e.g. %s) and %d promoter hypomethylation (predicted increased transcription; e.g. %s). Basis: %s.",
                nrow(up), paste(utils::head(up$gene, 8), collapse = ", "), nrow(dn), paste(utils::head(dn$gene, 8), collapse = ", "), r$genes$basis[1]))
  }
  if (!is.null(r$expr_int)) add(sprintf("<b>Expression integration.</b> Across %d genes with promoter changes, Spearman &rho; between promoter &Delta;&beta; and expression = %.2f (P = %s); %.0f%% show the canonical inverse relationship.",
                                        r$expr_int$n_promoter_genes, r$expr_int$promoter_rho, fmtp(r$expr_int$promoter_p), 100 * r$expr_int$frac_inverse))
  if (!is.null(r$enrich)) {
    hits <- do.call(rbind, lapply(setdiff(names(r$enrich), "basis"), function(k) { x <- r$enrich[[k]]; if (is.null(x)) NULL else data.frame(k = k, term = x$term, fdr = x$FDR, p = x$P.DE) }))
    if (!is.null(hits) && any(hits$fdr < 0.05)) { h <- utils::head(hits[hits$fdr < 0.05, ][order(hits$p[hits$fdr < 0.05]), ], 8)
      add("<b>Pathways.</b> Significant after probe-number bias correction: ", paste(sprintf("%s [%s]", h$term, h$k), collapse = "; "), ".")
    } else if (!is.null(hits)) add("<b>Pathways.</b> No gene set reaches FDR &lt; 0.05 after correcting for CpGs-per-gene bias (", r$enrich$basis, "); nominal leads: ",
                                  paste(utils::head(hits$term[order(hits$p)], 5), collapse = "; "), ".")
  }
  if (!is.null(r$motif) && nrow(r$motif)) {
    s <- r$motif[r$motif$neg_log10_padj > -log10(0.05) & r$motif$log2_enrichment > 0, ]
    if (nrow(s)) add("<b>Transcription factors.</b> Motifs enriched near DMPs (adj. P &lt; 0.05): ", paste(utils::head(sprintf("%s (%s)", s$tf, s$bin), 12), collapse = ", "),
                     ". Methylation at these TF sites can alter TF occupancy and so link the DMPs to transcriptional regulation.")
    else add("<b>Transcription factors.</b> No JASPAR motif is significantly enriched near DMPs relative to GC-matched background.")
  }
  if (!is.null(r$clocks)) {
    bm <- r$clocks$benchmark; ab <- bm[bm$type == "age estimate (years)" & !is.na(bm$MAE_years), ]
    if (nrow(ab)) { best <- ab[which.min(ab$MAE_years), ]; bestr <- ab[which.max(ab$pearson_r), ]
      add(sprintf("<b>Epigenetic age.</b> Best calibrated clock: %s (MAE %.1f y); best correlated: %s (r = %.2f).", best$clock, best$MAE_years, bestr$clock, bestr$pearson_r)) }
    gt <- r$clocks$group_tests
    if (!is.null(gt) && nrow(gt)) { s <- gt[gt$p_adjusted_model < 0.05, ]
      add(if (nrow(s)) paste0("Compared with controls, the ", lv[2], " group differs in covariate-adjusted acceleration on ", nrow(s), " of ", nrow(gt), " clocks: ",
                              paste(sprintf("%s (%+.2f SD, P = %s)", s$clock, s$effect_sd_units, fmtp(s$p_adjusted_model)), collapse = "; "),
                              ". Positive = older/faster in cases; for DunedinPACE the raw difference is ",
                              if ("DunedinPACE" %in% gt$clock) sprintf("%+.3f biological years per year", gt$effect_case_vs_control[gt$clock == "DunedinPACE"]) else "not available", ".")
          else "No clock shows a significant case-control difference in age acceleration after covariate adjustment.") }
  }
  if (!is.null(r$compare) && nrow(r$compare)) add("<b>Cross-dataset concordance.</b> ", paste(sprintf("%s: %.0f%% sign concordance of %d selected CpGs (binomial P = %s), %d replicate at P &lt; 0.05 in the same direction (%.0f expected by chance)",
                                                   r$compare$comparison, 100 * r$compare$sign_concordance, r$compare$n_selected, fmtp(r$compare$binom_p), r$compare$n_replicated_p05_same_dir, r$compare$expected_by_chance), collapse = "; "), ".")
  paste0("<ul>", paste0("<li>", b, "</li>", collapse = ""), "</ul>")
}

#' Run the downstream analysis on a methena dataset
#'
#' Writes to `<output_dir>/analysis/`: DMP/DMR/VMP tables, clock tables, enrichment,
#' regulatory context, TF motifs, gene-level predictions, cross-dataset comparison, figures,
#' analysis_report.html, results_bundle.rds (for the Shiny app) and run_manifest.json.
#' @param cfg config
#' @param dataset path to methena_dataset.rds (default: <output_dir>/preprocess/methena_dataset.rds)
#' @export
mf_analyze <- function(cfg, dataset = NULL) {
  mf_reset(); set.seed(cfg$seed)
  dataset <- dataset %||% file.path(cfg$output_dir, "preprocess", "methena_dataset.rds")
  od <- file.path(cfg$output_dir, "analysis"); fd <- file.path(od, "figures")
  dir.create(fd, recursive = TRUE, showWarnings = FALSE)
  mf_log_to(file.path(od, "analysis.log"))
  se <- mf_load_dataset(dataset)
  md <- S4Vectors::metadata(se); a <- cfg$analysis
  mf_log("Dataset: ", nrow(se), " CpGs x ", ncol(se), " samples (", md$platform, ")")
  r <- list(); figs <- list()
  r$dmp <- mf_step("differential methylation (DMP)", mf_dmp(se, cfg))
  mf_write_tsv(r$dmp$table, file.path(od, "dmp_results.tsv.gz"))
  r$pcassoc <- mf_try_step("PC-covariate diagnostics", mf_pc_association(r$dmp$beta, r$dmp$design$cd, r$dmp$design$covariates))
  if (!is.null(r$pcassoc)) { mf_write_tsv(r$pcassoc, file.path(od, "pc_covariate_association.tsv")); figs$pcassoc <- pcassoc_figure(r$pcassoc, fd) }
  figs$dmp <- dm_figures(r$dmp$table, fd, a$fdr, "DMP"); figs$box <- top_boxplots(r$dmp, fd)
  if (isTRUE(a$dmr)) { r$dmr <- mf_try_step("differential regions (DMR)", mf_dmr(r$dmp, cfg, md$genome))
    if (!is.null(r$dmr)) mf_write_tsv(r$dmr, file.path(od, "dmr_results.tsv")) }
  if (isTRUE(a$vmp)) { r$vmp <- mf_step("variably methylated positions (VMP)", tryCatch(mf_vmp(r$dmp, cfg), error = function(e) { mf_log("VMP failed: ", conditionMessage(e), level = "WARN"); NULL }))
    if (!is.null(r$vmp)) { mf_write_tsv(r$vmp, file.path(od, "vmp_results.tsv.gz")); figs$vmp <- dm_figures(r$vmp, fd, a$fdr, "VMP") } }
  r$genes <- mf_gene_level(r$dmp$table, a$fdr, a$enrichment_top_n)
  mf_write_tsv(r$genes, file.path(od, "gene_level_regulation.tsv"))
  if (!is.null(a$expression_results)) {
    ex <- as.data.frame(data.table::fread(file.path(cfg$config_dir %||% ".", a$expression_results)))
    r$expr_int <- mf_integrate_expression(r$genes, ex); mf_write_tsv(r$expr_int$table, file.path(od, "methylation_expression_integration.tsv"))
  }
  if (isTRUE(a$enrichment)) {
    r$enrich <- mf_try_step("gene-set enrichment", mf_enrichment(r$dmp$table, se, cfg, "DMP"))
    for (k in setdiff(names(r$enrich), "basis")) if (!is.null(r$enrich[[k]])) mf_write_tsv(r$enrich[[k]], file.path(od, paste0("enrichment_dmp_", k, ".tsv")))
    figs$enrich <- enrich_figure(r$enrich, fd, "DMP")
    if (!is.null(r$vmp)) { r$enrich_vmp <- mf_try_step("VMP gene-set enrichment", mf_enrichment(r$vmp, se, cfg, "VMP"))
      for (k in setdiff(names(r$enrich_vmp), "basis")) if (!is.null(r$enrich_vmp[[k]])) mf_write_tsv(r$enrich_vmp[[k]], file.path(od, paste0("enrichment_vmp_", k, ".tsv"))) }
  }
  r$reg <- mf_try_step("regulatory context", mf_regulatory_enrichment(r$dmp$table, cfg))
  if (!is.null(r$reg)) { mf_write_tsv(r$reg, file.path(od, "regulatory_context_enrichment.tsv")); figs$reg <- reg_figure(r$reg, fd) }
  if (isTRUE(a$motif)) { r$motif <- mf_try_step("TF motif enrichment", mf_motif_enrichment(r$dmp$table, cfg, md$genome))
    if (!is.null(r$motif)) { mf_write_tsv(r$motif, file.path(od, "tf_motif_enrichment.tsv")); figs$motif <- motif_figure(r$motif, fd) } }
  if (isTRUE(a$clocks)) {
    cl <- mf_try_step("epigenetic clocks", mf_clocks(se, cfg))
    if (!is.null(cl) && ncol(cl$values) > 1 && !is.null(r$dmp$design$cd$age)) {
      r$clocks <- mf_clock_benchmark(cl, r$dmp$design)
      mf_write_tsv(r$clocks$values, file.path(od, "clock_estimates.tsv")); mf_write_tsv(r$clocks$benchmark, file.path(od, "clock_benchmark.tsv"))
      if (!is.null(r$clocks$group_tests)) mf_write_tsv(r$clocks$group_tests, file.path(od, "clock_group_tests.tsv"))
      figs$clock <- clock_figures(r$clocks, fd)
    } else if (!is.null(cl) && ncol(cl$values) > 1) mf_write_tsv(cl$values, file.path(od, "clock_estimates.tsv"))
  }
  if (length(a$compare_with)) {
    r$compare <- do.call(rbind, lapply(a$compare_with, function(x) {
      if (!file.exists(x$path)) { mf_log("compare_with file missing: ", x$path, level = "WARN"); return(NULL) }
      b <- as.data.frame(data.table::fread(x$path)); mf_compare_results(r$dmp$table, b, x$label %||% basename(x$path), a$fdr, a$enrichment_top_n) }))
    if (!is.null(r$compare)) mf_write_tsv(r$compare, file.path(od, "cross_dataset_concordance.tsv"))
  }
  # bundle for Shiny: results + beta of top CpGs
  topc <- unique(c(utils::head(r$dmp$table$cpg, 2000), if (!is.null(r$vmp)) utils::head(r$vmp$cpg, 500)))
  bundle <- list(project = cfg$project, platform = md$platform, dmp = r$dmp$table, dmr = r$dmr,
                 vmp = if (!is.null(r$vmp)) utils::head(r$vmp, 50000) else NULL, genes = r$genes, reg = r$reg, motif = r$motif,
                 enrich = r$enrich[setdiff(names(r$enrich), "basis")], clocks = r$clocks, compare = r$compare,
                 beta_top = r$dmp$beta[topc, , drop = FALSE], pheno = r$dmp$design$cd, formula = r$dmp$design$formula,
                 lambda = r$dmp$lambda, fdr = a$fdr)
  saveRDS(bundle, file.path(od, "results_bundle.rds"))
  if (isTRUE(a$report)) mf_try_step("report", analysis_report(r, cfg, figs, file.path(od, "analysis_report.html")))
  outs <- list.files(od, full.names = TRUE, pattern = "\\.(tsv|gz|rds|html)$")
  mf_write_manifest(file.path(od, "run_manifest.json"), cfg, inputs = dataset, outputs = outs,
                    extra = list(n_dmp = sum(r$dmp$table$FDR < a$fdr), lambda = r$dmp$lambda,
                                 n_vmp = if (!is.null(r$vmp)) sum(r$vmp$significant) else NA, n_dmr = if (!is.null(r$dmr)) nrow(r$dmr) else 0))
  mf_log("Analysis finished: ", od)
  invisible(r)
}

analysis_report <- function(r, cfg, figs, path) {
  a <- cfg$analysis; t <- r$dmp$table
  cols <- c("cpg", "chr", "pos", "gene", "gene_group", "island", "delta_beta_adj", "mean_beta_control", "mean_beta_case", "t", "P.Value", "FDR")
  secs <- list(
    list(heading = "Key findings", html = paste0("<div class='box'>", insight_text(r, cfg), "</div>",
         "<div class='warn'>Interpretation notes: associations are cross-sectional and may reflect medication, smoking, cell composition or other exposures that differ between groups; promoter-methylation-based transcriptional predictions are hypotheses to be checked against expression data.</div>")),
    list(heading = "Differentially methylated positions", plots = unname(c(figs$dmp[c("manhattan", "volcano", "qq")], figs$box, figs$pcassoc)),
         tables = list("Top 50 DMPs" = utils::head(t[, cols], 50))),
    list(heading = "Differentially methylated regions", tables = list("Top DMRs" = if (!is.null(r$dmr)) utils::head(r$dmr[, intersect(c("region", "no.cpgs", "min_smoothed_fdr", "Stouffer", "meandiff", "maxdiff", "direction", "overlapping.genes"), names(r$dmr))], 50) else NULL)),
    list(heading = "Variably methylated positions", plots = unname(figs$vmp[c("qq", "manhattan")]),
         tables = list("Top 50 VMPs" = if (!is.null(r$vmp)) utils::head(r$vmp, 50) else NULL)),
    list(heading = "Genes and predicted transcriptional effects", tables = list("Gene-level summary" = utils::head(r$genes, 60))),
    list(heading = "Regulatory context and transcription-factor motifs", plots = c(figs$reg, figs$motif),
         tables = list("Context enrichment" = r$reg, "Top TF motifs" = if (!is.null(r$motif)) utils::head(r$motif, 40) else NULL)),
    list(heading = "Pathways", plots = figs$enrich, html = if (!is.null(r$enrich)) paste0("<p>Input: ", r$enrich$basis, ". Tests account for the number of CpGs per gene (missMethyl).</p>") else "",
         tables = lapply(r$enrich[setdiff(names(r$enrich), "basis")], function(x) if (!is.null(x)) utils::head(x[, intersect(c("term", "N", "DE", "P.DE", "FDR"), names(x))], 25)))
  )
  if (!is.null(r$clocks)) secs[[length(secs) + 1]] <- list(heading = "Epigenetic clock benchmark", plots = unname(figs$clock),
    html = "<p>Age estimates are compared with chronological age (Pearson r, mean/median absolute error, bias). Acceleration is the residual from regressing each clock on age; group effects are adjusted for the analysis covariates. DunedinPACE is a rate (years of biological ageing per calendar year) and DNAmTL a telomere-length estimate.</p>",
    tables = list("Benchmark" = r$clocks$benchmark, "Case-control acceleration tests" = r$clocks$group_tests))
  if (!is.null(r$compare)) secs[[length(secs) + 1]] <- list(heading = "Cross-dataset / cross-tissue concordance", tables = list(Concordance = r$compare))
  secs[[length(secs) + 1]] <- list(heading = "Methods", html = paste0(
    "<p>M-values (log2 ratio of methylated to unmethylated signal) were modelled with limma (robust empirical Bayes) using <code>", htmltools::htmlEscape(r$dmp$design$formula), "</code>; ",
    "adjusted &Delta;&beta; come from the same model fitted on beta values. P-values were corrected with Benjamini-Hochberg. DMRs: DMRcate kernel smoothing (lambda = ", a$dmr_lambda, ", C = ", a$dmr_C, "). ",
    "VMPs: missMethyl::varFit (Levene-type test on absolute deviations). Enrichment: missMethyl gometh/gsameth. Context: Fisher tests vs all tested CpGs. Motifs: monaLisa binned enrichment with JASPAR2024 CORE vertebrate PWMs (optional). Clocks: methylclock and DunedinPACE.</p>"))
  mf_html_report(path, paste("methena analysis report:", cfg$project), secs)
}
