subset_raw <- function(raw, keep_ids) {
  raw$ss <- raw$ss[raw$ss$sample_id %in% keep_ids, , drop = FALSE]
  k <- raw$ss$sample_id
  if (!is.null(raw$rgset)) raw$rgset <- raw$rgset[, k]
  if (!is.null(raw$mset)) raw$mset <- raw$mset[, k]
  if (!is.null(raw$beta)) raw$beta <- raw$beta[, k, drop = FALSE]
  if (!is.null(raw$detp)) raw$detp <- raw$detp[, k, drop = FALSE]
  if (!is.null(raw$snp_probes)) raw$snp_probes <- raw$snp_probes[, k, drop = FALSE]
  if (!is.null(raw$M)) { raw$M <- raw$M[, k, drop = FALSE]; raw$Cov <- raw$Cov[, k, drop = FALSE] }
  raw
}

#' Published cross-reactive probe IDs bundled with methylflow (Chen 2013; McCartney 2016; Pidsley 2016)
#' @export
mf_cross_reactive_probes <- function() {
  f <- system.file("extdata", "cross_reactive_probes.tsv.gz", package = "methylflow")
  data.table::fread(f, data.table = FALSE)$probe
}

snp_ch_filter <- function(probes, array, cfg) {
  q <- cfg$qc; drop <- list()
  ann <- mf_array_annotation(array)[probes, , drop = FALSE]
  if (isTRUE(q$remove_snp_probes)) {
    # SNPs (dbSNP common, minfi annotation) at the CpG site or single-base extension with MAF >= snp_maf
    snp <- (!is.na(ann$cpg_maf) & ann$cpg_maf >= q$snp_maf) | (!is.na(ann$sbe_maf) & ann$sbe_maf >= q$snp_maf)
    drop$snp <- probes[snp]
  }
  if (isTRUE(q$remove_cross_reactive)) {
    xr <- mf_cross_reactive_probes()
    drop$cross_reactive <- setdiff(probes[sub("_.*$", "", probes) %in% xr], unlist(drop))
  }
  if (isTRUE(q$remove_xy)) drop$sex_chromosome <- setdiff(probes[ann$chr %in% c("chrX", "chrY")], unlist(drop))
  drop
}

#' Collapse the two strands of each CpG (Bismark coverage) using a BSgenome
#' @keywords internal
collapse_strands <- function(raw, genome) {
  pkg <- paste0("BSgenome.Hsapiens.UCSC.", genome)
  if (!requireNamespace(pkg, quietly = TRUE)) { mf_log(pkg, " not installed: strands not collapsed", level = "WARN"); return(raw) }
  g <- get(pkg, envir = asNamespace(pkg))
  okchr <- raw$chr %in% GenomeInfoDb::seqnames(g)
  raw <- list(type = raw$type, ss = raw$ss, M = raw$M[okchr, , drop = FALSE], Cov = raw$Cov[okchr, , drop = FALSE],
              chr = raw$chr[okchr], pos = raw$pos[okchr])
  base <- as.character(Biostrings::getSeq(g, GenomicRanges::GRanges(raw$chr, IRanges::IRanges(raw$pos, width = 1))))
  newpos <- ifelse(base == "G", raw$pos - 1L, raw$pos)
  key <- paste(raw$chr, newpos, sep = ":")
  f <- factor(key, levels = unique(key))
  M <- rowsum(raw$M, f, reorder = FALSE); Cov <- rowsum(raw$Cov, f, reorder = FALSE)
  first <- !duplicated(f)
  mf_log("Strand collapse: ", length(key), " stranded sites -> ", nrow(M), " CpGs")
  list(type = raw$type, ss = raw$ss, M = M, Cov = Cov, chr = raw$chr[first], pos = newpos[first])
}

#' Annotate genomic CpGs (BS-seq) with nearest gene / promoter context and array probe IDs
#' @param chr chromosome vector
#' @param pos position vector (1-based C of the CpG on + strand)
#' @param genome hg19 or hg38
#' @export
mf_annotate_positions <- function(chr, pos, genome = "hg19") {
  out <- data.frame(probe = paste(chr, pos, sep = ":"), chr = chr, pos = pos, stringsAsFactors = FALSE)
  out$gene <- ""; out$gene_group <- ""; out$island <- NA_character_; out$enhancer <- FALSE; out$dhs <- FALSE
  out$array_probe <- NA_character_
  txpkg <- paste0("TxDb.Hsapiens.UCSC.", genome, ".knownGene")
  gr <- GenomicRanges::GRanges(chr, IRanges::IRanges(pos, width = 1))
  if (requireNamespace(txpkg, quietly = TRUE) && requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    tx <- get(txpkg, envir = asNamespace(txpkg))
    genes <- suppressMessages(GenomicFeatures::genes(tx, single.strand.genes.only = TRUE))
    sym <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, names(genes), "SYMBOL", "ENTREZID"))
    genes$symbol <- unname(sym)
    tss <- GenomicRanges::resize(genes, 1, fix = "start")
    nt <- GenomicRanges::distanceToNearest(gr, tss, ignore.strand = TRUE)
    d <- S4Vectors::mcols(nt)$distance
    sh <- S4Vectors::subjectHits(nt); qh <- S4Vectors::queryHits(nt)
    near <- genes$symbol[sh]
    grp <- ifelse(d <= 200, "TSS200", ifelse(d <= 1500, "TSS1500", ""))
    inbody <- IRanges::overlapsAny(gr, genes, ignore.strand = TRUE)
    out$gene[qh] <- ifelse(d <= 1500 | inbody[qh], ifelse(is.na(near), "", near), "")
    out$gene_group[qh] <- ifelse(grp != "", grp, ifelse(inbody[qh], "Body", ""))
    out$dist_to_tss <- NA_integer_; out$dist_to_tss[qh] <- d
  } else mf_log("TxDb/org.Hs.eg.db not installed: genes not annotated", level = "WARN")
  # map to array probe IDs for clocks / deconvolution
  arrays <- if (genome == "hg19") c("450k", "epic") else "epicv2"
  for (a in arrays) {
    ann <- tryCatch(mf_array_annotation(a), error = function(e) NULL)
    if (is.null(ann)) next
    key <- paste(ann$chr, ann$pos, sep = ":")
    hit <- match(out$probe, key)
    fill <- is.na(out$array_probe) & !is.na(hit)
    out$array_probe[fill] <- sub("_.*$", "", ann$probe[hit[fill]])
    isl <- ann$island[hit]; fill2 <- is.na(out$island) & !is.na(isl); out$island[fill2] <- isl[fill2]
  }
  rownames(out) <- out$probe
  out
}

#' Run the complete preprocessing workflow and write the standardised dataset
#'
#' Writes to `<output_dir>/preprocess/`: methylflow_dataset.rds (SummarizedExperiment),
#' sample_qc.tsv, probe_filter_summary.tsv, qc_report.html, figures/, run_manifest.json, log.
#' @param cfg config list (see [mf_read_config])
#' @return the SummarizedExperiment (invisibly)
#' @export
mf_preprocess <- function(cfg) {
  mf_reset(); set.seed(cfg$seed)
  od <- file.path(cfg$output_dir, "preprocess"); dir.create(file.path(od, "figures"), recursive = TRUE, showWarnings = FALSE)
  mf_log_to(file.path(od, "preprocess.log"))
  mf_log("methylflow ", as.character(utils::packageVersion("methylflow")), " preprocessing: ", cfg$project)
  ss <- mf_step("read sample sheet", mf_read_samples(cfg))
  mf_log(nrow(ss), " samples in sample sheet")
  is_array <- cfg$input$type != "bismark"
  filt <- list(); removed_beta <- NULL
  if (is_array) {
    raw <- mf_step("read array data", mf_read_array(cfg, ss))
    mf_log("Platform: ", raw$array, "; samples with data: ", nrow(raw$ss))
    qc <- mf_step("sample QC", mf_qc_array(raw, cfg))
    sq <- qc$sample
    reason <- character(nrow(sq))
    add <- function(flag, txt) { flag[is.na(flag)] <- FALSE; reason[flag] <<- paste0(reason[flag], ifelse(reason[flag] == "", "", ";"), txt) }
    add(sq$detection_fail, "detection_p")
    if (isTRUE(cfg$qc$sex_check) && isTRUE(cfg$qc$drop_sex_mismatch)) add(sq$sex_mismatch, "sex_mismatch")
    keep <- sq$sample_id[reason == ""]
    raw <- subset_raw(raw, keep)
    if (length(keep) < 3) mf_stop("Fewer than 3 samples passed QC")
    qc$mset <- NULL
    # keep only the pass/fail state of detection p-values (halves memory for large studies)
    if (!is.null(raw$detp)) { f <- raw$detp > cfg$qc$detp_threshold; f[is.na(f)] <- TRUE; raw$detp <- f; rm(f) }
    gc(verbose = FALSE)
    beta <- mf_step("normalisation", mf_normalize_array(raw, cfg))
    raw$mset <- NULL; raw$rgset <- if (cfg$input$type == "idat") raw$rgset else NULL; gc(verbose = FALSE)
    n0 <- nrow(beta); filt$start <- n0
    # probe filters
    if (!is.null(raw$detp)) {
      D <- raw$detp[intersect(rownames(beta), rownames(raw$detp)), colnames(beta), drop = FALSE]
      fail <- D
      bad <- rownames(D)[rowMeans(fail) > cfg$qc$probe_max_failed_frac]
      filt$detection_p <- length(bad)
      # mask individual failed measurements
      fm <- fail[setdiff(rownames(D), bad), , drop = FALSE]
      nmask <- sum(fm)
      idx <- which(fm, arr.ind = TRUE)
      if (nrow(idx)) beta[cbind(match(rownames(fm)[idx[, 1]], rownames(beta)), idx[, 2])] <- NA
      mf_log("Masked ", nmask, " individual failed measurements (detection p > ", cfg$qc$detp_threshold, ")")
      beta <- beta[setdiff(rownames(beta), bad), , drop = FALSE]
      filt$masked_values <- nmask
    }
    if (!is.null(qc$bead_fail)) {
      badb <- intersect(rownames(beta), names(qc$bead_fail)[qc$bead_fail])
      filt$low_beadcount <- length(badb); beta <- beta[setdiff(rownames(beta), badb), , drop = FALSE]
    }
    cd <- merge(raw$ss, sq, by = "sample_id", all.x = TRUE, suffixes = c("", "_qc")); rownames(cd) <- cd$sample_id
    cd <- cd[colnames(beta), , drop = FALSE]
    cd <- mf_step("estimate covariates", mf_add_covariates(beta, cd, cfg, raw$type))
    drops <- snp_ch_filter(rownames(beta), raw$array, cfg)
    for (n in names(drops)) filt[[n]] <- length(drops[[n]])
    rm_ids <- unique(unlist(drops))
    removed_beta <- beta[intersect(rm_ids, rownames(beta)), , drop = FALSE]
    beta <- beta[setdiff(rownames(beta), rm_ids), , drop = FALSE]
    ann <- mf_array_annotation(raw$array)
    ann <- ann[rownames(beta), , drop = FALSE]
    ok <- !is.na(ann$chr); beta <- beta[ok, , drop = FALSE]; ann <- ann[ok, , drop = FALSE]
    filt$no_annotation <- sum(!ok)
    platform <- raw$array; cov_assay <- NULL
  } else {
    raw <- mf_step("read bisulfite coverage", mf_read_bismark(cfg, ss))
    if (isTRUE(cfg$qc$bs_collapse_strands %||% TRUE)) raw <- mf_step("strand collapse", collapse_strands(raw, cfg$input$genome))
    qc <- mf_step("sample QC", mf_qc_bismark(raw, cfg))
    sq <- qc$sample
    reason <- ifelse(sq$low_cpg_count, "low_cpg_count", "")
    if (isTRUE(cfg$qc$sex_check) && isTRUE(cfg$qc$drop_sex_mismatch)) reason[sq$sex_mismatch %in% TRUE] <- paste0(reason[sq$sex_mismatch %in% TRUE], "sex_mismatch")
    keep <- sq$sample_id[reason == ""]
    raw <- subset_raw(raw, keep)
    if (length(keep) < 3) mf_stop("Fewer than 3 samples passed QC")
    Cov <- raw$Cov; Mc <- raw$M
    filt$start <- nrow(Cov)
    capq <- apply(Cov, 2, function(x) stats::quantile(x[x > 0], cfg$qc$bs_max_coverage_quantile))
    over <- sweep(Cov, 2, capq, ">")
    Cov[over] <- 0L; Mc[over] <- 0L
    okcov <- Cov >= cfg$qc$bs_min_coverage
    keep_rows <- rowMeans(okcov) >= cfg$qc$bs_min_sample_frac
    filt$low_coverage <- sum(!keep_rows)
    Cov <- Cov[keep_rows, , drop = FALSE]; Mc <- Mc[keep_rows, , drop = FALSE]; okcov <- okcov[keep_rows, , drop = FALSE]
    beta <- Mc / pmax(Cov, 1); beta[!okcov] <- NA
    chr <- raw$chr[keep_rows]; pos <- raw$pos[keep_rows]
    if (isTRUE(cfg$qc$remove_xy)) { xy <- chr %in% c("chrX", "chrY"); filt$sex_chromosome <- sum(xy)
      beta <- beta[!xy, , drop = FALSE]; Cov <- Cov[!xy, , drop = FALSE]; chr <- chr[!xy]; pos <- pos[!xy] }
    ann <- mf_step("annotate CpGs", mf_annotate_positions(chr, pos, cfg$input$genome))
    rownames(beta) <- rownames(Cov) <- ann$probe
    cd <- merge(raw$ss, sq, by = "sample_id", all.x = TRUE, suffixes = c("", "_qc")); rownames(cd) <- cd$sample_id
    cd <- cd[colnames(beta), , drop = FALSE]
    # covariates through array-probe mapping (works when coverage overlaps array CpGs)
    mapped <- !is.na(ann$array_probe)
    if (sum(mapped) > 1000) {
      bm <- beta[mapped, , drop = FALSE]; rownames(bm) <- ann$array_probe[mapped]
      bm <- bm[!duplicated(rownames(bm)), , drop = FALSE]
      cd <- tryCatch(mf_add_covariates(bm, cd, cfg, "mapped"), error = function(e) { mf_log("Covariate estimation failed: ", conditionMessage(e), level = "WARN"); cd })
    }
    platform <- paste0("bisulfite_seq_", cfg$input$genome); cov_assay <- Cov
  }
  grp_all <- if (!is.null(cd$group)) stats::setNames(as.character(cd$group), rownames(cd)) else NULL
  # outliers on autosomal normalised data
  out <- mf_step("outlier detection", mf_flag_outliers(beta, cfg$qc$outlier_sd))
  cd$pca_outlier <- out$outlier
  if (!is.null(out$cor_to_median)) cd$cor_to_median_profile <- out$cor_to_median
  if (!is.null(out$pca)) { pcs <- out$pca$pcs[, 1:min(5, ncol(out$pca$pcs)), drop = FALSE]; colnames(pcs) <- paste0("qc_", colnames(pcs)); cd <- cbind(cd, pcs[rownames(cd), , drop = FALSE]) }
  if (isTRUE(cfg$qc$drop_outliers) && any(out$outlier)) {
    # safeguard: automatic removal of many samples (or a large share of one group) can bias results
    frac <- mean(out$outlier)
    gshare <- if (!is.null(cd$group)) max(tapply(out$outlier, cd$group, mean), na.rm = TRUE) else frac
    if ((frac > (cfg$qc$max_outlier_frac %||% 0.1) || gshare > 0.25) && !isTRUE(cfg$qc$force_drop_outliers)) {
      mf_log(sum(out$outlier), " samples (", round(100 * frac), "% overall, up to ", round(100 * gshare),
             "% of one group) are PCA/profile outliers - too many to remove automatically. They are KEPT and flagged ",
             "(column pca_outlier); inspect the QC report and set qc$force_drop_outliers: true to remove them.", level = "WARN")
      out$outlier[] <- FALSE
    }
  }
  if (isTRUE(cfg$qc$drop_outliers) && any(out$outlier)) {
    mf_log("Dropping ", sum(out$outlier), " PCA/profile outlier samples: ", paste(colnames(beta)[out$outlier], collapse = ", "))
    sq$exclusion_reason <- reason
    sq$exclusion_reason[sq$sample_id %in% colnames(beta)[out$outlier]] <- "pca_outlier"
    keepc <- !out$outlier
    beta <- beta[, keepc, drop = FALSE]; cd <- cd[keepc, , drop = FALSE]
    if (!is.null(cov_assay)) cov_assay <- cov_assay[, keepc, drop = FALSE]
    if (!is.null(removed_beta)) removed_beta <- removed_beta[, keepc, drop = FALSE]
  } else sq$exclusion_reason <- reason
  sq$included <- sq$sample_id %in% colnames(beta)
  filt$final <- nrow(beta)
  mf_log("Final dataset: ", nrow(beta), " CpGs x ", ncol(beta), " samples")
  assays <- list(beta = beta); if (!is.null(cov_assay)) assays$coverage <- cov_assay
  se <- SummarizedExperiment::SummarizedExperiment(assays = assays, colData = S4Vectors::DataFrame(cd, check.names = FALSE),
                                                   rowData = S4Vectors::DataFrame(ann[rownames(beta), , drop = FALSE]))
  S4Vectors::metadata(se) <- list(methylflow_version = as.character(utils::packageVersion("methylflow")),
                                  input_type = cfg$input$type, platform = platform, genome = if (is_array) "hg19" else cfg$input$genome,
                                  tissue = cfg$input$tissue, probe_filters = filt, config = cfg,
                                  beta_filtered_probes = removed_beta, created = format(Sys.time()))
  if (is_array && platform == "epicv2") S4Vectors::metadata(se)$genome <- "hg38"
  # outputs
  ds <- file.path(od, "methylflow_dataset.rds"); saveRDS(se, ds)
  mf_write_tsv(sq, file.path(od, "sample_qc.tsv"))
  fsum <- data.frame(step = names(filt), n = unlist(filt)); mf_write_tsv(fsum, file.path(od, "probe_filter_summary.tsv"))
  mf_write_tsv(data.frame(sample_id = colnames(beta), as.data.frame(cd)[, setdiff(names(cd), "sample_id"), drop = FALSE], check.names = FALSE),
               file.path(od, "sample_annotation_final.tsv"))
  figs <- mf_step("QC figures", qc_figures(beta, cd, sq, od, out$pca, grp_all = grp_all))
  qc_report(cfg, se, sq, fsum, figs, file.path(od, "qc_report.html"))
  outs <- file.path(od, c("methylflow_dataset.rds", "sample_qc.tsv", "probe_filter_summary.tsv", "sample_annotation_final.tsv", "qc_report.html"))
  ins <- unlist(cfg$input[c("sample_sheet", "geo_series_matrix", "meth_file", "unmeth_file", "detp_file", "beta_file")])
  mf_write_manifest(file.path(od, "run_manifest.json"), cfg, inputs = ins, outputs = outs,
                    extra = list(n_samples = ncol(beta), n_cpgs = nrow(beta), probe_filters = filt))
  mf_log("Preprocessing finished. Dataset: ", ds)
  invisible(se)
}

qc_figures <- function(beta, cd, sq, od, pca, grp_all = NULL) {
  fd <- file.path(od, "figures"); f <- list()
  grp <- if (!is.null(cd$group)) as.character(cd$group) else rep("all", ncol(beta))
  set.seed(1); idx <- sample(nrow(beta), min(20000, nrow(beta)))
  dd <- data.frame(beta = as.vector(beta[idx, ]), sample = rep(colnames(beta), each = length(idx)),
                   group = rep(grp, each = length(idx)))
  dd <- dd[!is.na(dd$beta), ]
  f$density <- mf_save_plot(ggplot2::ggplot(dd, ggplot2::aes(beta, group = sample, colour = group)) +
    ggplot2::geom_density(linewidth = 0.2, alpha = 0.4) + ggplot2::theme_bw() +
    ggplot2::labs(title = "Beta-value density per sample (normalised)", x = "Beta", y = "Density", colour = "Group"),
    file.path(fd, "beta_density.png"), 7, 4)
  if (!is.null(sq$frac_failed_probes) && any(!is.na(sq$frac_failed_probes))) {
    s2 <- sq[order(sq$frac_failed_probes), ]; s2$rank <- seq_len(nrow(s2))
    f$detp <- mf_save_plot(ggplot2::ggplot(s2, ggplot2::aes(rank, frac_failed_probes, colour = included)) +
      ggplot2::geom_point(size = 1) + ggplot2::theme_bw() + ggplot2::scale_y_continuous(labels = scales::percent) +
      ggplot2::labs(title = "Fraction of probes failing detection p-value per sample", x = "Samples (ranked)", y = "Failed probes", colour = "Included"),
      file.path(fd, "detection_fail.png"), 7, 3.5)
  }
  if (!is.null(sq$xMed)) {
    f$sex <- mf_save_plot(ggplot2::ggplot(sq, ggplot2::aes(xMed, yMed, colour = predicted_sex, shape = sex_mismatch)) +
      ggplot2::geom_point() + ggplot2::theme_bw() + ggplot2::labs(title = "Sex check: X vs Y chromosome median intensity",
      colour = "Predicted sex", shape = "Mismatch with reported"), file.path(fd, "sex_check.png"), 6, 4)
  }
  if (!is.null(pca)) {
    g2 <- if (!is.null(grp_all)) unname(grp_all[rownames(pca$pcs)]) else rep("all", nrow(pca$pcs))
    g2[is.na(g2) | rownames(pca$pcs) %in% sq$sample_id[sq$exclusion_reason == "pca_outlier"]] <- "excluded (outlier)"
    p <- data.frame(pca$pcs[, 1:2], group = g2)
    ve <- round(100 * pca$var_explained[1:2], 1)
    f$pca <- mf_save_plot(ggplot2::ggplot(p, ggplot2::aes(PC1, PC2, colour = group)) + ggplot2::geom_point(size = 1.4, alpha = 0.8) +
      ggplot2::theme_bw() + ggplot2::labs(title = "PCA of the 20,000 most variable CpGs", x = paste0("PC1 (", ve[1], "%)"),
      y = paste0("PC2 (", ve[2], "%)"), colour = "Group"), file.path(fd, "pca.png"), 6, 4.5)
  }
  cc <- grep("^cell_", names(cd), value = TRUE)
  if (length(cc)) {
    m <- data.frame(group = rep(grp, length(cc)), cell = rep(sub("^cell_", "", cc), each = nrow(cd)), prop = unlist(as.data.frame(cd)[, cc]))
    f$cells <- mf_save_plot(ggplot2::ggplot(m, ggplot2::aes(cell, prop, fill = group)) + ggplot2::geom_boxplot(outlier.size = 0.4) +
      ggplot2::theme_bw() + ggplot2::labs(title = "Estimated cell-type proportions", x = NULL, y = "Proportion", fill = "Group"),
      file.path(fd, "cell_proportions.png"), 8, 4)
  }
  f
}

qc_report <- function(cfg, se, sq, fsum, figs, path) {
  md <- S4Vectors::metadata(se)
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  ex <- sq[!sq$included, c("sample_id", "exclusion_reason"), drop = FALSE]
  grp_tab <- if (!is.null(cd$group)) as.data.frame(table(group = cd$group, useNA = "ifany")) else NULL
  intro <- sprintf("<div class='box'><b>Input:</b> %s (%s). <b>Samples:</b> %d in sample sheet &rarr; %d retained. <b>CpGs:</b> %s retained.</div>",
                   md$input_type, md$platform, nrow(sq), ncol(se), format(nrow(se), big.mark = ","))
  secs <- list(
    list(heading = "Summary", html = intro, tables = list("Samples per group" = grp_tab, "Probe / CpG filtering" = fsum)),
    list(heading = "Sample QC", html = "<p>Samples are excluded for high detection-p failure, sex mismatch (reported vs. methylation-predicted) and PCA/profile outliers, according to the configured thresholds.</p>",
         plots = figs[intersect(c("detp", "sex"), names(figs))], tables = list("Excluded samples" = ex)),
    list(heading = "Normalised data", plots = figs[intersect(c("density", "pca"), names(figs))]),
    list(heading = "Estimated covariates", plots = figs[intersect("cells", names(figs))],
         html = "<p>Cell proportions (blood: EpiDISH 7-cell reference; brain: HiBED) and a methylation smoking score are added to the sample table for adjustment.</p>"))
  mf_html_report(path, paste("methylflow QC report:", cfg$project), secs)
}
