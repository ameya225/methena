select_sig <- function(tab, cfg, p_col = "P.Value", fdr_col = "FDR") {
  sig <- tab$cpg[tab[[fdr_col]] < cfg$analysis$fdr]
  basis <- sprintf("FDR < %s (n = %d)", cfg$analysis$fdr, length(sig))
  if (length(sig) < 50) {
    sig <- utils::head(tab$cpg[order(tab[[p_col]])], cfg$analysis$enrichment_top_n)
    basis <- sprintf("top %d CpGs by P (fewer than 50 at FDR < %s)", length(sig), cfg$analysis$fdr)
  }
  list(cpgs = sig, basis = basis)
}

array_type_missmethyl <- function(platform) switch(platform, "450k" = "450K", "epic" = "EPIC", "epicv2" = "EPIC_V2", NA)

#' Gene-set enrichment of DMPs corrected for the number of CpGs per gene
#'
#' Arrays: missMethyl::gometh (GO, KEGG) and gsameth (MSigDB collections from msigdbr).
#' Sequencing: limma::goana with the number of measured CpGs per gene as bias covariate.
#' @param tab DMP (or VMP) table with columns cpg, P.Value, FDR, gene
#' @param se dataset
#' @param cfg config
#' @param label name used in output
#' @export
mf_enrichment <- function(tab, se, cfg, label = "DMP") {
  plat <- S4Vectors::metadata(se)$platform
  s <- select_sig(tab, cfg)
  out <- list(basis = s$basis)
  mf_log(label, " enrichment input: ", s$basis)
  at <- array_type_missmethyl(plat)
  fmt <- function(x, coll) { if (is.null(x) || !nrow(x)) return(NULL)
    x$set_id <- rownames(x); x$collection <- coll; x$label <- label
    nm <- intersect(c("TERM", "Term", "Description", "ONTOLOGY"), names(x))
    x$term <- if (length(nm)) x[[nm[1]]] else x$set_id
    x <- x[order(x$P.DE), ]; utils::head(x, 500) }
  if (!is.na(at)) {
    need_pkg("missMethyl", "gene-set testing")
    invisible(mf_array_annotation(plat))   # attaches the annotation package that missMethyl looks up
    all <- tab$cpg
    out$GO <- fmt(tryCatch(suppressMessages(missMethyl::gometh(s$cpgs, all.cpg = all, collection = "GO", array.type = at, plot.bias = FALSE)),
                           error = function(e) { mf_log("gometh GO: ", conditionMessage(e), level = "WARN"); NULL }), "GO")
    out$KEGG <- fmt(tryCatch(suppressMessages(missMethyl::gometh(s$cpgs, all.cpg = all, collection = "KEGG", array.type = at, plot.bias = FALSE)),
                             error = function(e) { mf_log("KEGG unavailable (needs internet access to rest.kegg.jp): ", conditionMessage(e), level = "WARN"); NULL }), "KEGG")
    gmt <- gmt_sets(cfg$analysis$gene_set_files, cfg$config_dir)
    for (nm in names(gmt)) out[[nm]] <- fmt(tryCatch(suppressMessages(missMethyl::gsameth(s$cpgs, all.cpg = all, collection = gmt[[nm]], array.type = at)),
                                                     error = function(e) { mf_log("gsameth ", nm, ": ", conditionMessage(e), level = "WARN"); NULL }), nm)
    sets <- msigdb_sets(cfg$analysis$msigdb_collections)
    if (length(sets)) out$MSigDB <- fmt(tryCatch(suppressMessages(missMethyl::gsameth(s$cpgs, all.cpg = all, collection = sets, array.type = at)),
                                                 error = function(e) { mf_log("gsameth: ", conditionMessage(e), level = "WARN"); NULL }), "MSigDB")
  } else {
    ann <- as.data.frame(SummarizedExperiment::rowData(se))
    eg <- gene_entrez(ann$gene)
    names(eg) <- ann$probe
    univ <- eg[!is.na(eg)]
    covar <- table(univ)
    de <- unique(univ[intersect(names(univ), s$cpgs)])
    if (length(de) >= 5) {
      g <- tryCatch(limma::goana(de, universe = names(covar), covariate = as.numeric(covar), species = "Hs"), error = function(e) NULL)
      if (!is.null(g)) { g$P.DE <- g$P.DE; g$FDR <- stats::p.adjust(g$P.DE, "BH"); out$GO <- fmt(g, "GO") }
    }
  }
  out
}

gene_entrez <- function(sym) {
  s1 <- sub(";.*$", "", ifelse(is.na(sym), "", sym))
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) return(rep(NA_character_, length(s1)))
  m <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, unique(s1[s1 != ""]), "ENTREZID", "SYMBOL"))
  unname(m[s1])
}

#' Read gene sets from GMT files (gene symbols or Entrez IDs) as Entrez-ID lists
#' @param files named list/vector of GMT paths (names become collection labels)
#' @param base folder for relative paths
#' @export
gmt_sets <- function(files, base = ".") {
  if (!length(files)) return(list())
  if (is.null(names(files))) names(files) <- sub("\\.gmt$", "", basename(unlist(files)))
  out <- list()
  for (nm in names(files)) {
    f <- files[[nm]]; if (!grepl("^(/|~)", f)) f <- file.path(base, f)
    if (!file.exists(f)) { mf_log("Gene-set file not found: ", f, level = "WARN"); next }
    ln <- strsplit(readLines(f, warn = FALSE), "\t")
    sets <- lapply(ln, function(x) x[-(1:2)]); names(sets) <- vapply(ln, `[`, "", 1)
    allg <- unique(unlist(sets))
    if (mean(grepl("^[0-9]+$", allg)) < 0.5) {   # symbols -> Entrez
      need_pkg("org.Hs.eg.db", "converting gene symbols")
      map <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, allg, "ENTREZID", "SYMBOL"))
      sets <- lapply(sets, function(g) unique(stats::na.omit(unname(map[g]))))
    }
    sets <- sets[lengths(sets) >= 5 & lengths(sets) <= 1000]
    mf_log("Loaded ", length(sets), " gene sets from ", basename(f))
    out[[nm]] <- sets
  }
  out
}

msigdb_sets <- function(colls) {
  if (!length(colls) || !requireNamespace("msigdbr", quietly = TRUE)) return(list())
  sets <- list()
  for (cc in colls) {
    parts <- strsplit(cc, ":", fixed = TRUE)[[1]]
    sub <- if (length(parts) > 1) paste(parts[-1], collapse = ":") else NULL
    d <- tryCatch(suppressMessages(msigdbr::msigdbr(species = "Homo sapiens", collection = parts[1], subcollection = sub)),
                  error = function(e) tryCatch(suppressMessages(msigdbr::msigdbr(species = "Homo sapiens", category = parts[1], subcategory = sub)),
                                               error = function(e2) { mf_log("msigdbr ", cc, ": ", conditionMessage(e2), level = "WARN"); NULL }))
    if (is.null(d) || !nrow(d)) next
    idc <- intersect(c("ncbi_gene", "entrez_gene"), names(d))[1]
    sets <- c(sets, split(as.character(d[[idc]]), d$gs_name))
  }
  sets
}

#' Enrichment of DMPs in genomic / regulatory context (Fisher tests, hyper and hypo separately)
#' @param tab DMP table
#' @param cfg config
#' @export
mf_regulatory_enrichment <- function(tab, cfg) {
  s <- select_sig(tab, cfg)
  tab$sel <- tab$cpg %in% s$cpgs
  feat <- list()
  # CpG-island context is only known for array CpGs (sequencing CpGs inherit it only when they
  # overlap an array probe); skip it when most CpGs are unannotated
  if (mean(!is.na(tab$island)) >= 0.5) {
    isl <- ifelse(is.na(tab$island) | tab$island == "", "OpenSea", tab$island)
    for (lv in unique(isl)) feat[[paste0("CpG context: ", lv)]] <- isl == lv
  }
  g1 <- sub(";.*$", "", ifelse(is.na(tab$gene_group), "", tab$gene_group))
  prom <- g1 %in% c("TSS1500", "TSS200", "5'UTR", "1stExon")
  feat[["Gene region: promoter (TSS1500/TSS200/5'UTR/1stExon)"]] <- prom
  feat[["Gene region: gene body"]] <- g1 %in% c("Body", "ExonBnd")
  feat[["Gene region: 3'UTR"]] <- g1 == "3'UTR"
  feat[["Gene region: intergenic"]] <- g1 == ""
  if (!is.null(tab$enhancer) && any(tab$enhancer %in% TRUE)) feat[["Enhancer (array annotation)"]] <- tab$enhancer %in% TRUE
  res <- list()
  for (dirn in c("hyper", "hypo")) {
    sel <- tab$sel & tab$direction == dirn
    if (sum(sel) < 5) next
    for (f in names(feat)) {
      x <- feat[[f]]
      if (sum(x) < 10) next
      ft <- stats::fisher.test(matrix(c(sum(sel & x), sum(sel & !x), sum(!sel & x), sum(!sel & !x)), 2))
      res[[length(res) + 1]] <- data.frame(direction = dirn, feature = f, n_sig_in_feature = sum(sel & x), n_sig = sum(sel),
                                           frac_sig = mean(x[sel]), frac_background = mean(x[!sel]),
                                           odds_ratio = unname(ft$estimate), ci_low = ft$conf.int[1], ci_high = ft$conf.int[2], p = ft$p.value)
    }
  }
  if (!length(res)) return(NULL)
  out <- do.call(rbind, res); out$fdr <- stats::p.adjust(out$p, "BH"); out$basis <- s$basis
  out
}

#' Transcription-factor motif enrichment around DMPs (monaLisa + JASPAR2024)
#'
#' Sequences of +/- flank bp around hyper- and hypomethylated DMPs are compared with
#' GC-matched background CpGs from the same array/assay.
#' @param tab DMP table
#' @param cfg config
#' @param genome hg19 or hg38
#' @export
mf_motif_enrichment <- function(tab, cfg, genome = "hg19") {
  pk <- paste0("BSgenome.Hsapiens.UCSC.", genome)
  for (p in c("monaLisa", "TFBSTools", pk, "RSQLite")) if (!requireNamespace(p, quietly = TRUE)) {
    mf_log("Motif analysis skipped: package ", p, " missing", level = "WARN"); return(NULL) }
  set.seed(cfg$seed)
  tab <- tab[!is.na(tab$chr) & !is.na(tab$pos), ]
  top <- utils::head(tab[order(tab$P.Value), ], cfg$analysis$motif_top_n)
  if (sum(top$FDR < cfg$analysis$fdr) >= 50) top <- top[top$FDR < cfg$analysis$fdr, ]
  bg <- tab[!tab$cpg %in% top$cpg & tab$P.Value > 0.5, ]
  bg <- bg[sample(nrow(bg), min(nrow(bg), 3 * nrow(top))), ]
  d <- rbind(data.frame(top[, c("chr", "pos")], bin = ifelse(top$t > 0, "hyper", "hypo")), data.frame(bg[, c("chr", "pos")], bin = "background"))
  g <- get(pk, envir = asNamespace(pk))
  gr <- GenomicRanges::GRanges(d$chr, IRanges::IRanges(d$pos - cfg$analysis$motif_flank, d$pos + cfg$analysis$motif_flank))
  seqs <- Biostrings::getSeq(g, gr)
  names(seqs) <- paste0("s", seq_along(seqs))
  bins <- factor(d$bin, levels = c("background", "hypo", "hyper"))
  opts <- list(tax_group = "vertebrates", collection = "CORE", matrixtype = "PWM")
  # JASPAR2024 downloads its database on first use; fall back to the bundled JASPAR2020 offline
  pwms <- tryCatch({
    jdb <- suppressMessages(JASPAR2024::JASPAR2024())
    con <- RSQLite::dbConnect(RSQLite::SQLite(), JASPAR2024::db(jdb)); on.exit(RSQLite::dbDisconnect(con), add = TRUE)
    attr_src <- "JASPAR2024"; TFBSTools::getMatrixSet(con, opts) },
    error = function(e) {
      if (!requireNamespace("JASPAR2020", quietly = TRUE)) mf_stop("JASPAR2024 database unavailable and JASPAR2020 not installed")
      mf_log("JASPAR2024 database not reachable; using bundled JASPAR2020 motifs", level = "WARN")
      TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020, opts) })
  se <- tryCatch(monaLisa::calcBinnedMotifEnrR(seqs = seqs, bins = bins, pwmL = pwms, background = "otherBins",
                                               BPPARAM = BiocParallel::SerialParam(), verbose = FALSE),
                 error = function(e) { mf_log("monaLisa failed: ", conditionMessage(e), level = "WARN"); NULL })
  if (is.null(se)) return(NULL)
  rd <- as.data.frame(SummarizedExperiment::rowData(se))
  le <- SummarizedExperiment::assay(se, "log2enr"); pa <- SummarizedExperiment::assay(se, "negLog10Padj")
  out <- do.call(rbind, lapply(c("hyper", "hypo"), function(b) if (b %in% colnames(le))
    data.frame(bin = b, motif_id = rd$motif.id, tf = rd$motif.name, log2_enrichment = le[, b], neg_log10_padj = pa[, b])))
  out <- out[order(-out$neg_log10_padj), ]
  attr(out, "n") <- table(bins)
  out
}

#' Compare two DMP result tables (e.g. blood vs brain, discovery vs replication)
#' @param a,b DMP tables (columns cpg, t, P.Value, FDR)
#' @param label name of comparison
#' @param fdr threshold to select CpGs in a
#' @param top_n fallback if few significant
#' @export
mf_compare_results <- function(a, b, label, fdr = 0.05, top_n = 1000) {
  m <- merge(a[, c("cpg", "t", "P.Value", "FDR")], b[, c("cpg", "t", "P.Value")], by = "cpg", suffixes = c("_a", "_b"))
  if (!nrow(m)) return(NULL)
  selA <- m$FDR < fdr; basis <- "FDR"
  if (sum(selA) < 20) { selA <- rank(m$P.Value_a) <= top_n; basis <- paste0("top", top_n) }
  s <- m[selA, ]
  conc <- sum(sign(s$t_a) == sign(s$t_b))
  bt <- stats::binom.test(conc, nrow(s), 0.5)
  data.frame(comparison = label, selection = basis, n_cpgs_shared = nrow(m), n_selected = nrow(s),
             sign_concordance = conc / nrow(s), binom_p = bt$p.value,
             spearman_t_selected = suppressWarnings(stats::cor(s$t_a, s$t_b, method = "spearman")),
             spearman_t_all = suppressWarnings(stats::cor(m$t_a, m$t_b, method = "spearman")),
             n_replicated_p05_same_dir = sum(s$P.Value_b < 0.05 & sign(s$t_a) == sign(s$t_b)),
             expected_by_chance = round(0.025 * nrow(s), 1))
}

#' Integrate methylation with external differential expression results
#' @param gene_tab gene-level methylation table from [mf_gene_level]
#' @param expr data.frame with columns gene and logFC (or t)
#' @export
mf_integrate_expression <- function(gene_tab, expr) {
  names(expr) <- tolower(names(expr))
  stat <- intersect(c("logfc", "t", "stat"), names(expr))[1]
  if (is.na(stat) || !"gene" %in% names(expr)) mf_stop("Expression table needs columns 'gene' and 'logFC' (or 't')")
  m <- merge(gene_tab, expr[, c("gene", stat)], by = "gene")
  m$expr_stat <- m[[stat]]
  pr <- m[!is.na(m$promoter_mean_delta_beta), ]
  ct <- if (nrow(pr) >= 5) stats::cor.test(pr$promoter_mean_delta_beta, pr$expr_stat, method = "spearman") else NULL
  list(table = m, n_promoter_genes = nrow(pr),
       promoter_rho = if (!is.null(ct)) unname(ct$estimate) else NA, promoter_p = if (!is.null(ct)) ct$p.value else NA,
       frac_inverse = if (nrow(pr)) mean(sign(pr$promoter_mean_delta_beta) != sign(pr$expr_stat)) else NA)
}
