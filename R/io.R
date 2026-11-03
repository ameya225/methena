# ---- sample sheets ---------------------------------------------------------

clean_names <- function(x) {
  x <- tolower(gsub("[^A-Za-z0-9]+", "_", x))
  x <- gsub("^_+|_+$", "", x)
  make.unique(ifelse(x == "", "col", x), sep = "_")
}

#' Parse sample metadata from a GEO series-matrix file
#'
#' Every `!Sample_characteristics_ch1` line of the form `key: value` becomes a column
#' named after the (cleaned) key; title, accession and source are kept too.
#' @param path path to *_series_matrix.txt(.gz)
#' @return data.frame, one row per GSM
#' @export
mf_read_geo_pheno <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  lines <- character()
  repeat {
    chunk <- readLines(con, n = 5000, warn = FALSE)
    if (!length(chunk)) break
    keep <- startsWith(chunk, "!Sample_")
    lines <- c(lines, chunk[keep])
    if (any(startsWith(chunk, "!series_matrix_table_begin"))) break
  }
  if (!length(lines)) mf_stop("No !Sample_ lines found in ", path)
  split <- strsplit(lines, "\t", fixed = TRUE)
  key <- vapply(split, `[`, "", 1)
  vals <- lapply(split, function(v) gsub('^"|"$', "", v[-1]))
  n <- length(vals[[which(key == "!Sample_geo_accession")[1]]])
  out <- data.frame(geo_accession = vals[[which(key == "!Sample_geo_accession")[1]]], stringsAsFactors = FALSE)
  simple <- c(title = "!Sample_title", source_name = "!Sample_source_name_ch1",
              organism = "!Sample_organism_ch1", platform_id = "!Sample_platform_id",
              supplementary_file = "!Sample_supplementary_file")
  for (nm in names(simple)) { i <- which(key == simple[[nm]])[1]; if (!is.na(i)) out[[nm]] <- vals[[i]] }
  for (i in which(key == "!Sample_characteristics_ch1")) {
    v <- vals[[i]]
    has <- grepl(":", v, fixed = TRUE)
    k <- ifelse(has, sub(":.*$", "", v), NA)
    val <- ifelse(has, trimws(sub("^[^:]*:", "", v)), v)
    for (kk in unique(stats::na.omit(k))) {
      col <- clean_names(sub("\\s*\\(.*$", "", kk))
      if (is.null(out[[col]])) out[[col]] <- NA_character_
      sel <- which(k == kk)
      out[[col]][sel] <- val[sel]
    }
  }
  # convert numeric-looking columns
  for (c in names(out)) {
    v <- out[[c]]
    num <- suppressWarnings(as.numeric(v))
    if (all(is.na(v) | v %in% c("", "NA") | !is.na(num)) && any(!is.na(num))) out[[c]] <- num
  }
  stopifnot(nrow(out) == n)
  out
}

#' Read the sample sheet (CSV/TSV) and/or GEO series-matrix metadata
#' @param cfg config list
#' @export
mf_read_samples <- function(cfg) {
  inp <- cfg$input
  ss <- NULL
  if (!is.null(inp$sample_sheet)) {
    ss <- as.data.frame(data.table::fread(inp$sample_sheet, colClasses = "character", na.strings = c("", "NA")))
    for (c in names(ss)) { num <- suppressWarnings(as.numeric(ss[[c]]))
      if (all(is.na(ss[[c]]) | !is.na(num)) && any(!is.na(num))) ss[[c]] <- num }
  }
  if (!is.null(inp$geo_series_matrix)) {
    geo <- mf_read_geo_pheno(inp$geo_series_matrix)
    if (is.null(ss)) ss <- geo else {
      key <- intersect(c("geo_accession", "GSM", "gsm"), names(ss))
      if (!length(key)) mf_stop("To combine a sample sheet with GEO metadata, the sample sheet needs a 'geo_accession' column")
      names(ss)[names(ss) == key[1]] <- "geo_accession"
      ss <- merge(geo, ss, by = "geo_accession", all.y = TRUE, suffixes = c("_geo", ""))
    }
  }
  idc <- inp$sample_id_col %||% (if ("geo_accession" %in% names(ss)) "geo_accession" else names(ss)[1])
  if (!idc %in% names(ss)) mf_stop("sample_id_col '", idc, "' not found in sample sheet. Columns: ", paste(names(ss), collapse = ", "))
  if (anyDuplicated(ss[[idc]])) mf_stop("sample_id_col '", idc, "' has duplicated values")
  ss$sample_id <- as.character(ss[[idc]])
  rownames(ss) <- ss$sample_id
  mf_apply_phenotype(ss, cfg)
}

#' Recode/standardise phenotype columns (group, age, sex) per config
#' @param ss sample data.frame
#' @param cfg config list
#' @export
mf_apply_phenotype <- function(ss, cfg) {
  ph <- cfg$phenotype
  if (!is.null(ph$recode)) for (col in names(ph$recode)) {
    if (!col %in% names(ss)) mf_stop("phenotype$recode refers to missing column '", col, "'")
    map <- ph$recode[[col]]
    v <- as.character(ss[[col]])
    hit <- v %in% names(map)
    v[hit] <- unlist(map)[v[hit]]
    ss[[col]] <- v
  }
  if (!is.null(ph$keep_samples_query)) {
    keep <- tryCatch(with(ss, eval(parse(text = ph$keep_samples_query))),
                     error = function(e) mf_stop("keep_samples_query failed: ", conditionMessage(e)))
    keep[is.na(keep)] <- FALSE
    mf_log("keep_samples_query '", ph$keep_samples_query, "' keeps ", sum(keep), "/", nrow(ss), " samples")
    ss <- ss[keep, , drop = FALSE]
  }
  if (!is.null(ph$group_col)) {
    if (!ph$group_col %in% names(ss)) mf_stop("group_col '", ph$group_col, "' not in sample sheet. Columns: ", paste(names(ss), collapse = ", "))
    g <- as.character(ss[[ph$group_col]])
    lev <- c(ph$control_level, ph$case_level)
    if (length(lev) == 2) {
      bad <- !is.na(g) & !g %in% lev
      if (any(bad)) mf_log(sum(bad), " samples have group values outside {", paste(lev, collapse = ", "), "} and are set to NA", level = "WARN")
      g[bad] <- NA
      ss$group <- factor(g, levels = lev)
    } else ss$group <- factor(g)
  }
  if (!is.null(ph$age_col)) {
    if (!ph$age_col %in% names(ss)) mf_stop("age_col '", ph$age_col, "' not found")
    ss$age <- suppressWarnings(as.numeric(ss[[ph$age_col]]))
  }
  if (!is.null(ph$sex_col)) {
    if (!ph$sex_col %in% names(ss)) mf_stop("sex_col '", ph$sex_col, "' not found")
    s <- toupper(substr(trimws(as.character(ss[[ph$sex_col]])), 1, 1))
    s[!s %in% c("M", "F")] <- NA
    ss$sex <- factor(s, levels = c("F", "M"))
  }
  ss
}

# match data columns (matrix colnames / IDAT basenames) to sample ids
match_samples <- function(ss, data_ids, cfg) {
  dc <- cfg$input$data_id_col
  if (is.null(dc)) {
    ov <- vapply(names(ss), function(c) sum(as.character(ss[[c]]) %in% data_ids), 0)
    if (max(ov) == 0) {
      # fall back to substring matching (e.g. 'GSM123_9876_R01C01' contains 'GSM123')
      hit <- vapply(ss$sample_id, function(s) { m <- grep(s, data_ids, fixed = TRUE, value = TRUE); if (length(m) == 1) m else NA_character_ }, "")
      if (all(is.na(hit))) mf_stop("Could not match any data columns to sample sheet rows. Set input$data_id_col. Example data ids: ", paste(utils::head(data_ids, 3), collapse = ", "))
      ss$data_id <- hit
      mf_log("Matched ", sum(!is.na(hit)), " samples by substring of sample_id")
      return(add_sentrix(ss[!is.na(ss$data_id), , drop = FALSE]))
    }
    dc <- names(which.max(ov))
    mf_log("Auto-detected data_id_col = '", dc, "' (", max(ov), " matches)")
  }
  if (!dc %in% names(ss)) mf_stop("data_id_col '", dc, "' not found in sample sheet")
  ss$data_id <- as.character(ss[[dc]])
  miss <- !ss$data_id %in% data_ids
  if (any(miss)) mf_log(sum(miss), " sample-sheet rows have no data and are dropped (e.g. ", paste(utils::head(ss$sample_id[miss], 3), collapse = ", "), ")", level = "WARN")
  extra <- setdiff(data_ids, ss$data_id)
  if (length(extra)) mf_log(length(extra), " data columns have no sample-sheet row and are ignored", level = "WARN")
  add_sentrix(ss[!miss, , drop = FALSE])
}

# derive slide and array row from Illumina sentrix IDs ("200123456789_R01C01")
add_sentrix <- function(ss) {
  id <- basename(as.character(ss$data_id))
  m <- regmatches(id, regexec("([0-9]{9,})_(R[0-9]{2})(C[0-9]{2})", id))
  ok <- lengths(m) == 4
  if (mean(ok) > 0.9) {
    if (is.null(ss$slide)) ss$slide <- ifelse(ok, vapply(m, function(x) if (length(x) == 4) x[2] else NA_character_, ""), NA)
    if (is.null(ss$array_row)) ss$array_row <- ifelse(ok, vapply(m, function(x) if (length(x) == 4) x[3] else NA_character_, ""), NA)
  }
  ss
}

# ---- array platform ---------------------------------------------------------

#' Guess array type from probe IDs
#' @param probes character vector of probe IDs
#' @export
mf_guess_array <- function(probes) {
  n <- length(probes)
  if (any(grepl("^cg[0-9]{8}_(TC|BC|TO|BO)[0-9]+$", utils::head(probes, 5000)))) return("epicv2")
  if (n > 600000) return("epic")
  if (n > 300000) return("450k")
  # subsets: decide by EPIC-only probes
  if (any(probes %in% c("cg07881041", "cg23229610"))) "epic" else "450k"
}

array_annotation <- function(array) {
  switch(array,
         "450k" = c(array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19"),
         "epic" = c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b4.hg19"),
         "epicv2" = c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38"))
}

anno_pkg <- function(array) {
  switch(array, "450k" = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
         "epic" = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
         "epicv2" = "IlluminaHumanMethylationEPICv2anno.20a1.hg38")
}

#' Probe annotation table for an array type (chr, pos, genes, island context, regulatory flags)
#' @param array 450k, epic or epicv2
#' @export
mf_array_annotation <- function(array) {
  key <- paste0("anno_", array)
  if (!is.null(.mf_env[[key]])) return(.mf_env[[key]])
  pkg <- anno_pkg(array); need_pkg(pkg, "probe annotation"); need_pkg("minfi")
  # minfi annotation objects resolve their tables from the attached package environment
  if (!paste0("package:", pkg) %in% search()) suppressPackageStartupMessages(do.call("require", list(pkg, character.only = TRUE)))
  obj <- get(pkg, envir = as.environment(paste0("package:", pkg)))
  a <- as.data.frame(minfi::getAnnotation(obj))
  pick <- function(cands) { x <- intersect(cands, names(a)); if (length(x)) a[[x[1]]] else rep(NA, nrow(a)) }
  out <- data.frame(
    probe = rownames(a), chr = a$chr, pos = a$pos, strand = a$strand, type = a$Type,
    gene = pick(c("UCSC_RefGene_Name", "GencodeV41_Name")),
    gene_group = pick(c("UCSC_RefGene_Group", "GencodeV41_Group")),
    island = pick(c("Relation_to_Island", "Relation_to_UCSC_CpG_Island")),
    enhancer = pick(c("Enhancer", "X450k_Enhancer", "Phantom5_Enhancers", "Methyl450_Enhancer")),
    cpg_maf = suppressWarnings(as.numeric(pick("CpG_maf"))), sbe_maf = suppressWarnings(as.numeric(pick("SBE_maf"))),
    dhs = pick(c("DHS", "DNase_Hypersensitivity_NAME")),
    regulatory = pick(c("Regulatory_Feature_Group", "Regulatory_Feature_Name")),
    stringsAsFactors = FALSE)
  out$enhancer <- !is.na(out$enhancer) & out$enhancer != "" & out$enhancer != "FALSE"
  out$dhs <- !is.na(out$dhs) & out$dhs != "" & out$dhs != "FALSE"
  out$island[out$island %in% c("", NA)] <- "OpenSea"
  out$island <- sub("^[NS]_", "", out$island)
  rownames(out) <- out$probe
  .mf_env[[key]] <- out
  out
}

# ---- array readers ----------------------------------------------------------

read_matrix <- function(path, what) {
  mf_log("Reading ", what, " matrix: ", basename(path))
  dt <- data.table::fread(path, data.table = FALSE, showProgress = FALSE)
  idcol <- 1
  ids <- as.character(dt[[idcol]])
  m <- as.matrix(dt[, -idcol, drop = FALSE])
  storage.mode(m) <- "double"
  rownames(m) <- ids
  # some GEO matrices interleave 'Detection Pval' columns -- drop them
  dp <- grepl("detection|pval", colnames(m), ignore.case = TRUE)
  if (any(dp) && !all(dp)) { mf_log("Dropping ", sum(dp), " detection-p columns embedded in ", what); m <- m[, !dp, drop = FALSE] }
  mf_log(what, ": ", nrow(m), " probes x ", ncol(m), " samples")
  m
}

#' Read raw array data into a common list structure
#'
#' @param cfg config
#' @param ss sample data.frame from [mf_read_samples]
#' @return list(type, array, ss, rgset|mset|beta, detp, beadcount)
#' @export
mf_read_array <- function(cfg, ss) {
  inp <- cfg$input; need_pkg("minfi")
  if (inp$type == "idat") {
    files <- list.files(inp$idat_dir, pattern = "_Grn\\.idat(\\.gz)?$", recursive = TRUE, full.names = TRUE)
    if (!length(files)) mf_stop("No *_Grn.idat files in ", inp$idat_dir)
    bases <- sub("_Grn\\.idat(\\.gz)?$", "", files)
    if (!is.null(inp$basename_col)) {
      ss$data_id <- file.path(inp$idat_dir, ss[[inp$basename_col]])
      ok <- file.exists(paste0(ss$data_id, "_Grn.idat")) | file.exists(paste0(ss$data_id, "_Grn.idat.gz"))
      if (!all(ok)) mf_log(sum(!ok), " samples without IDATs dropped", level = "WARN")
      ss <- ss[ok, , drop = FALSE]
    } else {
      ss <- match_samples(ss, basename(bases), utils::modifyList(cfg, list(input = list(data_id_col = inp$data_id_col))))
      ss$data_id <- bases[match(ss$data_id, basename(bases))]
    }
    targets <- data.frame(Basename = ss$data_id, sample_id = ss$sample_id)
    rg <- minfi::read.metharray(targets$Basename, extended = TRUE, verbose = FALSE, force = TRUE)
    colnames(rg) <- ss$sample_id
    arr <- if (inp$array == "auto") switch(minfi::annotation(rg)[["array"]],
                IlluminaHumanMethylation450k = "450k", IlluminaHumanMethylationEPIC = "epic",
                IlluminaHumanMethylationEPICv2 = "epicv2", "450k") else inp$array
    if (arr == "epicv2") minfi::annotation(rg) <- array_annotation("epicv2")
    detp <- minfi::detectionP(rg)
    bc <- tryCatch(minfi::getNBeads(rg), error = function(e) NULL)
    return(list(type = "idat", array = arr, ss = ss, rgset = rg, detp = detp, beadcount_rg = bc))
  }
  if (inp$type == "intensity") {
    M <- read_matrix(inp$meth_file, "methylated")
    U <- read_matrix(inp$unmeth_file, "unmethylated")
    ss <- match_samples(ss, colnames(M), cfg)
    common <- intersect(rownames(M), rownames(U))
    M <- M[common, ss$data_id, drop = FALSE]; U <- U[common, ss$data_id, drop = FALSE]
    colnames(M) <- colnames(U) <- ss$sample_id
    detp <- NULL
    if (!is.null(inp$detp_file)) {
      P <- read_matrix(inp$detp_file, "detection p")
      if (!all(ss$data_id %in% colnames(P))) mf_stop("Detection-p file is missing samples present in the intensity files")
      detp <- P[intersect(common, rownames(P)), ss$data_id, drop = FALSE]; colnames(detp) <- ss$sample_id
    } else mf_log("No detection p-values supplied: probe/sample detection QC will use intensity thresholds only", level = "WARN")
    arr <- if (inp$array == "auto") mf_guess_array(rownames(M)) else inp$array
    keep <- !grepl("^rs", rownames(M))   # SNP control probes are not CpGs
    ms <- minfi::MethylSet(Meth = M[keep, ], Unmeth = U[keep, ], annotation = array_annotation(arr))
    return(list(type = "intensity", array = arr, ss = ss, mset = ms, detp = detp,
                snp_probes = if (any(!keep)) M[!keep, , drop = FALSE] / (M[!keep, , drop = FALSE] + U[!keep, , drop = FALSE] + 100) else NULL))
  }
  if (inp$type == "beta") {
    B <- read_matrix(inp$beta_file, "beta")
    ss <- match_samples(ss, colnames(B), cfg)
    B <- B[!grepl("^rs", rownames(B)), ss$data_id, drop = FALSE]; colnames(B) <- ss$sample_id
    if (max(B, na.rm = TRUE) > 1.5) {
      mf_log("Values > 1 detected: treating input as M-values and converting to beta", level = "WARN")
      B <- 2^B / (1 + 2^B)
    }
    detp <- NULL
    if (!is.null(inp$detp_file)) {
      P <- read_matrix(inp$detp_file, "detection p")
      detp <- P[intersect(rownames(B), rownames(P)), ss$data_id, drop = FALSE]; colnames(detp) <- ss$sample_id
    }
    arr <- if (inp$array == "auto") mf_guess_array(rownames(B)) else inp$array
    return(list(type = "beta", array = arr, ss = ss, beta = B, detp = detp))
  }
  mf_stop("mf_read_array called for non-array input")
}

# ---- bisulfite sequencing readers ------------------------------------------

#' Read Bismark coverage files (chr, start, end, percent methylated, count methylated, count unmethylated)
#' @param cfg config
#' @param ss sample data.frame
#' @return list(type = 'bismark', ss, M = methylated counts, Cov = coverage, gr = GRanges)
#' @export
mf_read_bismark <- function(cfg, ss) {
  inp <- cfg$input
  files <- list.files(inp$bismark_dir, pattern = inp$bismark_pattern, full.names = TRUE, recursive = TRUE)
  if (!length(files)) mf_stop("No coverage files matching '", inp$bismark_pattern, "' in ", inp$bismark_dir)
  if (!is.null(inp$bismark_file_col)) {
    ss$data_id <- file.path(inp$bismark_dir, ss[[inp$bismark_file_col]])
    ok <- file.exists(ss$data_id)
    if (!all(ok)) mf_log(sum(!ok), " samples have no coverage file and are dropped", level = "WARN")
    ss <- ss[ok, , drop = FALSE]
  } else {
    ss <- match_samples(ss, basename(files), cfg)
    ss$data_id <- files[match(ss$data_id, basename(files))]
  }
  if (nrow(ss) < 2) mf_stop("Fewer than 2 samples matched to coverage files")
  mf_log("Reading ", nrow(ss), " Bismark coverage files")
  per <- lapply(seq_len(nrow(ss)), function(i) {
    d <- data.table::fread(ss$data_id[i], header = FALSE, showProgress = FALSE,
                           col.names = c("chr", "start", "end", "pct", "m", "u"), select = 1:6)
    # sites below the analysis coverage threshold are treated as missing downstream, so drop them
    # early to keep memory bounded (half the threshold is kept so two strands can still be merged)
    pre <- max(1L, floor((cfg$qc$bs_min_coverage %||% 1L) / 2))
    d <- d[(d$m + d$u) >= pre, ]
    if (!grepl("^chr", d$chr[1])) d$chr <- paste0("chr", d$chr)
    d$key <- paste(d$chr, d$start, sep = ":")
    d
  })
  # keep only sites that can pass the 'present in >= bs_min_sample_frac of samples' rule
  allk <- unlist(lapply(per, `[[`, "key"), use.names = FALSE)
  uk <- unique(allk); cnt <- tabulate(match(allk, uk), nbins = length(uk)); rm(allk)
  need <- max(1, floor((cfg$qc$bs_min_sample_frac %||% 0) * nrow(ss) * 0.5))
  keys <- uk[cnt >= need]; rm(uk, cnt)
  mf_log(length(keys), " candidate CpG sites after per-sample coverage pre-filter")
  M <- Cov <- matrix(0L, nrow = length(keys), ncol = nrow(ss), dimnames = list(keys, ss$sample_id))
  for (i in seq_along(per)) {
    idx <- match(per[[i]]$key, keys); ok <- !is.na(idx)
    M[idx[ok], i] <- as.integer(per[[i]]$m[ok]); Cov[idx[ok], i] <- as.integer(per[[i]]$m[ok] + per[[i]]$u[ok])
  }
  chr <- sub(":.*", "", keys); pos <- as.integer(sub(".*:", "", keys))
  o <- order(chr, pos)
  list(type = "bismark", ss = ss, M = M[o, , drop = FALSE], Cov = Cov[o, , drop = FALSE],
       chr = chr[o], pos = pos[o])
}
