# ---- helpers ---------------------------------------------------------------

beta_from_mu <- function(M, U, offset = 100) M / (M + U + offset)

#' Convert beta values to M-values with clamping
#' @param b beta matrix
#' @param eps clamp
#' @export
mf_beta2m <- function(b, eps = 1e-4) { b <- pmin(pmax(b, eps), 1 - eps); log2(b / (1 - b)) }

row_impute <- function(x) {
  if (!anyNA(x)) return(x)
  med <- matrixStats::rowMedians(x, na.rm = TRUE)
  idx <- which(is.na(x), arr.ind = TRUE)
  x[idx] <- med[idx[, 1]]
  x
}

top_var_pca <- function(b, n = 10000, k = 10) {
  b <- row_impute(b)
  b <- b[stats::complete.cases(b), , drop = FALSE]
  v <- matrixStats::rowVars(b)
  sel <- order(v, decreasing = TRUE)[seq_len(min(n, nrow(b)))]
  x <- t(b[sel, , drop = FALSE]); x <- scale(x, scale = FALSE)
  k <- min(k, ncol(b) - 1, nrow(x) - 1)
  sv <- svd(x, nu = k, nv = 0)
  pcs <- sv$u[, seq_len(k), drop = FALSE] %*% diag(sv$d[seq_len(k)], k)
  colnames(pcs) <- paste0("PC", seq_len(k)); rownames(pcs) <- colnames(b)
  list(pcs = pcs, var_explained = sv$d[seq_len(k)]^2 / sum(sv$d^2))
}

# ---- array QC --------------------------------------------------------------

#' Sample- and probe-level QC metrics for array data (before filtering)
#' @param raw list from [mf_read_array]
#' @param cfg config
#' @return list(sample = data.frame, probe_fail = logical, probe_fail_frac = numeric, mset = MethylSet or NULL)
#' @export
mf_qc_array <- function(raw, cfg) {
  q <- cfg$qc; ss <- raw$ss
  mset <- switch(raw$type, idat = minfi::preprocessRaw(raw$rgset), intensity = raw$mset, beta = NULL)
  sq <- data.frame(sample_id = ss$sample_id, stringsAsFactors = FALSE)
  if (!is.null(mset)) {
    qc <- minfi::getQC(mset)
    sq$median_log2_meth <- qc$mMed; sq$median_log2_unmeth <- qc$uMed
    sq$low_intensity <- (qc$mMed + qc$uMed) / 2 < 10.5
  }
  probe_fail <- NULL; probe_fail_frac <- NULL
  if (!is.null(raw$detp)) {
    D <- raw$detp[, ss$sample_id, drop = FALSE]
    if (!is.null(mset)) D <- D[intersect(rownames(D), rownames(mset)), , drop = FALSE]
    failed <- D > q$detp_threshold
    failed[is.na(failed)] <- TRUE
    # chrY probes legitimately fail in females: exclude them from the per-sample failure rate
    ychr <- tryCatch(rownames(D) %in% mf_array_annotation(raw$array)$probe[mf_array_annotation(raw$array)$chr == "chrY"], error = function(e) rep(FALSE, nrow(D)))
    sq$frac_failed_probes <- colMeans(failed[!ychr, , drop = FALSE])
    probe_fail_frac <- rowMeans(failed)
    sq$detection_fail <- sq$frac_failed_probes > q$sample_max_failed_frac
  } else {
    sq$frac_failed_probes <- NA_real_; sq$detection_fail <- FALSE
  }
  # bead counts (IDAT only)
  bead_fail <- NULL
  if (raw$type == "idat" && requireNamespace("wateRmelon", quietly = TRUE)) {
    bc <- tryCatch(wateRmelon::beadcount(raw$rgset), error = function(e) NULL)
    if (!is.null(bc)) {
      low <- is.na(bc) | bc < q$min_beadcount
      bead_fail <- rowMeans(low) > q$probe_max_failed_frac
      names(bead_fail) <- rownames(bc)
      sq$frac_low_beadcount <- colMeans(low)
    }
  }
  # predicted sex from X/Y intensities
  if (!is.null(mset) && isTRUE(cfg$covariates$predict_sex)) {
    gm <- tryCatch(minfi::mapToGenome(mset), error = function(e) NULL)
    if (!is.null(gm)) {
      sx <- tryCatch(minfi::getSex(gm), error = function(e) NULL)
      if (!is.null(sx)) { sq$predicted_sex <- sx$predictedSex; sq$xMed <- sx$xMed; sq$yMed <- sx$yMed }
    }
  }
  if (is.null(sq$predicted_sex) && raw$type == "beta" && isTRUE(cfg$covariates$predict_sex)) {
    sq$predicted_sex <- predict_sex_beta(raw$beta, raw$array)
  }
  if (!is.null(sq$predicted_sex) && !is.null(ss$sex)) {
    sq$reported_sex <- as.character(ss$sex)
    sq$sex_mismatch <- !is.na(sq$reported_sex) & sq$predicted_sex != sq$reported_sex
  } else sq$sex_mismatch <- FALSE
  # genotype (rs) probes: unexpected duplicates
  snpb <- switch(raw$type, idat = tryCatch(minfi::getSnpBeta(raw$rgset), error = function(e) NULL),
                 intensity = raw$snp_probes, NULL)
  sq$possible_duplicate_of <- NA_character_
  if (!is.null(snpb) && nrow(snpb) >= 20 && ncol(snpb) > 1) {
    r <- suppressWarnings(stats::cor(snpb, use = "pairwise.complete.obs")); diag(r) <- NA
    for (i in seq_len(ncol(r))) { j <- which(r[i, ] > 0.95); if (length(j)) sq$possible_duplicate_of[i] <- paste(colnames(r)[j], collapse = ";") }
  }
  list(sample = sq, probe_fail_frac = probe_fail_frac, bead_fail = bead_fail, mset = mset)
}

predict_sex_beta <- function(beta, array) {
  ann <- tryCatch(mf_array_annotation(array), error = function(e) NULL)
  if (is.null(ann)) return(rep(NA_character_, ncol(beta)))
  x <- intersect(rownames(beta), ann$probe[ann$chr == "chrX"])
  if (length(x) < 100) return(rep(NA_character_, ncol(beta)))
  # females carry intermediate methylation at X-inactivated CpG islands
  isl <- intersect(x, ann$probe[ann$island == "Island"])
  m <- colMeans(beta[isl, , drop = FALSE] > 0.2 & beta[isl, , drop = FALSE] < 0.8, na.rm = TRUE)
  km <- stats::kmeans(m, centers = range(m))
  hi <- which.max(km$centers)
  ifelse(km$cluster == hi, "F", "M")
}

# ---- BS-seq QC -------------------------------------------------------------

#' QC for bisulfite-sequencing counts
#' @param raw list from [mf_read_bismark]
#' @param cfg config
#' @export
mf_qc_bismark <- function(raw, cfg) {
  q <- cfg$qc
  cov <- raw$Cov
  covered <- cov >= q$bs_min_coverage
  sq <- data.frame(sample_id = colnames(cov),
                   n_cpgs_any = colSums(cov > 0), n_cpgs_min_cov = colSums(covered),
                   median_coverage = apply(cov, 2, function(x) stats::median(x[x > 0])),
                   mean_methylation = colSums(raw$M) / pmax(1, colSums(cov)),
                   stringsAsFactors = FALSE)
  sq$low_cpg_count <- sq$n_cpgs_min_cov < q$bs_min_cpgs_per_sample
  ychr <- raw$chr %in% c("chrY")
  if (any(ychr)) {
    yfrac <- colSums(cov[ychr, , drop = FALSE]) / pmax(1, colSums(cov))
    sq$chrY_fraction <- yfrac
    if (isTRUE(cfg$covariates$predict_sex) && length(unique(yfrac > 1e-4)) > 1) {
      sq$predicted_sex <- ifelse(yfrac > 0.1 * max(yfrac), "M", "F")
      if (!is.null(raw$ss$sex)) {
        sq$reported_sex <- as.character(raw$ss$sex)
        sq$sex_mismatch <- !is.na(sq$reported_sex) & sq$predicted_sex != sq$reported_sex
      }
    }
  }
  if (is.null(sq$sex_mismatch)) sq$sex_mismatch <- FALSE
  list(sample = sq)
}

# ---- outliers on normalised data -------------------------------------------

#' Flag samples that are outliers on the top principal components
#' @param beta normalised beta matrix
#' @param sd_cut robust z threshold
#' @export
mf_flag_outliers <- function(beta, sd_cut = 4) {
  if (ncol(beta) < 5) return(list(outlier = rep(FALSE, ncol(beta)), pca = NULL))
  p <- top_var_pca(beta, n = 20000, k = 5)
  z <- apply(p$pcs[, 1:min(3, ncol(p$pcs)), drop = FALSE], 2, function(x) abs(x - stats::median(x)) / stats::mad(x))
  out <- apply(z, 1, function(r) any(r > sd_cut))
  # also a beta-distribution check: correlation to the median sample profile
  med <- matrixStats::rowMedians(row_impute(beta))
  r <- apply(row_impute(beta), 2, function(x) stats::cor(x, med))
  rz <- (r - stats::median(r)) / stats::mad(r)
  list(outlier = out | rz < -sd_cut * 1.5, pca = p, cor_to_median = r)
}
