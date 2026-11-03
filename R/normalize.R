#' Normalise array data
#'
#' Methods: noob, noob_bmiq, funnorm, quantile, swan (IDAT); dasen, quantile, bmiq, none
#' (intensity matrices); bmiq, none (beta matrices). 'auto' picks noob (IDAT), dasen
#' (intensities) or none (beta).
#' @param raw list from [mf_read_array] (already subset to passing samples)
#' @param cfg config
#' @return beta matrix (probes x samples)
#' @export
mf_normalize_array <- function(raw, cfg) {
  meth <- cfg$normalization$method
  if (meth == "auto") meth <- switch(raw$type, idat = "noob", intensity = "dasen", beta = "none")
  allowed <- switch(raw$type, idat = c("noob", "noob_bmiq", "funnorm", "quantile", "swan", "none"),
                    intensity = c("dasen", "quantile", "bmiq", "none"), beta = c("bmiq", "none"))
  if (!meth %in% allowed) mf_stop("Normalisation '", meth, "' is not available for input type '", raw$type,
                                  "'. Choose one of: ", paste(allowed, collapse = ", "))
  mf_log("Normalisation method: ", meth)
  design_v <- function(probes) {
    ann <- mf_array_annotation(raw$array)
    t <- ann[probes, "type"]; ifelse(t == "I", 1L, 2L)
  }
  bmiq <- function(b) {
    need_pkg("wateRmelon", "BMIQ")
    b <- b[stats::complete.cases(b), , drop = FALSE]
    dv <- design_v(rownames(b)); ok <- !is.na(dv); b <- b[ok, , drop = FALSE]; dv <- dv[ok]
    b <- pmin(pmax(b, 1e-6), 1 - 1e-6)
    out <- vapply(seq_len(ncol(b)), function(i)
      suppressMessages(wateRmelon::BMIQ(b[, i], dv, plots = FALSE, pri = FALSE)$nbeta), numeric(nrow(b)))
    dimnames(out) <- dimnames(b); out
  }
  if (raw$type == "idat") {
    rg <- raw$rgset
    b <- switch(meth,
      noob = minfi::getBeta(minfi::preprocessNoob(rg), offset = 100),
      noob_bmiq = bmiq(minfi::getBeta(minfi::preprocessNoob(rg), offset = 100)),
      funnorm = minfi::getBeta(minfi::preprocessFunnorm(rg, verbose = FALSE)),
      quantile = minfi::getBeta(minfi::preprocessQuantile(rg, verbose = FALSE)),
      swan = minfi::getBeta(minfi::preprocessSWAN(rg), offset = 100),
      none = minfi::getBeta(minfi::preprocessRaw(rg), offset = 100))
  } else if (raw$type == "intensity") {
    ms <- raw$mset
    b <- switch(meth,
      dasen = { need_pkg("wateRmelon", "dasen"); need_pkg("preprocessCore", "dasen")
                dv <- design_v(rownames(ms)); ok <- !is.na(dv)
                roco <- sentrix_roco(raw$ss$data_id)
                M <- dasen_channel(minfi::getMeth(ms)[ok, , drop = FALSE], dv[ok], roco)
                U <- dasen_channel(minfi::getUnmeth(ms)[ok, , drop = FALSE], dv[ok], NULL)
                rm(ms); raw$mset <- NULL; gc(verbose = FALSE)
                M / (M + U + 100) },
      quantile = minfi::getBeta(suppressWarnings(minfi::preprocessQuantile(ms, verbose = FALSE))),
      bmiq = bmiq(minfi::getBeta(ms, offset = 100)),
      none = minfi::getBeta(ms, offset = 100))
  } else {
    b <- switch(meth, bmiq = bmiq(raw$beta), none = raw$beta)
  }
  b <- as.matrix(b)
  b[b < 0] <- 0; b[b > 1] <- 1
  b
}

# Sentrix row/column (e.g. R01C02 -> row 1, col 2) used by dasen's positional background model
sentrix_roco <- function(ids) {
  m <- regmatches(ids, regexpr("R[0-9]{2}C[0-9]{2}", ids))
  if (length(m) != length(ids)) return(NULL)
  m
}

#' Memory-efficient re-implementation of wateRmelon::dasen for one channel
#'
#' Identical algorithm to dasen (Pidsley et al. 2013): Infinium I intensities are shifted by
#' the per-sample difference between Infinium I and II background peaks (wateRmelon's dfs2,
#' optionally smoothed by a linear model on array row/column), then type I and type II probes are
#' quantile-normalised separately. Works in place to stay within memory for ~1000-sample studies.
#' @param X intensity matrix (probes x samples)
#' @param dv design vector (1 = Infinium I, 2 = Infinium II)
#' @param roco optional sentrix positions (R01C01 ...) for the positional fit
#' @keywords internal
dasen_channel <- function(X, dv, roco = NULL) {
  onetwo <- ifelse(dv == 1, "I", "II")
  dfs2 <- utils::getFromNamespace("dfs2", "wateRmelon")
  mdf <- apply(X, 2, dfs2, onetwo)
  if (!is.null(roco)) {
    srow <- as.numeric(substr(roco, 2, 3)); scol <- as.numeric(substr(roco, 5, 6))
    fit <- try(stats::lm(mdf ~ srow + scol), silent = TRUE)
    if (!inherits(fit, "try-error")) mdf <- stats::fitted(fit)
  }
  i1 <- which(onetwo == "I"); i2 <- which(onetwo == "II")
  X[i1, ] <- sweep(X[i1, , drop = FALSE], 2, mdf)
  X[i1, ] <- preprocessCore::normalize.quantiles(X[i1, , drop = FALSE], copy = FALSE)
  X[i2, ] <- preprocessCore::normalize.quantiles(X[i2, , drop = FALSE], copy = FALSE)
  X
}
