#' Link GEO sample accessions to data-matrix columns by methylation profile
#'
#' Some GEO series provide raw intensity files whose columns are Illumina sentrix IDs
#' while the sample metadata only lists GSM accessions. The series-matrix file, however,
#' contains processed beta values per GSM. Because methylation profiles are
#' individual-specific, each GSM can be matched to its raw-data column by correlation.
#' A match is accepted only if it is reciprocal and the best correlation clearly exceeds
#' the second best.
#' @param series_matrix path to *_series_matrix.txt.gz containing a data table
#' @param meth_file,unmeth_file raw intensity files (or set beta_file)
#' @param beta_file alternative to meth/unmeth
#' @param n_probes number of high-variance probes used
#' @param min_gap minimum difference between best and second-best correlation
#' @return data.frame(geo_accession, data_id, r_best, r_second, reciprocal, accepted)
#' @export
mf_geo_link_by_profile <- function(series_matrix, meth_file = NULL, unmeth_file = NULL, beta_file = NULL,
                                   n_probes = 20000, min_gap = 0.02) {
  con <- gzfile(series_matrix, "rt"); on.exit(close(con))
  skip <- 0L
  repeat { l <- readLines(con, n = 1); skip <- skip + 1L; if (!length(l)) mf_stop("No data table in ", series_matrix)
           if (startsWith(l, "!series_matrix_table_begin")) break }
  close(con); on.exit()
  geo <- data.table::fread(series_matrix, skip = skip, fill = TRUE, data.table = FALSE, showProgress = FALSE)
  geo <- geo[!grepl("^!series_matrix_table_end", geo[[1]]), ]
  rownames(geo) <- geo[[1]]; G <- as.matrix(geo[, -1]); storage.mode(G) <- "double"
  if (!is.null(beta_file)) B <- read_matrix(beta_file, "beta") else {
    M <- read_matrix(meth_file, "methylated"); U <- read_matrix(unmeth_file, "unmethylated")
    common <- intersect(rownames(M), rownames(U)); B <- M[common, ] / (M[common, ] + U[common, colnames(M)] + 100)
  }
  pr <- intersect(rownames(G), rownames(B)); pr <- pr[!grepl("^rs", pr)]
  G <- G[pr, , drop = FALSE]; B <- B[pr, , drop = FALSE]
  ok <- stats::complete.cases(G) & stats::complete.cases(B)
  G <- G[ok, , drop = FALSE]; B <- B[ok, , drop = FALSE]
  v <- matrixStats::rowVars(G); sel <- order(v, decreasing = TRUE)[seq_len(min(n_probes, length(v)))]
  # centre per probe so correlation reflects individual-specific deviations
  Gc <- G[sel, ] - rowMeans(G[sel, ]); Bc <- B[sel, ] - rowMeans(B[sel, ])
  R <- stats::cor(Gc, Bc)
  best <- apply(R, 1, which.max)
  second <- apply(R, 1, function(r) sort(r, decreasing = TRUE)[2])
  back <- apply(R, 2, which.max)
  out <- data.frame(geo_accession = rownames(R), data_id = colnames(R)[best], r_best = R[cbind(seq_len(nrow(R)), best)],
                    r_second = second, stringsAsFactors = FALSE)
  out$reciprocal <- back[best] == seq_len(nrow(R))
  out$accepted <- out$reciprocal & (out$r_best - out$r_second) > min_gap
  mf_log("Profile matching: ", sum(out$accepted), "/", nrow(out), " GSMs linked unambiguously (median r = ",
         round(stats::median(out$r_best), 3), ", median second-best = ", round(stats::median(out$r_second), 3), ")")
  out
}
