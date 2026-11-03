#' Simulate a small case-control methylation dataset with planted effects (for tests/demo)
#'
#' Array mode writes methylated/unmethylated/detection-p CSVs using real 450K probe IDs plus
#' a sample sheet; sequencing mode writes Bismark coverage files. One sample is given poor
#' detection and (array mode) one a reported-sex mismatch, so QC exclusion can be verified.
#' @param dir output folder
#' @param type 'intensity' or 'bismark'
#' @param n_per_group samples per group
#' @param n_cpg number of CpGs
#' @param n_dmp planted differentially methylated CpGs
#' @param effect planted beta difference
#' @param seed RNG seed
#' @return list with file paths and planted CpG IDs
#' @export
mf_simulate <- function(dir, type = c("intensity", "bismark"), n_per_group = 20, n_cpg = 20000, n_dmp = 200,
                        effect = 0.12, seed = 1) {
  type <- match.arg(type); set.seed(seed)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  n <- 2 * n_per_group
  sid <- sprintf("S%03d", seq_len(n))
  group <- rep(c("control", "case"), each = n_per_group)
  sex <- rep(c("F", "M"), length.out = n)
  age <- round(stats::runif(n, 25, 70), 1)
  if (type == "intensity") {
    need_pkg("IlluminaHumanMethylation450kanno.ilmn12.hg19")
    ann <- mf_array_annotation("450k")
    auto <- ann$probe[!ann$chr %in% c("chrX", "chrY")]
    x <- ann$probe[ann$chr == "chrX"]; y <- ann$probe[ann$chr == "chrY"]
    x <- sample(x, 400)
    probes <- c(sample(auto, n_cpg), x, y)
    # include smoking CpGs so the score is exercised
    probes <- unique(c(probes, intersect(.smoking_cpgs, ann$probe)))
  } else {
    probes <- sprintf("chr%d:%d", sample(1:22, n_cpg, TRUE), sample(1e6:5e7, n_cpg))
    probes <- unique(probes)
  }
  p <- length(probes)
  base <- ifelse(stats::runif(p) < 0.5, stats::rbeta(p, 1.5, 12), stats::rbeta(p, 12, 1.5))
  B <- matrix(base, p, n) + matrix(stats::rnorm(p * n, 0, 0.025), p, n)
  auto_idx <- if (type == "intensity") which(probes %in% auto) else seq_len(p)
  mid <- auto_idx[base[auto_idx] > 0.15 & base[auto_idx] < 0.85]
  if (length(mid) < n_dmp) mid <- auto_idx
  planted <- sample(mid, n_dmp)
  sgn <- ifelse(seq_along(planted) %% 2 == 0, 1, -1)
  B[planted, group == "case"] <- B[planted, group == "case"] + sgn * effect
  agecpg <- sample(setdiff(auto_idx, planted), 200)
  B[agecpg, ] <- B[agecpg, ] + outer(rep(0.003, 200), age - 45)
  B <- pmin(pmax(B, 0.005), 0.995)
  dimnames(B) <- list(probes, sid)
  out <- list(planted = probes[planted], samples = NULL)
  ss <- data.frame(sample_id = sid, data_id = paste0("X", sid), group = group, age = age, sex = sex,
                   batch = rep(c("b1", "b2"), each = 2, length.out = n))
  if (type == "intensity") {
    B[x, sex == "F"] <- pmin(pmax(0.5 + stats::rnorm(length(x) * sum(sex == "F"), 0, 0.15), 0.01), 0.99)
    tot <- matrix(exp(stats::rnorm(p * n, 9, 0.35)), p, n, dimnames = dimnames(B))
    tot[y, sex == "F"] <- exp(stats::rnorm(length(y) * sum(sex == "F"), 5.5, 0.3))
    M <- round(B * tot); U <- round((1 - B) * tot)
    D <- matrix(0, p, n, dimnames = dimnames(B))
    D[cbind(sample(p * n, p * n * 0.001))] <- 0.5
    D[, 3] <- ifelse(stats::runif(p) < 0.15, 0.5, 0)       # bad sample
    D[y, sex == "F"] <- 0.5                                 # Y probes fail in females
    ss$sex[6] <- ifelse(ss$sex[6] == "F", "M", "F")          # sex mismatch (reported)
    wr <- function(m, f) { d <- data.frame(ID = rownames(m), m, check.names = FALSE); colnames(d)[-1] <- paste0("X", colnames(m))
                           data.table::fwrite(d, f) ; f }
    out$meth_file <- wr(M, file.path(dir, "meth.csv.gz")); out$unmeth_file <- wr(U, file.path(dir, "unmeth.csv.gz"))
    out$detp_file <- wr(D, file.path(dir, "detp.csv.gz"))
    out$bad_samples <- sid[c(3, 6)]
  } else {
    ss$data_id <- paste0(sid, ".bismark.cov.gz")
    for (i in seq_len(n)) {
      cv <- stats::rpois(p, 25); m <- stats::rbinom(p, cv, B[, i])
      ch <- sub(":.*", "", probes); ps <- as.integer(sub(".*:", "", probes))
      d <- data.frame(ch, ps, ps, round(100 * m / pmax(cv, 1), 2), m, cv - m)[cv > 0, ]
      d <- d[order(d[[1]], d[[2]]), ]
      data.table::fwrite(d, file.path(dir, ss$data_id[i]), sep = "\t", col.names = FALSE)
    }
    out$bismark_dir <- dir
  }
  out$sample_sheet <- file.path(dir, "samples.csv"); utils::write.csv(ss, out$sample_sheet, row.names = FALSE)
  out$beta <- B; out$ss <- ss
  out
}
