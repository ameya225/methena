skip_if_not_installed("minfi"); skip_if_not_installed("IlluminaHumanMethylation450kanno.ilmn12.hg19")

sim_dir <- withr::local_tempdir(.local_envir = teardown_env())
sim <- mf_simulate(file.path(sim_dir, "sim"), "intensity", n_per_group = 20, n_cpg = 15000, n_dmp = 150)

run_once <- function(out, norm = "quantile", extra = list()) {
  cfgf <- write_cfg(sim_dir, sim_cfg(sim, out, norm = norm, extra = extra))
  cfg <- mf_read_config(cfgf)
  suppressMessages({ mf_preprocess(cfg); r <- mf_analyze(cfg) })
  list(cfg = cfg, r = r, se = readRDS(file.path(out, "preprocess", "methena_dataset.rds")))
}
res1 <- run_once(file.path(sim_dir, "out1"))

test_that("QC removes the failing and sex-mismatched samples", {
  sq <- read.delim(file.path(sim_dir, "out1", "preprocess", "sample_qc.tsv"))
  expect_setequal(sq$sample_id[!sq$included], sim$bad_samples)
  expect_match(sq$exclusion_reason[sq$sample_id == sim$bad_samples[1]], "detection_p")
  expect_match(sq$exclusion_reason[sq$sample_id == sim$bad_samples[2]], "sex_mismatch")
  expect_equal(ncol(res1$se), 38)
  ann <- SummarizedExperiment::rowData(res1$se)
  expect_false(any(ann$chr %in% c("chrX", "chrY")))          # XY removed
})

test_that("planted DMPs are recovered with controlled false discoveries", {
  t <- res1$r$dmp$table
  sig <- t$cpg[t$FDR < 0.05]
  planted <- intersect(sim$planted, t$cpg)
  sens <- mean(planted %in% sig); fdp <- mean(!sig %in% planted)
  expect_gt(sens, 0.85); expect_lt(fdp, 0.1)
  # direction of effect is recovered
  tp <- t[t$cpg %in% planted, ]
  truth <- sign(rowMeans(sim$beta[tp$cpg, sim$ss$group == "case"]) - rowMeans(sim$beta[tp$cpg, sim$ss$group == "control"]))
  expect_gt(mean(sign(tp$delta_beta_adj) == truth), 0.95)
  expect_lt(abs(res1$r$dmp$lambda - 1), 0.3)
})

test_that("all expected output files exist and the report is self-contained", {
  a <- file.path(sim_dir, "out1", "analysis")
  for (f in c("dmp_results.tsv.gz", "vmp_results.tsv.gz", "gene_level_regulation.tsv", "results_bundle.rds",
              "analysis_report.html", "run_manifest.json")) expect_true(file.exists(file.path(a, f)), info = f)
  for (f in c("qc_report.html", "sample_qc.tsv", "probe_filter_summary.tsv", "run_manifest.json"))
    expect_true(file.exists(file.path(sim_dir, "out1", "preprocess", f)), info = f)
  html <- readLines(file.path(a, "analysis_report.html"), warn = FALSE)
  expect_true(any(grepl("data:image/png;base64", html)))
  man <- jsonlite::read_json(file.path(a, "run_manifest.json"))
  expect_true(nchar(man$parameter_digest) == 64)
})

test_that("re-running gives byte-identical result tables", {
  res2 <- run_once(file.path(sim_dir, "out2"))
  for (f in c("dmp_results.tsv.gz", "vmp_results.tsv.gz", "gene_level_regulation.tsv")) {
    a <- data.table::fread(file.path(sim_dir, "out1", "analysis", f)); b <- data.table::fread(file.path(sim_dir, "out2", "analysis", f))
    expect_identical(a, b)
  }
  expect_identical(digest::digest(file = file.path(sim_dir, "out1", "analysis", "dmp_results.tsv.gz")),
                   digest::digest(file = file.path(sim_dir, "out2", "analysis", "dmp_results.tsv.gz")))
})

test_that("alternative normalisations run", {
  skip_if_not_installed("wateRmelon")
  for (m in c("dasen", "none")) {
    r <- run_once(file.path(sim_dir, paste0("out_", m)), norm = m)
    expect_gt(sum(r$r$dmp$table$FDR < 0.05), 100)
  }
  expect_error(run_once(file.path(sim_dir, "out_bad"), norm = "funnorm"), "not available")
})

test_that("design handles missing covariates, NA values and bad groups", {
  cd <- as.data.frame(SummarizedExperiment::colData(res1$se))
  cfg <- res1$cfg
  cfg$analysis$covariates <- c("age", "sex", "not_a_column")
  d <- mf_design(cd, cfg); expect_false("not_a_column" %in% d$covariates)
  cd2 <- cd; cd2$age[1:3] <- NA
  expect_equal(nrow(mf_design(cd2, cfg)$cd), nrow(cd) - 3)
  cd3 <- cd; cd3$group <- factor("case")
  expect_error(mf_design(cd3, cfg), "exactly 2 groups")
  cd4 <- cd; cd4$group <- NULL
  expect_error(mf_design(cd4, cfg), "group_col")
})

test_that("beta-matrix and M-value inputs are accepted", {
  B <- res1$se@assays@data$beta
  bf <- file.path(sim_dir, "beta.csv"); data.table::fwrite(data.frame(ID = rownames(B), B, check.names = FALSE), bf)
  ss <- data.frame(sample_id = colnames(B), group = as.character(res1$se$group), age = res1$se$age, sex = as.character(res1$se$sex))
  sf <- file.path(sim_dir, "beta_samples.csv"); write.csv(ss, sf, row.names = FALSE)
  lst <- list(project = "b", output_dir = file.path(sim_dir, "out_beta"),
              input = list(type = "beta", beta_file = bf, sample_sheet = sf, sample_id_col = "sample_id", tissue = "other"),
              phenotype = list(group_col = "group", control_level = "control", case_level = "case", age_col = "age", sex_col = "sex"),
              qc = list(drop_outliers = FALSE), covariates = list(cell_type = "none"),
              analysis = list(covariates = c("age"), clocks = FALSE, enrichment = FALSE, vmp = FALSE, dmr = FALSE))
  cfg <- mf_read_config(write_cfg(sim_dir, lst))
  se <- suppressMessages(mf_preprocess(cfg)); expect_equal(ncol(se), ncol(B))
  M <- mf_beta2m(B); mf2 <- file.path(sim_dir, "mvals.csv"); data.table::fwrite(data.frame(ID = rownames(M), M, check.names = FALSE), mf2)
  lst$input$beta_file <- mf2; lst$output_dir <- file.path(sim_dir, "out_m")
  se2 <- suppressMessages(mf_preprocess(mf_read_config(write_cfg(sim_dir, lst))))
  expect_lt(max(abs(SummarizedExperiment::assay(se2) - SummarizedExperiment::assay(se)[rownames(se2), ]), na.rm = TRUE), 1e-3)
})

test_that("sample sheet mismatches are reported, not silently ignored", {
  lst <- sim_cfg(sim, file.path(sim_dir, "out_mis"))
  ss <- read.csv(sim$sample_sheet); ss$data_id[1:2] <- c("nope1", "nope2")
  f <- file.path(sim_dir, "ss_mis.csv"); write.csv(ss, f, row.names = FALSE); lst$input$sample_sheet <- f
  cfg <- mf_read_config(write_cfg(sim_dir, lst))
  raw <- suppressMessages(mf_read_array(cfg, mf_read_samples(cfg)))
  expect_equal(nrow(raw$ss), nrow(ss) - 2)
  lst$input$data_id_col <- "no_such_col"
  cfg <- mf_read_config(write_cfg(sim_dir, lst))
  expect_error(suppressMessages(mf_read_array(cfg, mf_read_samples(cfg))), "not found")
})

test_that("Shiny app builds from the results bundle and its server runs", {
  skip_if_not_installed("shiny"); skip_if_not_installed("plotly"); skip_if_not_installed("DT")
  b <- readRDS(file.path(sim_dir, "out1", "analysis", "results_bundle.rds"))
  app <- mf_shiny_app(b)
  expect_s3_class(app, "shiny.appobj")
  shiny::testServer(app, {
    session$setInputs(fdr = 0.05, db = 0, chr = "All", gene = "", cpg = b$dmp$cpg[1], coll = "none", clock = "none")
    expect_true(sum(filt()$sig) > 50)
    session$setInputs(db = 0.5); expect_equal(sum(filt()$sig), 0)
  })
})

test_that("outlier removal safeguard keeps samples when too many are flagged", {
  r <- run_once(file.path(sim_dir, "out_outl"), extra = list(qc = list(drop_outliers = TRUE, outlier_sd = 0.3)))
  expect_equal(ncol(r$se), 38)
  expect_true(any(r$se$pca_outlier))
})

test_that("bacon correction can be switched off and both settings recover the signal", {
  r <- run_once(file.path(sim_dir, "out_nobacon"), extra = list(analysis = list(inflation_correction = "none")))
  t <- r$r$dmp$table
  expect_null(t$P.uncorrected)
  expect_gt(mean(intersect(sim$planted, t$cpg) %in% t$cpg[t$FDR < 0.05]), 0.85)
  t1 <- res1$r$dmp$table
  if (!is.null(res1$r$dmp$bacon)) {           # correction applied only for inflation > 1, never increases hits
    expect_gt(res1$r$dmp$bacon[["inflation"]], 1)
    expect_lte(sum(t1$FDR < 0.05), sum(t1$FDR.uncorrected < 0.05))
  } else expect_null(t1$P.uncorrected)
})

test_that("age_or_dnam keyword fills missing ages", {
  cd <- as.data.frame(SummarizedExperiment::colData(res1$se)); cd$dnam_age <- cd$age + 1; cd$age[1:4] <- NA
  cfg <- res1$cfg; cfg$analysis$covariates <- c("age_or_dnam", "sex")
  d <- mf_design(cd, cfg)
  expect_equal(nrow(d$cd), nrow(cd)); expect_true("age_or_dnam" %in% d$covariates)
})
