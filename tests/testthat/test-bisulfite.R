skip_if_not_installed("GenomicRanges")
test_that("Bismark coverage input runs end to end and recovers planted DMPs", {
  d <- withr::local_tempdir()
  sim <- mf_simulate(file.path(d, "bs"), "bismark", n_per_group = 8, n_cpg = 6000, n_dmp = 100, effect = 0.25)
  cfg <- mf_read_config(write_cfg(d, sim_cfg(sim, file.path(d, "out"), type = "bismark",
                                             extra = list(qc = list(bs_min_coverage = 5L, bs_min_sample_frac = 0.8, drop_outliers = FALSE, bs_collapse_strands = FALSE)))))
  se <- suppressMessages(mf_preprocess(cfg))
  expect_true("coverage" %in% SummarizedExperiment::assayNames(se))
  expect_equal(ncol(se), 16)
  r <- suppressMessages(mf_analyze(cfg))
  t <- r$dmp$table; sig <- t$cpg[t$FDR < 0.05]
  expect_gt(mean(intersect(sim$planted, t$cpg) %in% sig), 0.7)
  expect_lt(mean(!sig %in% sim$planted), 0.15)
})

test_that("Bismark reader errors clearly when no files match", {
  d <- withr::local_tempdir(); write.csv(data.frame(sample_id = c("a", "b"), group = c("x", "y")), file.path(d, "s.csv"), row.names = FALSE)
  dir.create(file.path(d, "cov"))
  cfg <- mf_read_config(write_cfg(d, list(input = list(type = "bismark", bismark_dir = "cov", sample_sheet = "s.csv"))))
  expect_error(mf_read_bismark(cfg, mf_read_samples(cfg)), "No coverage files")
})
