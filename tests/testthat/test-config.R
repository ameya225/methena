test_that("config validation gives readable errors", {
  d <- withr::local_tempdir()
  expect_error(mf_read_config(file.path(d, "nope.yaml")), "not found")
  writeLines("input: [unclosed", file.path(d, "bad.yaml"))
  expect_error(mf_read_config(file.path(d, "bad.yaml")), "not valid YAML")
  f <- write_cfg(d, list(input = list(type = "foo")))
  expect_error(mf_read_config(f), "input\\$type must be one of")
  f <- write_cfg(d, list(input = list(type = "beta", beta_file = "missing.csv", sample_sheet = "missing.csv")))
  expect_error(mf_read_config(f), "does not exist")
  file.create(file.path(d, "s.csv")); file.create(file.path(d, "b.csv"))
  f <- write_cfg(d, list(input = list(type = "beta", beta_file = "b.csv", sample_sheet = "s.csv"), analysis = list(fdr = 2)))
  expect_error(mf_read_config(f), "fdr")
  f <- write_cfg(d, list(input = list(type = "beta", beta_file = "b.csv", sample_sheet = "s.csv"), normalization = list(method = "magic")))
  expect_error(mf_read_config(f), "normalization")
  f <- write_cfg(d, list(input = list(type = "beta", beta_file = "b.csv", sample_sheet = "s.csv")))
  cfg <- mf_read_config(f)
  expect_equal(normalizePath(cfg$input$beta_file), normalizePath(file.path(d, "b.csv")))   # relative paths resolved
  expect_equal(cfg$qc$detp_threshold, 0.01)                                                # defaults filled
})

test_that("GEO series matrix metadata is parsed", {
  d <- withr::local_tempdir(); f <- file.path(d, "GSE1_series_matrix.txt")
  writeLines(c('!Series_title\t"x"', '!Sample_title\t"a"\t"b"', '!Sample_geo_accession\t"GSM1"\t"GSM2"',
               '!Sample_characteristics_ch1\t"Sex: M"\t"Sex: F"', '!Sample_characteristics_ch1\t"age (years): 40.5"\t"age (years): 51"',
               '!Sample_characteristics_ch1\t"disease_status: 1"\t"disease_status: 2"', '!series_matrix_table_begin', '"ID_REF"\t"GSM1"\t"GSM2"'), f)
  ph <- mf_read_geo_pheno(f)
  expect_equal(ph$geo_accession, c("GSM1", "GSM2"))
  expect_equal(ph$sex, c("M", "F")); expect_equal(ph$age, c(40.5, 51)); expect_true(is.numeric(ph$disease_status))
})

test_that("array type is guessed from probe IDs", {
  expect_equal(mf_guess_array(c("cg00000029_TC21", "cg1")), "epicv2")
  expect_equal(mf_guess_array(paste0("cg", 1:10)), "450k")
})

test_that("beta/M conversion is bounded and invertible", {
  m <- mf_beta2m(c(0, 0.5, 1)); expect_true(all(is.finite(m))); expect_equal(m[2], 0)
})
