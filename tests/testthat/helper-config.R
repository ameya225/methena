write_cfg <- function(dir, lst) {
  f <- file.path(dir, "config.yaml"); yaml::write_yaml(lst, f); f
}
sim_cfg <- function(sim, out, type = "intensity", norm = "quantile", extra = list()) {
  base <- list(project = "sim", output_dir = out, seed = 7L,
    input = list(type = type, sample_sheet = sim$sample_sheet, sample_id_col = "sample_id", data_id_col = "data_id",
                 tissue = "other"),
    phenotype = list(group_col = "group", control_level = "control", case_level = "case", age_col = "age", sex_col = "sex"),
    normalization = list(method = norm),
    qc = list(drop_outliers = FALSE, bs_collapse_strands = FALSE, remove_cross_reactive = TRUE),
    covariates = list(cell_type = "none"),
    analysis = list(covariates = c("age", "sex"), clocks = FALSE, enrichment = FALSE, motif = FALSE, dmr = TRUE, vmp = TRUE))
  if (type == "intensity") base$input[c("meth_file", "unmeth_file", "detp_file")] <- sim[c("meth_file", "unmeth_file", "detp_file")]
  if (type == "bismark") { base$input$bismark_dir <- sim$bismark_dir; base$input$genome <- "hg19"; base$input$data_id_col <- "data_id" }
  utils::modifyList(base, extra)
}
