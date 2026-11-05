#' Launch the interactive results browser
#' @param results_dir folder containing results_bundle.rds (the 'analysis' output folder)
#' @param port port number
#' @param launch.browser open a browser window
#' @export
mf_shiny <- function(results_dir, port = 3838, launch.browser = interactive()) {
  for (p in c("shiny", "DT", "plotly")) need_pkg(p, "the Shiny app")
  f <- file.path(results_dir, "results_bundle.rds")
  if (!file.exists(f)) mf_stop("results_bundle.rds not found in ", results_dir, ". Run the analysis step first.")
  app <- mf_shiny_app(readRDS(f))
  shiny::runApp(app, port = port, launch.browser = launch.browser)
}

#' Build the Shiny app object from a results bundle (useful for testing / deployment)
#' @param bundle list read from results_bundle.rds
#' @export
mf_shiny_app <- function(bundle) {
  env <- new.env()
  sys.source(system.file("shiny", "app_core.R", package = "methena"), envir = env)
  env$build_app(bundle)
}
