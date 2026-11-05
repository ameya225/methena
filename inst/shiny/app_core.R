# Shiny UI/server for methena results. `build_app(bundle)` returns a shiny.appobj.
build_app <- function(bundle) {
  library(shiny)
  dmp <- bundle$dmp; fdr <- bundle$fdr %||% 0.05
  `%||%` <- function(a, b) if (is.null(a)) b else a
  chrs <- unique(dmp$chr[order(suppressWarnings(as.integer(sub("chr", "", sub("X", "23", sub("Y", "24", dmp$chr))))))])
  ui <- fluidPage(
    titlePanel(paste("methena results:", bundle$project)),
    sidebarLayout(
      sidebarPanel(width = 3,
        sliderInput("fdr", "FDR threshold", min = 0.001, max = 0.25, value = fdr, step = 0.001),
        sliderInput("db", "Minimum |adjusted delta beta|", min = 0, max = 0.2, value = 0, step = 0.005),
        selectInput("chr", "Chromosome", c("All", chrs)),
        textInput("gene", "Gene filter (symbol, regex)", ""),
        hr(), helpText(sprintf("Platform: %s | model: %s | lambda = %.3f", bundle$platform, bundle$formula, bundle$lambda)),
        downloadButton("dl", "Download filtered DMPs")),
      mainPanel(width = 9, tabsetPanel(
        tabPanel("Volcano", plotly::plotlyOutput("volcano", height = "520px"), p("Click a point to show the CpG in the 'CpG viewer' tab.")),
        tabPanel("Manhattan", plotly::plotlyOutput("manhattan", height = "460px")),
        tabPanel("DMP table", DT::DTOutput("dmp_table")),
        tabPanel("CpG viewer", selectizeInput("cpg", "CpG", choices = NULL), plotOutput("cpg_plot", height = "420px"), tableOutput("cpg_info")),
        tabPanel("DMRs", DT::DTOutput("dmr_table")),
        tabPanel("VMPs", DT::DTOutput("vmp_table")),
        tabPanel("Genes", DT::DTOutput("gene_table")),
        tabPanel("Pathways", selectInput("coll", "Collection", names(bundle$enrich) %||% "none"), DT::DTOutput("enrich_table")),
        tabPanel("Regulatory / TF", DT::DTOutput("reg_table"), hr(), DT::DTOutput("motif_table")),
        tabPanel("Epigenetic clocks", selectInput("clock", "Clock", setdiff(names(bundle$clocks$values), c("sample_id", "age", "group")) %||% "none"),
                 plotly::plotlyOutput("clock_plot", height = "420px"), DT::DTOutput("clock_bench"), DT::DTOutput("clock_tests")),
        tabPanel("Concordance", DT::DTOutput("conc_table"))
      ))))
  server <- function(input, output, session) {
    filt <- reactive({
      d <- dmp
      if (input$chr != "All") d <- d[d$chr %in% input$chr, ]
      if (nzchar(input$gene)) d <- d[grepl(input$gene, d$gene, ignore.case = TRUE), ]
      d$sig <- d$FDR < input$fdr & abs(d$delta_beta_adj) >= input$db
      d
    })
    updateSelectizeInput(session, "cpg", choices = rownames(bundle$beta_top), server = TRUE)
    output$volcano <- plotly::renderPlotly({
      d <- filt(); keep <- d$P.Value < 0.05 | d$sig
      if (sum(!keep) > 0) { set.seed(1); keep[!keep][sample(sum(!keep), min(sum(!keep), 15000))] <- TRUE }
      d <- d[keep, ]
      d$col <- ifelse(d$sig, ifelse(d$t > 0, "hyper", "hypo"), "n.s.")
      plotly::plot_ly(d, x = ~delta_beta_adj, y = ~-log10(P.Value), color = ~col, colors = c(hyper = "#d7301f", hypo = "#2b8cbe", n.s. = "grey"),
                      text = ~paste(cpg, gene, "<br>FDR", signif(FDR, 3)), key = ~cpg, type = "scattergl", mode = "markers",
                      marker = list(size = 4), source = "v") |>
        plotly::layout(xaxis = list(title = "Adjusted delta beta"), yaxis = list(title = "-log10 P"))
    })
    observeEvent(plotly::event_data("plotly_click", source = "v"), {
      k <- plotly::event_data("plotly_click", source = "v")$key
      if (!is.null(k) && k %in% rownames(bundle$beta_top)) updateSelectizeInput(session, "cpg", choices = rownames(bundle$beta_top), selected = k, server = TRUE)
    })
    output$manhattan <- plotly::renderPlotly({
      d <- filt(); d <- d[d$P.Value < 0.01 | d$sig, ]
      d$cn <- suppressWarnings(as.integer(sub("chr", "", sub("X", "23", sub("Y", "24", d$chr)))))
      d <- d[order(d$cn, d$pos), ]
      off <- tapply(d$pos, d$cn, max); cum <- c(0, cumsum(as.numeric(off)))[seq_along(off)]; names(cum) <- names(off)
      d$x <- d$pos + cum[as.character(d$cn)]
      plotly::plot_ly(d, x = ~x, y = ~-log10(P.Value), color = ~factor(cn %% 2), colors = c("#08519c", "#6baed6"),
                      text = ~paste(cpg, gene, chr, pos), type = "scattergl", mode = "markers", marker = list(size = 4)) |>
        plotly::layout(showlegend = FALSE, xaxis = list(title = "Genome position (P < 0.01 shown)"))
    })
    output$dmp_table <- DT::renderDT({ d <- filt(); d <- d[d$sig, setdiff(names(d), "sig")]
      DT::datatable(utils::head(d, 20000), filter = "top", options = list(pageLength = 15, scrollX = TRUE)) |>
        DT::formatSignif(c("P.Value", "FDR", "delta_beta_adj", "t"), 3) })
    output$dl <- downloadHandler(filename = function() "dmp_filtered.tsv",
                                 content = function(f) utils::write.table(filt()[filt()$sig, ], f, sep = "\t", row.names = FALSE, quote = FALSE))
    output$cpg_plot <- renderPlot({
      req(input$cpg); b <- bundle$beta_top[input$cpg, rownames(bundle$pheno)]
      d <- data.frame(beta = b, group = bundle$pheno$group, age = bundle$pheno$age %||% NA)
      ggplot2::ggplot(d, ggplot2::aes(group, beta, fill = group)) + ggplot2::geom_violin(alpha = 0.4) + ggplot2::geom_boxplot(width = 0.2, outlier.shape = NA) +
        ggplot2::geom_jitter(width = 0.08, size = 0.7, alpha = 0.5) + ggplot2::theme_bw(base_size = 14) + ggplot2::labs(title = input$cpg, y = "Beta", x = NULL) +
        ggplot2::theme(legend.position = "none")
    })
    output$cpg_info <- renderTable(dmp[dmp$cpg == input$cpg, c("cpg", "chr", "pos", "gene", "gene_group", "island", "delta_beta_adj", "P.Value", "FDR")], digits = 5)
    dt <- function(x) DT::renderDT(if (is.null(x)) NULL else DT::datatable(x, filter = "top", options = list(pageLength = 15, scrollX = TRUE)))
    output$dmr_table <- dt(bundle$dmr); output$vmp_table <- dt(bundle$vmp); output$gene_table <- dt(bundle$genes)
    output$enrich_table <- DT::renderDT({ x <- bundle$enrich[[input$coll]]; if (is.null(x)) return(NULL)
      DT::datatable(x[, intersect(c("term", "collection", "N", "DE", "P.DE", "FDR"), names(x))], filter = "top", options = list(pageLength = 15)) })
    output$reg_table <- dt(bundle$reg); output$motif_table <- dt(bundle$motif); output$conc_table <- dt(bundle$compare)
    output$clock_plot <- plotly::renderPlotly({ req(bundle$clocks, input$clock %in% names(bundle$clocks$values)); v <- bundle$clocks$values
      plotly::plot_ly(v, x = ~age, y = v[[input$clock]], color = ~group, type = "scatter", mode = "markers", text = ~sample_id) |>
        plotly::layout(xaxis = list(title = "Chronological age"), yaxis = list(title = input$clock)) })
    output$clock_bench <- dt(bundle$clocks$benchmark); output$clock_tests <- dt(bundle$clocks$group_tests)
  }
  shinyApp(ui, server)
}
