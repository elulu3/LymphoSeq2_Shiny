library(shiny)
library(shinycssloaders)
library(LymphoSeqTest)
library(DT)
library(tidyverse)
library(ggplot2)
library(plotly)
library(heatmaply)
library(pheatmap)
library(ggalluvial)
library(htmltools)
library(sp)
library(chorddiag)
library(writexl)
library(shinyjs)
library(wordcloud2)
library(webshot)
library(htmlwidgets)

options(shiny.maxRequestSize = 50 * 1024^2)

ui <- fluidPage(
    shinyjs::useShinyjs(),
    titlePanel("LymphoSeqTest Shiny"),
    tags$head(
        tags$style(
            HTML(".shiny-notification {
                position:fixed;
                top: calc(50%);
                left: calc(50%);
                }")
        )
    ),
    sidebarLayout(
        sidebarPanel(
            fileInput("airr_files", label = "Upload Files", multiple = TRUE, accept = ".tsv"),

            conditionalPanel(condition = "input.tabselected == 'chord_diagram'",
                selectizeInput("vdj_association", label = "Select VDJ Association", choices = c("", "VJ", "DJ")),
                actionButton("chord_button", label = "Create Chord Diagram")),

            conditionalPanel(condition = "input.tabselected == 'common_panel' && input.common_sub_tab == 'common_bar'",
                selectizeInput("bar_id", label = "Select at least 2 repertoire ids", choices = NULL, multiple = TRUE),
                selectizeInput("color_rep_id", label = "Select 1 repertoire id to color", choices = NULL, selected = NULL),
                selectizeInput("color_intersect", label = "Select at least 2 repertoire ids to color intersections", 
                                choices = NULL, multiple = TRUE, selected = NULL),
                actionButton("bar_button", label = "Create Bar Chart")),

            conditionalPanel(condition = "input.tabselected == 'common_panel' && input.common_sub_tab == 'common_plot'",
                selectizeInput("plot_id", label = "Select 2 repertoire ids", choices = NULL, multiple = TRUE),
                # radioButtons("show_type", "Show Sequences Type", choices = c("common" = "common", "all" = "all"), inline = TRUE),
                actionButton("plot_button", label = "Create Plot")),
            
            conditionalPanel(condition = "input.tabselected == 'common_panel' && input.common_sub_tab == 'common_venn'",
                selectizeInput("venn_id", label = "Select 2 or 3 repertoire ids", choices = NULL, multiple = TRUE),
                actionButton("venn_button", label = "Create Venn Diagram")),

            conditionalPanel(condition = "input.tabselected == 'gene_panel'",
                selectizeInput("locus", label = "Select which VDJ genes to include",
                                choices = c("VDJ", "DJ", "VJ", "DJ", "V", "D", "J")),
                radioButtons("family_bool", "Select which names to use",
                            choices = c("family" = TRUE, "gene" = FALSE), inline = TRUE)),
            
            conditionalPanel(condition = "input.tabselected == 'gene_panel' && input.gene_sub_tab == 'gene_freq_word'",
                selectizeInput("word_id", label = "Select repertoire id", choices = NULL),
                actionButton("word_button", label = "Create Word Cloud")),
            
            conditionalPanel(condition = "input.tabselected == 'clonal_relate'",
                radioButtons("merge_results", "Merge results with clonality table?",
                            choices = c("no", "yes"), inline = TRUE)),
            
            conditionalPanel(condition = "input.tabselected == 'diff_abundance'",
                selectizeInput("diff_id", label = "Select 2 repertoire ids to be compared",
                    choices = NULL, multiple = TRUE),
                numericInput("q_val", "q Value (False Discovery Rate):", value = 1, min = 0, max = 1),
                numericInput("zero_val", "Zero Value:", value = 1),
                actionButton("diff_button", "Create Table")),
            
            conditionalPanel(condition = "input.tabselected == 'count_kmers'",
                numericInput("k_val", "Length of kmer:", value = 2),
                actionButton("kmer_button", "Count Kmers")),

            conditionalPanel(condition = "input.tabselected == 'clone_track'",
                selectizeInput("track_id", label = "Select repertoire ids to track",
                    choices = NULL, multiple = TRUE),
                radioButtons("top_clones", "Plot top clones?",
                    choices = c("no", "yes"), inline = TRUE),
                numericInput("top_num", "Number of top clones to map:", value = NULL, min = 1),
                actionButton("track_button", "Create Plot")),
            
            conditionalPanel(condition = "input.tabselected == 'pairwise_sim'",
                radioButtons("mode", "Select mode to use for calculation:",
                    choices = c("Bhattacharyya", "Similarity", "Sorensen", "PSI"))),
            
            tags$br(),
            fluidRow (style = "border: 1px solid #d3d3d3; border-radius: 5px;",
                column(width = 6,
                    tags$br(),
                    radioButtons("download_type", "Download Type",
                        choices = c("PDF" = ".pdf", "TSV" = ".tsv")),
                    downloadButton("download"),
                    tags$br(),
                    tags$br()
                )
            )
        ),
        mainPanel(
            tabsetPanel(
                tabPanel("AIRR Data Table", value = "airr_table",
                        DT::dataTableOutput("table") %>% withSpinner()),
                tabPanel("Explore Common Sequences", value = "common_panel",
                            tabsetPanel(
                                tabPanel("Common Sequences Table", value = "common_seq_table",
                                    DT::dataTableOutput("common_seq_table") %>% withSpinner()),
                                tabPanel("Common Sequences Bar Chart", value = "common_bar",
                                    plotOutput("commonSeqs_bar") %>% withSpinner()),
                                tabPanel("Common Sequences Plot", value = "common_plot",
                                    plotlyOutput("commonSeqs_plot") %>% withSpinner()),
                                tabPanel("Common Sequences Venn Diagram", value = "common_venn",
                                    plotOutput("commonSeqs_venn") %>% withSpinner()),
                                id = "common_sub_tab"
                            )),
                tabPanel("Explore Productive Sequences", value = "prod_seq_panel",
                            tabsetPanel(
                                tabPanel("Unique Productive Sequences Plot", value = "produtive_seq_plot",
                                    plotlyOutput("productive_plot") %>% withSpinner()),
                                tabPanel("Top Productive Sequences Table", value = "top_10_seq_table",
                                    DT::dataTableOutput("top_productive_seq") %>% withSpinner()),
                                id = "prod_seq_sub_tab"
                            )),
                tabPanel("Explore Gene Frequencies", value = "gene_panel",
                            tabsetPanel(
                                tabPanel("Gene Frequencies Data Table", value = "gene_freq_table",
                                    DT::dataTableOutput("gene_freq_table") %>% withSpinner()),
                                tabPanel("Gene Frequencies Heat Map", value = "gene_freq_heat",
                                    plotlyOutput("gene_freq_heat") %>% withSpinner()),
                                tabPanel("Gene Frequencies Bar Chart", value = "gene_freq_bar",
                                    plotlyOutput("gene_freq_bar") %>% withSpinner()),
                                tabPanel("Word Cloud", value = "gene_freq_word",
                                    wordcloud2Output("gene_freq_word") %>% withSpinner()),
                                id = "gene_sub_tab"
                            )),
                tabPanel("Explore Clonality", value = "clonality_panel",
                            tabsetPanel(
                                tabPanel("Clonality", value = "clonality",
                                    DT::dataTableOutput("clonality") %>% withSpinner()),
                                tabPanel("Clonal Relatedness", value = "clonal_relate",
                                    DT::dataTableOutput("clonal_relate") %>% withSpinner()),
                                tabPanel("Clonality Statistics", value = "clonality_stats",
                                    DT::dataTableOutput("clonality_stats") %>% withSpinner()),
                                tabPanel("Clonality Plot", value = "clonality_plot",
                                    plotlyOutput("clonality_plot") %>% withSpinner()),
                                id = "clonal_sub_tab"
                            )),
                tabPanel("Chord Diagram VDJ", value = "chord_diagram",
                        chorddiagOutput("chord_diagram", width = '100%', height = '600px') %>% withSpinner()),
                tabPanel("Lorenz Curve", value = "lorenz_curve",
                        plotlyOutput("lorenz") %>% withSpinner()),
                tabPanel("Sequencing Counts", value = "seq_counts",
                        DT::dataTableOutput("seq_count") %>% withSpinner()),
                tabPanel("Count Statistics", value = "count_stats",
                        DT::dataTableOutput("stats_count") %>% withSpinner()),
                tabPanel("Top 10 Sequences Plot", value = "top_10_seq_plot",
                        plotlyOutput("top_10_seq") %>% withSpinner()),
                tabPanel("Clone Tracking", tags$div(
                        style = "position: relative;",
                        plotOutput("clone_track",
                        hover = hoverOpts(id = "clone_track_hover")) %>% withSpinner(),
                        htmlOutput("clone_track_tooltip")),
                        value = "clone_track"),
                tabPanel("Pairwise Similarity", value = "pairwise_sim",
                        plotlyOutput("pairwise_sim") %>% withSpinner()),
                tabPanel("Public TCRB Sequences", value = "public_tcrb_seq",
                        DT::dataTableOutput("public_tcrb") %>% withSpinner()),
                tabPanel("Differential Abundance", value = "diff_abundance",
                        DT::dataTableOutput("diff_abundance") %>% withSpinner()),
                tabPanel("Count Kmers", value = "count_kmers",
                        DT::dataTableOutput("count_kmers") %>% withSpinner()),
                id = "tabselected"
            )
        )
    )
)

'%then%' <- function(a, b) {
    if (is.null(a)) b else a
}

server <- function(input, output, session) {

    shinyjs::disable("download")
    shinyjs::disable("top_num")

    airr_data <- reactive({
        validate(
            need(!is.null(input$airr_files),
                    "Please select files to upload to render output.")
        )
        table <- LymphoSeqTest::readImmunoSeq(input$airr_files$datapath)
        # files <- c(input$airr_files$name)
        in_files <- lapply(c(input$airr_files$name), function(i) substr(i, 1, stringr::str_length(i) - 4))
        in_files <- unlist(in_files)
        table <- table %>%
                mutate(repertoire_id = if_else(as.integer(repertoire_id) <= length(in_files),
                        in_files[as.integer(repertoire_id) + 1], "unknown"))
        table
    })

    output$table <- DT::renderDataTable({
        airr_data()
    })

    productive_aa <- reactive({
        LymphoSeqTest::productiveSeq(study_table = airr_data(), aggregate = "junction_aa")
    })

    productive_nt <- reactive({
        LymphoSeqTest::productiveSeq(study_table = airr_data(), aggregate = "junction")
    })

    unique_prod_rep <- reactive({
        unique(productive_aa()[, "repertoire_id"])
    })

    observeEvent(input$airr_files, {
        if (input$chord_button || input$venn_button || input$bar_button ||
                input$plot_button || input$diff_button) {
            shinyjs::hide("chord_diagram")
            shinyjs::hide("commonSeqs_venn")
            shinyjs::hide("commonSeqs_bar")
            shinyjs::hide("commonSeqs_plot")
            shinyjs::hide("diff_abundance")
            shinyjs::hide("pairwise_sim")
            shinyjs::hide("gene_freq_word")
        }
        if (input$tabselected == "chord_diagram" || input$common_sub_tab == "common_venn" ||
                input$common_sub_tab == "common_bar" || input$common_sub_tab == "common_plot" ||
                input$tabselected == "diff_abundance" || input$tabselected == "pairwise_sim") {
            shinyjs::disable("download")
        } else {
            shinyjs::enable("download")
        }
        if (input$tabselected == "common_panel" || input$tabselected == "diff_abundance" ||
                input$tabselected == "clone_track") {
            shiny::updateSelectizeInput(session, "common_table_id", choices = unique_prod_rep())
            shiny::updateSelectizeInput(session, "bar_id", choices = unique_prod_rep())
            shiny::updateSelectizeInput(session, "color_rep_id", choices = c("none", unique_prod_rep()))
            shiny::updateSelectizeInput(session, "color_intersect", choices = c(""))
            shiny::updateSelectizeInput(session, "plot_id", choices = unique_prod_rep())
            shiny::updateSelectizeInput(session, "venn_id", choices = unique_prod_rep())
            shiny::updateSelectizeInput(session, "diff_id", choices = unique_prod_rep())
            shiny::updateSelectizeInput(session, "track_id", choices = unique_prod_rep())
        }
        if (input$tabselected == "gene_panel") {
            nt_rep_id <- unique(productive_nt()[, "repertoire_id"])
            shiny::updateSelectizeInput(session, "word_id", choices = nt_rep_id)
        }
    })

    observeEvent(input$bar_id, {
        shiny::updateSelectizeInput(session, "color_intersect", choices = input$bar_id)
    })

    observeEvent(input$tabselected, {
        data_frame_tabs = c("airr_table", "seq_counts", "count_stats", "clonality_panel",
                            "prod_seq_panel", "common_panel", "public_tcrb_seq",
                            "gene_panel", "diff_abundance", "count_kmers")
        if (input$tabselected %in% data_frame_tabs) {
            download_choices <- c("TSV" = ".tsv", "EXCEL" = ".xlsx", "Rda" = ".RData")
        } else {
            download_choices <- c("PDF" = ".pdf")
        }
        shiny::updateRadioButtons(session, "download_type", choices = download_choices)
        if (input$tabselected == "diff_abundance") {
            shiny::updateSelectizeInput(session, "diff_id", choices = unique_prod_rep())
        } else if (input$tabselected == "clone_track") {
            shiny::updateSelectizeInput(session, "track_id", choices = unique_prod_rep())
        }
    })

    update_common_tabs <- reactive({
        shiny::updateSelectizeInput(session, "common_table_id", choices = unique_prod_rep())
        shiny::updateSelectizeInput(session, "bar_id", choices = unique_prod_rep())
        shiny::updateSelectizeInput(session, "color_rep_id", choices = c("none", unique_prod_rep()))
        shiny::updateSelectizeInput(session, "color_intersect", choices = c(""))
        shiny::updateSelectizeInput(session, "plot_id", choices = unique_prod_rep())
        shiny::updateSelectizeInput(session, "venn_id", choices = unique_prod_rep())
    })

    observeEvent(input$common_sub_tab, {
        if (input$common_sub_tab == "common_seq_table") {
            shiny::updateRadioButtons(session, "download_type",
                choices = c("TSV" = ".tsv", "EXCEL" = ".xlsx", "Rda" = ".RData"))
        } else {
            shiny::updateRadioButtons(session, "download_type", choices = c("PDF" = ".pdf"))
        }
        update_common_tabs()
    })

    update_gene_tabs <- reactive({
        nt_rep_id <- unique(productive_nt()[, "repertoire_id"])
        shiny::updateSelectizeInput(session, "word_id", choices = nt_rep_id)
    })

    observeEvent(input$gene_sub_tab, {
        if (input$gene_sub_tab == "gene_freq_table") {
            shiny::updateRadioButtons(session, "download_type",
                choices = c("TSV" = ".tsv", "EXCEL" = ".xlsx", "Rda" = ".RData"))
        } else {
            shiny::updateRadioButtons(session, "download_type", choices = c("PDF" = ".pdf"))
        }
        update_gene_tabs()
    })

    observeEvent(input$clonal_sub_tab, {
        if (input$clonal_sub_tab == "clonality_plot") {
            shiny::updateRadioButtons(session, "download_type", choices = c("PDF" = ".pdf"))
        } else {
            shiny::updateRadioButtons(session, "download_type",
                choices = c("TSV" = ".tsv", "EXCEL" = ".xlsx", "Rda" = ".RData"))
        }
    })

    observeEvent(input$prod_seq_sub_tab, {
        if (input$prod_seq_sub_tab == "top_10_seq_table") {
            shiny::updateRadioButtons(session, "download_type",
                choices = c("TSV" = ".tsv", "EXCEL" = ".xlsx", "Rda" = ".RData"))
        }
    })

    observeEvent(input$chord_button, {
        shinyjs::enable("download")
        shinyjs::show("chord_diagram")
    })

    chord_data <- reactive({
        top_seqs <- LymphoSeqTest::topSeqs(productive_table = productive_nt(), top = 1)
        validate(
            need(!is.null(input$vdj_association), "Please select VDJ Association")
        )
        if (input$vdj_association == 'VJ') {
            vj <- top_seqs %>% 
                select(v_family, j_family) %>% 
                mutate(v_family = replace_na(v_family, "Unresolved"), j_family = replace_na(j_family, "Unresolved"))

            vj <- vj %>% 
                group_by(v_family, j_family) %>% 
                summarize(duplicate_count = n(), .groups = 'drop') %>% 
                pivot_wider(id_cols=v_family, names_from = j_family, values_from = duplicate_count) 
            row_names <- vj$v_family 
            vj <- vj %>% dplyr::select(-v_family)
            vj[is.na(vj)] <- 0
            vj <- as.matrix(vj)
            rownames(vj) <- row_names
            vj

        } else if (input$vdj_association == 'DJ') {
            dj <- top_seqs %>% 
                select(d_family, j_family) %>% 
                mutate(d_family = replace_na(d_family, "Unresolved"), j_family = replace_na(j_family, "Unresolved"))
            dj <- dj %>% 
                group_by(d_family, j_family) %>% 
                summarize(duplicate_count = n(), .groups = 'drop') %>% 
                pivot_wider(id_cols=d_family, names_from = j_family, values_from = duplicate_count)
            row_names <- dj$d_family
            dj <- dj %>% dplyr::select(-d_family)
            dj[is.na(dj)] <- 0
            dj <- as.matrix(dj)
            rownames(dj) <- row_names
            dj
        }
    })

    output$chord_diagram <- renderChorddiag({
        input$chord_button
        isolate ({
            chorddiag::chorddiag(chord_data(), type = "bipartite", groupnameFontsize = 15)
        })
    })

    observeEvent(input$bar_button, {
        shinyjs::enable("download")
        shinyjs::show("commonSeqs_bar")
    })

    common_bar_data <- reactive({
        color_rep_id <- input$color_rep_id
        color_intersect <- input$color_intersect
        if (input$color_rep_id == "none") {
            color_rep_id <- NULL
        }
        if (!is.null(input$color_intersect) & "none" %in% input$color_intersect) {
            color_intersect <- NULL
        }
        data_output <- LymphoSeqTest::commonSeqsBar(productive_aa(),
                    input$bar_id, color_rep_id, color_intersect,
                    labels = "yes")
        data_output
    })

    output$commonSeqs_bar <- renderPlot({
        input$bar_button
        isolate({
            validate(
                need(length(input$bar_id) > 1,
                        "Please select at least 2 repertoire ids")
                %then%
                need(length(input$color_rep_id) == 1,
                    "Please select only 1 repertoire id")
                %then%
                need(is.null(input$color_intersect) || input$color_intersect == "none" 
                    || length(input$color_intersect) > 1,
                    "Please select at least 2 repertoire ids")
            )
            tryCatch({
                common_bar_data()
                },
                error = function(e) {
                    shinyjs::disable("download")
                    shiny::showNotification("Cannot render plot", "", type = "error")
                    return()
                }
            )
        })
    })

    observeEvent(input$plot_button, {
        shinyjs::enable("download")
        shinyjs::show("commonSeqs_plot")
    })

    common_seqs_plot <- reactive({
        print("common")
        LymphoSeqTest::commonSeqsPlot(input$plot_id[1], input$plot_id[2], productive_aa())
    })

    output$commonSeqs_plot <- renderPlotly({
        input$plot_button
        isolate({
            validate(
                need(length(input$plot_id) == 2,
                    "Please select exactly 2 repertoire ids")
            )
            if (nrow(common_seqs_plot()$data) > 0) {
                common_seqs_plot()
            } else {
                shinyjs::disable("download")
                shiny::showNotification("No common sequences", "", type = "error")
                return()
            }
        })
    })

    observeEvent(input$venn_button, {
        shinyjs::show("commonSeqs_venn")
        shinyjs::enable("download")
    })

    common_venn_data <- reactive({
        if (length(input$venn_id) == 2) {
            a <- productive_aa() %>% 
                dplyr::filter(repertoire_id == input$venn_id[[1]])
            b <- productive_aa() %>% 
                dplyr::filter(repertoire_id == input$venn_id[[2]])
            grid::grid.newpage()
            venn <- VennDiagram::draw.pairwise.venn(area1 = length(a$junction_aa), 
                                                    area2 = length(b$junction_aa), 
                                                    cross.area = length(intersect(a$junction_aa, 
                                                                                b$junction_aa)), 
                                                    category = c(input$venn_id[1], 
                                                                input$venn_id[2]), 
                                                    cat.fontfamily = rep("sans", 2), 
                                                    fontfamily = rep("sans", 3), 
                                                    fill = c("#3288bd", "#d53e4f"), 
                                                    cat.pos = c(0, 0),
                                                    cat.dist = rep(0.025, 2),
                                                    cex = 1, 
                                                    cat.cex = 0.7,
                                                    lwd = rep(2, 2))
            data_output <- venn
        }
        if (length(input$venn_id) == 3) {
            a <- productive_aa() %>% 
                dplyr::filter(repertoire_id == input$venn_id[[1]])
            b <- productive_aa() %>% 
                dplyr::filter(repertoire_id == input$venn_id[[2]])
            c <- productive_aa() %>% 
                dplyr::filter(repertoire_id == input$venn_id[[3]])
            grid::grid.newpage()
            venn <- VennDiagram::draw.triple.venn(area1 = length(a$junction_aa), 
                                                area2 = length(b$junction_aa), 
                                                area3 = length(c$junction_aa), 
                                                n12 = length(intersect(a$junction_aa, 
                                                                        b$junction_aa)), 
                                                n23 = length(intersect(b$junction_aa, 
                                                                        c$junction_aa)), 
                                                n13 = length(intersect(a$junction_aa, 
                                                                        c$junction_aa)), 
                                                n123 = length(Reduce(intersect, 
                                                                    list(a$junction_aa, 
                                                                            b$junction_aa, 
                                                                            c$junction_aa))), 
                                                category = c(input$venn_id[1], 
                                                            input$venn_id[2], 
                                                            input$venn_id[3]), 
                                                cat.fontfamily = rep("sans", 3), 
                                                fontfamily = rep("sans", 7), 
                                                fill = c("#3288bd", "#abdda4", "#d53e4f"), 
                                                cat.pos = c(0, 0, 180), 
                                                cat.dist = rep(0.025, 3),
                                                cex = 1, 
                                                cat.cex = 0.7,
                                                lwd = rep(2, 3))
            data_output <- venn
        }
        data_output
    })

    output$commonSeqs_venn <- renderPlot({
        input$venn_button
        isolate({
            validate(
                need(length(input$venn_id) == 2 | length(input$venn_id) == 3,
                    "Please select only 2 or 3 repertoire ids")
            )
            common_venn_data()
        })
    })

    lorenz_data <- reactive({
        repertoire_ids <- productive_aa() %>% dplyr::pull(repertoire_id) %>% unique()
        data_output <- LymphoSeqTest::lorenzCurve(repertoire_ids, airr_data()) + coord_fixed(1/2)
        data_output
    })

    output$lorenz <- renderPlotly({
        lorenz_data()
    })

    itable <- reactive({
        LymphoSeqTest::clonality(airr_data())
    })

    seq_count_data <- reactive({
        data_output <- itable() %>%
            select(repertoire_id, total_sequences,
                    unique_productive_sequences, total_count)
        data_output
    })

    output$seq_count <- DT::renderDataTable({
        seq_count_data() %>%
            datatable(rownames = FALSE, colnames = c("Repertoire ID",
                        "Total Sequences", "Unique Productive Sequences",
                        "Total Count"), filter = "top")
    })

    stats_count_data <- reactive({
        stats <- itable() %>%
            select(repertoire_id, total_sequences,
                    unique_productive_sequences, total_count) %>%
            pivot_longer(!repertoire_id,
                        names_to = "metric",
                        values_to = "values") %>%
            group_by(metric) %>%
            summarize(Total = sum(values),
                    Mean = mean(values),
                    Minimum = min(values),
                    Maximum = max(values)) %>%
            drop_na() %>%
            ungroup() %>%
            mutate(metric = str_to_title(str_replace_all(metric, "_", " ")))
        stats
    })

    output$stats_count <- DT::renderDataTable({
        datatable(stats_count_data(), rownames = FALSE, colnames = c("Metric", "Total",
                    "Mean", "Minimum", "Maximum"), filter = "top")
    })

    prod_plot_data <- reactive({
        data_output <- ggplot(data = itable(), aes(x = repertoire_id, y = unique_productive_sequences)) +
            geom_bar(stat = "identity", position=position_dodge(), fill = "#3182bd") +
            theme_minimal() +
            theme(axis.text.x  = element_text(angle=90, vjust=0.5, hjust = 1), text = element_text(size = 10)) +
            scale_x_discrete(limits = itable()$repertoire_id) +
            labs(x = "", y = "Unique productive sequences")
        data_output
    })

    output$productive_plot <- renderPlotly({
        plotly::ggplotly(prod_plot_data())
    })

    top_prod_seq_data <- reactive({
        data_output <- LymphoSeqTest::topSeqs(productive_aa()) %>%
            select(repertoire_id, junction_aa, duplicate_count, duplicate_frequency,
                        v_call, d_call, j_call)
        data_output
    })

    output$top_productive_seq <- DT::renderDataTable({
        top_prod_seq_data() %>%
            datatable(rownames = FALSE, filter = "top")
    })

    top_seq_data <- reactive({
        data_output <- topSeqsPlot(productive_aa(), top = 10) +
            theme(text = element_text(size = 10)) +
            scale_x_discrete(limits = itable()$repertoire_id)
        data_output
    })

    output$top_10_seq <- renderPlotly({
        plotly::ggplotly(top_seq_data())
    })

    common_table_data <- reactive({
        LymphoSeqTest::commonSeqs(study_table = productive_aa())
    })

    output$common_seq_table <- DT::renderDataTable({
        common_table_data()
    })

    clone_table_data <- reactive({
        data_output <- itable() %>%
            select(repertoire_id, clonality, gini_coefficient, top_productive_sequence)
    })

    output$clonality <- DT::renderDataTable({
        clone_table_data() %>%
            datatable(rownames = FALSE, colnames = c("Repertoire ID", "Clonality", "Gini Coefficient", "Top Productive Sequence"), filter = "top")
    })

    clone_relate_data <- reactive({
        data_output <- LymphoSeqTest::clonalRelatedness(airr_data())
        if (input$merge_results == "yes") {
            clonality <- LymphoSeqTest::clonality(airr_data())
            data_output <- merge(clonality, data_output)
        }
        data_output
    })

    output$clonal_relate <- DT::renderDataTable({
        clone_relate_data()
    })

    clone_stats_data <- reactive({
        data_output <- itable() %>%
            select(repertoire_id, clonality, gini_coefficient, top_productive_sequence) %>% 
            pivot_longer(!repertoire_id,
                        names_to = "metric",
                        values_to = "values") %>%
            group_by(metric) %>%
            summarize(Mean = mean(values),
                    Minimum = min(values),
                    Maximum = max(values)) %>%
            drop_na() %>%
            ungroup() %>%
            mutate(metric = str_to_title(str_replace_all(metric, "_", " ")))
        data_output
    })

    output$clonality_stats <- DT::renderDataTable({
        clone_stats_data() %>%
            datatable(rownames = FALSE, colnames = c("Metric", "Mean", "Minimum", "Maximum"), filter = "top")
    })

    clone_plot_data <- reactive({
        data_output <- ggplot(data = itable(), aes(x = repertoire_id, y = clonality)) +
            geom_bar(stat = "identity", position=position_dodge(), fill = "#b2182b") +
            theme_minimal() +
            theme(axis.text.x  = element_text(angle=90, vjust=0.5, hjust = 1), text = element_text(size = 10)) +
            scale_x_discrete(limits = itable()$repertoire_id) +
            labs(x = "", y = "Clonality")
        data_output
    })

    output$clonality_plot <- renderPlotly({
        plotly::ggplotly(clone_plot_data())
    })

    interactive_alluvial <- function(p) {
        node_width <- 1/3
        alluvium_width <- 1/3
        pbuilt <<- ggplot_build(p)

        data_draw <- transform(pbuilt$data[[3]], width = 1/3)
        groups_to_draw <<- split(data_draw, data_draw$group)
        group_xsplines <- lapply(groups_to_draw,
                                data_to_alluvium)

        xspline_coords <- lapply(
            group_xsplines,
            function(coords) grid::xsplineGrob(x = coords$x,
                                            y = coords$y,
                                            shape = coords$shape,
                                            open = FALSE)
        )

        xspline_points <- lapply(xspline_coords, grid::xsplinePoints)

        xrange_old <- range(unlist(lapply(
            xspline_points,
            function(pts) as.numeric(pts$x)
        )))
        yrange_old <- range(unlist(lapply(
            xspline_points,
            function(pts) as.numeric(pts$y)
        )))
        xrange_new <- c(1 - alluvium_width/2, max(pbuilt$data[[3]]$x)
                        + alluvium_width/2)
        yrange_new <- c(0, sum(pbuilt$data[[3]]$count[pbuilt$data[[3]]$x == 1]))

        new_range_transform <- function(x_old, range_old, range_new) {
            (x_old - range_old[1])/(range_old[2] - range_old[1]) *
            (range_new[2] - range_new[1]) + range_new[1]
        }

        polygon_coords <<- lapply(xspline_points, function(pts) {
            x_trans <- new_range_transform(x_old = as.numeric(pts$x),
                                        range_old = xrange_old,
                                        range_new = xrange_new)
            y_trans <- new_range_transform(x_old = as.numeric(pts$y),
                                        range_old = yrange_old,
                                        range_new = yrange_new)
            list(x = x_trans, y = y_trans)
        })
    }

    alluvial_tooltip <- function(plot_hover) {
        node_width <- 1/3
        if (!is.null(plot_hover)) {
            hover <- plot_hover
            x_coord <- round(hover$x)

            if (abs(hover$x - x_coord) < (node_width / 2)) {
                node_row <- pbuilt$data[[3]]$x == x_coord &
                            hover$y > pbuilt$data[[3]]$ymin &
                            hover$y < pbuilt$data[[3]]$ymax
                node_label <- pbuilt$data[[3]]$stratum[node_row]
                node_n <- pbuilt$data[[3]]$n[node_row]
                offset <- 5

                renderTags(
                    tags$div(
                        node_label, #tags$br(),
                        #"n =", node_n,
                        style = paste0(
                        "position: absolute; ",
                        "top: ", hover$coords_css$y + offset, "px; ",
                        "left: ", hover$coords_css$x + offset, "px; ",
                        "background: gray; ",
                        "padding: 3px; ",
                        "color: white; "
                        )
                    )
                )$html
            } else {
                hover_within_flow <- sapply(
                    polygon_coords,
                    function(pol) point.in.polygon(point.x = hover$x,
                                                    point.y = hover$y,
                                                    pol.x = pol$x,
                                                    pol.y = pol$y)
                )
                if (any(hover_within_flow)) {
                    coord_id <- rev(which(hover_within_flow == 1))[1]
                    flow_label <- paste(groups_to_draw[[coord_id]]$stratum, collapse = ' -> ')
                    flow_n <- groups_to_draw[[coord_id]]$count[1]
                    offset <- 5
                    
                    renderTags(
                        tags$div(
                            flow_label, #tags$br(),
                            #"n =", flow_n,
                            style = paste0(
                                "position: absolute; ",
                                "top: ", hover$coords_css$y + offset, "px; ",
                                "left: ", hover$coords_css$x + offset, "px; ",
                                "background: gray; ",
                                "padding: 3px; ",
                                "color: white; "
                            )
                        )
                    )$html
                }
            }
        }
    }

    observeEvent(input$top_clones, {
        if (input$top_clones == "no") {
            shinyjs::disable("top_num")
        } else {
            shinyjs::enable("top_num")
        }
    })

    clone_track_data <- reactive({
        if (input$top_clones == "yes") {
            validate(
                need(input$top_num > 0, "Please enter a number for tracking top clones")
            )
            ttable <- LymphoSeqTest::topSeqs(productive_aa(), top = input$top_num)
            ctable <- LymphoSeqTest::cloneTrack(ttable, sample_list = input$track_id)

        } else {
            ctable <- LymphoSeqTest::cloneTrack(productive_aa(), sample_list = input$track_id)
        }
        # max_seen <- ctable %>%
        #     pull(seen) %>%
        #     max()
        # max_amino <- ctable %>%
        #      filter(seen == max_seen) %>%
        #      pull(junction_aa) %>%
        #      unique()
        # top_seen <- ctable %>%
        #     slice_max(seen, n = 100)
        # p <- LymphoSeqTest::plotTrack(clone_table = top_seen, alist = max_amino)
        LymphoSeqTest::plotTrack(clone_table = ctable)
    })

    output$clone_track <- renderPlot({
        input$track_button
        isolate({
            validate(
                need(length(input$track_id) > 1, "Please select repertoire ids to track")
            )
            interactive_alluvial(clone_track_data())
            clone_track_data() + theme(legend.position = "none")
        })
    })

    output$clone_track_tooltip <- renderText(
        alluvial_tooltip(input$clone_track_hover)
    )

    pairwise_sim_data <- reactive({
        similarity <- LymphoSeqTest::scoringMatrix(productive_aa(), mode = input$mode)
        data_output <- pairwisePlot(matrix = similarity) +
            labs(fill = paste(input$mode, "\ncoefficient")) +
            theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8)) +
            scale_fill_gradient(low = "#deebf7", high = "#3182bd") +
            theme(legend.position = c(0.8, 0.8))
        data_output
    })

    output$pairwise_sim <- renderPlotly({
        plotly::ggplotly(pairwise_sim_data())
    })

    public_table_data <- reactive({
        published <- LymphoSeqTest::searchPublished(productive_aa())
        data_output <- published %>%
            filter(!is.na(PMID)) %>%
            select(repertoire_id, junction_aa, duplicate_count, antigen, prevalence)
        data_output
    })

    output$public_tcrb <- DT::renderDataTable({
        public_table_data() %>%
                datatable(colnames = c("Repertoire ID", "Junction AA",
                            "Duplicate Count", "Antigen", "Prevalence"),
                            filter = "top")
    })

    gene_heatmap_data <- reactive({
        genes <- LymphoSeqTest::geneFreq(productive_nt(), locus = input$locus, family = input$family_bool) %>% 
            pivot_wider(id_cols = gene_name,
                        names_from = repertoire_id,
                        values_from = gene_frequency,
                        values_fn = sum,
                        values_fill = 0)
        gene_names <- genes %>%
            pull(gene_name)
        genes <- genes %>%
                select(-gene_name) %>%
                as.matrix()
        rownames(genes) <- gene_names
        genes
    })

    output$gene_freq_heat <- renderPlotly({
        heatmaply::heatmaply(gene_heatmap_data(), colors = "RdBu", scale = "row")
    })

    gene_bar_data <- reactive({
        genes <- LymphoSeqTest::geneFreq(productive_nt(), input$locus, input$family_bool)
        data_output <<- ggplot(genes, aes(x = repertoire_id, y = gene_frequency, fill = gene_name)) +
            geom_bar(stat = "identity") +
            theme_minimal() + 
            scale_y_continuous(expand = c(0, 0)) + 
            guides(fill = guide_legend(ncol = 3)) +
            labs(y = "Frequency (%)", x = "", fill = "") +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
        data_output
    })
    
    output$gene_freq_bar <- renderPlotly({
        plotly::ggplotly(gene_bar_data())
    })

    observeEvent(input$word_button, {
        shinyjs::enable("download")
        shinyjs::show("gene_freq_word")
    })

    gene_word_data <- reactive({
        genes <- LymphoSeqTest::geneFreq(productive_nt(), input$locus, input$family_bool)
        gene_data <- genes %>%
                        filter(repertoire_id == input$word_id)
        gene_data <- data.frame(gene_data$gene_name, gene_data$gene_frequency)
        data_output <- wordcloud2::wordcloud2(gene_data, color = 'random-dark', minSize = 10)
        data_output
    })

    output$gene_freq_word <- renderWordcloud2({
        input$word_button
        isolate ({
            gene_word_data()
        })
    })

    gene_table_data <- reactive({
        data_output <- LymphoSeqTest::geneFreq(productive_nt(), input$locus, input$family_bool)
        data_output
    })

    output$gene_freq_table <- DT::renderDataTable({
        gene_table_data()
    })

    observeEvent(input$diff_button, {
        shinyjs::enable("download")
        shinyjs::show("diff_abundance")
    })

    diff_table_data <- reactive({
        LymphoSeqTest::differentialAbundance(productive_aa(),
                        input$diff_id, q = input$q_val, zero = input$zero_val)
    })

    output$diff_abundance <- DT::renderDataTable({
        input$diff_button
        isolate({
            validate(
                need(length(input$diff_id) == 2, "Please select 2 repertoire ids")
            )
            diff_table_data()
        })
    })

    kmer_table_data <- reactive({
        LymphoSeqTest::countKmer(airr_data(), input$k_val)
    })

    output$count_kmers <- DT::renderDataTable({
        input$kmer_button
        isolate({
            validate(
                need(input$k_val > 2, "Please input length of at least 2")
            )
            kmer_table_data()
        })
    })

    output$download <- downloadHandler(
        filename <- function() {
            paste0(input$tabselected, input$download_type)
        },
        content <- function(file) {
            if (input$download_type == ".pdf") {
                pdf(file)
                if (input$tabselected == "gene_panel") {
                    if (input$gene_sub_tab == "gene_freq_heat") {
                        RedBlue <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(256)
                        data_output <- list(gene_heatmap_data(), RedBlue)
                        pheatmap(data_output[[1]], color = data_output[[2]], scale = "row")
                    } else if (input$gene_sub_tab == "gene_freq_bar") {
                        plot(gene_bar_data())
                    }
                } else if (input$tabselected == "chord_diagram") {
                    circlize::chordDiagram(chord_data(), annotationTrack = c("grid", "name"))
                } else if (input$tabselected == "common_panel") {
                    if (input$common_sub_tab == "common_bar") {
                        print(common_bar_data())
                    } else if (input$common_sub_tab == "common_venn") {
                        grid::grid.newpage()
                        grid::grid.draw(common_venn_data())
                    } else if (input$common_sub_tab == "common_plot") {
                        plot(common_seqs_plot())
                    }
                } else if (input$tabselected == "clonality_plot") {
                    plot(clone_plot_data())
                } else if (input$tabselected == "lorenz_curve") {
                    plot(lorenz_data())
                } else if (input$tabselected == "top_10_seq_plot") {
                    plot(top_seq_data())
                } else if (input$tabselected == "pairwise_sim") {
                    plot(pairwise_sim_data())
                } else if (input$tabselected == "prod_seq_panel" && input$prod_seq_sub_tab == "produtive_seq_plot") {
                    plot(prod_plot_data())
                } else if (input$tabselected == "clone_track") {
                    plot(clone_track_data())
                }
                dev.off()
            } else {
                if (input$tabselected == "airr_table") {
                    data_output <- airr_data()
                } else if (input$tabselected == "common_panel" && input$common_sub_tab == "common_seq_table") {
                    data_output <- common_table_data()
                } else if (input$tabselected == "prod_seq_panel" && input$prod_seq_sub_tab == "top_10_seq_table") {
                    data_output <- top_prod_seq_data()
                } else if (input$tabselected == "gene_panel" && input$gene_sub_tab == "gene_freq_table") {
                    data_output <- gene_table_data()
                } else if (input$tabselected == "clonality_panel") {
                    if (input$clonal_sub_tab == "clonality") {
                        data_output <- clone_table_data()
                    } else if (input$clonal_sub_tab == "clonal_relate") {
                        data_output <- clone_relate_data()
                    } else if (input$clonal_sub_tab == "clonality_stats") {
                        data_output <- clone_stats_data()
                    }
                } else if (input$tabselected == "seq_counts") {
                    data_output <- seq_count_data()
                } else if (input$tabselected == "count_stats") {
                    data_output <- stats_count_data()
                } else if (input$tabselected == "public_tcrb_seq") {
                    data_output <- public_table_data()
                } else if (input$tabselected == "count_kmers") {
                    data_output <- kmer_table_data()
                } else if (input$tabselected == "diff_abundance") {
                    data_output <- diff_table_data()
                }

                if (input$download_type == ".tsv") {
                    write.table(data_output, file, quote = FALSE, sep='\t', row.names = FALSE)
                } else if (input$download_type == ".xlsx") {
                    write_xlsx(data_output, path = file)
                } else if (input$download_type == ".RData") {
                    save(data_output, file = file)
                }
            }
        }
    )

}

shinyApp(ui = ui, server = server)
