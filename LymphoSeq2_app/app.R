library(shiny)
library(shinycssloaders)
library(LymphoSeq2)
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

library(shinyalert)
library(shinyscreenshot)
library(capture)

options(shiny.maxRequestSize = 500 * 1024^2)


ui <- 
navbarPage("LymphoSeq2 Application", theme = shinythemes::shinytheme("cerulean"),
    tabPanel("Explore Data",
        fluidPage(
            shinyjs::useShinyjs(),
            # tags$head(
            #     tags$style(
            #         HTML(".shiny-notification {
            #             position:fixed;
            #             top: calc(50%);
            #             left: calc(50%);
            #             }")
            #     )
            # ),
    sidebarLayout(
        sidebarPanel(
            fileInput("airr_files", label = "Upload Files", multiple = TRUE, accept = c(".tsv", ".rda", ".RData")),

            conditionalPanel(condition = "input.tabselected == 'chord_diagram'",
                selectizeInput("vdj_association", label = "Select VDJ Association",
                                choices = c("", "VJ", "DJ")),
                radioButtons("top_seq_chord", "Plot top sequences? (Recommended)",
                    choices = c("yes", "no"), inline = TRUE),
                numericInput("top_chord_num", "Number of top sequences to map:", value = 1, min = 1),
                actionButton("chord_button", label = "Create Chord Diagram")),

            conditionalPanel(condition = "input.tabselected == 'common_panel' && input.common_sub_tab == 'common_bar'",
                selectizeInput("bar_id", label = "Select at least 2 repertoire ids", choices = NULL, multiple = TRUE),
                h1("Options to color bar chart:", style = "font-size:15px"),
                selectizeInput("color_rep_id", label = "Select 1 repertoire id to color", choices = NULL, selected = NULL),
                h5("OR", align = "center"),
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
            
            conditionalPanel(condition = "input.tabselected == 'clonality_panel' && input.clonal_sub_tab == 'clonal_relate'",
                numericInput("edit_dis", label = "Minimum edit distance the sequence must be less than or equal to:",
                                value = 10, min = 0),
                radioButtons("merge_results", "Merge results with clonality table?",
                            choices = c("no", "yes"), inline = TRUE),
                actionButton("clonal_relate_button", "Create Table")),
            
            conditionalPanel(condition = "input.tabselected == 'diff_abundance'",
                selectizeInput("diff_id", label = "Select 2 repertoire ids to be compared",
                    choices = c(""), multiple = TRUE),
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
                    choices = c("yes", "no"), inline = TRUE),
                numericInput("top_num", "Number of top clones to map:", value = 50, min = 1),
                actionButton("track_button", "Create Plot")),
            
            conditionalPanel(condition = "input.tabselected == 'pairwise_sim'",
                radioButtons("mode", "Select mode to use for calculation:",
                    choices = c("Bhattacharyya", "Similarity", "Sorensen", "PSI"))),

            conditionalPanel(condition = "input.tabselected == 'prod_seq_panel' && 
                                (input.prod_seq_sub_tab == 'top_seq_table' || input.prod_seq_sub_tab == 'top_seq_plot')",
                numericInput("top_seq_num", "Number of top sequences:", value = NULL, min = 0),
                actionButton("top_seq_button", "Calculate Top Sequences")),
            
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
                tabPanel("Explore Common Productive Sequences", value = "common_panel",
                            tabsetPanel(
                                tabPanel("Data Table", value = "common_seq_table",
                                    DT::dataTableOutput("common_seq_table") %>% withSpinner()),
                                tabPanel("Bar Chart", value = "common_bar",
                                    plotOutput("commonSeqs_bar") %>% withSpinner()),
                                tabPanel("Plot", value = "common_plot",
                                    plotlyOutput("commonSeqs_plot") %>% withSpinner()),
                                tabPanel("Venn Diagram", value = "common_venn",
                                    plotOutput("commonSeqs_venn") %>% withSpinner()),
                                id = "common_sub_tab"
                            )),
                tabPanel("Explore Productive Sequences", value = "prod_seq_panel",
                            tabsetPanel(
                                tabPanel("Top Productive Sequences Table", value = "top_seq_table",
                                    DT::dataTableOutput("top_seq_table") %>% withSpinner()),
                                tabPanel("Top Productive Sequences Plot", value = "top_seq_plot",
                                    plotlyOutput("top_seq_plot") %>% withSpinner()),                                    
                                tabPanel("Unique Productive Sequences Plot", value = "produtive_seq_plot",
                                    plotlyOutput("productive_plot") %>% withSpinner()),
                                tabPanel("Unique Productive Sequences Table", value = "produtive_seq_table",
                                    DT::dataTableOutput("produtive_seq_table") %>% withSpinner()),
                                id = "prod_seq_sub_tab"
                            )),
                tabPanel("Explore Gene Frequencies", value = "gene_panel",
                            tabsetPanel(
                                tabPanel("Data Table", value = "gene_freq_table",
                                    DT::dataTableOutput("gene_freq_table") %>% withSpinner()),
                                tabPanel("Heat Map", value = "gene_freq_heat",
                                    plotlyOutput("gene_freq_heat") %>% withSpinner()),
                                tabPanel("Bar Chart", value = "gene_freq_bar",
                                    plotlyOutput("gene_freq_bar") %>% withSpinner()),
                                # tabPanel("Word Cloud", value = "gene_freq_word",
                                #     wordcloud2Output("gene_freq_word") %>% withSpinner()),
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
                                tabPanel("Sequencing Counts", value = "seq_counts",
                                    DT::dataTableOutput("seq_count") %>% withSpinner()),
                                tabPanel("Count Statistics", value = "count_stats",
                                    DT::dataTableOutput("stats_count") %>% withSpinner()),               
                                id = "clonal_sub_tab"
                            )),
                tabPanel("Chord Diagram VDJ", value = "chord_diagram",
                        chorddiagOutput("chord_diagram", width = '100%', height = '600px') %>% withSpinner()),
                tabPanel("Lorenz Curve", value = "lorenz_curve",
                        plotlyOutput("lorenz") %>% withSpinner()),
                tabPanel("Clone Tracking", value = "clone_track", tags$div(
                        style = "position: relative;",
                        plotOutput("clone_track",
                        hover = hoverOpts(id = "clone_track_hover")) %>% withSpinner(),
                        htmlOutput("clone_track_tooltip"))),
                tabPanel("Pairwise Similarity", value = "pairwise_sim",
                        plotlyOutput("pairwise_sim") %>% withSpinner()),
                tabPanel("Public TCRB Sequences", value = "public_tcrb_seq",
                        DT::dataTableOutput("public_tcrb") %>% withSpinner()),
                tabPanel("Differential Abundance", value = "diff_abundance",
                        DT::dataTableOutput("diff_abundance") %>% withSpinner()),
                tabPanel("Kmer Counts", value = "count_kmers",
                        DT::dataTableOutput("count_kmers") %>% withSpinner()),
                id = "tabselected"
            )
        )
    )
)
),
    tabPanel("About",
        fluidPage(
            tags$h2("More about the LymphoSeq2 package here:", style = "font-size:30px"),
            tags$br(),
            tags$a(href = "https://shashidhar22.github.io/LymphoSeq2/", "LymphoSeq2 Package Website", style = "font-size:20px"),
            tags$br(), tags$br(),
            tags$a(href = "https://github.com/shashidhar22/LymphoSeq2", "LymphoSeq2 on Github", style = "font-size:20px")
        ))
)

'%then%' <- function(a, b) {
    if (is.null(a)) b else a
}

server <- function(input, output, session) {

    shinyjs::disable("top_num")
    shinyjs::disable("top_chord_num")
    rda_envir <<- NULL

    airr_data <- reactive({
        validate(
            need(!is.null(input$airr_files),
                    "Please select files to upload to render output.")
        )
        tryCatch({
            if (tools::file_ext(input$airr_files$name) == "tsv") {
                rda_envir <<- NULL
                table <- LymphoSeq2::readImmunoSeq(input$airr_files$datapath)
                in_files <- lapply(c(input$airr_files$name), function(i) substr(i, 1, stringr::str_length(i) - 4))
                in_files <- unlist(in_files)
                table <- table %>%
                        mutate(repertoire_id = if_else(as.integer(repertoire_id) <= length(in_files),
                                in_files[as.integer(repertoire_id) + 1], "unknown"))
                table
            } else if (tools::file_ext(input$airr_files$name) == "rda") {
                rda_envir <<- new.env()
                name <- load(input$airr_files$datapath, envir = rda_envir)
                rda_envir$study_table
            }
        },
            error = function(e) {
                shinyalert::shinyalert("No Common Sequences!",
                    "Select different repertoire ids.", type = "error")
                # shiny::showNotification("Cannot convert data. Please choose other files.", "", type = "error")
                return()
            }
        )
    })

    output$table <- DT::renderDataTable({
        airr_table <- airr_data() %>%
                        purrr::discard(~all(is.na(.) | . == "")) %>%
                        DT::datatable(filter = "top", options = list(scrollX = TRUE))

    })

    productive_aa <- reactive({
        if (is.null(rda_envir)) {
            LymphoSeq2::productiveSeq(study_table = airr_data(), aggregate = "junction_aa")
        } else {
            rda_envir$amino_table
        }
    })

    productive_nt <- reactive({
        if (is.null(rda_envir)) {
            LymphoSeq2::productiveSeq(study_table = airr_data(), aggregate = "junction")
        } else {
            rda_envir$nucleotide_table
        }
    })

    unique_prod_rep <- reactive({
        unique(productive_aa()[, "repertoire_id"])
    })

    clonality_data <- reactive({
        if (is.null(rda_envir)) {
            LymphoSeq2::clonality(airr_data())
        } else {
            rda_envir$summary_table
        }
    })

    observeEvent(input$airr_files, {
        airr_data()
        if (input$chord_button || input$venn_button || input$bar_button ||
                input$plot_button || input$diff_button || input$track_button ||
                input$clonal_relate_button || input$word_button || input$kmer_button ||
                input$top_seq_button) {
            shinyjs::hide("chord_diagram")
            shinyjs::hide("commonSeqs_venn")
            shinyjs::hide("commonSeqs_bar")
            shinyjs::hide("commonSeqs_plot")
            shinyjs::hide("diff_abundance")
            shinyjs::hide("gene_freq_word")
            shinyjs::hide("clonal_relate")
            shinyjs::hide("clone_track")
            shinyjs::hide("count_kmers")
            shinyjs::hide("top_seq_table")
            shinyjs::hide("top_seq_plot")
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

    data_frame_tabs <- c("airr_table", "seq_counts", "count_stats",
                         "public_tcrb_seq", "diff_abundance", "count_kmers")

    observeEvent(input$tabselected, {
        if (input$tabselected %in% data_frame_tabs ||
                input$tabselected == "common_panel" && input$common_sub_tab == "common_seq_table" ||
                input$tabselected == "prod_seq_panel" && (input$prod_seq_sub_tab == "top_seq_table" || input$prod_seq_sub_tab == "produtive_seq_table") ||
                input$tabselected == "clonality_panel" && input$clonal_sub_tab != "clonality_plot" ||
                input$tabselected == "gene_panel" && input$gene_sub_tab == "gene_freq_table") {
            download_choices <- c("TSV" = ".tsv", "EXCEL" = ".xlsx", "RData" = ".rda")
            if (length(input$airr_files) > 0) {
                shinyjs::enable("download")
            }
        } else {
            download_choices <- c("PDF" = ".pdf", "RData" = ".rda")
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
                choices = c("TSV" = ".tsv", "EXCEL" = ".xlsx", "RData" = ".rda"))
        } else {
            shiny::updateRadioButtons(session, "download_type", choices = c("PDF" = ".pdf", "RData" = ".rda"))
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
                choices = c("TSV" = ".tsv", "EXCEL" = ".xlsx", "RData" = ".rda"))
        } else {
            shiny::updateRadioButtons(session, "download_type", choices = c("PDF" = ".pdf", "RData" = ".rda"))
        }
        update_gene_tabs()
    })

    observeEvent(input$clonal_sub_tab, {
        if (input$clonal_sub_tab == "clonality_plot") {
            shiny::updateRadioButtons(session, "download_type", choices = c("PDF" = ".pdf", "RData" = ".rda"))
        } else {
            shiny::updateRadioButtons(session, "download_type",
                choices = c("TSV" = ".tsv", "EXCEL" = ".xlsx", "RData" = ".rda"))
        }
    })

    observeEvent(input$prod_seq_sub_tab, {
        if (input$prod_seq_sub_tab == "top_seq_table" || input$prod_seq_sub_tab == "produtive_seq_table") {
            shiny::updateRadioButtons(session, "download_type",
                choices = c("TSV" = ".tsv", "EXCEL" = ".xlsx", "RData" = ".rda"))
        } else {
            shiny::updateRadioButtons(session, "download_type", choices = c("PDF" = ".pdf", "RData" = ".rda"))
        }
    })

    observeEvent(input$chord_button, {
        shinyjs::show("chord_diagram")
    })

    observeEvent(input$top_seq_chord, {
        if (input$top_seq_chord == "no") {
            shinyjs::disable("top_chord_num")
        } else {
            shinyjs::enable("top_chord_num")
        }
    })

    chord_data <- reactive({
        if (input$top_seq_chord == "yes") {
            validate(
                need(input$top_chord_num > 0, "Please enter number of top sequences.")
            )
            chord_table <- LymphoSeq2::topSeqs(productive_table = productive_nt(), top = input$top_chord_num)
        } else {
            chord_table <- productive_nt()
        }
        if (input$vdj_association == 'VJ') {
            vj <- chord_table %>% 
                select(v_family, j_family) %>% 
                mutate(v_family = replace_na(v_family, "Unresolved"),
                       j_family = replace_na(j_family, "Unresolved"))

            vj <- vj %>% 
                group_by(v_family, j_family) %>%
                summarize(duplicate_count = n(), .groups = 'drop') %>% 
                pivot_wider(id_cols=v_family, names_from = j_family, values_from = duplicate_count) 
            row_names <- vj$v_family 
            vj <- vj %>%
                    dplyr::select(-v_family)
            vj[is.na(vj)] <- 0
            vj <- as.matrix(vj)
            rownames(vj) <- row_names
            vj

        } else if (input$vdj_association == 'DJ') {
            dj <- chord_table %>% 
                select(d_family, j_family) %>% 
                mutate(d_family = replace_na(d_family, "Unresolved"), j_family = replace_na(j_family, "Unresolved"))
            dj <- dj %>% 
                group_by(d_family, j_family) %>% 
                summarize(duplicate_count = n(), .groups = 'drop') %>% 
                pivot_wider(id_cols=d_family, names_from = j_family, values_from = duplicate_count)
            row_names <- dj$d_family
            dj <- dj %>%
                    dplyr::select(-d_family)
            dj[is.na(dj)] <- 0
            dj <- as.matrix(dj)
            rownames(dj) <- row_names
            dj
        }
    })

    output$chord_diagram <- renderChorddiag({
        input$chord_button
        isolate({
            validate(
                need(input$vdj_association != "", "Please select VDJ Association")
            )
            print(chord_data())
            chorddiag::chorddiag(chord_data(), type = "bipartite",
                        showTicks = FALSE,
                        groupnameFontsize = 15,
                        groupnamePadding = 10, margin = 90)
        })
    }) %>%
    bindCache(chord_data()) %>%
    bindEvent(input$chord_button)

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
        print(c(color_intersect))
        LymphoSeq2::commonSeqsBar(productive_aa(),
                    input$bar_id, color_rep_id, c(color_intersect),
                    labels = "yes")
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
                    shinyalert::shinyalert("No Common Sequences!",
                        "Select different repertoire ids.", type = "error")
                    # shiny::showNotification("Cannot render plot", "", type = "error")
                    return()
                }
            )
        })
    })

    observeEvent(input$plot_button, {
        # shinyjs::enable("download")
        shinyjs::show("commonSeqs_plot")
    })

    common_seqs_plot <- reactive({
        LymphoSeq2::commonSeqsPlot(input$plot_id[1], input$plot_id[2], productive_aa())
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
                shinyalert::shinyalert("No Common Sequences!",
                    "Select different repertoire ids.", type = "error")
                # shiny::showNotification("No common sequences", "", type = "error")
                return()
            }
        })
    })

    observeEvent(input$venn_button, {
        shinyjs::show("commonSeqs_venn")
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
        LymphoSeq2::lorenzCurve(repertoire_ids, airr_data()) + ggplot2::coord_fixed(1/2)
    })

    output$lorenz <- renderPlotly({
        lorenz_data()
    })

    seq_count_data <- reactive({
        data_output <- clonality_data() %>%
            select(repertoire_id, total_sequences,
                    unique_productive_sequences, total_count)
        data_output
    })

    output$seq_count <- DT::renderDataTable({
        seq_count_data() %>%
            DT::datatable(rownames = FALSE,
                        colnames = c("Repertoire ID",
                            "Total Sequences", "Unique Productive Sequences",
                            "Total Count"),
                        filter = "top",
                        options = list(scrollX = TRUE))
    })

    stats_count_data <- reactive({
        stats <- clonality_data() %>%
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
        stats_count_data() %>%
            DT::datatable(rownames = FALSE,
                            colnames = c("Metric", "Total",
                                "Mean", "Minimum", "Maximum"),
                            filter = "top",
                            options = list(scrollX = TRUE))
    })

    prod_plot_data <- reactive({
        ggplot(data = clonality_data(), aes(x = repertoire_id, y = unique_productive_sequences)) +
            geom_bar(stat = "identity", position=position_dodge(), fill = "#3182bd") +
            theme_minimal() +
            theme(axis.text.x  = element_text(angle=90, vjust=0.5, hjust = 1), text = element_text(size = 10)) +
            scale_x_discrete(limits = clonality_data()$repertoire_id) +
            labs(x = "", y = "Unique productive sequences")
    }) 

    output$productive_plot <- renderPlotly({
        plotly::ggplotly(prod_plot_data())
    })

    unique_prod_data <- reactive({
        LymphoSeq2::uniqueSeqs(productive_table = productive_aa())
    })

    output$produtive_seq_table <- DT::renderDataTable({
        unique_prod_data() %>%
            DT::datatable(filter = "top", options = list(scrollX = TRUE))
    })

    observeEvent(input$top_seq_button, {
        shinyjs::show("top_seq_table")
        shinyjs::show("top_seq_plot")
    })

    top_prod_seq_data <- reactive({
        LymphoSeq2::topSeqs(productive_aa(), input$top_seq_num) %>%
            select(repertoire_id, junction_aa, duplicate_count, duplicate_frequency,
                        v_call, d_call, j_call)
    })

    output$top_seq_table <- DT::renderDataTable({
        input$top_seq_button
        isolate({
            validate(
                need(input$top_seq_num > 0, "Please choose number of top sequences.")
            )
            top_prod_seq_data() %>%
                DT::datatable(rownames = FALSE, filter = "top",
                            options = list(scrollX = TRUE))
        })
    })

    top_seq_data <- reactive({
        LymphoSeq2::topSeqsPlot(productive_aa(), top = input$top_seq_num)
    })

    output$top_seq_plot <- renderPlotly({
        input$top_seq_button
        isolate({
            validate(
                need(input$top_seq_num > 0, "Please choose number of top sequences.")
            )
            plotly::ggplotly(top_seq_data(), tooltip = c("x", "y", "fill", "text"))
        })
    })

    common_table_data <- reactive({
        LymphoSeq2::commonSeqs(study_table = productive_aa())

    })

    output$common_seq_table <- DT::renderDataTable({
        common_table_data() %>%
            DT::datatable(filter = "top", options = list(scrollX = TRUE))
    })

    output$clonality <- DT::renderDataTable({
        clonality_data() %>%
            DT::datatable(rownames = FALSE, filter = "top",
                            options = list(scrollX = TRUE))
    })

    clone_relate_data <- reactive({
        data_output <- LymphoSeq2::clonalRelatedness(airr_data(), editDistance = input$edit_dis)
        if (input$merge_results == "yes") {
            data_output <- merge(clonality_data(), data_output)
        }
        data_output
    })

    observeEvent(input$clonal_relate_button, {
        shinyjs::show("clonal_relate")
    })

    output$clonal_relate <- DT::renderDataTable({
        input$clonal_relate_button
        isolate({
            clone_relate_data() %>%
                DT::datatable(filter = "top", options = list(scrollX = TRUE))
        })
    })

    clone_stats_data <- reactive({
        clonality_data() %>%
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
    })

    output$clonality_stats <- DT::renderDataTable({
        clone_stats_data() %>%
            DT::datatable(rownames = FALSE,
                            colnames = c("Metric", "Mean", "Minimum", "Maximum"),
                            filter = "top",
                            options = list(scrollX = TRUE))
    })

    clone_plot_data <- reactive({
        ggplot(data = clonality_data(), aes(x = repertoire_id, y = clonality)) +
            geom_bar(stat = "identity", position=position_dodge(), fill = "#b2182b") +
            theme_minimal() +
            theme(axis.text.x  = element_text(angle=90, vjust=0.5, hjust = 1), text = element_text(size = 10)) +
            scale_x_discrete(limits = clonality_data()$repertoire_id) +
            labs(x = "", y = "Clonality")
    })

    output$clonality_plot <- renderPlotly({
        plotly::ggplotly(clone_plot_data())
    })

    interactive_alluvial <- function(p) {
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

    observeEvent(input$track_button, {
        shinyjs::show("clone_track")
    })

    clone_track_data <- reactive({
        if (input$top_clones == "yes") {
            validate(
                need(input$top_num > 0, "Please enter a number for tracking top clones")
            )
            ttable <- LymphoSeq2::topSeqs(productive_aa(), top = input$top_num)
            ctable <- LymphoSeq2::cloneTrack(ttable, sample_list = input$track_id)

        } else {
            ctable <- LymphoSeq2::cloneTrack(productive_aa(), sample_list = input$track_id)
        }

        ctable <- dplyr::left_join(data.frame(repertoire_id = input$track_id), ctable, by = "repertoire_id") 
        LymphoSeq2::plotTrack(clone_table = ctable)
    })

    output$clone_track <- renderPlot({
        input$track_button
        isolate({
            validate(
                need(length(input$track_id) > 0, "Please select repertoire ids to track")
            )
            interactive_alluvial(clone_track_data())
            clone_track_data() #+ 
                # geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
                # theme(legend.position = "right")
        })
    })

    output$clone_track_tooltip <- renderText(
        alluvial_tooltip(input$clone_track_hover)
    )

    pairwise_sim_data <- reactive({
        similarity <- LymphoSeq2::scoringMatrix(productive_aa(), mode = input$mode)
        LymphoSeq2::pairwisePlot(matrix = similarity) +
            labs(fill = paste(input$mode, "\ncoefficient")) +
            theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8)) +
            scale_fill_gradient(low = "#deebf7", high = "#3182bd") +
            theme(legend.position = c(0.8, 0.8))
    })

    output$pairwise_sim <- renderPlotly({
        plotly::ggplotly(pairwise_sim_data())
    }) %>%
    bindCache(productive_aa(), input$mode)

    public_table_data <- reactive({
        published <- LymphoSeq2::searchPublished(productive_aa())
        published %>%
            filter(!is.na(PMID))
    })

    output$public_tcrb <- DT::renderDataTable({
        public_table_data() %>%
                DT::datatable(filter = "top", options = list(scrollX = TRUE))
    })

    gene_heatmap_data <- reactive({
        genes <- LymphoSeq2::geneFreq(productive_nt(), locus = input$locus, family = input$family_bool) %>% 
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
        heatmaply::heatmaply(gene_heatmap_data(), scale = "row", fontsize_row = 6)
    }) %>%
    bindCache(gene_heatmap_data())

    gene_bar_data <- reactive({
        genes <- LymphoSeq2::geneFreq(productive_nt(), input$locus, input$family_bool)
        ggplot(genes, aes(x = repertoire_id, y = gene_frequency, fill = gene_name)) +
            geom_bar(stat = "identity") +
            theme_minimal() + 
            scale_y_continuous(expand = c(0, 0)) + 
            guides(fill = guide_legend(ncol = 3)) +
            labs(y = "Frequency (%)", x = "", fill = "") +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    })

    output$gene_freq_bar <- renderPlotly({
        plotly::ggplotly(gene_bar_data())
    }) %>%
    bindCache(gene_bar_data())

    observeEvent(input$word_button, {
        # shinyjs::enable("download")
        shinyjs::show("gene_freq_word")
    })

    gene_word_data <- reactive({
        genes <- LymphoSeq2::geneFreq(productive_nt(), input$locus, input$family_bool)
        gene_data <- genes %>%
                        filter(repertoire_id == input$word_id)
        gene_data <- data.frame(gene_data$gene_name, gene_data$gene_frequency)
        wordcloud2::wordcloud2(gene_data, color = 'random-dark', minSize = 10)
    })

    output$gene_freq_word <- renderWordcloud2({
        input$word_button
        isolate ({
            gene_word_data()
        })
    }) %>%
    bindCache(gene_word_data()) %>%
    bindEvent(input$word_button)

    gene_table_data <- reactive({
        LymphoSeq2::geneFreq(productive_nt(), input$locus, input$family_bool)
    })

    output$gene_freq_table <- DT::renderDataTable({
        gene_table_data() %>%
            DT::datatable(filter = "top", options = list(scrollX = TRUE))
    })

    observeEvent(input$diff_button, {
        shinyjs::show("diff_abundance")
    })

    diff_table_data <- reactive({
        LymphoSeq2::differentialAbundance(productive_aa(),
                        input$diff_id, q = input$q_val, zero = input$zero_val)
    })

    output$diff_abundance <- DT::renderDataTable({
        input$diff_button
        validate(
            need(length(input$diff_id) == 2, "Please select 2 repertoire ids")
        )
        isolate({
            diff_table_data() %>%
                DT::datatable(filter = "top", options = list(scrollX = TRUE))
        })
    })

    kmer_table_data <- reactive({
        LymphoSeq2::countKmer(airr_data(), input$k_val)
    })

    output$count_kmers <- DT::renderDataTable({
        input$kmer_button
        isolate({
            validate(
                need(input$k_val > 1, "Please input length of at least 2")
            )
            kmer_table_data() %>%
                DT::datatable(filter = "top", options = list(scrollX = TRUE))
        })
    })

    output$download <- downloadHandler(
        filename <- function() {
            paste0(input$tabselected, input$download_type)
        },
        content <- function(file) {
            if (input$tabselected == "gene_panel") {
                if (input$gene_sub_tab == "gene_freq_heat") {
                    RedBlue <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(256)
                    data <- list(gene_heatmap_data(), RedBlue)
                    data_output <- pheatmap::pheatmap(data[[1]], color = data[[2]], scale = "row")
                } else if (input$gene_sub_tab == "gene_freq_bar") {
                    data_output <- plot(gene_bar_data())
                } else {
                    data_output <- gene_table_data()
                }

            } else if (input$tabselected == "chord_diagram") {
                data_output <- circlize::chordDiagram(chord_data(), annotationTrack = c("grid", "name"))

            } else if (input$tabselected == "common_panel") {
                if (input$common_sub_tab == "common_bar") {
                    data_output <- common_bar_data()
                } else if (input$common_sub_tab == "common_venn") {
                    data_output <- common_venn_data()
                } else if (input$common_sub_tab == "common_plot") {
                    data_output <- common_seqs_plot()
                } else {
                    data_output <- common_table_data()
                }
            } else if (input$tabselected == "prod_seq_panel") {
                if (input$prod_seq_sub_tab == "produtive_seq_plot") {
                    data_output <- prod_plot_data()
                } else if (input$prod_seq_sub_tab == "top_seq_plot") {
                    data_output <- top_seq_data()
                } else if (input$prod_seq_sub_tab == "produtive_seq_table") {
                    data_output <- unique_prod_data()
                } else {
                    data_output <- top_prod_seq_data()
                }

            } else if (input$tabselected == "clonality_plot") {
                data_output <- clone_plot_data()
            } else if (input$tabselected == "lorenz_curve") {
                data_output <- lorenz_data()
            } else if (input$tabselected == "pairwise_sim") {
                data_output <- pairwise_sim_data()
            } else if (input$tabselected == "clone_track") {
                data_output <- clone_track_data()
            }

            else if (input$tabselected == "airr_table") {
                data_output <- airr_data()
            } else if (input$tabselected == "clonality_panel") {
                if (input$clonal_sub_tab == "clonality") {
                    data_output <- clonality_data()
                } else if (input$clonal_sub_tab == "clonal_relate") {
                    data_output <- clone_relate_data()
                } else if (input$clonal_sub_tab == "clonality_stats") {
                    data_output <- clone_stats_data()
                } else if (input$tabselected == "seq_counts") {
                    data_output <- seq_count_data()
                } else if (input$tabselected == "count_stats") {
                    data_output <- stats_count_data()
                }
            } else if (input$tabselected == "public_tcrb_seq") {
                data_output <- public_table_data()
            } else if (input$tabselected == "count_kmers") {
                data_output <- kmer_table_data()
            } else if (input$tabselected == "diff_abundance") {
                data_output <- diff_table_data()
            }

            if (input$download_type == ".rda") {
                shiny::withProgress(
                    message = "This make take a few minutes...",
                    {
                        airr_table <- airr_data()
                        prod_aa <- productive_aa()
                        prod_nt <- productive_nt()
                        clone_data <- clonality_data()
                        save(data_output, airr_table, prod_aa, prod_nt, clone_data, file = file)
                    }
                )
                # airr_table <- airr_data()
                # prod_aa <- productive_aa()
                # prod_nt <- productive_nt()
                # clone_data <- clonality_data()
                # save(data_output, airr_table, prod_aa, prod_nt, clone_data, file = file)

            } else if (input$download_type == ".pdf") {
                pdf(file, width = 11, height = 8.5)

                if (input$tabselected == "gene_panel") {
                    print(data_output)
                } else if (input$tabselected == "chord_diagram") {
                    data_output <- circlize::chordDiagram(chord_data(), annotationTrack = c("grid", "name"))
                } else if (input$tabselected == "common_panel") {
                    if (input$common_sub_tab == "common_bar") {
                        print(data_output)
                    } else if (input$common_sub_tab == "common_venn") {
                        grid::grid.newpage()
                        grid::grid.draw(data_output)
                    } else {
                        plot(data_output)
                    }
                } else {
                    plot(data_output)
                }
                dev.off()
                # shinyscreenshot::screenshot(id = "chord_diagram")

            } else if (input$download_type == ".tsv") {
                    write.table(data_output, file, quote = FALSE, sep='\t', row.names = FALSE)
            } else if (input$download_type == ".xlsx") {
                    writexl::write_xlsx(data_output, path = file)
            }
        }
    )
}

shinyApp(ui = ui, server = server)
