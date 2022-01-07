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

options(shiny.maxRequestSize = 30 * 1024^2)

ui <- fluidPage(

    useShinyjs(),

    titlePanel("LymphoSeq2 Shiny"),

    sidebarLayout(
        sidebarPanel(
            fileInput("airr_files", label = "Upload Files", multiple = TRUE, accept = ".tsv"),

            conditionalPanel(condition = "input.tabselected == 'chord_diagram'",
                selectizeInput("vdj_association", label = "Select VDJ Association", choices = c("", "VJ", "DJ"), selected = NULL),
                actionButton("chord_button", label = "Create Chord Diagram")),

            conditionalPanel(condition = "input.tabselected == 'common_bar'",
                selectizeInput("bar_id", label = "Select at least 2 repertoire ids", choices = NULL, multiple = TRUE),
                selectizeInput("color_rep_id", label = "Select 1 repertoire id to color", choices = NULL, selected = NULL),
                selectizeInput("color_intersect", label = "Select at least 2 repertoire ids to color intersections", 
                                choices = NULL, multiple = TRUE, selected = NULL),
                actionButton("bar_button", label = "Create Bar Chart")),

            conditionalPanel(condition = "input.tabselected == 'common_plot'",
                selectizeInput("plot_id", label = "Select 2 repertoire ids", choices = NULL, multiple = TRUE),
                radioButtons("show_type", "Show Sequences Type", choices = c("common", "all"), inline = TRUE),
                actionButton("plot_button", label = "Create Plot")),
            
            conditionalPanel(condition = "input.tabselected == 'common_venn'",
                selectizeInput("venn_id", label = "Select 2 or 3 repertoire ids", choices = NULL, multiple = TRUE),
                actionButton("venn_button", label = "Create Venn Diagram")),
            
            conditionalPanel(condition = "input.tabselected == 'gene_freq'",
                selectizeInput("locus", label = "Select which VDJ genes to include",
                                choices = c("VDJ", "DJ", "VJ", "DJ", "V", "D", "J")),
                radioButtons("family_bool", "Select which names to use",
                            choices = c("family" = TRUE, "gene" = FALSE), inline = TRUE)),
            
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
                tabPanel("Data Table", value = "airr_table",
                        DT::dataTableOutput("table") %>% withSpinner()),
                tabPanel("Chord Diagram VDJ", value = "chord_diagram",
                        chorddiagOutput("chord_diagram") %>% withSpinner()),
                tabPanel("Common Sequences Bar Chart", value = "common_bar",
                        plotOutput("commonSeqs_bar") %>% withSpinner()),
                tabPanel("Common Sequences Plot", value = "common_plot",
                        plotlyOutput("commonSeqs_plot") %>% withSpinner()),
                tabPanel("Common Sequences Venn Diagram", value = "common_venn",
                        plotOutput("commonSeqs_venn") %>% withSpinner()),
                tabPanel("Lorenz Curve", value = "lorenz_curve",
                        plotlyOutput("lorenz") %>% withSpinner()),
                tabPanel("Sequencing Counts", value = "seq_counts",
                        DT::dataTableOutput("seq_count") %>% withSpinner()),
                tabPanel("Count Statistics", value = "count_stats",
                        DT::dataTableOutput("stats_count") %>% withSpinner()),
                tabPanel("Unique Productive Sequences Plot", value = "produtive_seq_plot",
                        plotlyOutput("productive_plot") %>% withSpinner()),
                tabPanel("Clonality", value = "clonality",
                        DT::dataTableOutput("clonality") %>% withSpinner()),
                tabPanel("Clonality Statistics", value = "clonality_stats",
                        DT::dataTableOutput("clonality_stats") %>% withSpinner()),
                tabPanel("Clonality Plot", value = "clonality_plot",
                        plotlyOutput("clonality_plot") %>% withSpinner()),
                tabPanel("Top 10 Sequences Plot", value = "top_10_seq_plot",
                        plotlyOutput("top_10_seq") %>% withSpinner()),
                tabPanel("Top Productive Sequences Table", value = "top_10_seq_table",
                        DT::dataTableOutput("top_productive_seq") %>% withSpinner()),
                tabPanel("Common Sequences Table", value = "common_seq_table",
                        DT::dataTableOutput("common_seq_table") %>% withSpinner()),
                tabPanel("Clone Tracking", tags$div(
                        style = "position: relative;",
                        plotOutput("clone_track",
                        hover = hoverOpts(id = "clone_track_hover")) %>% withSpinner(),
                        htmlOutput("clone_track_tooltip")),
                        value = "clone_tacking"),
                tabPanel("Common Sequences Alluvial Plot", tags$div(
                        style = "position: relative;",
                        plotOutput("common_seq_plot",
                        hover = hoverOpts(id = "common_seq_hover")) %>% withSpinner(),
                        htmlOutput("common_seq_tooltip")),
                        value = "common_seq_alluvial"),
                tabPanel("Pairwise Similarity", value = "pairwise_sim",
                        plotlyOutput("pairwise_sim") %>% withSpinner()),
                tabPanel("Public TCRB Sequences", value = "public_tcrb_seq",
                        DT::dataTableOutput("public_tcrb") %>% withSpinner()),
                tabPanel("V Gene Frequency", value = "v_gene_freq",
                        plotlyOutput("v_gene_freq") %>% withSpinner()),
                tabPanel("Gene Frequencies", value = "gene_freq",
                        DT::dataTableOutput("gene_freq") %>% withSpinner()),
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

    airr_data <- reactive({
        validate(
            need(!is.null(input$airr_files), "Please select files to upload to render output.")
                # %then%
                # need(try(LymphoSeq2::readImmunoSeq(input$airr_files$datapath)),
                #      "Invalid File")
        )
        table <- LymphoSeq2::readImmunoSeq(input$airr_files$datapath)
        files <- c(input$airr_files$name)
        in_files <- c()
        for (i in files) {
            in_files <- append(in_files, substr(i, 1, str_length(i) - 4))
        }
        table <- table %>%
                mutate(repertoire_id = if_else(as.integer(repertoire_id) <= length(in_files),
                        in_files[as.integer(repertoire_id) + 1], "5"))
        table
    })

    output$table <- DT::renderDataTable({
        data_output <<- airr_data()
        DT::datatable(data_output)
    })

    productive_aa <- reactive({
        LymphoSeq2::productiveSeq(airr_data(), "junction_aa")
    })

    observeEvent(input$airr_files, {
        if (input$tabselected == "chord_diagram" || input$tabselected == "common_venn" ||
                input$tabselected == "common_bar" || input$tabselected == "common_plot") {
            shinyjs::disable("download")
        } else {
            shinyjs::enable("download")
        }
        productive_rep_id <- unique(productive_aa()[, "repertoire_id"])
        rep_with_null <- c("none", productive_rep_id)
        shiny::updateSelectizeInput(session, "bar_id", choices = productive_rep_id)
        shiny::updateSelectizeInput(session, "color_rep_id", choices = rep_with_null)
        shiny::updateSelectizeInput(session, "color_intersect", choices = rep_with_null)
        shiny::updateSelectizeInput(session, "plot_id", choices = productive_rep_id)
        shiny::updateSelectizeInput(session, "venn_id", choices = productive_rep_id)
    })

    observeEvent(input$tabselected, {
        data_frame_tabs = c("airr_table", "seq_counts", "count_stats", "clonality",
                            "clonality_stats", "top_10_seq_table", "common_seq_table", 
                            "public_tcrb_seq", "gene_freq")
        if (input$tabselected %in% data_frame_tabs) {
            download_choices <- c("TSV" = ".tsv", "EXCEL" = ".xlsx")
        } else {
            download_choices <- c("PDF" = ".pdf")
        }
        updateRadioButtons(session, "download_type", choices = download_choices)
        if (input$tabselected == "common_venn") {
            if (length(input$venn_id) == 0) {
                shinyjs::disable("download")
                hide("commonSeqs_venn")
            } else {
                shinyjs::enable("download")
                show("commonSeqs_venn")
            }
        } else if (input$tabselected == "chord_diagram") {
            if (input$vdj_association == "") {
                shinyjs::disable("download")
                hide("chord_diagram")
            } else {
                shinyjs::enable("download")
                show("chord_diagram")
            }
        } else if (input$tabselected == "common_bar") {
            if (length(input$bar_id) == 0) {
                shinyjs::disable("download")
                hide("commonSeqs_bar")
            } else {
                shinyjs::enable("download")
                show("commonSeqs_bar")
            }
        } else if (input$tabselected == "common_plot") {
            if (length(input$plot_id) == 0) {
                shinyjs::disable("download")
                hide("commonSeqs_plot")
            } else {
                shinyjs::enable("download")
                show("commonSeqs_plot")
            }
        } else if (length(input$airr_files) != 0) {
            shinyjs::enable("download")
        }
    })

    observeEvent(input$chord_button, {
        shinyjs::enable("download")
        show("chord_diagram")
    })

    output$chord_diagram <- renderChorddiag({
        input$chord_button
        isolate ({
            productive_nt <- LymphoSeq2::productiveSeq(study_table = airr_data(), aggregate = "junction")
            top_seqs <- LymphoSeq2::topSeqs(productive_table = productive_nt, top = 1)
            if (input$vdj_association == 'VJ') {
                vj <- top_seqs %>% 
                    select(v_family, j_family) %>% 
                    mutate(v_family = replace_na(v_family, "Unresolved"), j_family = replace_na(j_family, "Unresolved"))

                vj <- vj %>% 
                    group_by(v_family, j_family) %>% 
                    summarize(duplicate_count = n(), .groups = 'drop') %>% 
                    pivot_wider(id_cols=v_family, names_from = j_family, values_from = duplicate_count) 
                row_names <- vj$v_family 
                vj <- vj %>% select(-v_family)
                vj[is.na(vj)] <- 0
                vj <- as.matrix(vj)
                rownames(vj) <- row_names
                data_output <<- vj
                chorddiag(vj, type = "bipartite")

            } else if (input$vdj_association == 'DJ') {
                dj <- top_seqs %>% 
                    select(d_family, j_family) %>% 
                    mutate(d_family = replace_na(d_family, "Unresolved"), j_family = replace_na(j_family, "Unresolved"))
                dj <- dj %>% 
                    group_by(d_family, j_family) %>% 
                    summarize(duplicate_count = n(), .groups = 'drop') %>% 
                    pivot_wider(id_cols=d_family, names_from = j_family, values_from = duplicate_count)
                row_names <- dj$d_family
                dj <- dj %>% select(-d_family)
                dj[is.na(dj)] <- 0
                dj <- as.matrix(dj)
                rownames(dj) <- row_names
                data_output <<- dj
                chorddiag(dj, type = "bipartite")
            }
        })

    })

    observeEvent(input$bar_button, {
        show("commonSeqs_bar")
        shinyjs::enable("download")
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
            color_rep_id <- input$color_rep_id
            color_intersect <- input$color_intersect
            if (input$color_rep_id == "none") {
                color_rep_id <- NULL
            }
            if (!is.null(input$color_intersect) & "none" %in% input$color_intersect) {
                color_intersect <- NULL
            }
            data_output <<- LymphoSeq2::commonSeqsBar(productive_aa(),
                input$bar_id, color_rep_id, color_intersect,
                labels = "yes")
            data_output
        })
    })

    observeEvent(input$plot_button, {
        show("commonSeqs_plot")
        shinyjs::enable("download")
    })

    output$commonSeqs_plot <- renderPlotly({
        input$plot_button
        isolate({
            validate(
                need(length(input$plot_id) == 2,
                    "Please select exactly 2 repertoire ids")
            )
            rep_id_legal <- c()
            aa <- productive_aa()
            digit_vector <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")
            replaced <- c()

            for (i in input$plot_id) {
                if (startsWith(substr(i, 1, 1), digit_vector)) {
                    rep_id_legal <- append(rep_id_legal, paste("plot_", i, sep=''))
                    replaced <- append(replaced, i)
                } else {
                    rep_id_legal <- append(rep_id_legal, i)
                }
            }
            if (length(replaced) != 0) {
                aa <- aa %>% mutate(repertoire_id = dplyr::if_else(repertoire_id %in% replaced, 
                    paste0("plot_", repertoire_id), repertoire_id))
            }
            data_output <<- LymphoSeq2::commonSeqsPlot(rep_id_legal[1], rep_id_legal[2], aa)
            data_output
        })
    })

    observeEvent(input$venn_button, {
        show("commonSeqs_venn")
        shinyjs::enable("download")
    })

    output$commonSeqs_venn <- renderPlot({
        input$venn_button
        isolate({
            validate(
                need(length(input$venn_id) == 2 | length(input$venn_id) == 3,
                    "Please select only 2 or 3 repertoire ids")
            )
            repertoire_ids <- input$venn_id
            if (length(repertoire_ids) == 2) {
                a <- productive_aa() %>% 
                    filter(repertoire_id == repertoire_ids[[1]])
                b <- productive_aa() %>% 
                    filter(repertoire_id == repertoire_ids[[2]])
                grid::grid.newpage()
                venn <- VennDiagram::draw.pairwise.venn(area1 = length(a$junction_aa), 
                                                        area2 = length(b$junction_aa), 
                                                        cross.area = length(intersect(a$junction_aa, 
                                                                                    b$junction_aa)), 
                                                        category = c(repertoire_ids[1], 
                                                                    repertoire_ids[2]), 
                                                        cat.fontfamily = rep("sans", 2), 
                                                        fontfamily = rep("sans", 3), 
                                                        fill = c("#3288bd", "#d53e4f"), 
                                                        cat.pos = c(0, 0),
                                                        cat.dist = rep(0.025, 2),
                                                        cex = 1, 
                                                        cat.cex = 0.7,
                                                        lwd = rep(2, 2))
                data_output <<- venn
                grid::grid.draw(venn)
            }
            if (length(repertoire_ids) == 3) {
                a <- productive_aa %>% 
                    filter(repertoire_id == repertoire_ids[[1]])
                b <- productive_aa %>% 
                    filter(repertoire_id == repertoire_ids[[2]])
                c <- productive_aa %>% 
                    filter(repertoire_id == repertoire_ids[[3]])
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
                                                    category = c(repertoire_ids[1], 
                                                                repertoire_ids[2], 
                                                                repertoire_ids[3]), 
                                                    cat.fontfamily = rep("sans", 3), 
                                                    fontfamily = rep("sans", 7), 
                                                    fill = c("#3288bd", "#abdda4", "#d53e4f"), 
                                                    cat.pos = c(0, 0, 180), 
                                                    cat.dist = rep(0.025, 3),
                                                    cex = 1, 
                                                    cat.cex = 0.7,
                                                    lwd = rep(2, 3))
                data_output <<- venn
                grid::grid.draw(venn)
            }
        })
    })

    output$lorenz <- renderPlotly({
        repertoire_ids <- productive_aa() %>% pull(repertoire_id) %>% unique()
        data_output <<- LymphoSeq2::lorenzCurve(repertoire_ids, airr_data()) + coord_fixed(1/2)
        data_output
    })

    itable <- reactive({
        LymphoSeq2::clonality(airr_data())
    })

    output$seq_count <- DT::renderDataTable({
        data_output <<- itable() %>%
            select(repertoire_id, total_sequences,
                    unique_productive_sequences, total_count)
        data_output %>%
            datatable(rownames = FALSE, colnames = c("Repertoire ID",
                        "Total Sequences", "Unique Productive Sequences",
                        "Total Count"), filter = "top")
    })

    output$stats_count <- DT::renderDataTable({
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
        data_output <<- stats
        datatable(stats, rownames = FALSE, colnames = c("Metric", "Total",
                    "Mean", "Minimum", "Maximum"), filter = "top")
    })

    output$productive_plot <- renderPlotly({
        data_output <<- ggplot(data = itable(), aes(x = repertoire_id, y = unique_productive_sequences)) +
            geom_bar(stat = "identity", position=position_dodge(), fill = "#3182bd") +
            theme_minimal() +
            theme(axis.text.x  = element_text(angle=90, vjust=0.5, hjust = 1), text = element_text(size = 10)) +
            scale_x_discrete(limits = itable()$repertoire_id) +
            labs(x = "", y = "Unique productive sequences")
        ggplotly(data_output)
    })

    output$top_productive_seq <- DT::renderDataTable({
        data_output <<- LymphoSeq2::topSeqs(productive_aa()) %>%
            select(repertoire_id, junction_aa, duplicate_count, duplicate_frequency, v_call, d_call, j_call)
        data_output %>%
            datatable(rownames = FALSE, colnames = c("Repertoire ID", "Junction AA", "Duplicate Count", "Duplicate Frequency", "V Call", "D Call", "J Call"), filter = "top") 

    })

    output$top_10_seq <- renderPlotly({
        data_output <<- topSeqsPlot(productive_aa(), top = 10) +
            theme(text = element_text(size = 10)) +
            scale_x_discrete(limits = itable()$repertoire_id)
        ggplotly(data_output)
    })

    output$common_seq_table <- DT::renderDataTable({
        atable <- LymphoSeq2::productiveSeq(airr_data())
        top_freq <- LymphoSeq2::topFreq(atable, frequency = 0.001)
        clone_table <- atable %>%
                        select(repertoire_id, junction_aa)
        top_freq <- left_join(top_freq, atable, by = "junction_aa") %>%
                        select(repertoire_id, junction_aa, duplicate_frequency,
                            numberSamples, antigen, prevalence)
        data_output <<- top_freq
        datatable(top_freq, colnames = c("Repertoire ID", "Junction AA",
                "Duplicate Frequency", "Number of Samples", "Antigen",
                "Prevalence"), filter = "top")
    })

    output$clonality <- DT::renderDataTable({
        data_output <<- itable() %>%
            select(repertoire_id, clonality, gini_coefficient, top_productive_sequence)
        data_output %>%
            datatable(rownames = FALSE, colnames = c("Repertoire ID", "Clonality", "Gini Coefficient", "Top Productive Sequence"), filter = "top")
    })

    output$clonality_stats <- DT::renderDataTable({
        data_output <<- itable() %>%
            select(repertoire_id, clonality, gini_coefficient, top_productive_sequence) %>% 
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
        data_output %>%
            datatable(rownames = FALSE, colnames = c("Metric", "Total", "Mean", "Minimum", "Maximum"), filter = "top")
    })

    output$clonality_plot <- renderPlotly({
        data_output <<- ggplot(data = itable(), aes(x = repertoire_id, y = clonality)) +
            geom_bar(stat = "identity", position=position_dodge(), fill = "#b2182b") +
            theme_minimal() +
            theme(axis.text.x  = element_text(angle=90, vjust=0.5, hjust = 1), text = element_text(size = 10)) +
            scale_x_discrete(limits = itable()$repertoire_id) +
            labs(x = "", y = "Clonality")
        ggplotly(data_output)
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
        xrange_new <- c(1 - alluvium_width/2, max(pbuilt$data[[1]]$x)
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

    output$clone_track <- renderPlot({
        ctable <- LymphoSeq2::cloneTrack(productive_aa())
        max_seen <- ctable %>%
            pull(seen) %>%
            max()
        max_amino <- ctable %>%
             filter(seen == max_seen) %>%
             pull(junction_aa) %>%
             unique()
        top_seen <- ctable %>%
            slice_max(seen, n = 100)
        p <- LymphoSeq2::plotTrack(clone_table = top_seen, alist = max_amino)
        interactive_alluvial(p)
        data_output <<- p
        p + theme(legend.position = "none")
    })

    output$clone_track_tooltip <- renderText(
        alluvial_tooltip(input$clone_track_hover)
    )

    output$common_seq_plot <- renderPlot({
        ttable <- LymphoSeq2::topSeqs(productive_aa(), top = 50)
        ctable <- LymphoSeq2::cloneTrack(ttable)
        p <- LymphoSeq2::plotTrack(clone_table = ctable)
        interactive_alluvial(p)
        data_output <<- p
        p
    })

    output$common_seq_tooltip <- renderText (
        alluvial_tooltip(input$common_seq_hover)
    )

    output$pairwise_sim <- renderPlotly({
        similarity <- LymphoSeq2::scoringMatrix(productive_aa())
        data_output <<- pairwisePlot(matrix = similarity) +
            labs(fill = "Bhattacharyya\n coefficient") +
            theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8)) +
            scale_fill_gradient(low = "#deebf7", high = "#3182bd") +
            theme(legend.position = c(0.8, 0.8))
        ggplotly(data_output)
    })

    output$public_tcrb <- DT::renderDataTable({
        published <- LymphoSeq2::searchPublished(productive_aa())
        data_output <<- published %>%
            filter(!is.na(PMID)) %>%
            select(repertoire_id, junction_aa, duplicate_count, antigen, prevalence) #%>%
        data_output %>%
            datatable(colnames = c("Repertoire ID", "Junction AA", "Duplicate Count", "Antigen", "Prevalence"), filter = "top")
    })

    output$v_gene_freq <- renderPlotly({
        vGenes <- LymphoSeq2::geneFreq(productive_aa(), locus = "V", family = TRUE) %>% 
            pivot_wider(id_cols = gene_name,
                        names_from = repertoire_id,
                        values_from = gene_frequency,
                        values_fn = sum,
                        values_fill = 0)
        gene_names <- vGenes %>%
              pull(gene_name)
        vGenes <- vGenes %>%
                  select(-gene_name) %>%
                  as.matrix()
        rownames(vGenes) <- gene_names
        RedBlue = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(256)
        data_output <<- list(vGenes, RedBlue) # pheatmap(vGenes, color = RedBlue, scale = "row")
        heatmaply(vGenes, colors = "RdBu", scale = "row")
    })

    output$gene_freq <- DT::renderDataTable({
        productive_nt <- LymphoSeq2::productiveSeq(study_table = airr_data(), aggregate = "nucleotide")
        data_output <<- LymphoSeq2::geneFreq(productive_nt, input$locus, input$family_bool)
        data_output
    })

    output$download <- downloadHandler(
        filename <- function() {
            paste0(input$tabselected, input$download_type)
        },
        content <- function(file) {
            if (input$download_type == ".pdf") {
                pdf(file)
                if (input$tabselected == "v_gene_freq") {
                    pheatmap(data_output[[1]], color = data_output[[2]], scale = "row")
                } else if (input$tabselected == "chord_diagram") {
                    circlize::chordDiagram(data_output, annotationTrack = c("grid", "name"))
                } else if (input$tabselected == "common_venn") {
                    grid::grid.newpage()
                    grid::grid.draw(data_output)
                } else if (input$tabselected == "common_bar") {
                    print(data_output)
                } else {
                    plot(data_output)
                }
                dev.off()
            } else if (input$download_type == ".tsv") {
                write.table(data_output, file, quote = FALSE, sep='\t', row.names = FALSE)
            } else if (input$download_type == ".xlsx") {
                write_xlsx(data_output, path = file)
            }
        }
    )

}

shinyApp(ui = ui, server = server)
