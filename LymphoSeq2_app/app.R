library(shiny)
library(shinycssloaders)
library(LymphoSeq2)
library(DT)
library(tidyverse)
library(ggplot2)
library(plotly)
library(heatmaply)
library(ggalluvial)
library(htmltools)
library(sp)
library(chorddiag)


ui <- fluidPage(
    sidebarLayout(
        sidebarPanel(
            fileInput("airr_files", label = "Upload Files", multiple = TRUE, accept = ".tsv"),

            conditionalPanel(condition = "input.tabselected == 'chord_diagram'",
                selectizeInput("vdj_association", label = "Select VDJ Association", choices = c("VJ", "DJ")),
                actionButton("chord_button", label = "Create Chord Diagram")),

            conditionalPanel(condition = "input.tabselected == 'common_bar'",
                selectizeInput("rep_id", label = "Select repertoire ids", choices = NULL, multiple = TRUE),
                actionButton("bar_button", label = "Create Bar Chart")),

            conditionalPanel(condition = "input.tabselected == 'common_plot'",
                selectizeInput("plot_id", label = "Select 2 repertoire ids", choices = NULL, multiple = TRUE),
                actionButton("plot_button", label = "Create Plot")),
            
            conditionalPanel(condition = "input.tabselected == 'common_venn'",
                selectizeInput("venn_id", label = "Select 2 or 3 repertoire ids", choices = NULL, multiple = TRUE),
                actionButton("venn_button", label = "Create Venn Diagram")),
                
            conditionalPanel(condition = "input.tabselected == 'phylo_tree'",
                selectizeInput("phylo_id", label = "Select one repertoire id", choices = NULL, multiple = TRUE),
                actionButton("phylo_button", label = "Create Phylogenetic Tree")),
        ),

        mainPanel(
            tabsetPanel(
                tabPanel("Data Table",
                        DT::dataTableOutput("table") %>% withSpinner()),
                tabPanel("Chord Diagram VDJ", value = "chord_diagram",
                        chorddiagOutput("chord_diagram") %>% withSpinner()),
                tabPanel("Common Sequences Bar Chart", value = "common_bar",
                        plotOutput("commonSeqs_bar") %>% withSpinner()),
                tabPanel("Common Sequences Plot", value = "common_plot",
                        plotlyOutput("commonSeqs_plot") %>% withSpinner()),
                tabPanel("Common Seqs Venn Diagram", value = "common_venn",
                        plotOutput("commonSeqs_venn") %>% withSpinner()),
                tabPanel("Lorenz Curve",
                        plotlyOutput("lorenz") %>% withSpinner()),
                tabPanel("Phylogenetic Tree", value = "phylo_tree",
                        plotOutput("phylogenetic_tree") %>% withSpinner()),
                tabPanel("Sequencing Counts",
                        DT::dataTableOutput("seq_count") %>% withSpinner()),
                tabPanel("Count Statistics",
                        DT::dataTableOutput("stats_count") %>% withSpinner()),
                tabPanel("Unique Productive Sequences Plot",
                        plotlyOutput("productive_plot") %>% withSpinner()),
                tabPanel("Clonality",
                        DT::dataTableOutput("clonality") %>% withSpinner()),
                tabPanel("Clonality Statistics",
                        DT::dataTableOutput("clonality_stats") %>% withSpinner()),
                tabPanel("Clonality Plot",
                        plotlyOutput("clonality_plot") %>% withSpinner()),
                tabPanel("Top 10 Sequences Plot",
                        plotlyOutput("top_10_seq") %>% withSpinner()),
                tabPanel("Top Productive Sequences Table",
                        DT::dataTableOutput("top_productive_seq") %>% withSpinner()),
                tabPanel("Common Sequences Table",
                        DT::dataTableOutput("common_seq_table") %>% withSpinner()),
                tabPanel("Clone Tracking", tags$div(
                        style = "position: relative;",
                        plotOutput("clone_track",
                        hover = hoverOpts(id = "clone_track_hover")) %>% withSpinner(),
                        htmlOutput("clone_track_tooltip"))),
                tabPanel("Common Sequences Alluvial Plot", tags$div(
                        style = "position: relative;",
                        plotOutput("common_seq_plot",
                        hover = hoverOpts(id = "common_seq_hover")) %>% withSpinner(),
                        htmlOutput("common_seq_tooltip"))),
                tabPanel("Pairwise Similarity",
                        plotlyOutput("pairwise_sim") %>% withSpinner()),
                tabPanel("Public TCRB Sequences",
                        DT::dataTableOutput("public_tcrb") %>% withSpinner()),
                tabPanel("V Gene Frequency",
                        plotlyOutput("v_gene_freq") %>% withSpinner()),
                tabPanel("Gene Frequencies",
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

    airr_data <- reactive({
        validate(
            need(!is.null(input$airr_files), "No file uploaded to render any output")
                %then%
                need(try(LymphoSeq2::readImmunoSeq(input$airr_files$datapath)),
                     "Invalid File")
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

    atable <- reactive({
        LymphoSeq2::productiveSeq(airr_data())
    })

    productive_aa <- reactive({
        LymphoSeq2::productiveSeq(airr_data(), "junction_aa", TRUE)
    })


    output$table <- DT::renderDataTable({
        DT::datatable(airr_data())
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
                chorddiag(vj, type = "bipartite")
            }
            if (input$vdj_association == 'DJ') {
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
                chorddiag(dj, type = "bipartite")
            }
        })

    })


    observeEvent(input$airr_files, {
        rep_id_selection <- unique(airr_data()[, "repertoire_id"])
        productive_rep_id <- unique(productive_aa()[, "repertoire_id"])
        shiny::updateSelectizeInput(session, "rep_id", choices = productive_rep_id)
        shiny::updateSelectizeInput(session, "plot_id", choices = productive_rep_id)
        shiny::updateSelectizeInput(session, "venn_id", choices = productive_rep_id)
        shiny::updateSelectizeInput(session, "phylo_id", choices = rep_id_selection)
    })

    output$commonSeqs_bar <- renderPlot({
        input$bar_button
        isolate({
            validate(
                need(length(input$rep_id) > 1,
                        "You need to select at least 2 repertoire ids")
            )
            LymphoSeq2::commonSeqsBar(productive_aa(), input$rep_id)
        })
    })

    output$commonSeqs_plot <- renderPlotly({
        input$plot_button
        isolate({
            validate(
                need(length(input$plot_id) == 2,
                    "You can only select 2 repertoire ids")
            )
            LymphoSeq2::commonSeqsPlot(input$plot_id[1], input$plot_id[2], productive_aa())
        })
    })

    output$commonSeqs_venn <- renderPlot({
        input$venn_button
        isolate({
            validate(
                need(length(input$venn_id) == 2 | length(input$venn_id) == 3,
                    "You can only select 2 or 3 repertoire ids")
            )
            LymphoSeq2::commonSeqsVenn(input$venn_id, productive_aa())
        })
    })

    output$lorenz <- renderPlotly({
        repertoire_ids <- productive_aa() %>% pull(repertoire_id) %>% unique()
        LymphoSeq2::lorenzCurve(repertoire_ids, airr_data())
    })

    output$phylogenetic_tree <- renderPlot({
        input$phylo_button
        isolate({
            validate(
                need(length(input$phylo_id) == 1,
                "You can only select one repertoire id")
            )
            LymphoSeq2::phyloTree(airr_data(), input$phylo_id,
                        type = "junction", layout = "rectangular", label = FALSE) +
                        ggplot2::theme(legend.position = "none")
        })
    })

    itable <- reactive({
        LymphoSeq2::clonality(airr_data())
    })

    output$seq_count <- DT::renderDataTable({
        itable() %>%
            select(repertoire_id, total_sequences,
                    unique_productive_sequences, total_count) %>%
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
        datatable(stats, rownames = FALSE, colnames = c("Metric", "Total",
                    "Mean", "Minimum", "Maximum"), filter = "top")
    })

    output$productive_plot <- renderPlotly({
        ggplotly(ggplot(data = itable(), aes(x = repertoire_id, y = unique_productive_sequences)) +
            geom_bar(stat = "identity", position=position_dodge(), fill = "#3182bd") +
            theme_minimal() +
            theme(axis.text.x  = element_text(angle=90, vjust=0.5, hjust = 1), text = element_text(size = 10)) +
            scale_x_discrete(limits = itable()$repertoire_id) +
            labs(x = "", y = "Unique productive sequences"))
    })

    output$top_productive_seq <- DT::renderDataTable({
       LymphoSeq2::topSeqs(atable()) %>%
                    select(repertoire_id, junction_aa, duplicate_count, duplicate_frequency, v_call, d_call, j_call) %>%
                    datatable(rownames = FALSE, colnames = c("Repertoire ID", "Junction AA", "Duplicate Count", "Duplicate Frequency", "V Call", "D Call", "J Call"), filter = "top") 
    })

    output$top_10_seq <- renderPlotly({
        ggplotly(topSeqsPlot(atable(), top = 10) +
            theme(text = element_text(size = 10)) +
            scale_x_discrete(limits = itable()$repertoire_id))
    })

    output$common_seq_table <- DT::renderDataTable({
        top_freq <- LymphoSeq2::topFreq(atable(), frequency = 0.001)
        clone_table <- atable() %>%
                        select(repertoire_id, junction_aa)
        top_freq <- left_join(top_freq, atable(), by = "junction_aa") %>%
                        select(repertoire_id, junction_aa, duplicate_frequency,
                            numberSamples, antigen, prevalence)
        datatable(top_freq, colnames = c("Repertoire ID", "Junction AA",
                "Duplicate Frequency", "Number of Samples", "Antigen",
                "Prevalence"), filter = "top")
    })

    output$clonality <- DT::renderDataTable({
        itable() %>%
            select(repertoire_id, clonality, gini_coefficient, top_productive_sequence) %>% 
            datatable(rownames = FALSE, colnames = c("Repertoire ID", "Clonality", "Gini Coefficient", "Top Productive Sequence"), filter = "top")
    })

    output$clonality_stats <- DT::renderDataTable({
        itable() %>%
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
            mutate(metric = str_to_title(str_replace_all(metric, "_", " "))) %>%
            datatable(rownames = FALSE, colnames = c("Metric", "Total", "Mean", "Minimum", "Maximum"), filter = "top")
    })

    output$clonality_plot <- renderPlotly({
        ggplotly(ggplot(data = itable(), aes(x = repertoire_id, y = clonality)) +
            geom_bar(stat = "identity", position=position_dodge(), fill = "#b2182b") +
            theme_minimal() +
            theme(axis.text.x  = element_text(angle=90, vjust=0.5, hjust = 1), text = element_text(size = 10)) +
            scale_x_discrete(limits = itable()$repertoire_id) +
            labs(x = "", y = "Clonality"))
    })

    sequence_matrix <- reactive({
        atable()
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
        ctable <- LymphoSeq2::cloneTrack(atable())
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
        p + theme(legend.position = "none")
    })

    output$clone_track_tooltip <- renderText(
        alluvial_tooltip(input$clone_track_hover)
    )

    output$common_seq_plot <- renderPlot({
        ttable <- LymphoSeq2::topSeqs(atable(), top = 50)
        ctable <- LymphoSeq2::cloneTrack(ttable)
        p <- LymphoSeq2::plotTrack(clone_table = ctable)
        interactive_alluvial(p)
        p
    })

    output$common_seq_tooltip <- renderText (
        alluvial_tooltip(input$common_seq_hover)
    )

    output$pairwise_sim <- renderPlotly({
        similarity <- LymphoSeq2::scoringMatrix(atable())
        ggplotly(pairwisePlot(matrix = similarity) +
            labs(fill = "Bhattacharyya\n coefficient") +
            theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8)) +
            scale_fill_gradient(low = "#deebf7", high = "#3182bd") +
            theme(legend.position = c(0.8, 0.8)))
    })

    output$public_tcrb <- DT::renderDataTable({
        published <- LymphoSeq2::searchPublished(atable())
        published %>%
            filter(!is.na(PMID)) %>%
            select(repertoire_id, junction_aa, duplicate_count, antigen, prevalence) %>%
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
        heatmaply(vGenes, colors = "RdBu", scale = "row")
    })

    output$gene_freq <- DT::renderDataTable({
        productive_nt <- LymphoSeq2::productiveSeq(study_table = airr_data(), aggregate = "nucleotide")
        LymphoSeq2::geneFreq(productive_nt, locus = "VDJ", family = FALSE)
    })

}

shinyApp(ui = ui, server = server)
