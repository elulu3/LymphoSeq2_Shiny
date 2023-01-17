# Install dependencies if missing
required_packages <- c(
  "shinyjs", "shinycssloaders", "shinythemes", "shinyalert",
  "ggplot2", "plotly", "DT", "heatmaply", "pheatmap",
  "tidyverse", "htmltools", "sp", "writexl", "wordcloud2",
  "remotes"
)
lapply(
  required_packages,
  function(x) if (!require(x, character.only = TRUE)) install.packages(x)
)

if (!require("ggalluvial", character.only = TRUE)) {
  remotes::install_github("corybrunson/ggalluvial@main", build_vignettes = FALSE)
}
if (!require("chorddiag", character.only = TRUE)) {
  remotets::install_github("mattflor/chorddiag", build_vignettes = FALSE)
}
if (!require("LymphoSeq2", character.only = TRUE)) {
    devtools::install_github("WarrenLabFH/LymphoSeq2", ref="v1", build_vignette=FALSE)
}
library(LymphoSeq2)
# library(BiocManager)
# options(repos = BiocManager::repositories())

# ------------------------------------------------------------------------------------------------------------------ # # nolint

# max upload size (current: 2 GB)
options(shiny.maxRequestSize = 2000 * 1024^2)

ui <-
  navbarPage("LymphoSeq2",
    theme = shinythemes::shinytheme("flatly"),
    id = "tabselected",
    tabPanel(
      "Input Data",
      value = "airr_table",
      sidebarLayout(
        sidebarPanel(
          fileInput("airr_files", label = "Upload Files",
                    multiple = TRUE, accept = c(".tsv", ".rda", ".RData")),
          radioButtons("analysis_table",
                        "Select table to view and perform analysis on",
            choices = c("AIRR", "Productive Amino Acid Sequences",
                        "Productive Nucleotide Sequences")
          ),
          checkboxInput("hide_null", label = "Hide empty columns",
                        value = FALSE)
        ),
        mainPanel(
          tabsetPanel(
            tabPanel("AIRR Data Table",
              DT::dataTableOutput("table") %>% withSpinner()
            )
          )
        )
      )
    ),
    navbarMenu(
      "Filtering",
      tabPanel("Top Sequences", value = "top_seq_panel",
        sidebarLayout(
          sidebarPanel(
              numericInput("top_seq_num", "Number of top sequences:", value = NULL, min = 0),
              actionButton("top_seq_button", "Calculate Top Sequences")
          ),
          mainPanel(
            tabsetPanel(
              tabPanel("Top Productive Sequences Table",
                value = "top_seq_table",
                DT::dataTableOutput("top_seq_table") %>% withSpinner()
              ),
              tabPanel("Top Productive Sequences Plot",
                value = "top_seq_plot",
                plotlyOutput("top_seq_plot", height = "600px") %>% withSpinner()
              ),
              id = "top_seq_sub_tab"
            )
          )
        )
      ),
      tabPanel("Unique Sequences", value = "unique_seq_panel",
        # sidebarLayout(
          # sidebarPanel(
          # ),
          # mainPanel(
            tabsetPanel(
              tabPanel("Unique Productive Sequences Table",
                value = "produtive_seq_table",
                DT::dataTableOutput("produtive_seq_table") %>% withSpinner()
              ),
              tabPanel("Unique Productive Sequences Plot",
                value = "produtive_seq_plot",
                plotlyOutput("productive_plot", height = "600px") %>% withSpinner()
              ),
              id = "unique_seq_sub_tab"
            )
          # )
        # )
      )
    ),
    navbarMenu(
      "Repertoire Overlap",
      tabPanel("Common Sequences", value = "common_panel",
        sidebarLayout(
          sidebarPanel(
            conditionalPanel(
              condition = "input.tabselected == 'common_panel' &&
                          input.common_sub_tab == 'common_bar'",
              selectizeInput("bar_id",
                             label = "Select at least 2 repertoire ids",
                             choices = NULL, multiple = TRUE),
              h1("Options to color bar chart:", style = "font-size:15px"),
              selectizeInput("color_rep_id",
                             label = "Select 1 repertoire id to color",
                             choices = NULL, selected = NULL),
              h5("OR", align = "center"),
              selectizeInput("color_intersect",
                label = "Select at least 2 repertoire ids to color intersections",
                choices = NULL, multiple = TRUE, selected = NULL
              ),
              actionButton("bar_button", label = "Create Bar Chart")
            ),
            conditionalPanel(
              condition = "input.tabselected == 'common_panel' &&
                          input.common_sub_tab == 'common_plot'",
              selectizeInput("plot_id",
                label = "Select 2 repertoire ids",
                choices = NULL, multiple = TRUE
              ),
              # radioButtons("show_type", "Show Sequences Type", choices = c("common" = "common", "all" = "all"), inline = TRUE),
              actionButton("plot_button", label = "Create Plot")
            ),
            conditionalPanel(
              condition = "input.tabselected == 'common_panel' && input.common_sub_tab == 'common_venn'",
              selectizeInput("venn_id",
                label = "Select 2 or 3 repertoire ids",
                choices = NULL, multiple = TRUE
              ),
              actionButton("venn_button", label = "Create Venn Diagram")
            )
          ),
          mainPanel(
            tabsetPanel(
              tabPanel("Data Table",
                value = "common_seqs_table",
                DT::dataTableOutput("common_seqs_table") %>% withSpinner()
              ),
              tabPanel("Bar Chart",
                value = "common_bar",
                plotOutput("common_seqs_bar") %>% withSpinner()
              ),
              tabPanel("Plot",
                value = "common_plot",
                plotlyOutput("common_seqs_plot") %>% withSpinner()
              ),
              tabPanel("Venn Diagram",
                value = "common_venn",
                plotOutput("common_seqs_venn") %>% withSpinner()
              ),
              id = "common_sub_tab"
            )
          )
        )
      ),
      tabPanel("Pairwise Similarity", value = "pairwise_sim",
        sidebarLayout(
          sidebarPanel(
            radioButtons("mode", "Select mode to use for calculation:",
              choices = c("Bhattacharyya", "Similarity", "Sorensen", "PSI")
            )
          ),
          mainPanel(
            tabsetPanel(
              tabPanel("Pairwise Similarity",
                value = "pairwise_sim",
                plotlyOutput("pairwise_sim", height = "600px") %>% withSpinner()
              )
            )
          )
        )
      )
    ),
    navbarMenu(
      "Repertoire Diversity Curves",
      tabPanel("Rarefaction Curve", value = "rarefaction_curve",
        tabsetPanel(
          tabPanel("Rarefaction Curve",
            plotlyOutput("rarefaction_curve") %>% withSpinner()
          )
        )
      ),
      tabPanel("Lorenz Curve", value = "lorenz_curve",
        tabsetPanel(
          tabPanel("Lorenz Curve",
            plotlyOutput("lorenz") %>% withSpinner()
          )
        )
      )
    ),
    navbarMenu(
      "Gene Frequency",
      tabPanel("Gene Frequency", value = "gene_panel",
        sidebarLayout(
          sidebarPanel(
            selectizeInput("locus",
              label = "Select which VDJ genes to include",
              choices = c("VDJ", "DJ", "VJ", "DJ", "V", "D", "J")
            ),
            radioButtons("family_bool", "Select which names to use",
              choices = c("family" = TRUE, "gene" = FALSE), inline = TRUE
            ),
            conditionalPanel(
              condition = "input.gene_sub_tab == 'gene_freq_word'",
              selectizeInput("word_id", label = "Select repertoire id", choices = NULL),
              actionButton("word_button", label = "Create Word Cloud")
            )
          ),
          mainPanel(
            tabsetPanel(
              tabPanel("Data Table",
                value = "gene_freq_table",
                DT::dataTableOutput("gene_freq_table") %>% withSpinner()
              ),
              tabPanel("Heat Map",
                value = "gene_freq_heat",
                plotlyOutput("gene_freq_heat", height = "800px") %>% withSpinner()
              ),
              tabPanel("Bar Chart",
                value = "gene_freq_bar",
                plotlyOutput("gene_freq_bar", height = "600px") %>% withSpinner()
              ),
              tabPanel("Word Cloud",
                value = "gene_freq_word",
                wordcloud2Output("gene_freq_word") %>% withSpinner()
              ),
              id = "gene_sub_tab"
            )
          )
        )
      ),
      tabPanel("VDJ Chord Diagram", value = "chord_diagram",
        sidebarLayout(
          sidebarPanel(
            selectizeInput("vdj_association",
              label = "Select VDJ Association",
              choices = c("", "VJ", "DJ")
            ),
            radioButtons("top_seq_chord", "Plot top sequences? (Recommended)",
              choices = c("yes", "no"), inline = TRUE
            ),
            numericInput("top_chord_num", "Number of top sequences to map:", value = 1, min = 1),
            actionButton("chord_button", label = "Create Chord Diagram")
          ),
          mainPanel(
            tabsetPanel(
              tabPanel("VDJ Chord Diagram",
                chorddiagOutput("chord_diagram", width = "100%", height = "600px") %>% withSpinner()
              )
            )
          )
        )
      )
    ),
    tabPanel(
      "Clonality Statistics", value = "clonality_panel",
      sidebarLayout(
        sidebarPanel(
          conditionalPanel(
            condition = "input.clonal_sub_tab == 'clonality'",
            radioButtons("calc_relate", "Calculate clonal relatedness?",
              choices = c("yes", "no"), inline = TRUE
            ),
            numericInput("edit_dis",
              label = "Minimum edit distance the sequence must be less than or equal to:",
              value = 10, min = 0
            ),
            actionButton("clonal_relate_button", "Create Table")
          ),
          conditionalPanel(
            condition = "input.clonal_sub_tab == 'clonality_plot'",
            selectInput("clonal_x",
                        label = "Measurement to graph",
                        choices = c("total_sequences", "unique_productive_sequences",
                                   "total_count", "clonality", "gini_coefficient"),
                        multiple = FALSE)
          )
        ),
        mainPanel(
          tabPanel("Explore Clonality",
            tabsetPanel(
              tabPanel("Clonality",
                value = "clonality",
                DT::dataTableOutput("clonality") %>% withSpinner()
              ),
              tabPanel("Clonality Statistics",
                value = "clonality_stats",
                DT::dataTableOutput("clonality_stats") %>% withSpinner()
              ),
              tabPanel("Clonality Plot",
                value = "clonality_plot",
                plotlyOutput("clonality_plot") %>% withSpinner()
              ),
              tabPanel("Count Statistics",
                value = "count_stats",
                DT::dataTableOutput("stats_count") %>% withSpinner()
              ),
              id = "clonal_sub_tab"
            )
          )
        )
      )
    ),
    tabPanel(
      "Clone Tracking",
      value = "clone_track",
      sidebarLayout(
        sidebarPanel(
          selectizeInput("track_id",
            label = "Select repertoire ids to track",
            choices = NULL, multiple = TRUE
          ),
          radioButtons("top_clones", "Plot top clones?",
            choices = c("yes", "no"), inline = TRUE
          ),
          numericInput("top_num", "Number of top clones to map:",
                       value = 50, min = 1),
          radioButtons("track_aa_choice",
            "Track specific amino acid sequences?",
            choices = c("yes", "no"), selected = "no", inline = TRUE
          ),
          selectizeInput("track_aa",
            label = "Select amino acid sequences to track",
            choices = c(""), multiple = TRUE
          ),
          actionButton("track_button", "Create Plot")
        ),
        mainPanel(
          tabsetPanel(
            tabPanel("Clone Tracking", tags$div(
                style = "position: relative;",
                plotOutput("clone_track",
                  height = "600px",
                  hover = hoverOpts(id = "clone_track_hover")
                ) %>% withSpinner(),
                htmlOutput("clone_track_tooltip")
              )
            )
          )
        )
      )
    ),
    tabPanel(
      "K-mer Distribution",
      value = "kmer_panel",
      sidebarLayout(
        sidebarPanel(
          conditionalPanel(
            condition = "input.kmer_sub_tab == 'count_kmers'",
            numericInput("k_val", "Length of kmer:", value = 5),
            radioButtons("separate_by_rep", "View counts by repertoire?",
              choices = c("yes" = TRUE, "no" = FALSE), inline = TRUE
            ),
            actionButton("kmer_button", "Count Kmers")
          ),
          conditionalPanel(
            condition = "input.kmer_sub_tab == 'kmer_distrib'",
            numericInput("k_val", "Length of kmer:", value = 5),
            numericInput("k_top", "Number of top kmers", value = 10),
            actionButton("kmer_distrib_button", "Plot Distribution")
          )
        ),
        mainPanel(
          tabsetPanel(
            tabPanel("Kmer Counts",
              value = "count_kmers",
              DT::dataTableOutput("count_kmers") %>% withSpinner()
            ),
            tabPanel("Kmer Distribution",
              value = "kmer_distrib",
              plotlyOutput("kmer_distrib", height = "600px") %>% withSpinner()
            ),
            id = "kmer_sub_tab"
          )
        )
      )
    ),
    tabPanel(
      "Differential Abundance", value = "diff_abundance",
      sidebarLayout(
        sidebarPanel(
          selectizeInput("diff_id",
            label = "Select 2 repertoire ids to be compared",
            choices = c(""), multiple = TRUE
          ),
          # radioButtons("diff_junction_type", "Select sequences to use",
          #     choices = c("junction_aa", "junction"), inline = TRUE
          # ),
          numericInput("q_val", "q Value (False Discovery Rate):",
            value = 1, min = 0, max = 1
          ),
          numericInput("zero_val", "Zero Value:", value = 1),
          actionButton("diff_button", "Create Table")
        ),
        mainPanel(
          tabsetPanel(
            tabPanel("Differential Abundance",
              DT::dataTableOutput("diff_abundance") %>% withSpinner())
          )
        )
      )
    ),
    navbarMenu(
      "Public Databases",
      tabPanel("LymphoSeqDB", value = "search_lymphoseqdb",
        sidebarLayout(
          sidebarPanel(
          ),
          mainPanel(
            tabsetPanel(
              tabPanel("Search LymphoSeqDB",
                DT::dataTableOutput("lymphoseqdb_match") %>% withSpinner()
              )
            )
          )
        )
      ),
      tabPanel("iReceptor", value = "search_ireceptor",
        sidebarLayout(
          sidebarPanel(
            radioButtons("search_table", "Select table to search database with",
              choices = c("Published Sequences from LymphoSeqDB", "Top Sequences")
            ),
            numericInput("search_top", "Top Number of Sequences:", value = 10, min = 1),
            actionButton("search_button", label = "Search Database")
          ),
          mainPanel(
            tabsetPanel(
              tabPanel("Search iReceptor",
                DT::dataTableOutput("ireceptor_match") %>% withSpinner()
              )
            )
          )
        )
      )
    ),
    tabPanel(
      "About",
      fluidPage(
        tags$h2("Full Documentation for LymphoSeq2_Shiny app here:",
                style = "font-size:30px"),
        tags$br(),
        tags$a(href = "https://elulu3.github.io/LymphoSeq2_Shiny/",
               "LymphoSeq2_Shiny Documentation", style = "font-size:20px"),
        tags$h2("More about the LymphoSeq2 package here:", style = "font-size:30px"),
        tags$br(),
        tags$a(href = "https://shashidhar22.github.io/LymphoSeq2/",
               "LymphoSeq2 Package Website", style = "font-size:20px"),
        tags$br(), tags$br(),
        tags$a(href = "https://github.com/shashidhar22/LymphoSeq2",
               "LymphoSeq2 on Github", style = "font-size:20px")
      )
    ),
    navbarMenu(
      "Download",
      tabPanel(list(
        downloadButton("rda", label = "rda", class = "butt"),
        downloadButton("tsv", label = "tsv", class = "butt"),
        downloadButton("excel", label = "excel", class = "butt"),
        downloadButton("pdf", label = "pdf", class = "butt"),
        tags$head(tags$style(".butt{background:#add8e6;}
                              .butt{color: #337ab7;}"))

      ))
    ),
    shinyjs::useShinyjs()
  )

"%then%" <- function(a, b) {
  if (is.null(a)) b else a
}

# Code for alluvial plots interactive
interactive_alluvial <- function(p) {
  alluvium_width <- 1 / 3
  pbuilt <<- ggplot_build(p)
  data_draw <- transform(pbuilt$data[[3]], width = 1 / 3)
  groups_to_draw <<- split(data_draw, data_draw$group)
  group_xsplines <- lapply(
    groups_to_draw,
    data_to_alluvium
  )
  xspline_coords <- lapply(
    group_xsplines,
    function(coords) {
      grid::xsplineGrob(
        x = coords$x,
        y = coords$y,
        shape = coords$shape,
        open = FALSE
      )
    }
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
  xrange_new <- c(1 - alluvium_width / 2, max(pbuilt$data[[3]]$x)
  + alluvium_width / 2)
  yrange_new <- c(0, sum(pbuilt$data[[3]]$count[pbuilt$data[[3]]$x == 1]))
  new_range_transform <- function(x_old, range_old, range_new) {
    (x_old - range_old[1]) / (range_old[2] - range_old[1]) *
      (range_new[2] - range_new[1]) + range_new[1]
  }
  polygon_coords <<- lapply(xspline_points, function(pts) {
    x_trans <- new_range_transform(
      x_old = as.numeric(pts$x),
      range_old = xrange_old,
      range_new = xrange_new
    )
    y_trans <- new_range_transform(
      x_old = as.numeric(pts$y),
      range_old = yrange_old,
      range_new = yrange_new
    )
    list(x = x_trans, y = y_trans)
  })
}

# Generates tooltip for alluvial plots
alluvial_tooltip <- function(plot_hover) {
  node_width <- 1 / 3
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
          node_label, # tags$br(),
          # "n =", node_n,
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
        function(pol) {
          point.in.polygon(
            point.x = hover$x,
            point.y = hover$y,
            pol.x = pol$x,
            pol.y = pol$y
          )
        }
      )
      if (any(hover_within_flow)) {
        coord_id <- rev(which(hover_within_flow == 1))[1]
        flow_label <- paste(groups_to_draw[[coord_id]]$stratum,
                            collapse = " -> ")
        flow_n <- groups_to_draw[[coord_id]]$count[1]
        offset <- 5

        renderTags(
          tags$div(
            flow_label, # tags$br(),
            # "n =", flow_n,
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

show_data_formats <- function() {
  shinyjs::show("rda")
  shinyjs::show("tsv")
  shinyjs::show("excel")
  shinyjs::hide("pdf")
}

show_img_formats <- function() {
  shinyjs::hide("tsv")
  shinyjs::hide("excel")
  shinyjs::show("rda")
  shinyjs::show("pdf")
}

check_file_ext <- function(files, ext) {
  check <- lapply(files, function(x) tools::file_ext(x) == ext)
  check <- unlist(check)
  if (all(check)) {
    return(TRUE)
  }
  return(FALSE)
}

# Color blind friendly palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# ------------------------------------------------------------------------------------------------------------------ # # nolint

server <- function(input, output, session) {
  shinyjs::disable("top_num")
  shinyjs::disable("track_aa")
  shinyjs::disable("top_chord_num")
  shinyjs::disable("edit_dis")
  shinyjs::disable("search_top")
  shinyjs::hide("pdf")
  rda_envir <<- NULL

  # Standardizes input files into AIRR-compliant format
  airr_data <- reactive({
    validate(
      need(
        !is.null(input$airr_files),
        "Please select files to upload to render output."
      )
    )
    tryCatch(
      {
        if (check_file_ext(input$airr_files$name, "tsv")) {
          rda_envir <<- NULL
          table <- LymphoSeq2::readImmunoSeq(input$airr_files$datapath)
          in_files <- lapply(c(input$airr_files$name),
                             function(i) substr(i, 1, stringr::str_length(i) - 4))
          in_files <- unlist(in_files)
          table <- table %>%
            dplyr::mutate(repertoire_id =
                          if_else(as.integer(repertoire_id) <= length(in_files),
              in_files[as.integer(repertoire_id) + 1], "unknown"
            ))
          table
        } else if (check_file_ext(input$airr_files$name, "rda")) {
          rda_envir <<- new.env()
          name <- load(input$airr_files$datapath, envir = rda_envir)
          rda_envir$study_table
        }
      },
      error = function(e) {
        shinyalert::shinyalert("Cannot read data format.",
          "Please upload other files.",
          type = "error"
        )
      }
    )
  })

  # Filters for productive amino acid sequences
  productive_aa <- reactive({
    validate(
      need(
        !is.null(input$airr_files),
        "Please select files to upload to render output."
      )
    )
    if (check_file_ext(input$airr_files$name, "tsv")) {
      LymphoSeq2::productiveSeq(
        study_table = airr_data(),
        aggregate = "junction_aa"
      )
    } else if (input$airr_files$name == "amino_table.rda") {
      rda_envir <<- new.env()
      name <- load(input$airr_files$datapath, envir = rda_envir)
      rda_envir$amino_table
    } else {
      rda_envir$amino_table
    }
  })

  # Filters for productive nucleotide sequences
  productive_nt <- reactive({
    validate(
      need(
        !is.null(input$airr_files),
        "Please select files to upload to render output."
      )
    )
    if (check_file_ext(input$airr_files$name, "tsv")) {
      LymphoSeq2::productiveSeq(
        study_table = airr_data(),
        aggregate = "junction"
      )
    } else if (check_file_ext(input$airr_files$name, "rda")) {
      if (input$airr_files$name == "nucleotide_table.rda") {
        rda_envir <<- new.env()
        name <- load(input$airr_files$datapath, envir = rda_envir)
        rda_envir$nucleotide_table
      } else {
        rda_envir$nucleotide_table
      }
    } else {
      NULL
    }
  })

  # Computes the summary statistics table
  clonality_data <- reactive({
    validate(
      need(
        !is.null(input$airr_files),
        "Please select files to upload to render output."
      )
    )
    if (check_file_ext(input$airr_files$name, "tsv")) {
      LymphoSeq2::clonality(study_table = airr_data())
    } else if (input$airr_files$name == "summary_table.rda") {
      rda_envir <<- new.env()
      name <- load(input$airr_files$datapath, envir = rda_envir)
      rda_envir$summary_table
    } else {
      rda_envir$summary_table
    }
  })

  stable <- reactive({
    validate(
      need(
        !is.null(input$airr_files),
        "Please select files to upload to render output."
      )
    )
    # if (check_file_ext(input$airr_files$name, "tsv")) {
    #   if (input$analysis_table == "AIRR") {
    #     airr_data()
    #   } else if (input$analysis_table == "Productive Amino Acid Sequences") {
    #     productive_aa()
    #   } else if (input$analysis_table == "Productive Nucleotide Sequences") {
    #     productive_nt()
    #   }
    # } else 
    # if (check_file_ext(input$airr_files$name, "rda")) {
    #   if (input$airr_files$name == "amino_table.rda") {
    #     productive_aa()
    #   } else if (input$airr_files$name == "nucleotide_table.rda") {
    #     productive_nt()
    #   } else if (input$airr_files$name == "summary_table.rda") {
    #     clonality_data()
    #   } else  if (input$airr_files$name == "study_table.rda"){
    #     airr_data()
    #   }
    # } else 
    if (input$analysis_table == "AIRR") {
        airr_data()
    } else if (input$analysis_table == "Productive Amino Acid Sequences") {
      productive_aa()
    } else if (input$analysis_table == "Productive Nucleotide Sequences") {
      productive_nt()
    }
  })

  # ------------------------------------------------------------------------------------------------------------------ # # nolint

  # Grabs all unqiue repertoire ids.
  # Intended to be used for drop down selections.
  unique_prod_rep <- reactive({
    unique(productive_aa()[, "repertoire_id"])
  })

  # Grabs list of unique productive amino acid sequences.
  # Intended to be used for drop down selections.
  unique_prod_aa <- reactive({
    unique_aa <- unique(productive_aa()[, "junction_aa"])
    df <- data.frame(value = unique_aa, label = unique_aa)
    colnames(df) <- c("value", "label")
    df
  })

  # When files are uploaded, the following should occur:
  #   - all plots should be cleared
  #   - all drop down selections should be updated to reflect uploaded data
  observeEvent(input$airr_files, {
    # if (check_file_ext(input$airr_files$name, "tsv")) {
    #   shinyjs::enable("analysis_table")
    # } 
    # else if (check_file_ext(input$airr_files$name, "rda")) {
    #   shinyjs::disable("analysis_table")
    # }
    stable()
    if (input$chord_button || input$venn_button || input$bar_button ||
      input$plot_button || input$diff_button || input$track_button ||
      input$clonal_relate_button || input$word_button || input$kmer_button ||
      input$top_seq_button) {
      purrr::map(
        c("chord_diagram", "common_seqs_venn", "common_seqs_bar",
          "common_seqs_plot", "diff_abundance", "gene_freq_word",
          "clonality", "clone_track",
          "top_seq_table", "top_seq_plot",
          "count_kmers", "kmer_distrib"),
        function(x) shinyjs::hide(x)
      )
    }
  })

  observeEvent(input$bar_id, {
    shiny::updateSelectizeInput(session, "color_intersect",
      choices = input$bar_id
    )
  })

  # Update download choices appropriately when tabs are clicked
  observeEvent(input$tabselected, {
    if (input$tabselected %in% c("airr_table", "diff_abundance",
                                 "search_lymphoseqdb", "search_ireceptor") ||
        input$tabselected == "common_panel" && input$common_sub_tab == "common_seqs_table" ||
        input$tabselected == "top_seq_panel" && input$top_seq_sub_tab == "top_seq_table" ||
        input$tabselected == "unique_seq_panel" && input$unique_seq_sub_tab == "produtive_seq_table" ||
        input$tabselected == "clonality_panel" && input$clonal_sub_tab != "clonality_plot" ||
        input$tabselected == "gene_panel" && input$gene_sub_tab == "gene_freq_table" ||
        input$tabselected == "kmer_panel" && input$kmer_sub_tab == "count_kmers") {
      show_data_formats()
    } else if (input$tabselected == "About") {
      shinyjs::hide("rda")
      shinyjs::hide("tsv")
      shinyjs::hide("excel")
      shinyjs::hide("pdf")
    } else {
      show_img_formats()
    }

    if (input$tabselected %in% c("common_panel", "diff_abundance", "clone_track")) {
      purrr::map(
        c(
          "common_table_id", "bar_id", "plot_id",
          "venn_id", "diff_id", "track_id"
        ),
        function(x) {
          shiny::updateSelectizeInput(session,
            inputId = x,
            choices = unique_prod_rep(),
            selected = ""
          )
        }
      )
      shiny::updateSelectizeInput(session,
        inputId = "color_rep_id",
        choices = c("none", unique_prod_rep())
      )
      shiny::updateSelectizeInput(session,
        inputId = "color_intersect",
        choices = c("")
      )
    }
    if (input$tabselected == "gene_panel") {
      nt_rep_id <- unique(productive_nt()[, "repertoire_id"])
      shiny::updateSelectizeInput(session, "word_id", choices = nt_rep_id)
    }
  })

  # Logic for updating repertoire_id selections for "Explore Repertoire Overlap"
  update_common_tabs <- reactive({
    purrr::map(
      c("common_table_id", "bar_id", "color_intersect", "plot_id", "venn_id"),
      function(x) {
        shiny::updateSelectizeInput(session, x,
          choices = unique_prod_rep()
        )
      }
    )
    shiny::updateSelectizeInput(session, "color_rep_id",
      choices = c("none", unique_prod_rep())
    )
  })

  # When "Explore Repertoire Overlap" tab is clicked, the following occurs:
  #   - download options are updated for each sub tabs
  #   - repertoire_id selections are updated
  observeEvent(input$common_sub_tab, {
    if (input$common_sub_tab == "common_seqs_table") {
      show_data_formats()
    } else {
      show_img_formats()
    }
    update_common_tabs()
  })

  # Logic for updating repertoire_id selections for "Explore Gene Frequencies"
  update_gene_tabs <- reactive({
    shiny::updateSelectizeInput(session, "word_id",
      choices = unique(productive_nt()[, "repertoire_id"])
    )
  })

  observeEvent(input$gene_sub_tab, {
    if (input$gene_sub_tab == "gene_freq_table") {
      show_data_formats()
    } else {
      show_img_formats()
    }
    update_gene_tabs()
  })

  observeEvent(input$clonal_sub_tab, {
    if (input$clonal_sub_tab == "clonality_plot") {
      show_img_formats()
    } else {
      show_data_formats()
    }
  })

  observeEvent(input$top_seq_sub_tab, {
    if (input$top_seq_sub_tab == "top_seq_table") {
      show_data_formats()
    } else {
      show_img_formats()
    }
  })

  observeEvent(input$unique_seq_sub_tab, {
    if (input$unique_seq_sub_tab == "produtive_seq_table") {
      show_data_formats()
    } else {
      show_img_formats()
    }
  })

  observeEvent(input$kmer_sub_tab, {
    if (input$kmer_sub_tab == "count_kmers") {
      show_data_formats()
    } else {
      show_img_formats()
    }
  })

  # ------------------------------------------------------------------------------------------------------------------ # # nolint

  # Shows uploaded data in "Input Table" tab
  output$table <- DT::renderDataTable({
    if (input$hide_null) {
      stable() %>%
        purrr::discard(~ all(is.na(.) | . == "")) %>%
        DT::datatable(filter = "top", options = list(scrollX = TRUE))
    } else {
      stable() %>%
        DT::datatable(filter = "top", options = list(scrollX = TRUE))
    }
  })

  observeEvent(input$chord_button, {
    shinyjs::show("chord_diagram")
  })

  # Disable input selection if user selects "no"
  observeEvent(input$top_seq_chord, {
    if (input$top_seq_chord == "no") {
      shinyjs::disable("top_chord_num")
    } else {
      shinyjs::enable("top_chord_num")
    }
  })

  # Compute data table to create VDJ chord diagram
  chord_data <- reactive({
    if (input$top_seq_chord == "yes") {
      validate(
        need(input$top_chord_num > 0, "Please enter number of top sequences.")
      )
      chord_table <- LymphoSeq2::topSeqs(stable(), input$top_chord_num)
    } else {
      chord_table <- stable()
    }
    chord_table <- as.data.frame(apply(chord_table, 2, function(x) gsub("TCRB", "", x)))
    if (input$vdj_association == "VJ") {
      vj <- chord_table %>%
        dplyr::select(v_family, j_family) %>%
        dplyr::mutate(
          v_family = tidyr::replace_na(v_family, "Unresolved"),
          j_family = tidyr::replace_na(j_family, "Unresolved")
        )
      vj <- vj %>%
        dplyr::group_by(v_family, j_family) %>%
        dplyr::summarize(duplicate_count = n(), .groups = "drop") %>%
        tidyr::pivot_wider(
          id_cols = v_family,
          names_from = j_family,
          values_from = duplicate_count
        )
      row_names <- vj$v_family
      vj <- vj %>%
        dplyr::select(-v_family)
      vj[is.na(vj)] <- 0
      vj <- as.matrix(vj)
      rownames(vj) <- row_names
      vj
    } else if (input$vdj_association == "DJ") {
      dj <- chord_table %>%
        dplyr::select(d_family, j_family) %>%
        dplyr::mutate(
          d_family = replace_na(d_family, "Unresolved"),
          j_family = replace_na(j_family, "Unresolved")
        )
      dj <- dj %>%
        dplyr::group_by(d_family, j_family) %>%
        dplyr::summarize(duplicate_count = n(), .groups = "drop") %>%
        tidyr::pivot_wider(
          id_cols = d_family,
          names_from = j_family,
          values_from = duplicate_count
        )
      row_names <- dj$d_family
      dj <- dj %>%
        dplyr::select(-d_family)
      dj[is.na(dj)] <- 0
      dj <- as.matrix(dj)
      rownames(dj) <- row_names
      dj
    }
  })

  # Outputs an interactive chord diagram with caching
  output$chord_diagram <- renderChorddiag({
    input$chord_button
    isolate({
      validate(
        need(input$vdj_association != "", "Please select VDJ Association")
      )
      chorddiag::chorddiag(chord_data(),
        type = "bipartite",
        showTicks = FALSE,
        groupColors = cbPalette,
        groupnameFontsize = 15,
        groupnamePadding = 10, margin = 90
      )
    })
  }) %>%
    shiny::bindCache(chord_data()) %>%
    shiny::bindEvent(input$chord_button)

  # Computes data table for common sequences
  common_table_data <- reactive({
    if (is.null(productive_aa())) {
      shinyalert::shinyalert("Incorrect input data.",
        "Can only perform function on productive amino acid sequences.",
        type = "error"
      )
      validate(
        need(
          !is.null(productive_aa()),
          "Can only perform function on productive amino acid sequences."
        )
      )
      return()
    }
    LymphoSeq2::commonSeqs(study_table = productive_aa())
  })

  # Outputs common sequences data table
  output$common_seqs_table <- DT::renderDataTable({
    common_table_data() %>%
      DT::datatable(filter = "top", options = list(scrollX = TRUE))
  })

  observeEvent(input$bar_button, {
    shinyjs::show("common_seqs_bar")
  })

  # Computes data table for UpSetR bar plot of common sequences
  common_bar_data <- reactive({
    color_rep_id <- input$color_rep_id
    color_intersect <- input$color_intersect
    if (input$color_rep_id == "none") {
      color_rep_id <- NULL
    }
    if (!is.null(input$color_intersect) & "none" %in% input$color_intersect) {
      color_intersect <- NULL
    }
    LymphoSeq2::commonSeqsBar(productive_aa(),
      input$bar_id, color_rep_id, c(color_intersect),
      labels = "yes"
    )
  })

  # Outputs UpSetR bar plot of common sequences
  # Error message will show when there are no common sequences
  output$common_seqs_bar <- renderPlot({
    input$bar_button
    isolate({
      validate(
        need(
          !is.null(productive_aa()),
          "Can only perform function on productive amino acid sequences."
        )
        %then%
          need(
            length(input$bar_id) > 1,
            "Please select at least 2 repertoire ids"
          )
          %then%
          need(
            length(input$color_rep_id) == 1,
            "Please select only 1 repertoire id"
          )
          %then%
          need(
            is.null(input$color_intersect) ||
              input$color_intersect == "none" ||
              length(input$color_intersect) > 1,
            "Please select at least 2 repertoire ids"
          )
      )
      tryCatch(
        {
          common_bar_data()
        },
        error = function(e) {
          shinyalert::shinyalert("No Common Sequences.",
            "Select different repertoire ids.",
            type = "error"
          )
        }
      )
    })
  })

  observeEvent(input$plot_button, {
    shinyjs::show("common_seqs_plot")
  })

  # Creates scatterplot of common sequences
  common_seqs_plot <- reactive({
    LymphoSeq2::commonSeqsPlot(
      input$plot_id[1],
      input$plot_id[2],
      productive_aa()
    )
  })

  # Outputs interactive scatterplot of common sequences
  # Error message will show when there are no common sequences
  output$common_seqs_plot <- renderPlotly({
    input$plot_button
    isolate({
      validate(
        need(
          !is.null(productive_aa()),
          "Can only perform function on productive amino acid sequences."
        )
        %then%
          need(
            length(input$plot_id) == 2,
            "Please select exactly 2 repertoire ids"
          )
      )
      if (nrow(common_seqs_plot()$data) > 0) {
        common_seqs_plot()
      } else {
        shinyalert::shinyalert("No Common Sequences!",
          "Select different repertoire ids.",
          type = "error"
        )
        return()
      }
    })
  })

  observeEvent(input$venn_button, {
    shinyjs::show("common_seqs_venn")
  })

  # Creates Venn diagram of common sequences
  common_venn_data <- reactive({
    if (length(input$venn_id) == 2) {
      a <- productive_aa() %>%
        dplyr::filter(repertoire_id == input$venn_id[[1]])
      b <- productive_aa() %>%
        dplyr::filter(repertoire_id == input$venn_id[[2]])
      grid::grid.newpage()
      venn <- VennDiagram::draw.pairwise.venn(
        area1 = length(a$junction_aa),
        area2 = length(b$junction_aa),
        cross.area = length(intersect(a$junction_aa, b$junction_aa)),
        category = c(input$venn_id[1], input$venn_id[2]),
        cat.fontfamily = rep("sans", 2),
        fontfamily = rep("sans", 3),
        fill = c("#3288bd", "#d53e4f"),
        cat.pos = c(0, 0),
        cat.dist = rep(0.025, 2),
        cex = 1,
        cat.cex = 0.7,
        lwd = rep(2, 2)
      )
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
      venn <- VennDiagram::draw.triple.venn(
        area1 = length(a$junction_aa),
        area2 = length(b$junction_aa),
        area3 = length(c$junction_aa),
        n12 = length(intersect(a$junction_aa, b$junction_aa)),
        n23 = length(intersect(b$junction_aa, c$junction_aa)),
        n13 = length(intersect(a$junction_aa, c$junction_aa)),
        n123 = length(Reduce(
          intersect,
          list(
            a$junction_aa,
            b$junction_aa,
            c$junction_aa
          )
        )),
        category = c(
          input$venn_id[1], input$venn_id[2],
          input$venn_id[3]
        ),
        cat.fontfamily = rep("sans", 3),
        fontfamily = rep("sans", 7),
        fill = c("#3288bd", "#abdda4", "#d53e4f"),
        cat.pos = c(0, 0, 180),
        cat.dist = rep(0.025, 3),
        cex = 1,
        cat.cex = 0.7,
        lwd = rep(2, 3)
      )
      data_output <- venn
    }
    data_output
  })

  # Outputs a Venn diagram of common sequences
  output$common_seqs_venn <- renderPlot({
    validate(
      need(
        !is.null(productive_aa()),
        "Can only perform function on productive amino acid sequences."
      )
    )
    input$venn_button
    isolate({
      validate(
        need(
          !is.null(productive_aa()),
          "Can only perform function on productive amino acid sequences."
        )
        %then%
          need(
            length(input$venn_id) == 2 | length(input$venn_id) == 3,
            "Please select only 2 or 3 repertoire ids"
          )
      )
      grid::grid.draw(common_venn_data())
    })
  })

  # Creates a pairwise comparison heat map of repertoires
  pairwise_sim_data <- reactive({
    validate(
      need(
        !is.null(productive_aa()),
        "Can only perform function on productive amino acid sequences."
      )
    )
    similarity <- LymphoSeq2::scoringMatrix(productive_aa(), mode = input$mode)
    LymphoSeq2::pairwisePlot(matrix = similarity) +
      labs(fill = paste(input$mode, "\ncoefficient")) +
      theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8)) +
      scale_fill_gradient(low = "#deebf7", high = "#3182bd") +
      theme(legend.position = c(0.8, 0.8))
  })

  # Ouputs an interactive pairwise comparison heat map with caching
  output$pairwise_sim <- renderPlotly({
    plotly::ggplotly(pairwise_sim_data())
  }) %>%
    shiny::bindCache(productive_aa(), input$mode)

  # Creates a Lorenz curve
  lorenz_data <- reactive({
    repertoire_ids <- stable() %>%
      dplyr::pull(repertoire_id) %>%
      unique()
    LymphoSeq2::lorenzCurve(repertoire_ids, stable()) +
      ggplot2::coord_fixed(1 / 2)
  })

  # Ouputs an interactive line plot of the Lorenz curve
  output$lorenz <- renderPlotly({
    lorenz_data()
  })

  rarefaction_data <- reactive({
    LymphoSeq2::plotRarefactionCurve(airr_data())
  })

  # Outputs an interative line plot of rarefactionan curves for repertoires
  output$rarefaction_curve <- renderPlotly({
    rarefaction_data()
  })

  stats_count_data <- reactive({
    clonality_data() %>%
      dplyr::select(
        repertoire_id, total_sequences,
        unique_productive_sequences, total_count
      ) %>%
      tidyr::pivot_longer(!repertoire_id,
        names_to = "metric",
        values_to = "values"
      ) %>%
      dplyr::group_by(metric) %>%
      dplyr::summarize(
        Total = sum(values),
        Mean = mean(values),
        Minimum = min(values),
        Maximum = max(values)
      ) %>%
      tidyr::drop_na() %>%
      dplyr::ungroup() %>%
      dplyr::mutate(metric = str_to_title(str_replace_all(metric, "_", " ")))
  })

  output$stats_count <- DT::renderDataTable({
    stats_count_data() %>%
      DT::datatable(
        rownames = FALSE,
        colnames = c(
          "Metric", "Total",
          "Mean", "Minimum", "Maximum"
        ),
        filter = "top",
        options = list(scrollX = TRUE)
      )
  })

  prod_plot_data <- reactive({
    ggplot2::ggplot(data = clonality_data(), aes(x = repertoire_id, y = unique_productive_sequences)) +
      geom_bar(stat = "identity", position = position_dodge(), fill = "#3182bd") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), text = element_text(size = 10)) +
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
    # TODO: check stable is not complete data (must be from productive sequences)
    LymphoSeq2::topSeqs(stable(), input$top_seq_num) %>%
      dplyr::select(
        repertoire_id, junction_aa, duplicate_count, duplicate_frequency,
        v_call, d_call, j_call
      )
  })

  output$top_seq_table <- DT::renderDataTable({
    input$top_seq_button
    isolate({
      validate(
        need(input$top_seq_num > 0, "Please choose number of top sequences.")
      )
      top_prod_seq_data() %>%
        DT::datatable(
          rownames = FALSE, filter = "top",
          options = list(scrollX = TRUE)
        )
    })
  })

  top_seq_data <- reactive({
    LymphoSeq2::topSeqsPlot(stable(), top = input$top_seq_num)
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

  observeEvent(input$calc_relate, {
    if (input$calc_relate == "no") {
      shinyjs::disable("edit_dis")
    } else {
      shinyjs::enable("edit_dis")
    }
  })

  clone_data <- reactive({
    clone_data_copy <- clonality_data()
    if (input$calc_relate == "yes") {
      clone_relate_data <- LymphoSeq2::clonalRelatedness(stable(), editDistance = input$edit_dis) %>%
        dplyr::select(clonalRelatedness)
      clone_data_copy$clonal_relatedness <- clone_relate_data
    }
    clone_data_copy
  })

  observeEvent(input$clonal_relate_button, {
    shinyjs::show("clonality")
  })

  output$clonality <- DT::renderDataTable({
    input$clonal_relate_button
    isolate({
      clone_data() %>%
        DT::datatable(
          rownames = FALSE, filter = "top",
          options = list(scrollX = TRUE)
        )
    })
  })

  clone_stats_data <- reactive({
    clonality_data() %>%
      dplyr::select(repertoire_id, clonality, gini_coefficient, top_productive_sequence) %>%
      tidyr::pivot_longer(!repertoire_id,
        names_to = "metric",
        values_to = "values"
      ) %>%
      dplyr::group_by(metric) %>%
      dplyr::summarize(
        Mean = mean(values),
        Minimum = min(values),
        Maximum = max(values)
      ) %>%
      tidyr::drop_na() %>%
      dplyr::ungroup() %>%
      dplyr::mutate(metric = str_to_title(str_replace_all(metric, "_", " ")))
  })

  output$clonality_stats <- DT::renderDataTable({
    clone_stats_data() %>%
      DT::datatable(
        rownames = FALSE,
        colnames = c("Metric", "Mean", "Minimum", "Maximum"),
        filter = "top",
        options = list(scrollX = TRUE)
      )
  })

  clone_plot_data <- reactive({
    ggplot2::ggplot(data = clonality_data(),
                    aes(x = repertoire_id, y = .data[[input$clonal_x]],
                        fill = repertoire_id)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        text = element_text(size = 10)
      ) +
      scale_x_discrete(limits = clonality_data()$repertoire_id) +
      labs(x = "", y = input$clonal_x)
    # data <- clonality_data() %>% dplyr::mutate(repertoire_id = "TRB")
    # ggplot2::ggplot(data = data,
    #                 aes(x = repertoire_id, y = clonality)) +
    #   geom_boxplot() +
    #   scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
    #   geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
    #   theme(axis.text.x = element_text(angle = 90))
  })

  output$clonality_plot <- renderPlotly({
    plotly::ggplotly(clone_plot_data())
  })

  observeEvent(input$top_clones, {
    if (input$top_clones == "no") {
      shinyjs::disable("top_num")
    } else {
      shinyjs::enable("top_num")
    }
  })

  observeEvent(input$track_aa_choice, {
    if (input$track_aa_choice == "no") {
      shinyjs::disable("track_aa")
    } else {
      shiny::updateSelectizeInput(session,
        inputId = "track_aa",
        choices = unique_prod_aa(), server = TRUE
      )
      shinyjs::enable("track_aa")
    }
  })

  observeEvent(input$track_button, {
    shinyjs::show("clone_track")
  })

  clone_track_data <- reactive({
    clone_track_table <- productive_aa()
    aa_list <- NULL
    if (input$top_clones == "yes") {
      validate(
        need(input$top_num > 0, "Please enter a number for tracking top clones")
      )
      clone_track_table <- LymphoSeq2::topSeqs(productive_aa(), top = input$top_num)
    }
    if (input$track_aa_choice == "yes") {
      validate(
        need(length(input$track_aa) > 0, "Please select amino acid sequences to track")
      )
      aa_list <- input$track_aa
    }
    ctable <- LymphoSeq2::cloneTrack(clone_track_table,
      sample_list = input$track_id,
      sequence_track = aa_list
    )
    LymphoSeq2::plotTrack(clone_table = ctable)
    # + geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
    # + theme(legend.position = "bottom")
  })

  output$clone_track <- renderPlot({
    input$track_button
    isolate({
      validate(
        need(length(input$track_id) > 0, "Please select repertoire ids to track")
      )
      tryCatch(
        {
          interactive_alluvial(clone_track_data())
          clone_track_data()
        },
        error = function(e) {
          shinyalert::shinyalert("Amino acid sequences not found",
            "Please select other amino acid sequences
                                            or increase number of top clones to plot.",
            type = "error"
          )
        }
      )
    })
  })

  output$clone_track_tooltip <- renderText(
    alluvial_tooltip(input$clone_track_hover)
  )

  lymphoseqdb_data <- reactive({
    published <- LymphoSeq2::searchPublished(stable())
    published %>%
      dplyr::filter(!is.na(PMID))
  })

  output$lymphoseqdb_match <- DT::renderDataTable({
    lymphoseqdb_data() %>%
      DT::datatable(filter = "top", options = list(scrollX = TRUE))
  })

  observeEvent(input$search_table, {
    if (input$search_table == "Top Sequences") {
      shinyjs::enable("search_top")
    } else {
      shinyjs::disable("search_top")
    }
  })

  ireceptor_data <- reactive({
    if (input$search_table == "Published Sequences from LymphoSeqDB") {
      LymphoSeq2::searchDB(lymphoseqdb_data())
    } else {
      top_seqs <- LymphoSeq2::topSeqs(productive_table = stable(),
                                      top = input$search_top)
      LymphoSeq2::searchDB(top_seqs)
    }
  })

  output$ireceptor_match <- DT::renderDataTable({
    input$search_button
    isolate({
      ireceptor_data() %>%
        DT::datatable(filter = "top", options = list(scrollX = TRUE))
    })
  })

  gene_heatmap_data <- reactive({
    genes <- LymphoSeq2::geneFreq(productive_nt(),
                                  locus = input$locus,
                                  family = input$family_bool) %>%
      tidyr::pivot_wider(
        id_cols = gene_name,
        names_from = repertoire_id,
        values_from = gene_frequency,
        values_fn = sum,
        values_fill = 0
      )
    gene_names <- genes %>%
      dplyr::pull(gene_name)
    genes <- genes %>%
      dplyr::select(-gene_name) %>%
      as.matrix()
    rownames(genes) <- gene_names
    genes
  })

  output$gene_freq_heat <- renderPlotly({
    heatmaply::heatmaply(gene_heatmap_data(), scale = "row", fontsize_row = 6)
  }) %>%
    shiny::bindCache(gene_heatmap_data())

  gene_bar_data <- reactive({
    genes <- LymphoSeq2::geneFreq(productive_nt(), input$locus, input$family_bool)
    ggplot2::ggplot(genes, aes(x = repertoire_id, y = gene_frequency, fill = gene_name)) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      scale_y_continuous(expand = c(0, 0)) +
      guides(fill = guide_legend(ncol = 3)) +
      labs(y = "Frequency (%)", x = "repertoire_id", fill = "") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  })

  output$gene_freq_bar <- renderPlotly({
    plotly::ggplotly(gene_bar_data())
  }) %>%
    shiny::bindCache(gene_bar_data())

  observeEvent(input$word_button, {
    shinyjs::show("gene_freq_word")
  })

  gene_word_data <- reactive({
    genes <- LymphoSeq2::geneFreq(productive_nt(), input$locus, input$family_bool)
    gene_data <- genes %>%
      dplyr::filter(repertoire_id == input$word_id)
    gene_data <- data.frame(gene_data$gene_name, gene_data$gene_frequency)
    wordcloud2::wordcloud2(gene_data, color = "random-dark", minSize = 10)
  })

  output$gene_freq_word <- renderWordcloud2({
    input$word_button
    isolate({
      gene_word_data()
    })
  }) %>%
    shiny::bindCache(gene_word_data()) %>%
    shiny::bindEvent(input$word_button)

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
    # if (input$diff_junction_type == "junction_aa") {
    LymphoSeq2::differentialAbundance(study_table = productive_aa(),
                                      repertoire_ids = input$diff_id,
                                      type = "junction_aa",
                                      q = input$q_val,
                                      zero = input$zero_val
    )
    # } else if (input$diff_junction_type == "junction") {
    #   LymphoSeq2::differentialAbundance(study_table = productive_nt(),
    #                                     repertoire_ids = input$diff_id,
    #                                     type = "junction",
    #                                     q = input$q_val,
    #                                     zero = input$zero_val
    #   )
    # }

  })

  output$diff_abundance <- DT::renderDataTable({
    input$diff_button
    isolate({
      validate(
        need(length(input$diff_id) == 2, "Please select 2 repertoire ids")
      )
      diff_table_data() %>%
        DT::datatable(filter = "top", options = list(scrollX = TRUE))
    })
  })

  kmer_table_data <- reactive({
    LymphoSeq2::countKmer(productive_nt(), input$k_val, input$separate_by_rep)
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

  output$kmer_distrib <- renderPlotly({
    input$kmer_distrib_button
    isolate({
      validate(
        need(input$k_val > 1, "Please input length of at least 2")
        %then%
          need(input$k_top > 0, "Please enter how many repertoires to visualize")
      )
    })
    shiny::updateRadioButtons(session, "separate_by_rep", selected = "yes")
    LymphoSeq2::kmerPlot(kmer_table_data(), input$k_top)
  }) %>%
    shiny::bindCache(kmer_table_data()) %>%
    shiny::bindEvent(input$kmer_distrib_button)

  output$seq_align_plot <- renderPlot({
    msa <- LymphoSeq2::alignSeq(stable(), repertoire_id = "IGH_MVQ92552A_BL", type = "junction", method = "ClustalW")
    LymphoSeq2::plotAlignment(msa)
  })

  # ------------------------------------------------------------------------------------------------------------------ # # nolint

  data_output <- reactive({
    if (input$tabselected == "gene_panel") {
        if (input$gene_sub_tab == "gene_freq_heat") {
            RedBlue <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(256)
            data <- list(gene_heatmap_data(), RedBlue)
            data_output <- pheatmap::pheatmap(data[[1]],
                                              color = data[[2]],
                                              scale = "row")
        } else if (input$gene_sub_tab == "gene_freq_bar") {
            data_output <- plot(gene_bar_data())
        } else {
            data_output <- gene_table_data()
        }
      } else if (input$tabselected == "chord_diagram") {
        # data_output <- circlize::chordDiagram(chord_data(), annotationTrack = c("grid", "name"))
        data_output <- chorddiag::chorddiag(chord_data(),
                                            type = "bipartite",
                                            showTicks = FALSE,
                                            groupColors = cbPalette,
                                            groupnameFontsize = 15,
                                            groupnamePadding = 10,
                                            margin = 90
                                        )
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
      } else if (input$tabselected == "pairwise_sim") {
        data_output <- pairwise_sim_data()
      } else if (input$tabselected == "top_seq_panel") {
        if (input$top_seq_sub_tab == "top_seq_plot") {
          data_output <- top_seq_data()
        } else {
          data_output <- top_prod_seq_data()
        }
      } else if (input$tabselected == "unique_seq_panel") {
        if (input$unique_seq_sub_tab == "produtive_seq_plot") {
          data_output <- prod_plot_data()
        } else if (input$unique_seq_sub_tab == "produtive_seq_table") {
          data_output <- unique_prod_data()
        }
      } else if (input$tabselected == "clonality_plot") {
        data_output <- clone_plot_data()
      } else if (input$tabselected == "lorenz_curve") {
        data_output <- lorenz_data()
      } else if (input$tabselected == "clone_track") {
        data_output <- clone_track_data()
      } else if (input$tabselected == "airr_table") {
        data_output <- airr_data()
      } else if (input$tabselected == "clonality_panel") {
        if (input$clonal_sub_tab == "clonality") {
          data_output <- clone_data()
        } else if (input$clonal_sub_tab == "clonality_stats") {
          data_output <- clone_stats_data()
        } else if (input$clonal_sub_tab == "count_stats") {
          data_output <- stats_count_data()
        }
      } else if (input$tabselected == "kmer_panel") {
        if (input$kmer_sub_tab == "count_kmers") {
          data_output <- kmer_table_data()
        } else if (input$kmer_sub_tab == "kmer_distrib") {
          data_output <- LymphoSeq2::kmerPlot(kmer_table_data(), input$k_top)
        }
      } else if (input$tabselected == "search_lymphoseqdb") {
        data_output <- lymphoseqdb_data()
      } else if (input$tabselected == "search_ireceptor") {
        data_output <- ireceptor_data()
      } else if (input$tabselected == "diff_abundance") {
        data_output <- diff_table_data()
      } else if (input$tabselected == "rarefaction_curve") {
        data_output <- rarefaction_data()
      }
  })

  output$rda <- downloadHandler(
    filename <- function() {
      paste0(input$tabselected, ".rda")
    },
    content <- function(file) {
      shiny::withProgress(
        message = "This make take a few minutes...",
        {
          study_table <- airr_data()
          amino_table <- productive_aa()
          nucleotide_table <- productive_nt()
          summary_table <- clonality_data()
          download_data <- data_output()
          save(download_data, study_table, amino_table,
                nucleotide_table, summary_table, file = file)
        }
      )
    }
  )

  output$tsv <- downloadHandler(
    filename <- function() {
      paste0(input$tabselected, ".tsv")
    },
    content <- function(file) {
      write.table(data_output(), file,
                  quote = FALSE, sep = "\t",
                  row.names = FALSE)
    }
  )

  output$excel <- downloadHandler(
    filename <- function() {
      paste0(input$tabselected, ".xlsx")
    },
    content <- function(file) {
      writexl::write_xlsx(data_output(), path = file)
    }
  )

  output$pdf <- downloadHandler(
    filename <- function() {
      paste0(input$tabselected, ".pdf")
    },
    content <- function(file) {
      pdf(file, width = 11, height = 8.5)

      if (input$tabselected == "gene_panel" || input$tabselected == "kmer_panel") {
        if(input$gene_sub_tab == "gene_freq_word") {
          htmlwidgets::saveWidget(gene_word_data(), "gene_wordcloud.html", selfcontained = FALSE)
          webshot2::webshot("gene_wordcloud.html", file, delay = 60)
        } else {
          print(data_output())
        }
      } else if (input$tabselected == "chord_diagram") {
          # data_output <- circlize::chordDiagram(chord_data(), annotationTrack = c("grid", "name"))
          chord_download <- chorddiag::chorddiag(chord_data(),
                                              type = "bipartite",
                                              showTicks = FALSE,
                                              groupColors = cbPalette,
                                              groupnameFontsize = 15,
                                              groupnamePadding = 10,
                                              margin = 90
                                          )
          htmlwidgets::saveWidget(chord_download, "temp.html", selfcontained = FALSE)
          webshot2::webshot("temp.html", file, cliprect = "viewport")
          # file.copy("chord_diagram.pdf", file)
      } else if (input$tabselected == "common_panel") {
        if (input$common_sub_tab == "common_bar") {
          print(data_output())
        } else if (input$common_sub_tab == "common_venn") {
          grid::grid.newpage()
          grid::grid.draw(data_output())
        } else {
          plot(data_output())
        }
      } else {
        plot(data_output())
      }
      dev.off()
    }
  )

}

shinyApp(ui = ui, server = server)
