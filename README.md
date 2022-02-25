# LymphoSeq2_Shiny

LymphoSeq2_Shiny is an interactive web app version of [LymphoSeq2](https://github.com/shashidhar22/LymphoSeq2/tree/v0.0.0.9000) built with the Shiny R package. This Shiny application will allow you to upload Adaptive Immune Receptor Repertoire Sequencing (AIRR-seq) data to display and perform various analytical visualizations. All graphical visualizations will be available to be downloaded as PDF or RData files. Data tables will be available to be downloaded as TSV, Excel, and RData files.

## Required Packages

The following packages are used by LymphoSeq2_Shiny and will need to be installed (if not already) before running the application: 

1. shiny
2. shinyjs
3. shinycssloaders
4. ggplot2
5. plotly
6. DT
7. heatmaply
8. pheatmap
9. wordcloud2
10. tidyverse
11. htmltools
12. sp
13. writexl

These packages will also need to be installed under special instructions:

1. LymphoSeq2
   
   The **LymphoSeq2** package can be installed from Github:
   ```
   # install the devtools package
   install.packages("devtools")

   devtools::install_github("shashidhar22/LymphoSeq2", build_vignette=TRUE)
   ```

2. ggalluvial
  
   This application uses the development version of **ggalluvial** and needs to be installed from Github using the following code:
   ```
   # install the remotes package
   install.packages("remotes)
  
   remotes::install_github("corybrunson/ggalluvial@main", build_vignettes = FALSE)
   ```
3. chorddiag
   
   The **chordiagg** package may be incomaptible with newer versions of R. If this is the case, the package will need to be installed directly from Github  
   the `devtools` package:
   ```
   devtools::install_github("mattflor/chorddiag", build_vignettes = TRUE)
   ```

## Running the Application

The **Shiny** package must first be loaded into the R session. Then, to start the application, simply call Shiny's `runApp` function.

```
# Load the Shiny package
library(shiny)

# Run the application
runApp("LymphoSeq2_app")

```
