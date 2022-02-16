# LymphoSeq2_Shiny

LymphoSeq2_Shiny is an interactive web app version of [LymphoSeq2](https://github.com/shashidhar22/LymphoSeq2/tree/v0.0.0.9000) built with the Shiny R package. This Shiny application will allow you to upload Adaptive Immune Receptor Repertoire Sequencing (AIRR-seq) data to display and perform various analytical visualizations. All graphical visualizations will also be available to be downloaded as PDF files. Data tables will be available to be downloaded as TSV, Excel, and RData files.

## Required Packages

The following packages are used by LymphoSeq2_Shiny and will need to be installed (if not already) before running the application (The * will indicate special instructions for installation. See below the list for more details.): 

1. shiny
2. shinycssloaders
3. ggplot2
4. plotly
5. DT
6. heatmaply
7. chorddiag
8. tidyverse
9. htmltools
10. sp
11. ggalluvial*

  ***\* Special Installation Instuctions***
  
  This application uses the development version of **ggalluvial** and needs to be installed using the following code:
  ```
  remotes::install_github("corybrunson/ggalluvial@main", build_vignettes = FALSE)
  ```

## Running the Application

The **Shiny** package must first be loaded into the R session and then call Shiny's `runApp` function to start the application.

```
# Load the Shiny package
library(shiny)

# Run the application
runApp("LymphoSeq2_app")

```
