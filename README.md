# LymphoSeq2_Shiny

LymphoSeq2_Shiny is an interactive web app version of [LymphoSeq2](https://github.com/shashidhar22/LymphoSeq2/tree/v0.0.0.9000) built with the Shiny R package. This Shiny application will allow you to upload Adaptive Immune Receptor Repertoire Sequencing (AIRR-seq) data to display and download as PDF various analytical visualizations. 

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

* Special Installation Instuctions
This application uses the development version of **ggalluvial** and needs to be installed using the following code
```
remotes::install_github("corybrunson/ggalluvial@main", build_vignettes = FALSE)
```



To run the application, the **Shiny** package must be installed and loaded in an R session.
```
install.packages("shiny")
library(shiny)
```
