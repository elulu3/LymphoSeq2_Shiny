# LymphoSeq2_Shiny

LymphoSeq2_Shiny is an interactive web application of [LymphoSeq2](https://github.com/shashidhar22/LymphoSeq2/tree/v0.0.0.9000) built with the R/Shiny package. This Shiny application will allow you to upload Adaptive Immune Receptor Repertoire Sequencing (AIRR-seq) data to display and perform various analytical visualizations. All graphical visualizations are available to be downloaded as PDF or RData files. All Data tables will be available to be downloaded as TSV, Excel, or RData files.

## Dependencies

All of the required packages will be automatically installed when the application is run.

**Note:** LymphoSeq2 is currently undergoing updates and so some functions in the application may not work as intended.

## Running the Application Locally (from source)

First, clone the repository in the terminal.

```
# retrieve a local copy of the code
git clone https://github.com/elulu3/LymphoSeq2_Shiny.git
```

Start an R session and in that session, run the application from the copy of the code.

```
# install shiny package
install.packages("shiny")

# load the shiny package
library(shiny)

# run the application by providing runApp the path to app.R
runApp("LymphoSeq2_Shiny/LymphoSeq2_app/")
```

## Using the Application

A comprehensive User Guide is available [here](https://elulu3.github.io/LymphoSeq2_Shiny/) [update in progress].

## Adjusting for Data Size and Memory Usage

The Shiny application enforces a file upload limit which can be adjusted in the source code in the following line:

```
options(shiny.maxRequestSize = [SIZE] * 1024^2)
```

Should you run into an error message that states the following:

```
simpleError: vector memory exhausted (limit reached?)
```
The following commands entered into a terminal followed by restarting the R session should help resolve the error:

```
cd ~ 
touch .Renviron 
open .Renviron 
R_MAX_VSIZE=100Gb
```


