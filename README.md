# LymphoSeq2_Shiny

LymphoSeq2_Shiny is an interactive web application of [LymphoSeq2](https://github.com/shashidhar22/LymphoSeq2/tree/v0.0.0.9000) built with the R/Shiny package. This Shiny application will allow you to upload Adaptive Immune Receptor Repertoire Sequencing (AIRR-seq) data to display and perform various analytical visualizations. All graphical visualizations are available to be downloaded as PDF or RData files. All Data tables will be available to be downloaded as TSV, Excel, or RData files.

## Dependencies

All of the required packages will be automtically installed when the application is run. However, the user will need to manually install **LymphoSeq2**. Instructions are below.

### Installing LymphoSeq2

LymphoSeq2 can be installed from an R session by running the commands below:

   ```
   # install the devtools package first (skip if already installed)
   install.packages("devtools")
   
   # load devtools package
   library(devtools)

   # install LymphoSeq2
   devtools::install_github("WarrenLabFH/LymphoSeq2", ref="v1", build_vignette=FALSE)
   ```

## Running the Application Locally (from source)

The **Shiny** package must first be loaded into the R session. Then, to start the application, simply execute Shiny's `runApp` function.

```
# load the Shiny package
library(shiny)

# run the application (assuming current directory is LymphoSeq2_Shiny)
runApp("LymphoSeq2_app/")
```

## Using the Application

A comprehensive User Guide is available [here](https://elulu3.github.io/LymphoSeq2_Shiny/).

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


