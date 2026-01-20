
############################################################
# Auto-detect path - works in both Docker and local environments
if (file.exists("/app/Scripts/setup.R")) {
  # Running in Docker
  path_to_duplica <- "/app"
} else if (file.exists("/Users/grantmerrigan/Documents/DuplicA/app/Scripts/setup.R")) {
  # Running locally on Mac
  path_to_duplica <- "/Users/grantmerrigan/Documents/DuplicA/app"
} else {
  # Fallback - try to detect from current file location
  path_to_duplica <- dirname(dirname(sys.frame(1)$ofile))
}
############################################################


#library(shiny)
#library(shinyFiles)
library(bslib)
library(htmltools)
#library(shinyalert)
library(fs)
library(Biostrings) 
library(tools)
library(testthat) 
library(tidyverse)
library(readxl)
#library(shinyjs)
library(ape)
library(R.utils)
library(seqinr)
library(evemodel)
library(igraph)
library(biomartr)
library(data.table)
library(biomaRt)
library(RColorBrewer)
library(jsonlite)
library(future)
library(furrr)

# Load clusterProfiler (used by model_GO.R) - optional to avoid startup errors
tryCatch({
  suppressPackageStartupMessages(library(clusterProfiler))
}, error = function(e) {
  warning("clusterProfiler not available - GO analysis will not work. Install with: BiocManager::install('clusterProfiler')")
})

#library(r3dmol)
#library(httr)
#library(bio3d)




here_duplica <- path_to_duplica

source(paste0(here_duplica, '/Scripts/multispecies_functions.R'))
source(paste0(here_duplica, '/Scripts/model_OrthoFinder.R'))
source(paste0(here_duplica, '/Scripts/model_DnDs.R'))
source(paste0(here_duplica, '/Scripts/model_EVE.R'))


here_temp <- paste0(here_duplica, '/Temp')
here_results <- paste0(here_duplica, '/Results')
here_static_images <- paste0(here_duplica, '/Static/Images')


#here_duplica_linux <- '/mnt/c/Users/17735/Downloads/DuplicA'
#here_linux_temp <- paste0(here_duplica_linux, '/app/Temp')
#here_linux_dep <- paste0(here_duplica_linux, '/app/Dependencies')
#here_linux_results <- paste0(here_duplica_linux, '/app/Results')

# wsl is required 


