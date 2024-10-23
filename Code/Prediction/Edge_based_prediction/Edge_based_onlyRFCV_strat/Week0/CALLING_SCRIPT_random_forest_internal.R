library(lionessR)
library(igraph)
library(reshape2)
# BiocManager::install("limma")
library(limma)
library(SummarizedExperiment)
library(ROCR)
library(e1071)
library(randomForest)
library(dplyr)
library(readr)
library(AUC)
library(data.table)

# Parameters --------------------------------------------------------------
library("randomForest")
eliminate = T
diagnosi = c("CD", "UC") ## diagnosi = c("CD","UC")
terapia  =c("VDZ")# c("VDZ") # terapia = c("TNF", "UST", "VDZ")
setwd("---")
path_base = "---"
getwd()


files <- list.files(path=".", pattern="group*", 
                    full.names=TRUE, recursive=FALSE)
curr_wd = getwd()

files
el
args = list()
args$direttoria = NA
outcomes = rep(c("Clinical_outcome_combined", "Biomarker_outcome_combined","Endoscopic_outcome_combined" ), length(files) )
files_per_outcome = rep(files, each = 3 )
data_frame_loop  = data.frame(files = files_per_outcome, outcome = outcomes)

## NEW 19 --> until
for (riga_data in (1:nrow(data_frame_loop)) ){
  setwd(curr_wd)
  el = data_frame_loop[riga_data,]
  input_file = el$files
  outcome_file = el$outcome
  args = list()
  args$direttoria = gsub( "./","",input_file)
  args$outcome = outcome_file
  source ("---/script_computing_Random_forest_internal_v2.R")
  print(args$direttoria)
}

