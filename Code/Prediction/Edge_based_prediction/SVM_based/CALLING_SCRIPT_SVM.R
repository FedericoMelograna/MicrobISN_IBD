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
library("randomForest")

# Parameters --------------------------------------------------------------
eliminate = T
diagnosi = c("CD", "UC") 
terapia  =c("VDZ")# terapia = c("TNF", "UST", "VDZ")

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



# JUST for the SVM --------------------------------------------------------
for (riga_data in (1:nrow(data_frame_loop)) ){
  setwd(curr_wd)
  el = data_frame_loop[riga_data,]
  input_file = el$files
  outcome_file = el$outcome
  args = list()
  args$direttoria = gsub( "./","",input_file)
  args$outcome = outcome_file
  source ("Code/Prediction/script_computing_svm_v2.R")
  print(args$direttoria)
}

