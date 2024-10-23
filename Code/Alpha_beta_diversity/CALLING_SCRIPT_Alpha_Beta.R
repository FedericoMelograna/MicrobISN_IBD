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

path_base = "C:/Users/fmelo/Desktop/Microbiome_drug_repurpose_project/"
path_outcome = "C:/Users/fmelo/Desktop/Microbiome_drug_repurpose_project/Alpha_Beta_diversity_results/"

setwd(path_base)


files <- list.files(path=".", pattern="group*", 
                    full.names=TRUE, recursive=FALSE)
curr_wd = getwd()



args = list()
args$direttoria = NA



outcomes = c("Clinical_outcome_combined", "Biomarker_outcome_combined","Endoscopic_outcome_combined" )
data_frame_loop = expand.grid(files, outcomes)

colnames(data_frame_loop) = c("files", "outcomes")


# Canonical Cohort_based_Alpha_Beta_Diversity.R ---------------------------


for (riga_data in (1:nrow(data_frame_loop)) ){ ## TOTAL of 27
  setwd(curr_wd)
  el = data_frame_loop[riga_data,]
  input_file = el$files
  outcome_file = el$outcome

  args = list()
  args$direttoria = gsub( "./","",input_file)
  args$outcome = outcome_file

  start = Sys.time()
  source ("Cohort_based_Alpha_Beta_Diversity.R")
  stop = Sys.time()
  stop-start
  print(args$direttoria)
}




