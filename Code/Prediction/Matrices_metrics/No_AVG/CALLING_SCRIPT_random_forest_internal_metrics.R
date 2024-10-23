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
terapia  =c("VDZ") # terapia = c("TNF", "UST", "VDZ")

path_base = "---/metricmatrices/"
path_outcome = "---"

setwd(path_base)

getwd()


files <- list.files(path=".", pattern="group*", 
                    full.names=TRUE, recursive=FALSE)
curr_wd = getwd()

files
el
args = list()
args$direttoria = NA
outcomes = rep(c("Clinical_outcome_combined", "Biomarker_outcome_combined","Endoscopic_outcome_combined" ), 2*length(files) )



files_per_outcome = rep(files, each = 6 )
data_frame_loop  = data.frame(files = files_per_outcome, outcome = outcomes)

outcomes = c("Clinical_outcome_combined", "Biomarker_outcome_combined","Endoscopic_outcome_combined" )
na_treatment = c("zero", "elim")
feat_sel = c("AverageShortestPathLength","ClusteringCoefficient",    "ClosenessCentrality",       "Eccentricity","Stress","Degree",
      "BetweennessCentrality","NeighborhoodConnectivity","Radiality","TopologicalCoefficient") 
data_frame_loop = expand.grid(files, feat_sel, outcomes, na_treatment)

colnames(data_frame_loop) = c("files","feat_sel", "outcomes", "na_treatment")

for (riga_data in (1:nrow(data_frame_loop)) ){
  setwd(curr_wd)
  el = data_frame_loop[riga_data,]
  input_file = el$files
  outcome_file = el$outcome
  if (outcome_file == "Clinical_outcome_combined"){
    next;
  }
  args = list()
  args$direttoria = gsub( "./","",input_file)
  args$outcome = outcome_file
  args$feat_sel = el$feat_sel
  args$na_treatment = el$na_treatment
  start = Sys.time()
  source (paste0(path_outcome, "Code/Prediction/Matrices_metrics/No_AVG/script_computing_Random_forest_internal_v2_metrics.R"))
  stop = Sys.time()
  stop-start
  print(args$direttoria)
}


