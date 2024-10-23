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
library(RColorBrewer)


# Parameters --------------------------------------------------------------
library("randomForest")

path_base = "..."
path_outcome = "..."

setwd(path_base)

getwd()


files <- list.files(path=".", pattern="group*", 
                    full.names=TRUE, recursive=FALSE)
curr_wd = getwd()

## Choose your hierarchy
component_tackled = "Phylum" # "Family" "Class
args = list()
args$compon = component_tackled
args$direttoria = NA



outcomes = c("Clinical_outcome_combined", "Biomarker_outcome_combined","Endoscopic_outcome_combined" )
data_frame_loop = expand.grid(files)
colnames(data_frame_loop) = c("files")



if (component_tackled == "Phylum"){
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  coloring_df = data.frame("Phylums" = c("Actinobacteria",   "Bacteroidetes",      "Firmicutes",  "Proteobacteria", "Verrucomicrobia", NA), 
                           "Colours"= cbPalette[1:6])

} else if (component_tackled == "Family"){
  gl  =c("Desulfovibrionaceae"    ,   "Enterobacteriaceae",        "Sutterellaceae",           "Erysipelotrichaceae"   ,
         "Lachnospiraceae" ,      "Bifidobacteriaceae",  "Ruminococcaceae",    
         "Rikenellaceae",       "Veillonellaceae",     "Bacillales_Incertae_Sedis_XI"  ,    "Enterococcaceae",    
         "Lactobacillaceae",    "Verrucomicrobiaceae", "Pasteurellaceae",     "Coriobacteriaceae",  
         "Streptococcaceae",    "Bacteroidaceae",      "Porphyromonadaceae",  "Acidaminococcaceae", 
         "Actinomycetaceae",    "Clostridiaceae_1",    "Prevotellaceae",      "Peptostreptococcaceae",         
         "Clostridiales_Incertae_Sedis_XIII" ,"Eubacteriaceae",      "Clostridiales_Incertae_Sedis_XI"  )
  gl = gl[order(gl)]
  n <- length(gl)
  set.seed(111)
  colrs <- brewer.pal.info[brewer.pal.info$colorblind == TRUE, ]
  col_vec = unlist(mapply(brewer.pal, colrs$maxcolors, rownames(colrs)))
  col <- sample(col_vec, length(gl))
  col <- col_vec[1:length(gl)]
  coloring_df = data.frame("Phylums" = gl, "Colours" = col)
  
} else if(component_tackled == "Class"){
  gl = c("Deltaproteobacteria","Gammaproteobacteria","Betaproteobacteria","Erysipelotrichia","Clostridia",          "Actinobacteria" ,     "Bacteroidia"    ,    
         "Negativicutes"  ,     "Bacilli"       ,      "Verrucomicrobiae"    )
  gl = gl[order(gl)]
  n <- length(gl)
  set.seed(111)
  colrs <- brewer.pal.info[brewer.pal.info$colorblind == TRUE, ]
  col_vec = unlist(mapply(brewer.pal, colrs$maxcolors, rownames(colrs)))
  col <- sample(col_vec, length(gl))
  col <- col_vec[1:length(gl)]

  coloring_df = data.frame("Phylums" = gl, "Colours" = col)
  
}



for (riga_data in (1:nrow(data_frame_loop)) ){
  setwd(curr_wd)
  el = data_frame_loop[riga_data,,drop =F]
  input_file = el$files
  

  args$direttoria = gsub( "./","",input_file)
  
  start = Sys.time()
  source ("Fractional_abundances_Pie_chart/PIE_chart_maker.R")
  stop = Sys.time()
  stop-start
  print(args$direttoria)
}
