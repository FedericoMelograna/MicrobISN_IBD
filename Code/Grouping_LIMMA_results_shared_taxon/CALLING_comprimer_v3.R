

## Directory combination, the base directory is simply the path to the LIMMA 05 and 025 results
base_directory = c(rep("---",2))
method_directory = c("LIMMA_05_result_correct","LIMMA_025_result_correct")
imp_dir = data.frame(cbind(base_directory,method_directory))

outcome = c("Endoscopic", "Biomarker","Clinical")
level = c("Taxa", "Genus", "Family")
Combination = expand.grid("outcome" = outcome, "level" = level)

for (dirs in 1:nrow(imp_dir)){
  base_dir = imp_dir$base_directory[dirs]
  method_dir = imp_dir$method_directory[dirs]
  for (comb in 1:nrow(Combination)){
    outcome = Combination$outcome[comb]
    level = Combination$level[comb]
    source("---/comprimer_v3.R")
    
  }
}

for (dirs in 1:nrow(imp_dir)){
  base_dir = imp_dir$base_directory[dirs]
  method_dir = imp_dir$method_directory[dirs]
  for (comb in 1:nrow(Combination)){
    outcome = Combination$outcome[comb]
    level = Combination$level[comb]
    source("---/comprimer_v3_interactions.R")
    
  }
}
