
base_dir ="C:/Users/fmelo/Desktop/MAGMA_result_CORRECT_CONF_preprocessing_per_group_filtered_IND/.../"
method_dir = "LASSO"

base_directory = c(rep("C:/Users/fmelo/Desktop/MAGMA_result_CORRECT_CONF_preprocessing_per_group_filtered_IND/RANDOM_FOREST_results/",3),rep("C:/Users/fmelo/Desktop/MAGMA_result_CORRECT_CONF_preprocessing_per_group_filtered_IND/",2))
method_directory = c("LASSO","RFCV","SVM","LIMMA_05_result_correct","LIMMA_025_result_correct")
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
    source("C:/Users/fmelo/Documents/Github/MicroISN_IBD/Code/Grouping_LIMMA_results_shared_taxon/comprimer_v3.R")
    
  }
}

for (dirs in 1:nrow(imp_dir)){
  base_dir = imp_dir$base_directory[dirs]
  method_dir = imp_dir$method_directory[dirs]
  for (comb in 1:nrow(Combination)){
    outcome = Combination$outcome[comb]
    level = Combination$level[comb]
    source("C:/Users/fmelo/Documents/Github/MicroISN_IBD/Code/Grouping_LIMMA_results_shared_taxon/comprimer_v3_interactions.R")
    
  }
}
