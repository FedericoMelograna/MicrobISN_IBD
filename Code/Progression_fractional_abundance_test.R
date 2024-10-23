library(data.table)
library(dplyr)
library(tidyr)
library(xlsx)
library(vegan)
# Parameters --------------------------------------------------------------library("randomForest")

path_base = "---"
path_outcome = "---/Alpha_Beta_diversity_results/"

setwd(path_base)
files <- list.files(path=".", pattern="group*", 
                    full.names=TRUE, recursive=FALSE)
curr_wd = getwd()



data_frame_loop = expand.grid(files)

colnames(data_frame_loop) = c("files")

only_re  = data.frame(files = data_frame_loop[c(2,10,13,14,15),])
                

for (riga_data in (1:nrow(only_re)) ){
  setwd(curr_wd)
  el = only_re[riga_data,,drop =F]
  input_file = el$files
  direttoria = gsub( "./","",input_file)

  data_path   = "---/Data/"
  result_path = "---/Significant_diff_frac_abundance"
  
  setwd(path_base)
  

# MATCHING ----------------------------------------------------------------

  w0_metadata <- read.csv("./Data/w0_n335_metadata.txt", sep="")
  w0_metadata$VLECC.ID = paste0("A", w0_metadata$VLECC.ID)
  last_two_characters <- substr(direttoria, nchar(direttoria) - 1, nchar(direttoria))
  
  if (last_two_characters == "14"){
    wafter_meta <- read.csv("./Data/w14_n188_metadata.txt", sep="") 
  }   else {
    wafter_meta <- read.csv("./Data/w24_n160_metadata.txt", sep="")
  }
  wafter_meta$VLECC.ID = paste0("A", wafter_meta$VLECC.ID)
  otu_table = read.table(paste0("./", direttoria, "/Data/otu_table_selected_",direttoria, ".tsv" ))
  modified_string <- sub("_[^_]*$", "", direttoria)
  otu_table_wo = read.table(paste0("./Analysis/", modified_string, "/Data/otu_table_selected_",modified_string, ".tsv" ))
  
  new_rownames <- wafter_meta$VLECC.ID[match(rownames(otu_table), wafter_meta$FC.nummer)]
  rownames(otu_table) <- new_rownames
  
  new_rownames <- w0_metadata$VLECC.ID[match(rownames(otu_table_wo), w0_metadata$FC.nummer)]
  rownames(otu_table_wo) <- new_rownames
  
  
  match_otu_table = otu_table[rownames(otu_table) %in% rownames(otu_table_wo),]
  match_otu_table_wo = otu_table_wo[rownames(otu_table_wo) %in% rownames(otu_table),]
  
  ASV_taxonomic_annotation <- read.delim("---/Data/ASV_taxonomic_annotation.txt")
  ASV_reduced = ASV_taxonomic_annotation[ASV_taxonomic_annotation$ASV %in% colnames(match_otu_table),]
  df = (t(ASV_reduced))
  dim(df)
  colnames(df) = df[1,]
  df = df[-1,]
  family_info <- df[5, ]  # This is a named vector where names are OTUs and values are families
  names(family_info)
  otu_table_rel<-as.data.frame(t(apply(match_otu_table,1, function(x) x/sum(x))))
  
  ddf = data.frame(family = family_info, OTU = names(family_info))
  otu_long <- otu_table_rel %>%
    tibble::rownames_to_column("Individual") %>%
    tidyr::pivot_longer(-Individual, names_to = "OTU", values_to = "Abundance")
  
  otu_with_family <- merge(otu_long, ddf, by = "OTU")
  
  
  # Sum abundances by Individual and Family
  final_dataset <- otu_with_family %>%
    group_by(Individual, family) %>%
    summarise(Total_Abundance = sum(Abundance), .groups = 'drop')
  
  # Reshape the final dataset into wide format (optional)
  final_wide <- pivot_wider(final_dataset, names_from = family, values_from = Total_Abundance, values_fill = 0)
  

  # Same for w0 -------------------------------------------------------------
  setwd(path_base)
  
  
  ASV_reduced_wo = ASV_taxonomic_annotation[ASV_taxonomic_annotation$ASV %in% colnames(match_otu_table_wo),]
  df_wo = (t(ASV_reduced_wo))
  dim(df_wo)
  colnames(df_wo) = df_wo[1,]
  df_wo = df_wo[-1,]
  family_info_wo <- df_wo[5, ]  # This is a named vector where names are OTUs and values are families
  names(family_info_wo)
  otu_table_rel_wo<-as.data.frame(t(apply(match_otu_table_wo,1, function(x) x/sum(x))))
  
  ddf_wo = data.frame(family = family_info_wo, OTU = names(family_info_wo))
  otu_long_wo <- otu_table_rel_wo %>%
    tibble::rownames_to_column("Individual") %>%
    tidyr::pivot_longer(-Individual, names_to = "OTU", values_to = "Abundance")
  
  otu_with_family_wo <- merge(otu_long_wo, ddf_wo, by = "OTU")
  
  
  # Sum abundances by Individual and Family
  final_dataset_wo <- otu_with_family_wo %>%
    group_by(Individual, family) %>%
    summarise(Total_Abundance = sum(Abundance), .groups = 'drop')
  
  # Reshape the final dataset into wide format
  final_wide_wo <- pivot_wider(final_dataset_wo, names_from = family, values_from = Total_Abundance, values_fill = 0)
  
  merged_final_dataset = merge(final_dataset, final_dataset_wo, by = c("Individual", "family"), all = T)
  colnames(merged_final_dataset)[c(3,4)] = c("After", "Before")
  
  nomiss = merged_final_dataset[!is.na(merged_final_dataset$family),]
  nomiss[is.na(nomiss)] = 0
  
  # Perform paired t-test for each family
  results <- nomiss %>%
    group_by(family) %>%
    summarise(
      p_value = t.test(Before, After, paired = TRUE)$p.value, 
      
    )
  df_summary <- nomiss %>%
    group_by(family) %>%
    summarise(
      Avg_After = mean(After, na.rm = TRUE),
      Avg_Before = mean(Before, na.rm = TRUE)
    )
  More_10pc = df_summary[df_summary$Avg_After > 0.1 & df_summary$Avg_Before > 0.1,]
  results <- results %>%
    mutate(p_adj = p.adjust(p_value, method = "bonferroni")) 
  results10pc = results[results$family %in% More_10pc$family,]
  results10pc <- results10pc %>%
    mutate(p_adj = p.adjust(p_value, method = "bonferroni")) 
  
  setwd(result_path)
  write.xlsx(results, file = paste0("Adj_p_value_frac_abund", modified_string, ".xlsx"))
  write.xlsx(results10pc, file = paste0("TenPC_Adj_p_value_frac_abund", modified_string, ".xlsx"))
  saveRDS(object = final_dataset_wo, file = paste0("w0_", modified_string, ".rds"))
  saveRDS(object = final_dataset, file = paste0("wafter_", direttoria, ".rds"))
  

  # paired test -------------------------------------------------------------

  df_long <- nomiss %>%
    pivot_longer(cols = c("Before", "After"), names_to = "TimePoint", values_to = "FractionalAbundance") %>%
    mutate(
      IndividualID = as.numeric(factor(Individual, levels = unique(Individual))) # Create numeric IndividualID
    ) %>%
    select(IndividualID, TimePoint, family, FractionalAbundance)

    df_matrix <- df_long %>%
    pivot_wider(names_from = family, values_from = FractionalAbundance) %>%
    select(-IndividualID)

  # MANOVA
  adonis_result <- adonis(df_matrix[, -1] ~ df_matrix$TimePoint, method = "euclidean")
  print(adonis_result$aov.tab)
  
  saveRDS(adonis_result, file = paste0("Adonis_result", modified_string, ".rds"))
  
}
                