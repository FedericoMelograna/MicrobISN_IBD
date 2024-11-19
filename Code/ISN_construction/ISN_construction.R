# if (!exists("args")) {
#   suppressPackageStartupMessages(library("argparse"))
#   parser <- ArgumentParser()
#   parser$add_argument("-a", "--arg1", type="character", defalt="a",
#                       help="First parameter [default %(defult)s]")
#   parser$add_argument("-b", "--arg2", type="character", defalt="b",
#                       help="Second parameter [default %(defult)s]")
#   args <- parser$parse_args()
# }

print (args)

print(args$eliminate)
print(args$diagnosi)
print(args$terapia)


set.seed(123)


# Parameters --------------------------------------------------------------
eliminate = args$eliminate
diagnosi = args$diagnosi ## diagnosi = c("CD","UC")
terapia  =  args$terapia# c("VDZ") # terapia = c("TNF", "UST", "VDZ")
setwd(base_dir)

source (paste0(base_dir, "/Code/ISN_construction/",function_script ))




# Data importing ----------------------------------------------------------


#### ONLY FOR THE w0--> 
if (weeks == "0"){
mapping_nameIndividual_SPARCC <- read.csv(paste0(base_dir, "/Data/",individual_map ), sep="")
}
w0_n335_OTU_table = data.frame(fread(file = paste0(base_dir, "/Data/",OTU_dataset ) ))
rownames(w0_n335_OTU_table) = w0_n335_OTU_table[,1] # First row the identifier
w0_n335_OTU_table = w0_n335_OTU_table %>% dplyr::select(-c("SampleList_Metadata_Prediction"))
if ("Sequencing.QC" %in% colnames(w0_n335_OTU_table)){
  w0_n335_OTU_table = w0_n335_OTU_table %>% dplyr::select(-c("Sequencing.QC"))
  
}

in_cov = read.delim(paste0(base_dir, "/Data/",OTU_metadata ))
# ELIMINATE 0 and NA on the active disease!
in_cov = in_cov[in_cov$active_disease_at_baseline != "0",]
# Reg the re-runs, you are absolutely right - we need need to exclude the samples (rows) 
# who are identified with 0 in column U ("active_disease_at_baseline") as indicated in "w0_n335_metadata" file. 


# covariates_kept = c("Age", "Gender","FC.nummer", "Disease_duration", "Disease_location_baseline")
import_cov = in_cov %>% dplyr::select(all_of(covariates_kept))
# import_cov =in_cov %>% dplyr::select(c(Age, Gender,FC.nummer, Disease_duration, Disease_location_baseline))


# identifier_disease = c("FC.nummer","Diagnosis", "Therapy_2" )
filtering_disease = in_cov %>% dplyr::select(all_of(identifier_disease))
head(filtering_disease); dim(filtering_disease)
# check 
sum(is.na(import_cov)) ; head(import_cov)


ASV_taxonomic_annotation <- read.delim(paste0(base_dir, "/Data/",ASV_annotation ) )



# PREPROCESSING covariates ------------------------------------------------


# IMPUTE the mode for the MISSING VALUES

Mode(import_cov$Disease_location_baseline)
import_cov$Disease_location_baseline[is.na(import_cov$Disease_location_baseline)] = Mode(import_cov$Disease_location_baseline)
table(import_cov$Disease_location_baseline)
sum(is.na(import_cov)) ; head(import_cov)


# CREATING DUMMY BY HAND --------------------------------------------------

DummyL1 = ifelse(import_cov$Disease_location_baseline == "L1",1,0)
DummyL2 = ifelse(import_cov$Disease_location_baseline == "L2",1,0)
DummyL3 = ifelse(import_cov$Disease_location_baseline == "L3",1,0)
sum(c(DummyL1,DummyL2,DummyL3))

import_cov$Gender = ifelse(import_cov$Gender == "F",1,0)
import_cov$Age = as.numeric(gsub(",",".",import_cov$Age))
import_cov$Disease_duration = as.numeric(gsub(",",".",import_cov$Disease_duration))

import_cov$DummyL1 = DummyL1
import_cov$DummyL2 = DummyL2
import_cov$DummyL3 = DummyL3

# The important confounders are age, gender, Disease_duration and Disease_location_baseline.

import_cov_final = import_cov %>% dplyr::select(Age,Gender,FC.nummer,Disease_duration,DummyL1,DummyL2,DummyL3)


backup_w0_n335_OTU_table = w0_n335_OTU_table
backup_import_cov_final = import_cov_final

covs_filt = data.frame(Age = import_cov_final$Age,Gender = import_cov_final$Gender,Disease_duration = import_cov_final$Disease_duration,
                       DummyL1 = import_cov_final$DummyL1, DummyL2 = import_cov_final$DummyL2, DummyL3= import_cov_final$DummyL3)
rownames(covs_filt)= import_cov_final$FC.nummer


group_selection = selecting_only_group(w0_n335_OTU_table, covs_filt, select = filtering_disease, diagnosi =diagnosi , terapia = terapia)
otu_table_selected = group_selection[[1]] ; covariates_selected = group_selection[[2]] ; groupz = group_selection[[3]]
if (dim(otu_table_selected)[1] != 0) {
  
  final = "group"
  for (el in diagnosi){
    print(el)
    final = paste(final,el,sep = "_")
  }
  for(el in terapia){
    final = paste(final, el,sep = "_")
  }
  
  if (eliminate)
  {
    final = paste(final, "eliminate",sep = "_")
  } else {
    final = paste(final, "Not_eliminate",sep = "_")
    
  }
  final = paste(final, weeks, sep = "_")
  otu_table_selected = group_selection[[1]] ; covariates_selected = group_selection[[2]] ; groupz = group_selection[[3]]
  otu_table_selected[1:5,1:5]; dim(otu_table_selected)
  
  merged = merge(import_cov, filtering_disease, by = "FC.nummer")
  table(merged$Diagnosis, merged$Disease_location_baseline)
  
  nrow(covariates_selected)
  table(import_cov$Disease_location_baseline)
  res = importing_base_data_and_select_cropping(dati = otu_table_selected, covariate = covariates_selected,crop = NULL)
  dati = res[[1]]
  covariate = res[[2]]
  covariates = cbind(covariate$Age,covariate$Gender,covariate$Disease_duration, covariate$DummyL1
                     , covariate$DummyL2, covariate$DummyL3)
  diagnosi
  
  final_numb_cov = 6
  if (diagnosi == "UC" && length(diagnosi) == 1){
    covariates = covariates[,1:3]
    # covariates_selected = covariates_selected[,1:3]
    final_numb_cov = 3
    ## no need for the dummy, all of them are 0
  } else if(sum(apply(covariates[,4:6],1,sum) == 0) == 0){
    covariates = covariates[,1:5]
    # covariates_selected = covariates_selected[,1:5]
    final_numb_cov = 5
  }
  
  
  
  rownames(covariates)= covariate$FC.nummer
  
  magma_Stool_AllCov <- magma(data = dati,X = covariates)
  ### FIND THE ones that are not intersecting ! AND ELIMINATE THEM
  tt = tryCatch( magma(data = dati,X = covariates) ,error=function(e) e, warning=function(w) w)
  
  
  if (eliminate == T & is(tt,"warning") ){
    ## if warning means that raise a warning, that there are sample to eliminate
    # first_warn = names(last.warning)[1]
    
    to_eliminate_10 = stri_extract_all(tt[1],regex="FC-(\\d)+")[[1]]
    
    otu_table_eliminated =  dati[!(rownames(dati) %in%to_eliminate_10) ,]
    rownames(dati) == rownames(covariates_selected)
    covariates_eliminated = covariates_selected[!(rownames(covariates_selected) %in% to_eliminate_10) ,]
    covariates = cbind(covariates_eliminated$Age,covariates_eliminated$Gender,covariates_eliminated$Disease_duration, covariates_eliminated$DummyL1
                       , covariates_eliminated$DummyL2, covariates_eliminated$DummyL3)
    covariates = covariates[,1:final_numb_cov]
    
    
  } else {
    otu_table_eliminated = dati
    covariates_eliminated = covariates_selected
  }
  
  
  if(final_numb_cov == 5){
    cov_finals = covariates_eliminated %>% dplyr::select(-c(DummyL3))
  } else if (final_numb_cov == 3){
    cov_finals = covariates_eliminated %>% dplyr::select(-c(DummyL3,DummyL2,DummyL1))
    
  } else {
    cov_finals = covariates_eliminated
  }
  
  
  curr_wd = getwd()
  dir.create(file.path(curr_wd, final))
  dir.create(file.path(curr_wd, final,"Data"))
  setwd(file.path(curr_wd, final))
  getwd()
  write.table(otu_table_eliminated, file = paste0("Data/otu_table_selected_",final,  ".tsv" ), quote = FALSE
              , row.names = TRUE, col.names = TRUE,sep = "\t")
  write.table(covariates_eliminated, file = paste0("Data/covariates_total_", final, ".tsv" ), quote = FALSE
              , row.names = TRUE, col.names = TRUE,sep = "\t")
  write.table(groupz, file = paste0("Data/Group_selection_",final, ".tsv" ), quote = FALSE, row.names = TRUE
              , col.names = FALSE,sep = "\t")
  write.table(cov_finals, file = paste0("Data/covariates_selected_",final, ".tsv" ), quote = FALSE, row.names = TRUE
              , col.names = FALSE,sep = "\t")
  if (eliminate == T & is(tt,"warning") ){
    write(to_eliminate_10, paste0("Data/ELIMINATED_samples_",final, ".txt" ))
  }
  
  curr_wd = getwd()
  # dir.create(file.path(curr_wd, "LooNet"))
  # try(magma(data = otu_table_eliminated,X = covariates ) )
  magma_Stool_AllCov <- magma(data = otu_table_eliminated,X = covariates) #,distrib = "ZIP")
  if (weeks == "0"){
  build_LOO_net(otu_table_eliminated, covariates, mapping_nameIndividual_SPARCC, weeks)
  } else {
    build_LOO_net(otu_table_eliminated, covariates, NA, weeks)
  }
  
  GLOB_RESULT = build_GLOB_net(otu_table_eliminated, covariates)
  
  glob_model = GLOB_RESULT[[1]]
  binary_result = GLOB_RESULT[[2]]
  continuous_result = GLOB_RESULT[[3]]
  
  write.table(binary_result, file = paste0("Data/binary_result_",final, ".tsv" ), quote = FALSE, row.names = TRUE
              , col.names = TRUE,sep = "\t")
  write.table(continuous_result, file = paste0("Data/continuous_result_",final, ".tsv" ), quote = FALSE, row.names = TRUE
              , col.names = TRUE,sep = "\t")
  
  # GRAPHS ------------------------------------------------------------------
  
  
  ASV_final = ASV_taxonomic_annotation[ASV_taxonomic_annotation$ASV %in% colnames(otu_table_eliminated),]
  df = (t(ASV_final))
  dim(df)
  colnames(df) = df[1,]
  df = df[-1,]
  
  
  
  ll <- igraph::layout_with_fr(igraph::graph_from_adjacency_matrix(
    magma_Stool_AllCov$refit, mode="undirected"))
  
  write.table(ll,file = "LAYOUT.txt")
  ll2 = as.matrix(read.table("LAYOUT.txt"))
  
  dir.create(file.path(curr_wd, "Graphs"))
  setwd(file.path(curr_wd, "Graphs"))
  
  
  
  graphical_printing(magma_Stool_AllCov, ll2, level = 2, dataset_match_taxa = df)
  graphical_printing(magma_Stool_AllCov, ll2, level = 3, dataset_match_taxa = df)
  graphical_printing(magma_Stool_AllCov, ll2, level = 4, dataset_match_taxa = df)
  graphical_printing(magma_Stool_AllCov, ll2, level = 5, dataset_match_taxa = df)
  graphical_printing(magma_Stool_AllCov, ll2, level = 6, dataset_match_taxa = df)
  
  # LAYOUT = "C:/Users/fmelo/Desktop/Backup_Federico/work_microbiome/SPARCC_and_MAGMA/MAGMA_result_CORRECT_CONF_preprocessing_per_group_filtered_IND/LAYOUT_PHYLUM__graph.txt"
  # Do something, or tell me why it failed
  LAYOUT = paste0(base_dir, "/Graphs/","LAYOUT_PHYLUM__graph_4nodes.txt" )
  # IF POSSIBLE the layout will be the same as the one in the first group_CD_TNF
  
  # IN ANY CASE layout binary and continuous will be the same 
  # N.b. they can differ casue a lot of continuous are not significant OR cause 
  # the binary exist but the continuous counterpart is too low (e.g. less than a 0.1 threshold, with round(,1) )
  # and so it is not detected
  tryCatch(
    # This is what I want to do...
    {
      plot_graphical_phylum(
        OTU_RESULT = binary_result,
        label_name_final = "_graphs_Covariates_MAGMA_Phylum_shrinking",
        otu_table = otu_table_eliminated,
        Taxonomy_ASVs = ASV_final,
        continuous_or_binary = "BINARY",
        month = "",
        import_layout = LAYOUT
      )
      plot_graphical_phylum(
        OTU_RESULT = continuous_result,
        label_name_final = "_graphs_Covariates_MAGMA_Phylum_shrinking",
        otu_table = otu_table_eliminated,
        Taxonomy_ASVs = ASV_final,
        continuous_or_binary = "CONTINUOUS",
        month = "",
        import_layout = LAYOUT
        # substitute to the place where the "ORIGINAL" one is saved
      )
    
    },
    # ... but if an error occurs, tell me what happened: 
    error=function(error_message) {
      message("This is my custom message.")
      
      # The problem is that it breaks after having opened the png
      dev.off()
      plot_graphical_phylum(
        OTU_RESULT = binary_result,
        label_name_final = "_graphs_Covariates_MAGMA_Phylum_shrinking",
        otu_table = otu_table_eliminated,
        Taxonomy_ASVs = ASV_final,
        continuous_or_binary = "BINARY",
        month = "",
        import_layout = F
      )
      plot_graphical_phylum(
        OTU_RESULT = continuous_result,
        label_name_final = "_graphs_Covariates_MAGMA_Phylum_shrinking",
        otu_table = otu_table_eliminated,
        Taxonomy_ASVs = ASV_final,
        continuous_or_binary = "CONTINUOUS",
        month = "",
        import_layout = "LAYOUT_PHYLUM__graph.txt"
        # substitute to the place where the "ORIGINAL" one is saved
      )
  
    }
  )
  
  
  # ISN computing -----------------------------------------------------------
  
  
  
  
  
  setwd(file.path(curr_wd, "LooNet"))
  getwd()
  
  
  
  
  res_files = matching_list_creator("MAGMA_continuous_Data")
  files = res_files[[1]] ; files_name_matching = res_files[[2]]; names_file = res_files[[3]]
  
  global_net = read.table(file= "../MAGMA_continuous_Data_GLOBAL.tsv", sep = '\t', header = TRUE, row.names = 1)
  
  
  
  result_prep = preprocessing_global_net(global_net)
  
  global_net = result_prep[[1]]; global_net_vect= result_prep[[2]]; lionessOutput = result_prep[[3]];
  Sequences_nodes = result_prep[[4]]
  
  
  
  
  lionessOutput_fin =  ISN_computation(files = files, global_net = global_net, global_net_vect = global_net_vect,
                                       lionessOutput = lionessOutput, Sequences_nodes = Sequences_nodes,matching_string= "MAGMA_continuous_Data")
  # Calculation ISNs and average strenght -----------------------------------
  
  
  
  
  # WRITING -----------------------------------------------------------------
  
  dir.create(file.path(curr_wd,"ISNs"))
  setwd(file.path(curr_wd,"ISNs"))
  
  fwrite(lionessOutput_fin[,3:ncol(lionessOutput_fin)],file="Resulting_net_from_corr_MAGMACONF.txt",
         sep = " ", row.names = TRUE ) # keeps the rownames
  # 
  
  colnames(Sequences_nodes) = c("Node_list","NAME")
  write.table(Sequences_nodes,file = "Sequences_nodes.txt", quote = FALSE,row.names = F, col.names = T)
  # 
  write.table(files_name_matching,file = "files_name_matching.txt")
  Resulting_net = lionessOutput_fin[lionessOutput_fin$reg >= lionessOutput_fin$tar,] ## eliminated double rows
  
  fwrite(Resulting_net[,-c(1,2)],file="Resulting_net_from_corr_eliminate_MAGMACONF.txt",
         sep = " ", row.names = TRUE ) # keeps the rownames
  
  # Since is full of rows with just 0 --> we create a network with j --------
  
  
  aa = apply(Resulting_net[,-c(1,2)], 1, function(x) sum(abs(x))) > 0# )
  not_null = Resulting_net[aa,]
  dim(not_null)
  # [1] 585 81
  
  
  fwrite(not_null[,-c(1,2)],file="Resulting_net_notNULL_MAGMACONF.txt",
         sep = " ", row.names = TRUE ) # keeps the rownames
  
  # arg1 = 1
  # arg2 = 2
  # system(paste("cmd.exe",arg1,arg2), input = paste('"C:\\Program Files\\R\\R-4.0.3/bin/Rscript.exe" "C:/Users/fmelo/Desktop/Backup_Federico/work_microbiome/SPARCC_and_MAGMA/MAGMA_result_CORRECT_CONF_preprocessing_per_group/Pipeline_and_functions/script.R'))
  # args$arg1 = "c" 
  # args$arg2 = 2

  not_null_no_node = not_null[not_null$reg != not_null$tar,]
  fwrite(not_null_no_node[,-c(1,2)],file="Resulting_net_notNULL_NONODE_MAGMACONF.txt",
         sep = " ", row.names = TRUE ) # keeps the rownames
  getwd()


} else {print("NO ind")}
