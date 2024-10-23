library(vegan)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(xlsx)
library(ggpubr)
library(data.table)

data_path   = "C:/Users/fmelo/Desktop/Microbiome_drug_repurpose_project/Data/"
result_path = "C:/Users/fmelo/Desktop/Microbiome_drug_repurpose_project/Alpha_Beta_diversity_results/Significant/"

## NB to do all, not only the significants, took out ../Significant in result_path and the test for p.value < 0.05



# FUNCTIONS ---------------------------------------------------------------
checkStrict <- function(f, silent=FALSE) {
  vars <- codetools::findGlobals(f)
  found <- !vapply(vars, exists, logical(1), envir=as.environment(2))
  if (!silent && any(found)) {
    warning("global variables used: ", paste(names(found)[found], collapse=', '))
    return(invisible(FALSE))
  }
  
  !any(found)
}

calculate_violin_plot_test_del = function(dataset, matching_data, ext_pheno = "Delivery", name_g,dataset2 = NULL, low_l = 0.3, high_l = 3){
  if (is.null(dataset2)){ ## ONLY 1 dataset --> 1 Responder
    vSHAN = vegan::diversity(dataset, index = "shannon", MARGIN = 1, base = exp(1))
    if (ext_pheno == "Delivery"){ ## DIFFERENT per phenotype
      is_variable = ifelse(names(vSHAN) %in% matching_data$Child[matching_data$delivery_type.x == "Vaginal"],"Vaginal", "C-section")
    } else if (ext_pheno == "Diet"){
      is_same = matching_data$Child[matching_data$Diet_TP.x == matching_data$Diet_TP.y] %>% na.omit()
      is_variable = ifelse(names(vSHAN) %in% is_same,"Persist.", "Non persist.")
    }
    
  } else { ## two datasts
    vSHAN1 = vegan::diversity(dataset, index = "shannon", MARGIN = 1, base = exp(1))
    vSHAN2 = vegan::diversity(dataset2, index = "shannon", MARGIN = 1, base = exp(1))
    if (ext_pheno == "Delivery"){
      is_variable1 = ifelse(names(vSHAN1) %in% matching_data$Child[matching_data$delivery_type.x == "Vaginal"],"Vaginal", "C-section")
      is_variable2 = ifelse(names(vSHAN2) %in% matching_data$Child[matching_data$delivery_type.x == "Vaginal"],"Vaginal", "C-section")
    } else if (ext_pheno == "Diet"){
      is_same = matching_data$Child[matching_data$Diet_TP.x == matching_data$Diet_TP.y] %>% na.omit()
      is_variable1 = ifelse(names(vSHAN1) %in% is_same,"Persist.", "Non persist.")
      is_variable2 = ifelse(names(vSHAN2) %in% is_same,"Persist.", "Non persist.")
    }
    vSHAN = c(vSHAN1, vSHAN2) ; is_variable = c(is_variable1, is_variable2)
  }
 
  vplotdf = data.frame("Alpha_div" = vSHAN, "Type_of_Delivery" = is_variable) 
  
  graph_name = ifelse(ext_pheno == "Delivery", "Mode of Delivery", "Diet")
  
  g1 = ggplot(data.frame(Trait=vplotdf$Alpha_div, Cluster=as.factor(vplotdf$Type_of_Delivery)), aes(x=Cluster, y=Trait, fill = Cluster)) + 
    geom_violin() + ggtitle(paste0("Violin plot ",name_g, " Shannon diversity")) +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1), geom="pointrange", size=2, color="gray") +
    # DOT IS THE MEAN ; the extremities are + and - SD
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme_minimal(base_size = 17, base_line_size = 1.1)+
    labs(y= expression(alpha*" diversity Shannon"), x = graph_name, fill = "")+ #, x = "x axis name")
    coord_cartesian(ylim = c(low_l, high_l)) +
    theme(axis.text=element_text(face="bold"),
          axis.title=element_text(face="bold"), title = element_text(face = "bold",),text = element_text(size = 35),legend.position="none")
  
  if (ext_pheno == "Delivery"){
    (w4 = wilcox.test(vSHAN[is_variable == "Vaginal"], vSHAN[is_variable != "Vaginal"], paired = F))
  } else if (ext_pheno == "Diet"){
    (w4 = wilcox.test(vSHAN[is_variable == "Persist."], vSHAN[is_variable != "Persist."], paired = F))
    
  }
  return(list(graph = g1, test = w4, summ =vplotdf))
}
checkStrict(calculate_violin_plot_test_del)

# ALPHA and BETA diversity  -----------------------------------------------

setwd(data_path)


direttoria = args$direttoria
outcome_file = args$outcome
cancellation = args$cancelled
setwd(path_base)

### DATA 
otu_table = read.table(paste0("./", direttoria, "/Data/otu_table_selected_",direttoria, ".tsv" ))

if (grepl("14", direttoria)){
  metadata = read.delim("C:/Users/fmelo/Desktop/Microbiome_drug_repurpose_project/Data/w14w24_Data/w14_n188_metadata.txt")
  merged_outcome =  metadata[metadata$FC.nummer %in% rownames(otu_table), ]
  
  
} else if (grepl("24", direttoria)){
  metadata = read.delim("C:/Users/fmelo/Desktop/Microbiome_drug_repurpose_project/Data/w14w24_Data/w24_n160_metadata.txt")
  merged_outcome =  metadata[metadata$FC.nummer %in% rownames(otu_table), ]
  
  
  
} else {
  metadata <- read.delim("C:/Users/fmelo/Desktop/Microbiome_drug_repurpose_project/Data/w0_n335_metadata.txt")
  merged_outcome =  metadata[metadata$FC.nummer %in% rownames(otu_table), ]
  direttoria = paste0(direttoria,"_w0" )
}

if (sum(rownames(otu_table) %in% merged_outcome$FC.nummer) != nrow(otu_table)){
  abort("MISMATCH --> the metadata does not have all the individual ")
}


if (cancellation == F){ ## If we are using the cancellation or not
  
  direttoria = paste0(direttoria, "_selected_microbes_no_canc_un_name")
} else {
  direttoria = paste0(direttoria, "_selected_microbes_canc")
  
}

direttoria = paste0(direttoria, "_divided_responder")
# JUST different microbes -------------------------------------------------

## KNOWN key microbes
Selected_microbes = c("Fusicatenibacter", "Roseburia", "Dorea", "Blautia", "Intestinimonas", "Adlercreutzia", "Clostridium_IV", "Parabacteroides", "Collinsella", 
                      "Anaerostipes", "Subdoligranulum", "Butyricicoccus", "Oscillibacter", "Akkermansia", "Alistipes", "Ruminococcus2", "Ruminococcus", "Veillonella", 
                      "Barnesiella", "Bifidobacterium")
ASV_taxonomic_annotation <- read.delim("./ASV_taxonomic_annotation.txt")


setwd("C:/Users/fmelo/Desktop/Microbiome_drug_repurpose_project/")
files <- list.files(path=".", pattern="group*", 
                    full.names=TRUE, recursive=FALSE)
curr_wd = getwd()

ASV_reduced = ASV_taxonomic_annotation[ASV_taxonomic_annotation$ASV %in% colnames(otu_table),]
ASV_sel = ASV_reduced[ASV_reduced$Genus %in% Selected_microbes, ]
ASV_nonsel = ASV_reduced[!(ASV_reduced$Genus %in% Selected_microbes), ]

OTU_TABLE_MAGMA_selected = otu_table  %>% dplyr::select(ASV_sel$ASV)
OTU_TABLE_MAGMA_nonselected = otu_table %>% dplyr::select(ASV_nonsel$ASV)
if (cancellation == T){
  ASV_nonsel = ASV_reduced[!(ASV_reduced$Genus %in% Selected_microbes) & !(is.na(ASV_reduced$Genus)), ]
  OTU_TABLE_MAGMA_nonselected = otu_table %>% dplyr::select(ASV_nonsel$ASV)
}



# OUTCOME for selected and not --------------------------------------------


is_missing_ = is.na(merged_outcome %>% dplyr::select(outcome_file))
otu_table_no_miss_SELECTED = OTU_TABLE_MAGMA_selected %>% filter(!is_missing_)
merged_outcome_no_miss = merged_outcome %>% filter(!is_missing_)


outcome = rownames(otu_table_no_miss_SELECTED) %in% merged_outcome_no_miss$FC.nummer[merged_outcome_no_miss %>% dplyr::select(c(outcome_file)) == 1]


OTU_TABLE_MAGMA_no_remission_selected <- otu_table_no_miss_SELECTED[!outcome,] 
# 69 x 95
OTU_TABLE_MAGMA_remission_selected <- otu_table_no_miss_SELECTED[outcome,]

# NON selected ------------------------------------------------------------


is_missing_ = is.na(merged_outcome %>% dplyr::select(outcome_file))
otu_table_no_miss_NONSELECTED = OTU_TABLE_MAGMA_nonselected %>% filter(!is_missing_)
merged_outcome_no_miss = merged_outcome %>% filter(!is_missing_)


outcome = rownames(otu_table_no_miss_SELECTED) %in% merged_outcome_no_miss$FC.nummer[merged_outcome_no_miss %>% dplyr::select(c(outcome_file)) == 1]


OTU_TABLE_MAGMA_no_remission_nonselected <- otu_table_no_miss_NONSELECTED[!outcome,] 
OTU_TABLE_MAGMA_remission_nonselected <- otu_table_no_miss_NONSELECTED[outcome,]


# GRAPHS ------------------------------------------------------------------
# for selected ------------------------------------------------------------


vSHAN_DIV6M = vegan::diversity(OTU_TABLE_MAGMA_no_remission_selected, index = "shannon", MARGIN = 1, base = exp(1))
vSHAN_DIV9M = vegan::diversity(OTU_TABLE_MAGMA_remission_selected, index = "shannon", MARGIN = 1, base = exp(1))

vSHAN_DIV6M

# TEST --------------------------------------------------------------------


(w3 = wilcox.test(vSHAN_DIV6M,vSHAN_DIV9M, paired = F)  )

print(w3$p.value)
w3_selected = w3
if (w3$p.value <= 0.05 ) {
  dir.create(paste0(result_path,"/", direttoria) )
  setwd(paste0(result_path,"/",direttoria ))
  f_ex = file.exists("Statistic_alpha.xlsx")
  
  saveRDS(w3, paste0("TEST_WILCOXON_age_SELECTED_", outcome_file,".rds"))
  
  
  vplotdf = data.frame("Alpha_div" = c(vSHAN_DIV6M, vSHAN_DIV9M), "Responder" = c(rep("No remission", nrow(OTU_TABLE_MAGMA_no_remission_selected)),rep("Remission", nrow(OTU_TABLE_MAGMA_remission_selected)))) 
  
  
  # mult = 1 ----------------------------------------------------------------
  dir.create(paste0(result_path,"/", direttoria) )
  setwd(paste0(result_path,"/",direttoria ))
  f_ex = file.exists("Statistic_alpha.xlsx")
  
  vplotdf%>%
    group_by(Responder)%>% 
    summarise(Mean=mean(Alpha_div), Max=max(Alpha_div), Min=min(Alpha_div), Median=median(Alpha_div), Std=sd(Alpha_div), IQR = IQR(Alpha_div)) %>%
    write.xlsx(.,file = paste0("Statistic_alpha_Selected", ".xlsx"),
               sheetName = as.character(outcome_file), append = f_ex)
  
  
  outcome_name = gsub("_outcome_combined", "", outcome_file)
  low_l = 0.3 ; high_l = 5
  v_age = ggplot(data.frame(Trait=vplotdf$Alpha_div, Cluster=as.factor(vplotdf$Responder)), aes(x=Cluster, y=Trait, fill = Cluster)) + 
    geom_violin() + ggtitle(paste0("Shannon diversity Sel. ", outcome_name)) +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1), geom="pointrange", size=2, color="gray") +
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme_minimal(base_size = 17, base_line_size = 1.1)+
    coord_cartesian(ylim = c(low_l, high_l)) +
    labs(y= expression(alpha*" diversity Shannon"), x = "Responder", fill = "")+ #, x = "x axis name")
    theme(axis.text=element_text(face="bold"),
          axis.title=element_text(face="bold"), title = element_text(face = "bold",),text = element_text(size = 35),legend.position="none")
  png(paste0("Violin_Shannon_SELECTED_alphadiv_",outcome_name,".png"), width = 300, height = 300, units='mm', res = 300)
  print(v_age)
  dev.off()
  png(paste0("Violin_Shannon_SELECTED_alphadiv_ds_",outcome_name,".png"), width = 465, height = 225, units='mm', res = 300)
  print(v_age)
  dev.off()
  }



# Non selected ------------------------------------------------------------





vSHAN_DIV6M = vegan::diversity(OTU_TABLE_MAGMA_no_remission_nonselected, index = "shannon", MARGIN = 1, base = exp(1))
vSHAN_DIV9M = vegan::diversity(OTU_TABLE_MAGMA_remission_nonselected, index = "shannon", MARGIN = 1, base = exp(1))

vSHAN_DIV6M

w3 = wilcox.test(vSHAN_DIV6M,vSHAN_DIV9M, paired = F) 

print(w3$p.value)
w3_nonselected = w3

if (w3$p.value <= 0.05){
  dir.create(paste0(result_path,"/", direttoria) )
  setwd(paste0(result_path,"/",direttoria ))
  f_ex = file.exists("Statistic_alpha.xlsx")
  
  saveRDS(w3, paste0("TEST_WILCOXON_age_SELECTED_", outcome_file,".rds"))
  
  
  vplotdf = data.frame("Alpha_div" = c(vSHAN_DIV6M, vSHAN_DIV9M), "Responder" = c(rep("No remission", nrow(OTU_TABLE_MAGMA_no_remission_nonselected)),rep("Remission", nrow(OTU_TABLE_MAGMA_remission_nonselected)))) 
  
  
  vplotdf%>%
    group_by(Responder)%>% 
    summarise(Mean=mean(Alpha_div), Max=max(Alpha_div), Min=min(Alpha_div), Median=median(Alpha_div), Std=sd(Alpha_div), IQR = IQR(Alpha_div)) %>%
    write.xlsx(.,file = paste0("Statistic_alpha_NonSelected", ".xlsx"),
               sheetName = as.character(outcome_file), append = f_ex)
  
  
  outcome_name = gsub("_outcome_combined", "", outcome_file)
  low_l = 0.3 ; high_l = 5
  v_age = ggplot(data.frame(Trait=vplotdf$Alpha_div, Cluster=as.factor(vplotdf$Responder)), aes(x=Cluster, y=Trait, fill = Cluster)) + 
    geom_violin() + ggtitle(paste0("Shannon diversity Non Sel. ", outcome_name)) +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1), geom="pointrange", size=2, color="gray") +
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme_minimal(base_size = 17, base_line_size = 1.1)+
    coord_cartesian(ylim = c(low_l, high_l)) +
    labs(y= expression(alpha*" diversity Shannon"), x = "Responder", fill = "")+ #, x = "x axis name")
    theme(axis.text=element_text(face="bold"),
          axis.title=element_text(face="bold"), title = element_text(face = "bold",),text = element_text(size = 35),legend.position="none")
  png(paste0("Violin_Shannon_NonSELECTED_alphadiv_",outcome_name,".png"), width = 300, height = 300, units='mm', res = 300)
  print(v_age)
  dev.off()
  png(paste0("Violin_Shannon_NonSELECTED_alphadiv_ds_",outcome_name,".png"), width = 465, height = 225, units='mm', res = 300)
  print(v_age)
  dev.off()
  
  
}





# Beta diversity ----------------------------------------------------------

if (w3_nonselected$p.value > 0.05 && w3_selected$p.value > 0.05){
  next;
}
## SELECTED --------

library(compositions)
library(mixOmics)



vSHAN_DIV6M = vegan::diversity(OTU_TABLE_MAGMA_no_remission_selected, index = "shannon", MARGIN = 1, base = exp(1))
vSHAN_DIV9M = vegan::diversity(OTU_TABLE_MAGMA_remission_selected, index = "shannon", MARGIN = 1, base = exp(1))

OTU_TABLE_MAGMA_nonremission_CLR = data.frame( clr( OTU_TABLE_MAGMA_no_remission_selected ))
OTU_TABLE_MAGMA_remission_CLR = data.frame( clr( OTU_TABLE_MAGMA_remission_selected ))


PCA_non_remission_CLR = mixOmics::pca(OTU_TABLE_MAGMA_nonremission_CLR)

PCA_remission_CLR = mixOmics::pca(OTU_TABLE_MAGMA_remission_CLR)


# NOVEL analysis as per Pender's request  -----------------------------------------------------




OTUS = rbind(OTU_TABLE_MAGMA_nonremission_CLR,OTU_TABLE_MAGMA_remission_CLR)
PCA_CLR = mixOmics::pca(OTUS)

png(paste0("Beta_diversity_SELECTED_microbes__",outcome_file,".png"), width = 300, height = 300, units='mm', res = 300)
p_age = plotIndiv(PCA_CLR, ind.names = FALSE,
                  group = c(rep("No remission", nrow(OTU_TABLE_MAGMA_nonremission_CLR)),rep("Remission", nrow(OTU_TABLE_MAGMA_remission_CLR))),
                  pch = 16, #rep(as.numeric(as.factor(Children_6M$delivery_type.x))+15,2),
                  legend = TRUE, title = paste0("Sel. Microbes ", outcome_name),# 'Beta diversity: PCA comp 1 - 2',
                  legend.title = '', legend.title.pch = 'Diet', legend.position = "top",cex = 3
                  , size.title = 30, point.lwd = 2, size.legend = 25, size.legend.title = 30, size.xlabel = 25, size.ylabel = 25, size.axis = 18)

dev.off()




## NONSELECTED --------


OTU_TABLE_MAGMA_nonremission_CLR = data.frame( clr( OTU_TABLE_MAGMA_no_remission_nonselected ))
OTU_TABLE_MAGMA_remission_CLR = data.frame( clr( OTU_TABLE_MAGMA_remission_nonselected ))

PCA_non_remission_CLR = mixOmics::pca(OTU_TABLE_MAGMA_nonremission_CLR)
PCA_remission_CLR = mixOmics::pca(OTU_TABLE_MAGMA_remission_CLR)



OTUS = rbind(OTU_TABLE_MAGMA_nonremission_CLR,OTU_TABLE_MAGMA_remission_CLR)
PCA_CLR = mixOmics::pca(OTUS)

png(paste0("Beta_diversity_NONSELECTED_microbes__",outcome_file,".png"), width = 300, height = 300, units='mm', res = 300)
p_age = plotIndiv(PCA_CLR, ind.names = FALSE,
                  group = c(rep("No remission", nrow(OTU_TABLE_MAGMA_nonremission_CLR)),rep("Remission", nrow(OTU_TABLE_MAGMA_remission_CLR))),
                  pch = 16,#rep(as.numeric(as.factor(Children_6M$delivery_type.x))+15,2),
                  legend = TRUE, title = paste0("Non sel. Microbes ", outcome_name),# 'Beta diversity: PCA comp 1 - 2',
                  legend.title = '', legend.title.pch = 'Diet', legend.position = "top",cex = 3
                  , size.title = 30, point.lwd = 2, size.legend = 25, size.legend.title = 30, size.xlabel = 25, size.ylabel = 25, size.axis = 18)

dev.off()



