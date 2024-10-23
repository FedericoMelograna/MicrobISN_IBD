library(vegan)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(xlsx)
library(ggpubr)
library(data.table)

## NB to do all, not only the significant ones, took out ../Significant in result_path and the test for p.value < 0.05

data_path   = "C:/Users/fmelo/Desktop/Microbiome_drug_repurpose_project/Data/"
result_path = "C:/Users/fmelo/Desktop/Microbiome_drug_repurpose_project/Alpha_Beta_diversity_results/Significant/"


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


if (cancellation == F){
  
  direttoria = paste0(direttoria, "_selected_microbes_no_canc_un_name")
} else {
  direttoria = paste0(direttoria, "_selected_microbes_canc")
  
}
# JUST different microbes -------------------------------------------------

## Known key microbes
Selected_microbes = c("Fusicatenibacter", "Roseburia", "Dorea", "Blautia", "Intestinimonas", "Adlercreutzia", "Clostridium_IV", "Parabacteroides", "Collinsella", 
  "Anaerostipes", "Subdoligranulum", "Butyricicoccus", "Oscillibacter", "Akkermansia", "Alistipes", "Ruminococcus2", "Ruminococcus", "Veillonella", 
  "Barnesiella", "Bifidobacterium")
ASV_taxonomic_annotation <- read.delim("C:/Users/fmelo/Desktop/Backup_Federico/Paddy_work_microbiome/ASV_taxonomic_annotation.txt")


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


# GRAPHS ------------------------------------------------------------------


vSHAN_DIV6M = vegan::diversity(OTU_TABLE_MAGMA_nonselected, index = "shannon", MARGIN = 1, base = exp(1))
vSHAN_DIV9M = vegan::diversity(OTU_TABLE_MAGMA_selected, index = "shannon", MARGIN = 1, base = exp(1))


# TEST --------------------------------------------------------------------

sort(names(vSHAN_DIV6M)) == sort(names(vSHAN_DIV9M))
vSHAN_DIV6M[sort(names(vSHAN_DIV6M))]; vSHAN_DIV9M[sort(names(vSHAN_DIV9M))]
(w4 = wilcox.test(vSHAN_DIV6M[sort(names(vSHAN_DIV6M))], vSHAN_DIV9M[sort(names(vSHAN_DIV9M))], paired = T))

print(w4$p.value)
if (w4$p.value > 0.05 ) { next;}








dir.create(paste0(result_path,"/", direttoria) )
setwd(paste0(result_path,"/",direttoria ))
f_ex = file.exists("Statistic_alpha.xlsx")

saveRDS(w4,  paste0("TEST_WILCOXON_age_", outcome_file,".rds"))



vplotdf = data.frame("Alpha_div" = c(vSHAN_DIV6M, vSHAN_DIV9M), "Responder" = c(rep("Non selected", nrow(OTU_TABLE_MAGMA_nonselected)),rep("Selected", nrow(OTU_TABLE_MAGMA_selected)))) 
# Basic violin plot
p <- ggplot(vplotdf, aes(x=Responder, y=Alpha_div)) + 
  geom_violin()
p

ggplot(data.frame(Trait=vplotdf$Alpha_div, Cluster=as.factor(vplotdf$Responder)), aes(x=Cluster, y=Trait, fill = Cluster)) + 
  geom_violin() + ggtitle("Violin plot Shannon diversity") +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1), geom="pointrange",shape=3, size=3, color="gray") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_gray(base_size = 17, base_line_size = 1.1)+
  labs(y= expression(alpha*" diversity Shannon"), x = expression("Responder"), fill = "")+ #, x = "x axis name")
  theme(axis.text=element_text(face="bold"),
        axis.title=element_text(face="bold"), title = element_text(face = "bold",),text = element_text(size = 35),legend.position="none")


# Graphs  violin plot  ----------------------------------------------------------------
dir.create(paste0(result_path,"/", direttoria) )
setwd(paste0(result_path,"/",direttoria ))
f_ex = file.exists("Statistic_alpha.xlsx")

vplotdf%>%
  group_by(Responder)%>% 
  summarise(Mean=mean(Alpha_div), Max=max(Alpha_div), Min=min(Alpha_div), Median=median(Alpha_div), Std=sd(Alpha_div), IQR = IQR(Alpha_div)) %>%
  write.xlsx(.,file = paste0("Statistic_alpha", ".xlsx"),
             sheetName = as.character(outcome_file), append = f_ex)


outcome_name = gsub("_outcome_combined", "", outcome_file)
low_l = 0.3 ; high_l = 5
v_age = ggplot(data.frame(Trait=vplotdf$Alpha_div, Cluster=as.factor(vplotdf$Responder)), aes(x=Cluster, y=Trait, fill = Cluster)) + 
  geom_violin() + ggtitle(paste0("Shannon diversity ", outcome_name)) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1), geom="pointrange", size=2, color="gray") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_minimal(base_size = 17, base_line_size = 1.1)+
  coord_cartesian(ylim = c(low_l, high_l)) +
  labs(y= expression(alpha*" diversity Shannon"), x = "Responder", fill = "")+ #, x = "x axis name")
  theme(axis.text=element_text(face="bold"),
        axis.title=element_text(face="bold"), title = element_text(face = "bold",),text = element_text(size = 35),legend.position="none")
png(paste0("Violin_Shannon_alphadiv_",outcome_name,".png"), width = 300, height = 300, units='mm', res = 300)
print(v_age)
dev.off()
png(paste0("Violin_Shannon_alphadiv_ds_",outcome_name,".png"), width = 465, height = 225, units='mm', res = 300)
print(v_age)
dev.off()





# Beta diversity ----------------------------------------------------------



library(compositions)
library(mixOmics)


OTU_TABLE_MAGMA_nonselected_CLR = data.frame( clr( OTU_TABLE_MAGMA_nonselected ))
OTU_TABLE_MAGMA_selected_CLR = data.frame( clr( OTU_TABLE_MAGMA_selected ))

round(apply(OTU_TABLE_MAGMA_selected_CLR,1,sum),2) ; round(apply(OTU_TABLE_MAGMA_nonselected_CLR,1,sum),2)



PCA_non_selected_CLR = mixOmics::pca(OTU_TABLE_MAGMA_nonselected_CLR)

plotVar(PCA_non_selected_CLR)
plot(PCA_non_selected_CLR)

PCA_selected_CLR = mixOmics::pca(OTU_TABLE_MAGMA_selected_CLR)

plotVar(PCA_selected_CLR)
plot(PCA_selected_CLR)


# NOVEL Beta diversity analysis  -----------------------------------------------------



round(apply(OTU_TABLE_MAGMA_selected_CLR,1,sum),2) ; round(apply(OTU_TABLE_MAGMA_nonselected_CLR,1,sum),2)


OTU_TABLE_MAGMA_selected_CLR

OTUS = cbind(OTU_TABLE_MAGMA_nonselected_CLR,OTU_TABLE_MAGMA_selected_CLR)
PCA_CLR = mixOmics::pca(OTUS)

png(paste0("Beta_diversity_togheter_",outcome_file,".png"), width = 300, height = 300, units='mm', res = 300)
p_age = plotIndiv(PCA_CLR, ind.names = FALSE,
                  group = c(rep("Together",nrow(OTUS))),
                  pch = 16,
                  legend = TRUE, title = paste0("Together ", outcome_name),# 'Beta diversity: PCA comp 1 - 2',
                  legend.title = '', legend.title.pch = 'Diet', legend.position = "top",cex = 3
                  , size.title = 30, point.lwd = 2, size.legend = 25, size.legend.title = 30, size.xlabel = 25, size.ylabel = 25, size.axis = 18)

dev.off()

png(paste0("Beta_diversity_SELECTED_",outcome_file,".png"), width = 300, height = 300, units='mm', res = 300)
p_age = plotIndiv(PCA_selected_CLR, ind.names = FALSE,
                  group = c(rep("Selected",nrow(OTUS))), 
                  pch = 16,
                  legend = TRUE, title = paste0("Selected ", outcome_name),# 'Beta diversity: PCA comp 1 - 2',
                  legend.title = '', legend.title.pch = 'Diet', legend.position = "top",cex = 3
                  , size.title = 30, point.lwd = 2, size.legend = 25, size.legend.title = 30, size.xlabel = 25, size.ylabel = 25, size.axis = 18)

dev.off()

png(paste0("Beta_diversity_NON_SELECTED_",outcome_file,".png"), width = 300, height = 300, units='mm', res = 300)
p_age = plotIndiv(PCA_non_selected_CLR, ind.names = FALSE,
                  group = c(rep("Non selected",nrow(OTUS))), 
                  pch = 16,
                  legend = TRUE, title = paste0("Non selected ", outcome_name),# 'Beta diversity: PCA comp 1 - 2',
                  legend.title = '', legend.title.pch = 'Diet', legend.position = "top",cex = 3
                  , size.title = 30, point.lwd = 2, size.legend = 25, size.legend.title = 30, size.xlabel = 25, size.ylabel = 25, size.axis = 18)

dev.off()





