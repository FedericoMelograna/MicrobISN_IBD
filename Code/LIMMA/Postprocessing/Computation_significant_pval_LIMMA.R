
setwd("---")
files <- list.files(path=".", pattern="group*", 
                    full.names=TRUE, recursive=FALSE)
curr_wd = getwd()

files
el
p_val_thr = 0.05
outcomes = rep(c("Clinical_outcome_combined", "Biomarker_outcome_combined","Endoscopic_outcome_combined" ), length(files) )
files_per_outcome = rep(files, each = 3 )
data_frame_loop  = data.frame(files = files_per_outcome, outcome = outcomes)

args = list()
args$direttoria = NA
### LIMMA 0.25


setwd(curr_wd) 

dir.create("LIMMA_025_result_correct")
for (i in 1:nrow(data_frame_loop)){
  el = data_frame_loop[i,]
  setwd(curr_wd) 
  
  setwd(paste0(el,"/Results_LIMMA_CORRECT_025_",el$outcome, "/") );    
  if (!file.exists("MAGMA_CONF_all_edges_withCORRECT_TAXA_pval_005.txt")){next;}
  
  tab <- read.delim("MAGMA_CONF_all_edges_withCORRECT_TAXA_pval_005.txt")
  tab = tab[tab$P.Value< p_val_thr,]
  if (nrow(tab) == 0 ) {next;}
  ASV_taxonomic_annotation <- read.delim("C:/Users/fmelo/Desktop/.../ASV_taxonomic_annotation.txt")
  tab$Node2
  tab2 = tab
  # IT MEANS it is NOT ALREADY with FAMILY and GENUS
  if (!"Gene_n2" %in% colnames(tab2) | !"Gene_n1" %in% colnames(tab2)){
    
    tab2$Family_n2 = NA ; 
    tab2$Family_n1 = NA ; 
    tab2$Gene_n2 = NA ; 
    tab2$Gene_n1 = NA ; 
    for (j in 1: nrow(tab)){
      node1_f_g = ASV_taxonomic_annotation[ASV_taxonomic_annotation$ASV == tab[j,"Node1_Taxa"], c("Family", "Genus")]
      node2_f_g = ASV_taxonomic_annotation[ASV_taxonomic_annotation$ASV == tab[j,"Node2_Taxa"], c("Family", "Genus")]
      tab2$Family_n1[j] = node1_f_g$Family ; 
      tab2$Family_n2[j] = node2_f_g$Family ; 
      tab2$Gene_n1[j] = node1_f_g$Genus ; 
      tab2$Gene_n2[j] = node2_f_g$Genus; 
    }
    
  }
  final_el = gsub("./","",el$files)
  setwd(curr_wd) 
  write.table(tab2, file = paste0("LIMMA_025_result_correct/",final_el,"_",el$outcome,"_pval_005.txt"),quote = F, sep = "\t", col.names = T, row.names = F)
}


