setwd("---")

files <- list.files(path=".", pattern="group*", 
                    full.names=TRUE, recursive=FALSE)
curr_wd = getwd()

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
  tab = tab[tab$P.Value< 0.05,]
  if (nrow(tab) == 0 ) {next;}
  final_el = gsub("./","",el$files)
  setwd(curr_wd) 
  write.table(tab, file = paste0("LIMMA_025_result_correct/",final_el,"_",el$outcome,"_pval_005.txt"),quote = F, sep = "\t", col.names = T, row.names = F)
}


