

ASV_taxonomic_annotation <- read.delim(".../ASV_taxonomic_annotation.txt")



setwd(paste0(base_dir, method_dir))
cwd = getwd()



files <- list.files(path=".", pattern=paste0(outcome, "*"), 
                    full.names=TRUE, recursive=FALSE)

eliminated_files = files[!grepl("Not", files)]
eliminated_files_t1 = gsub(paste0("_eliminate_",outcome,".+$"), "", eliminated_files)

group_names = gsub("./group_", "", eliminated_files_t1)

df = data.frame(col.names = group_names ) 
df = matrix(0, ncol = length(group_names)+1)
colnames(df) 
colnames(df) = c("Taxa",group_names)
df = data.frame(df)

if (level == "Taxa"){
  v1 = "Node1_Taxa" ; v2 = "Node2_Taxa"
} else if (level == "Genus"){
  v1 = "Gene_n1" ; v2 = "Gene_n2"
} else if (level == "Family"){
  v1 = "Family_n1" ; v2 = "Family_n2"
} 

for (i in 1:length(group_names)){
  results = read.delim(eliminated_files[i])
  tt = table(c(results[,c(v1)],results[,c(v2)]))
  dd = data.frame(tt)
  colnames(dd) = c("Taxa", "Freq")
  dd$Taxa = as.character(dd$Taxa)
  for (j in 1:nrow(dd)){
    if (!(dd$Taxa[j] %in% df$Taxa)) {
      df[nrow(df) + 1,] = c(as.character(dd$Taxa[j]), rep(0,ncol(df)-1))
    }
    df[df$Taxa == dd$Taxa[j],group_names[i]] = as.numeric(df[df$Taxa == dd$Taxa[j],group_names[i]]) + dd$Freq[j]
  }
  
}


dir.create(file.path("Results"), showWarnings = FALSE)
setwd(file.path("Results"))

sum_taxas = apply(data.frame(lapply(df[,-1],as.numeric)), 1, sum)
appearance= apply(data.frame(lapply(df[,-1],as.numeric)), 1, function(c)sum(c!=0))
DF = cbind(df, sum_taxas, appearance)
if (level == "Taxa"){
DF_with_GENE = merge(DF, ASV_taxonomic_annotation[,c("ASV","Genus","Family")], by.x = "Taxa", by.y = "ASV")
} else if (level == "Genus"){
  ASV_taxonomic_annotation_dist =distinct(ASV_taxonomic_annotation, Genus, .keep_all = TRUE)
  DF_with_GENE = merge(DF, ASV_taxonomic_annotation_dist[,c("Genus","Family")], by.x = "Taxa", by.y = "Genus")
} else if (level == "Family"){
  DF_with_GENE = DF
} 

DF_with_GENE_ORDERED = DF_with_GENE[order(DF_with_GENE$sum_taxas,decreasing = T),]

if (level == "Genus" | level == "Family"){
  colnames(DF_with_GENE_ORDERED)[1] = as.character(level)
}

write.table(x = DF_with_GENE_ORDERED, file = paste0("Final_LIST_",level,"_",outcome,"_",method_dir,".txt" ),quote = F,sep = "\t",col.names = T, row.names = F)
library("xlsx")
write.xlsx(x = DF_with_GENE_ORDERED, file = paste0("Final_LIST_",level,"_",outcome,"_",method_dir,".xlsx" ), sheetName = "Sheet1", 
           col.names = TRUE, row.names = FALSE, append = FALSE)

setwd(cwd)
