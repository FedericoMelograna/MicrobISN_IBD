library(vegan)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(xlsx)
library(ggpubr)
library(data.table)
library(scales)

data_path   = "----"
result_path = "----"

library(RColorBrewer)


  
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

computing_appearing_rel = function(otu_table, Taxonomy_ASVs,  level_taxa_order, na_rmv = T){
  
  level_dataset = data.frame(name = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), rank = seq(1,7))
  
  otu_table_rel<-as.data.frame(t(apply(otu_table,1, function(x) x/sum(x))))
  rank_current_taxa = level_dataset[level_taxa_order == level_dataset$name, "rank"]
  summed_col = colMeans(otu_table_rel)
  
  ASV_final = Taxonomy_ASVs[Taxonomy_ASVs$ASV %in% colnames(otu_table_rel),]
  df = (t(ASV_final))
  dim(df)
  colnames(df) = df[1,]
  df = df[-1,]
  level_taxa = unlist(unique(ASV_final %>% dplyr::select (level_taxa_order)))
  if (na_rmv) {
    level_taxa = level_taxa[!is.na(level_taxa)]
  } else {
    level_taxa[is.na(level_taxa)] = "<Unnamed>"
  }
  Taxas_name = colnames(otu_table_rel)
  
  phy_sum = matrix(data = 0, nrow = 1, ncol = length(level_taxa))
  colnames(phy_sum) = level_taxa
  rownames(phy_sum) = c("count")
  phy_cont = phy_sum
  for (i in 1:ncol(otu_table)){
    if (na_rmv) {
      if (is.na(df[rank_current_taxa,Taxas_name[i]])) {next;}
      level_taxa = df[rank_current_taxa,Taxas_name[i]]
    } else {
      phyl_temp = ifelse(is.na(df[rank_current_taxa,Taxas_name[i]]), "<Unnamed>", df[rank_current_taxa,Taxas_name[i]]) 
    } 
    phyl_temp = ifelse(is.na(df[rank_current_taxa,Taxas_name[i]]), "<Unnamed>", df[rank_current_taxa,Taxas_name[i]]) 
    phy_sum[1,phyl_temp ] = phy_sum[1,phyl_temp ] + summed_col[Taxas_name[i]]
    phy_cont[1,phyl_temp ] = phy_cont[1,phyl_temp ] +1 
    
  }
  taxas_per_phylum = table(df[rank_current_taxa,])
  ordered_taxa_per_phylum = taxas_per_phylum[colnames(phy_cont)]
  return(list(phy_sum= phy_sum, phy_cont = phy_cont) ) 
}

checkStrict(computing_appearing_rel)


setwd(data_path)

direttoria = args$direttoria
hierarch_component = args$compon #  "Family" ## One of Class/Family/Phylum
base_size = ifelse(args$compon== "Phylum", 25, 15) ## If phylum, 25, otherwise 15

setwd(path_base)

### DATA 
otu_table = read.table(paste0("./", direttoria, "/Data/otu_table_selected_",direttoria, ".tsv" ))

ASV_taxonomic_annotation <- read.delim(".../ASV_taxonomic_annotation.txt")
ASV_reduced = ASV_taxonomic_annotation[ASV_taxonomic_annotation$ASV %in% colnames(otu_table),]

setwd(result_path)



dist_Phyl = ASV_reduced %>% dplyr::distinct(!!as.name(hierarch_component))
print(dist_Phyl)

colnames(dist_Phyl) = "Compon" 

subset_coloring_df = coloring_df %>% dplyr::inner_join(dist_Phyl, by = "Compon")




list_res = computing_appearing_rel(otu_table = otu_table, Taxonomy_ASVs = ASV_reduced, level_taxa_order = hierarch_component)
phy_cont = list_res$phy_cont ; phy_sum = list_res$phy_sum

# We use the frac abundance ------------------------------------


df_phy = data.frame("group" = colnames(phy_sum), "value" = t(phy_sum))
ord_df_phy = df_phy[order(df_phy$group),]
ord_df_phy
print(ord_df_phy)


Change_unnamed_NA = ifelse(ord_df_phy$group =="<Unnamed>", NA, ord_df_phy$group)
right_colours = subset_coloring_df$Colours[ match(subset_coloring_df$Compon, Change_unnamed_NA,nomatch = 0)]



gsub_dir = gsub("eliminate ", "",gsub("_"," ",gsub("group_","", direttoria)))
if (!grepl("14|24", gsub_dir)){
  gsub_dir = paste0(gsub("eliminate", "", gsub_dir), "w0")
  direttoria = paste0(direttoria, "_w0")
}


# GRAPHS ------------------------------------------------------------------


bp<- ggplot(ord_df_phy, aes(x="", y=count, fill=group))+
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values=right_colours, breaks=ord_df_phy$group)+
  labs(y= "Fractional abundances", x = "Cohort",fill='') +
  ggtitle(paste0("Barplot - ",gsub_dir ))+
  theme_classic(base_size = base_size, base_line_size = 1.1)+
  theme(axis.text=element_text(face="bold"),
        axis.title=element_text(face="bold"), title = element_text(face = "bold"),legend.text=element_text(size=base_size-5))
bp

bp2<- bp +   geom_text(aes(label=paste0(sprintf("%1.f", count*100),"%")),
                       position=position_stack(vjust=0.5), fontface = "bold")





pie <- ggplot(ord_df_phy, aes(x="", y=count, fill=group))+
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values=right_colours, breaks=ord_df_phy$group)+
  labs(y= "Fractional abundances", x = "",fill='') +
  ggtitle(paste0("Pie Chart - ",gsub_dir ))+
  theme_minimal(base_size = base_size, base_line_size = 1.1)+
  theme(axis.text=element_text(face="bold"),
        axis.title=element_text(face="bold"), title = element_text(face = "bold"),legend.text=element_text(size=base_size-5))+
  coord_polar("y", start=0)



png(paste0(direttoria,"_",hierarch_component,"_Fractional_abundance_barplot",".png"), width = 300, height = 300, units='mm', res = 300)
print(bp)
dev.off()
png(paste0(direttoria,"_",hierarch_component,"_Fractional_abundance_barplot_perc",".png"), width = 300, height = 300, units='mm', res = 300)
print(bp2)
dev.off()

png(paste0(direttoria,"_",hierarch_component, "_Fractional_abundance_piechart",".png"), width = 300, height = 300, units='mm', res = 300)
print(pie)
dev.off()
