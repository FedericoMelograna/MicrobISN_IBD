# Load necessary libraries
library(igraph)  # For graph creation and visualization
library(stringi) # For string manipulation
library(xlsx)    # For reading and writing Excel files

# Define the base directory and ASV annotation file
base_dir <- "C:/Users/fmelo/Desktop/Microbiome_drug_repurpose_project/"
ASV_annotation <- "ASV_taxonomic_annotation.txt"

# Switch to select the cohort for analysis
# Options: "CD_TNF", "UC_TNF", "CD_VDZ", "UC_VDZ", "CD_UST"
cohort <- "CD_TNF"

switch(cohort,
       "CD_TNF" = {
         after <- read.delim(paste0(base_dir, "group_CD_TNF_eliminate_24/Data/continuous_result_group_CD_TNF_eliminate_24.tsv"))
         before <- read.delim(paste0(base_dir, "Analysis/group_CD_TNF_eliminate/Data/continuous_result_group_CD_TNF_eliminate.tsv"))
         name_coh <- "CD_TNF"
       },
       "UC_TNF" = {
         after <- read.delim(paste0(base_dir, "group_UC_TNF_eliminate_14/Data/continuous_result_group_UC_TNF_eliminate_14.tsv"))
         before <- read.delim(paste0(base_dir, "Analysis/group_UC_TNF_eliminate/Data/continuous_result_group_UC_TNF_eliminate.tsv"))
         name_coh <- "UC_TNF"
       },
       "CD_VDZ" = {
         after <- read.delim(paste0(base_dir, "group_CD_VDZ_eliminate_24/Data/continuous_result_group_CD_VDZ_eliminate_24.tsv"))
         before <- read.delim(paste0(base_dir, "Analysis/group_CD_VDZ_eliminate/Data/continuous_result_group_CD_VDZ_eliminate.tsv"))
         name_coh <- "CD_VDZ"
       },
       "UC_VDZ" = {
         after <- read.delim(paste0(base_dir, "group_UC_VDZ_eliminate_14/Data/continuous_result_group_UC_VDZ_eliminate_14.tsv"))
         before <- read.delim(paste0(base_dir, "Analysis/group_UC_VDZ_eliminate/Data/continuous_result_group_UC_VDZ_eliminate.tsv"))
         name_coh <- "UC_VDZ"
       },
       "CD_UST" = {
         after <- read.delim(paste0(base_dir, "group_CD_UST_eliminate_24/Data/continuous_result_group_CD_UST_eliminate_24.tsv"))
         before <- read.delim(paste0(base_dir, "Analysis/group_CD_UST_eliminate/Data/continuous_result_group_CD_UST_eliminate.tsv"))
         name_coh <- "CD_UST"
       },
       stop("Invalid cohort selection"))

setwd(base_dir)

ASV_taxonomic_annotation <- read.delim(paste0(base_dir, "/Data/",ASV_annotation ) )

m1 = as.matrix(after)
m2 = as.matrix(before)
# Get the unique row and column names
all_row_names <- union(rownames(m1), rownames(m2))
all_col_names <- union(colnames(m1), colnames(m2))

# Create the result matrix with all zeros
result_matrix <- matrix(0, nrow = length(all_row_names), ncol = length(all_col_names),
                        dimnames = list(all_row_names, all_col_names))

# Fill in the values based on the conditions
result_matrix[rownames(m1), colnames(m1)] <- m1
result_matrix[rownames(m2), colnames(m2)] <- result_matrix[rownames(m2), colnames(m2)] - m2

result_matrix[result_matrix < 0.05 & result_matrix > - 0.05] = 0



# ASV notations -----------------------------------------------------------



ASV_taxonomic_annotation[2,]
ASV_final = ASV_taxonomic_annotation[ASV_taxonomic_annotation$ASV %in% colnames(result_matrix),]
df = (t(ASV_final))
dim(df)
colnames(df) = df[1,]
df = df[-1,]



# Create graph ------------------------------------------------------------



Taxas_name = colnames(result_matrix)


gg = graph_from_adjacency_matrix(
  result_matrix,
  mode = "lower", ## lower traingolar 
  weighted = TRUE,
  diag = FALSE,
  add.colnames = NULL,
  add.rownames = NA
)

## define gg structure; make it weight positive but with the sign saved
E(gg)$sign = ifelse(E(gg)$weight>0," ","-")
E(gg)$weight = abs(E(gg)$weight)

rank_current_taxa = 6
## Define the colors


# Genus level -------------------------------------------------------------



diff_ranks = df[rank_current_taxa,V(gg)$name]
diff_ranks[is.na(diff_ranks)] = "NoAv"
V(gg)$color = as.numeric(factor(((diff_ranks))) )
V.color.factor = levels(factor(diff_ranks))
V.color = NA
if (any(!is.na(V.color.factor))) {
  V.color.factor <- as.factor(V.color.factor)
  ll <- length(levels(V.color.factor))
  if (any(is.na(V.color))) {
    if (ll < 10) {
      V.color <- c("red", "green3", "cyan", "blue", 
                   "yellow", "magenta", "whitesmoke", "gray", 
                   "black")
    }
    else {
      V.color <- grDevices::rainbow(length(levels(V.color.factor)))
    }
    igraph::V(gg)$color <- V.color[V(gg)$color ]
  }
}

  set.seed(1234)
  LO = layout_nicely(gg)
  png(paste0("Comparison_",name_coh,".png") , width = 500, height = 300, units='mm', res = 300)
  plot(gg,
       vertex.size = ifelse(degree(gg)>=1, degree(gg)+2, 0.5), # 3, # (ordered_taxa_per_phylum^(7/10)+10),
       edge.width = (E(gg)$weight * 30)-0.5,
       V.color.factor = df[2,],
       layout = LO, 
       edge.label = NA ,# round(E(gg)$weight,1),
       edge.color = ifelse(E(gg)$sign ==" ", "blue", "tomato"), 
       vertex.label = NA ,#; ordered_taxa_per_phylum,
       vertex.label.font=2, edge.label.font =4, edge.label.cex = 1.2, 
       vertex.label.cex = 1 , edge.label.color = "black", vertex.label.color = "black")
  graphics::legend(x = +1.2, y = +1.,ifelse(stri_replace_all_regex(V.color.factor[order(V.color.factor)], "[a-z]__", "") == "", "<Unknown>",  stri_replace_all_regex(V.color.factor[order(V.color.factor)], "[a-z]__", "")), 
                   pch = 21, col = V.color, pt.bg = V.color, pt.cex = 0.5, 
                   cex = 0.8, bty = "n", ncol = 1)
dev.off()  






# Family level ------------------------------------------------------------
rank_current_taxa = 5

diff_ranks = df[rank_current_taxa,V(gg)$name]
diff_ranks[is.na(diff_ranks)] = "NoAv"
V(gg)$color = as.numeric(factor(((diff_ranks))) )
V.color.factor = levels(factor(diff_ranks))
V.color = NA
if (any(!is.na(V.color.factor))) {
  V.color.factor <- as.factor(V.color.factor)
  ll <- length(levels(V.color.factor))
  if (any(is.na(V.color))) {
    if (ll < 10) {
      V.color <- c("red", "green3", "cyan", "blue", 
                   "yellow", "magenta", "whitesmoke", "gray", 
                   "black")
    }
    else {
      V.color <- grDevices::rainbow(length(levels(V.color.factor)))
    }
    igraph::V(gg)$color <- V.color[V(gg)$color ]
  }
}

set.seed(1234)
LO = layout_nicely(gg)
png(paste0("Comparison_",name_coh,"_FAMILY.png") , width = 500, height = 300, units='mm', res = 300)
plot(gg,
     vertex.size = ifelse(degree(gg)>=1, degree(gg)+2, 0.5), # 3, # (ordered_taxa_per_phylum^(7/10)+10),
     edge.width = (E(gg)$weight * 30)-0.5,
     V.color.factor = df[2,],
     layout = LO, 
     edge.label = NA ,# round(E(gg)$weight,1),
     edge.color = ifelse(E(gg)$sign ==" ", "blue", "tomato"), 
     vertex.label = NA ,#; ordered_taxa_per_phylum,
     vertex.label.font=2, edge.label.font =4, edge.label.cex = 1.2, 
     vertex.label.cex = 1 , edge.label.color = "black", vertex.label.color = "black")
graphics::legend(x = +1.2, y = +1.1,ifelse(stri_replace_all_regex(V.color.factor[order(V.color.factor)], "[a-z]__", "") == "", "<Unknown>",  stri_replace_all_regex(V.color.factor[order(V.color.factor)], "[a-z]__", "")), 
                 pch = 21, col = V.color, pt.bg = V.color, pt.cex = 0.5, 
                 cex = 1, bty = "n", ncol = 1)
dev.off()  




# Genus level shrinked ----------------------------------------------------

# Identify nodes with at least one edge
nodes_with_edges <- which(degree(gg) > 0)

# Create a subgraph with nodes that have at least one edge
g_filtered <- induced_subgraph(gg, nodes_with_edges)



rank_current_taxa = 6
## Define the colours

diff_ranks = df[rank_current_taxa,V(g_filtered)$name]
diff_ranks[is.na(diff_ranks)] = "NoAv"
V(g_filtered)$color = as.numeric(factor(((diff_ranks))) )
V.color.factor = levels(factor(diff_ranks))
V.color = NA
if (any(!is.na(V.color.factor))) {
  V.color.factor <- as.factor(V.color.factor)
  ll <- length(levels(V.color.factor))
  if (any(is.na(V.color))) {
    if (ll < 10) {
      V.color <- c("red", "green3", "cyan", "blue", 
                   "yellow", "magenta", "whitesmoke", "gray", 
                   "black")
    }
    else {
      V.color <- grDevices::rainbow(length(levels(V.color.factor)))
    }
    igraph::V(g_filtered)$color <- V.color[V(g_filtered)$color ]
  }
}

set.seed(1234)
LO = layout_nicely(g_filtered)
png(paste0("Comparison_",name_coh,"_shrink.png") , width = 500, height = 300, units='mm', res = 300)
plot(g_filtered,
     vertex.size = ifelse(degree(g_filtered)>=1, degree(g_filtered)+2, 0.5), # 3, # (ordered_taxa_per_phylum^(7/10)+10),
     edge.width = (E(g_filtered)$weight * 30)-0.5,
     V.color.factor = df[2,],
     layout = LO, 
     edge.label = NA ,# round(E(g_filtered)$weight,1),
     edge.color = ifelse(E(g_filtered)$sign ==" ", "blue", "tomato"), 
     vertex.label = NA ,#; ordered_taxa_per_phylum,
     vertex.label.font=2, edge.label.font =4, edge.label.cex = 1.2, 
     vertex.label.cex = 1 , edge.label.color = "black", vertex.label.color = "black")
graphics::legend(x = +1.2, y = +1.,ifelse(stri_replace_all_regex(V.color.factor[order(V.color.factor)], "[a-z]__", "") == "", "<Unknown>",  stri_replace_all_regex(V.color.factor[order(V.color.factor)], "[a-z]__", "")), 
                 pch = 21, col = V.color, pt.bg = V.color, pt.cex = 0.5, 
                 cex = 0.8, bty = "n", ncol = 1)
dev.off()  

saveRDS(table(diff_ranks), file = paste0("Comparison_",name_coh,"_shrink.rds"))
write.xlsx(table(diff_ranks), file = paste0("Comparison_",name_coh,"_shrink.xlsx"))
# FAMILY shrinked ---------------------------------------------------------

rank_current_taxa = 5
## Define the colours


diff_ranks = df[rank_current_taxa,V(g_filtered)$name]
diff_ranks[is.na(diff_ranks)] = "NoAv"
V(g_filtered)$color = as.numeric(factor(((diff_ranks))) )
V.color.factor = levels(factor(diff_ranks))
V.color = NA
if (any(!is.na(V.color.factor))) {
  V.color.factor <- as.factor(V.color.factor)
  ll <- length(levels(V.color.factor))
  if (any(is.na(V.color))) {
    if (ll < 10) {
      V.color <- c("red", "green3", "cyan", "blue", 
                   "yellow", "magenta", "whitesmoke", "gray", 
                   "black")
    }
    else {
      V.color <- grDevices::rainbow(length(levels(V.color.factor)))
    }
    igraph::V(g_filtered)$color <- V.color[V(g_filtered)$color ]
  }
}

set.seed(1234)
LO = layout_nicely(g_filtered)
png(paste0("Comparison_",name_coh,"_FAMILY_shrink.png") , width = 500, height = 300, units='mm', res = 300)
plot(g_filtered,
     vertex.size = ifelse(degree(g_filtered)>=1, degree(g_filtered)+2, 0.5), # 3, # (ordered_taxa_per_phylum^(7/10)+10),
     edge.width = (E(g_filtered)$weight * 30)-0.5,
     V.color.factor = df[2,],
     layout = LO, 
     edge.label = NA ,# round(E(g_filtered)$weight,1),
     edge.color = ifelse(E(g_filtered)$sign ==" ", "blue", "tomato"), 
     vertex.label = NA ,#; ordered_taxa_per_phylum,
     vertex.label.font=2, edge.label.font =4, edge.label.cex = 1.2, 
     vertex.label.cex = 1 , edge.label.color = "black", vertex.label.color = "black")
graphics::legend(x = +1.2, y = +1.,ifelse(stri_replace_all_regex(V.color.factor[order(V.color.factor)], "[a-z]__", "") == "", "<Unknown>",  stri_replace_all_regex(V.color.factor[order(V.color.factor)], "[a-z]__", "")), 
                 pch = 21, col = V.color, pt.bg = V.color, pt.cex = 0.5, 
                 cex = 1, bty = "n", ncol = 1)
dev.off()  

saveRDS(table(diff_ranks), file = paste0("Comparison_",name_coh,"_FAMILY_shrink.rds"))
write.xlsx(table(diff_ranks), file = paste0("Comparison_",name_coh,"_FAMILY_shrink.xlsx"))





# Aggregate ---------------------------------------------------------------



small_result = result_matrix[(rowSums(abs(result_matrix)) !=0 ), (rowSums(abs(result_matrix)) !=0 )]
OTU_RESULT = small_result

level_taxa_order = "Family"
level_dataset = data.frame(name = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), rank = seq(1,7))
rank_current_taxa = level_dataset[level_taxa_order == level_dataset$name, "rank"]
  
  
Taxonomy_ASVs = ASV_taxonomic_annotation
ASV_final = Taxonomy_ASVs[Taxonomy_ASVs$ASV %in% colnames(OTU_RESULT),]
df = (t(ASV_final))

# # Use the taxa name as the colname
colnames(df) = df[1,]
df = df[-1,]
# ##
library(dplyr)
  
level_taxa = unlist(unique(ASV_final %>% dplyr::select (level_taxa_order)))
level_taxa = level_taxa[!is.na(level_taxa)]
# CREATE AGGREGATION
Taxas_name = colnames(OTU_RESULT)
# We already have taxa_names
phy_sum = matrix(data = 0, nrow = length(level_taxa), ncol = length(level_taxa))
colnames(phy_sum) = rownames(phy_sum) = level_taxa
phy_cont = phy_sum
  
  
diag(OTU_RESULT) = 0 ## We are not intereted to the diagonal 

for (i in 1:nrow(OTU_RESULT)){
  for (j in 1:ncol(OTU_RESULT)){
    if (is.na(df[rank_current_taxa,Taxas_name[i]]) || is.na(df[rank_current_taxa,Taxas_name[j]]) ) {next;}
    
    phy_sum[df[rank_current_taxa,Taxas_name[i]],df[rank_current_taxa,Taxas_name[j]] ] = phy_sum[df[rank_current_taxa,Taxas_name[i]],df[rank_current_taxa,Taxas_name[j]] ] + (OTU_RESULT[i,j])
    phy_cont[df[rank_current_taxa,Taxas_name[i]],df[rank_current_taxa,Taxas_name[j]] ] = phy_cont[df[rank_current_taxa,Taxas_name[i]],df[rank_current_taxa,Taxas_name[j]] ] +1 
  }
}

phy_small = phy_sum ; phy_small_cont = phy_cont
taxas_per_phylum = table(df[rank_current_taxa,])
ordered_taxa_per_phylum = taxas_per_phylum[colnames(phy_cont)]
phy_cont_less_diag = phy_cont  
diag(phy_cont_less_diag) = diag(phy_cont) - ordered_taxa_per_phylum 

  
phy_sum_norm = round((phy_sum / 2 ), 1)

# NA CREATED because 0 / 0 --> there is no signal, superimpose 0 ----------
phy_sum_norm[is.na(phy_sum_norm) | is.infinite(phy_sum_norm)] = 0

gg = graph_from_adjacency_matrix(
  (phy_sum_norm),
  mode = "lower", ## lower triangolar 
  weighted = TRUE,
  diag = TRUE,
  add.colnames = NULL,
  add.rownames = NA
)
  

## define gg structure; make it weight positive but with the sign saved
E(gg)$sign = ifelse(E(gg)$weight>0," ","-")
E(gg)$weight = abs(E(gg)$weight)

  
## Define the colors
uniq_ones = df[rank_current_taxa,] %>% unique()
V(gg)$color = uniq_ones[complete.cases(uniq_ones)] %>%  factor() %>% as.numeric() 

V.color.factor = uniq_ones[complete.cases(uniq_ones)] %>%  factor() 
V.color = NA
if (any(!is.na(V.color.factor))) {
  V.color.factor <- as.factor(V.color.factor)
  ll <- length(levels(V.color.factor))
  if (any(is.na(V.color))) {
    if (ll < 10) {
      V.color <- c("red", "green3", "cyan", "blue", 
                   "yellow", "magenta", "whitesmoke", "gray", 
                   "black")
    }
    else {
      V.color <- grDevices::rainbow(length(levels(V.color.factor)))
    }
    igraph::V(gg)$color <- V.color[V.color.factor]
  }
}


  set.seed(1234)
  LO = layout_nicely(gg)


# PLOTTING ----------------------------------------------------------------
  


png(paste0("Grouped_",name_coh,"_Family_","withEDGE_LABEL.png" ), width = 600, height = 300, units='mm', res = 300)
# if (rank_current_taxa >= 4){ ## IF YOU HAVE MANY different categories
  plot(gg,
       vertex.size = (ordered_taxa_per_phylum^(7/10)+10),
       edge.width = (E(gg)$weight / 2)+1,
       # V.color.factor = df[2,],
       layout = LO, 
       edge.label = paste0(E(gg)$sign,round(E(gg)$weight,1)),
       edge.color = ifelse(E(gg)$sign ==" ", "blue", "tomato"), 
       vertex.label = sapply(names(ordered_taxa_per_phylum), function(x) substr(x, 1, 4)), # names(ordered_taxa_per_phylum),
       vertex.label.font=2, edge.label.font = 2, edge.label.cex = 1.2, 
       vertex.label.cex = 1 , edge.label.color = "black",vertex.label.color = "black")
  graphics::legend(x = +1.2, y = +1,ifelse(stri_replace_all_regex(V.color.factor[order(V.color.factor)], "[a-z]__", "") == "", "<Unknown>",  stri_replace_all_regex(V.color.factor[order(V.color.factor)], "[a-z]__", "")), 
                   pch = 21, col = V.color, pt.bg = V.color, pt.cex = 2.5, 
                   cex = 1.3, bty = "n", ncol = 1)

dev.off()


# without edge weight label  
png(paste0("Grouped_",name_coh,"_Family_","withoutEDGE_LABEL.png" ), width = 600, height = 300, units='mm', res = 300)
  plot(gg,
       vertex.size = (ordered_taxa_per_phylum^(7/10)+10),
       edge.width = (E(gg)$weight / 2)+1,
       V.color.factor = df[2,],
       layout = LO, 
       edge.label = NA ,# round(E(gg)$weight,1),
       edge.color = ifelse(E(gg)$sign ==" ", "blue", "tomato"), 
       vertex.label = ordered_taxa_per_phylum,
       vertex.label.font=2, edge.label.font =4, edge.label.cex = 1.2, 
       vertex.label.cex = 1 , edge.label.color = "black", vertex.label.color = "black")
  graphics::legend(x = +1.2, y = +1.,ifelse(stri_replace_all_regex(V.color.factor[order(V.color.factor)], "[a-z]__", "") == "", "<Unknown>",  stri_replace_all_regex(V.color.factor[order(V.color.factor)], "[a-z]__", "")), 
                   pch = 21, col = V.color, pt.bg = V.color, pt.cex = 2.5, 
                   cex = 1.3, bty = "n", ncol = 1)

dev.off()

