

# Libraries ---------------------------------------------------------------




library(lionessR)
library(igraph)
library(reshape2)
# BiocManager::install("limma")
library(limma)
library(SummarizedExperiment)
library(data.table)
library(dplyr)

## references

# https://bioconductor.org/packages/release/bioc/vignettes/lionessR/inst/doc/lionessR.html




calculating_family_phylum = function(z_ll_merged3, ASV_taxonomic_annotation){
  if(dim(z_ll_merged3)[1] > 0){
  tab = z_ll_merged3
  tab2 = tab
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
  } else {tab2 = z_ll_merged3}
  return(tab2)
}


# Importing datasets ------------------------------------------------------
setwd("---")
getwd()


files <- list.files(path=".", pattern="group*", 
                    full.names=TRUE, recursive=FALSE)
curr_wd = getwd()

files
el
outcomes = rep(c("Clinical_outcome_combined", "Biomarker_outcome_combined","Endoscopic_outcome_combined" ), length(files) )
files_per_outcome = rep(files, each = 3 )
data_frame_loop  = data.frame(files = files_per_outcome, outcome = outcomes)

### LIMMA 0.25
for (i in (1:nrow(data_frame_loop)) ){
  el = data_frame_loop[i,]
  input_file = el$files
  outcome_file = el$outcome
  
  setwd(curr_wd)
  dir.create(file.path(curr_wd, input_file, paste0("/Results_LIMMA_CORRECT_025_",outcome_file) ) ) 
  setwd(file.path(curr_wd, input_file,paste0("/Results_LIMMA_CORRECT_025_",outcome_file)))
  
  
  NoLOOP_magma_net= data.frame(fread(file="../ISNs/Resulting_net_notNULL_NONODE_MAGMACONF.txt",
                                     sep = " " )) 
  
  # SPECIFIC metadata!
  if (grepl("14", el$files, fixed = TRUE)){
    metadata <- read.delim("C:/Users/fmelo/Desktop/Microbiome_drug_repurpose_project/Data/Follow_up_data/w14_n188_metadata.txt")
    metadata$FC.nummer = gsub("-",".",metadata$FC.nummer)
  } else if (grepl("24", el$files, fixed = TRUE)){
    metadata <- read.delim("C:/Users/fmelo/Desktop/Microbiome_drug_repurpose_project/Data/Follow_up_data/w24_n160_metadata.txt")
    metadata$FC.nummer = gsub("-",".",metadata$FC.nummer)
  } else {
    metadata <- read.delim("C:/Users/fmelo/Desktop/Microbiome_drug_repurpose_project/Data/w0_n335_metadata.txt")
  }

  
  metadata[1:5,1:5]
  NoLOOP_magma_net[1:5,1:5]
  
  ## FILTER MAPPING 
  

  (metadata$FC.nummer %in% colnames(NoLOOP_magma_net))
  ## MERGE OUTCOME
  merged_outcome =  metadata[metadata$FC.nummer %in% colnames(NoLOOP_magma_net), ]
  
  # Name change of the nodes ------------------------------------------------
  
  rownames(NoLOOP_magma_net) = NoLOOP_magma_net$V1
  name_vect = NoLOOP_magma_net$V1
  
  
  
  cc       <- strsplit(name_vect,'_')
  part1    <- unlist(cc)[2*(1:length(name_vect))-1]
  part2    <- unlist(cc)[2*(1:length(name_vect))  ]
  head(part1) ; head(part2) ; head(name_vect)
  
  
  
 
  cvar <- apply(NoLOOP_magma_net[,-1], 1, sd)
  
  dat = NoLOOP_magma_net %>% select (-c(merged_outcome$FC.nummer[is.na(merged_outcome %>% select(c(outcome_file))   )]) )
  dat = dat %>% select (-c(V1))
  is_responding = colnames(dat) %in% merged_outcome$FC.nummer[merged_outcome %>% select(c(outcome_file)) == 1]
  
  write(colnames(dat), "INDIVIDUAL_without_NA_outcome.txt")
  
  
  ## Average all the network of dummy = TRUE
  netyes <- apply(dat[,is_responding],1,mean)
  
  ## Average all the network of dummy = FALSE
  netno  <- apply(dat[,!is_responding ],1,mean)
  
  # calculate the difference
  netdiff <- netyes-netno
  
  
  
  # take the name of the selected edges -------------------------------------
  
  
  name_vect_500 = names(netdiff)
  cc       <- strsplit(name_vect_500,'_')
  part1    <- unlist(cc)[2*(1:length(name_vect_500))-1]
  part2    <- unlist(cc)[2*(1:length(name_vect_500))  ]
  head(part1) ; head(part2) ; head(name_vect_500)
  els <- data.frame(cbind(part1, part2, c(netdiff)))
  tosub <- els
  
  
  
  # Select only edges with an abs difference of 0.25 between 2 groups --------
  
  tosel <- row.names(tosub[which(abs(as.numeric(tosub[,3]))>0.25),])
  final_el = els[abs(as.numeric(els$V3)) > 0.25,]
  
  # ELIMINATE COLUMN (individual) where the given outcome is NA
  cut_complete_Inet = NoLOOP_magma_net[rownames(NoLOOP_magma_net) %in% rownames(final_el),]
  cut_complete_Inet = cut_complete_Inet %>% select (-c(merged_outcome$FC.nummer[is.na( merged_outcome %>% select(c(outcome_file)) )]) )
  cut_complete_Inet = cut_complete_Inet %>% select(-c(V1))
  
  
  if (dim(cut_complete_Inet)[1] == 0 ){
    ### NO NODE BEYOND THE THREHSOLD
    print("NO NODE BEYOND THE THREHSOLD")
    next
  }
  # Create the lm model and the fit -----------------------------------------
  
  
  
  ff = ifelse(is_responding, "yes", "no")  
  group <- factor(ff)
  design <- model.matrix(~0+group)
  cont.matrix <- makeContrasts(yesvsno = (groupyes - groupno), levels = design)  
  
  fit <- lmFit(cut_complete_Inet, design)
  # Fit linear model for each gene given a series of arrays
  
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2e <- eBayes(fit2)
  # Given a linear model fit from lmFit, compute moderated t-statistics, moderated F-statistic, 
  # and log-odds of differential expression by empirical Bayes moderation of the standard errors towards a global value.
  
  
  ## adjusting with fdr for multiple testing 
  toptable <- topTable(fit2e, number=nrow(cut_complete_Inet), adjust="fdr")
  
  
  if (nrow(toptable) > 1){
    toptable_edges <- t(matrix(unlist(c(strsplit(row.names(toptable), "_"))),2))
    
  } else{
    toptable_edges <- t(matrix(c(unlist(c(strsplit(rownames(cut_complete_Inet), "_"))))))
  }
  
  # Select top50 ------------------------------------------------------------
  
  # less than 50, top 50
  base = 50
  top = ifelse(nrow(toptable) > base, base, nrow(toptable))
  
  z <- cbind(toptable_edges[1:top,], toptable$logFC[1:top])
  if (nrow(toptable) > 1){
    z <- cbind(toptable_edges[1:top,], toptable[1:top,])
  } else {
    z <-  cbind(toptable_edges, toptable)
    
  }
  ## all
  z_all <- cbind(toptable_edges, toptable)
  
  
  topt_sign = toptable[toptable$P.Value < 0.05,]
  
  
  ## if no sign p value, just take minimum
  if (dim(topt_sign)[1] == 0){
    topt_sign = toptable[toptable$P.Value == min(toptable$P.Value),]
    z_sign <- cbind(toptable_edges[1,1],toptable_edges[1,2], topt_sign)
    
  } else if (dim(topt_sign)[1] == 1){
    topt_sign = toptable[toptable$P.Value == min(toptable$P.Value),]
    z_sign <- cbind(toptable_edges[1,1],toptable_edges[1,2], topt_sign)
  } else{
    z_sign <- cbind(toptable_edges[1:nrow(topt_sign),], topt_sign)
    
  }
  
  
  
  
  
  # Bonferroni  --------------------------------------------------------------------
  
  # In our context, Bonferroni ranking == FDR ranking
  toptable_b <- topTable(fit2e, number=nrow(cut_complete_Inet), adjust="bonferroni")
  
  if (nrow(toptable_b) > 1){
    toptable_edges_b <- t(matrix(unlist(c(strsplit(row.names(toptable_b), "_"))),2))
    
  } else{
    toptable_edges_b <- t(matrix(c(unlist(c(strsplit(rownames(cut_complete_Inet), "_"))))))
  }
  z_b <- cbind(toptable_edges_b[1:top,], toptable_b$logFC[1:top])
  sum(z_b != z)
  getwd()
  write.table(z, file = "MAGMA_CONF_all_edges_thre025.txt",quote = FALSE, sep = " ")
  write.table(z_all, file = "MAGMA_CONF_all_edges_thre025_top50.txt",quote = FALSE, sep = " ")
  write.table(z_sign, file = "MAGMA_CONF_all_edges_thre025_top50_SIGN.txt",quote = FALSE, sep = " ")
  
  
  
  # Visualization  ----------------------------------------------------------
  create_save_graph = function(edge_lista,name,layout, el = el)
  {
    z = cbind(edge_lista[,1], edge_lista[,2], edge_lista[,3])
    
    g <- graph.data.frame(z, directed=FALSE)
    V(g)$label.cex = 0.7
    V(g)$size = degree(g)*3
    E(g)$weight <- 3*as.numeric(edge_lista[,3]) 
    E(g)$color[E(g)$weight<0] <- "blue"
    E(g)$color[E(g)$weight>0] <- "red"
    E(g)$weight <- 1
    print(g)
    png(paste0(name,".png"), width = 465, height = 225, units='mm', res = 300)
    
    plot(g,
         layout = layout
    )
    dev.off()
    return(g)
  }
  g <- graph.data.frame(z_all, directed=FALSE)
  
  layout_f <- data.frame(layout.fruchterman.reingold(g))
  rownames(layout_f) = names(V(g))
  layout <-cbind(layout_f[,1],layout_f[,2])
  
  g_all = create_save_graph(z_all, "LIMMA_ALL_Nodes_unpaired_MAGMA_noCONF",layout)
  
  remove_node = function(g_full, g_low,layout_start)
  { ### we want to have the g_low with the same plotting setting as the g_full
    ###
    g_trunc = delete.edges(g_full, E(g_full)[!E(g_full) %in% E(g_low)] )
    remove = setdiff(names(V(g_trunc)), names(V(g_low)))
    retain = names(V(g_trunc))[!names(V(g_trunc)) %in% remove]
    s= (layout_start[rownames(layout_start) %in% retain, ])
    layout_2 = cbind(s[,1],s[,2])
    g_trunc2 <- induced_subgraph(g_trunc, V(g_trunc)[retain])
    
    return(list(g_trunc2, layout_2))
  }
  g_top = graph.data.frame(z, directed=FALSE)
  trunc_top = remove_node(g_all, g_top, layout_f)
  g_trunc = trunc_top[[1]]; layout_trunc = trunc_top[[2]]
  
  png(paste0("LIMMA_ALL_Nodes_unpaired_MAGMA_noCONF_top50BB",".png"), width = 465, height = 225, units='mm', res = 300)
  
  plot(g_trunc,
       layout = layout_trunc

  )
  dev.off()
  
  
  g_sig = graph.data.frame(z_sign, directed=FALSE)
  
  trunc_sig = remove_node(g_all, g_sig, layout_f)
  g_sig = trunc_sig[[1]]; layout_sig = trunc_sig[[2]]
  
  png(paste0("LIMMA_ALL_Nodes_unpaired_MAGMA_noCONF_SIGN",".png"), width = 465, height = 225, units='mm', res = 300)
  
  plot(g_sig,
       layout = layout_sig
  )
  dev.off()
  
  
  # THE list with TAXA ------------------------------------------------------
  Sequences_nodes = read.csv(file="../ISNs/Sequences_nodes.txt", sep = "" )

  ASV_taxonomic_annotation <- read.delim("---/ASV_taxonomic_annotation.txt")
  
  
  z_all_df = data.frame(z_all)
  z_all_df[1:5,1:3]
  
  colnames(z_all_df)[1:3] = c("Node1", "Node2", "Value")
  z_ll_merged1 = merge(z_all_df, Sequences_nodes, by.x = "Node1", by.y = "Node_list")
  colnames(z_ll_merged1)[length(colnames(z_ll_merged1))] = "Node1_Taxa"
  z_ll_merged2 = merge(z_ll_merged1, Sequences_nodes, by.x = "Node2", by.y = "Node_list")
  colnames(z_ll_merged2)[length(colnames(z_ll_merged2))] = "Node2_Taxa"
  z_ll_merged2$Node1_Taxa == z_ll_merged2$Node2_Taxa
  
  z_ll_merged3 = z_ll_merged2# %>% select(c(Node1, Node1_Taxa, Node2, Node2_Taxa, Value))
  z_ll_merged3$Node1_Taxa[2] == z_ll_merged3$Node1_Taxa[1]
  
  final_z_all = calculating_family_phylum(z_ll_merged3, ASV_taxonomic_annotation)
  write.table(final_z_all, file = "MAGMA_CONF_all_edges_withCORRECT_TAXA.txt", quote = F
              , sep = "\t", col.names = T, row.names = F)
  
  
  z_all_df = data.frame(z)
  z_all_df[1:5,1:3]
  
  colnames(z_all_df)[1:3] = c("Node1", "Node2", "Value")
  z_ll_merged1 = merge(z_all_df, Sequences_nodes, by.x = "Node1", by.y = "Node_list")
  colnames(z_ll_merged1)[length(colnames(z_ll_merged1))] = "Node1_Taxa"
  z_ll_merged2 = merge(z_ll_merged1, Sequences_nodes, by.x = "Node2", by.y = "Node_list")
  colnames(z_ll_merged2)[length(colnames(z_ll_merged2))] = "Node2_Taxa"
  z_ll_merged2$Node1_Taxa == z_ll_merged2$Node2_Taxa
  
  z_ll_merged3 = z_ll_merged2#  %>% select(c(Node1, Node1_Taxa, Node2, Node2_Taxa, Value))
  z_ll_merged3$Node1_Taxa[2] == z_ll_merged3$Node1_Taxa[1]
  final_z_all = calculating_family_phylum(z_ll_merged3, ASV_taxonomic_annotation)
  write.table(final_z_all, file = "MAGMA_CONF_all_edges_withCORRECT_TAXA_TOP50.txt", quote = F
              , sep = "\t", col.names = T, row.names = F)
  
  
  
  z_all_df = data.frame(z_sign)
  z_all_df[1:5,1:3]
  
  colnames(z_all_df)[1:3] = c("Node1", "Node2", "Value")
  z_ll_merged1 = merge(z_all_df, Sequences_nodes, by.x = "Node1", by.y = "Node_list")
  colnames(z_ll_merged1)[length(colnames(z_ll_merged1))] = "Node1_Taxa"
  z_ll_merged2 = merge(z_ll_merged1, Sequences_nodes, by.x = "Node2", by.y = "Node_list")
  colnames(z_ll_merged2)[length(colnames(z_ll_merged2))] = "Node2_Taxa"
  z_ll_merged2$Node1_Taxa == z_ll_merged2$Node2_Taxa
  
  z_ll_merged3 = z_ll_merged2#  %>% select(c(Node1, Node1_Taxa, Node2, Node2_Taxa, Value))
  z_ll_merged3$Node1_Taxa[2] == z_ll_merged3$Node1_Taxa[1]
  final_z_all = calculating_family_phylum(z_ll_merged3, ASV_taxonomic_annotation)
  write.table(final_z_all, file = "MAGMA_CONF_all_edges_withCORRECT_TAXA_pval_005.txt", quote = F
              , sep = "\t", col.names = T, row.names = F)
  
  
  
}


