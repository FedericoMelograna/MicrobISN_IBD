# Libraries ---------------------------------------------------------------
# table(merged$Diagnosis, merged$Gender)
# 
# L1  L2  L3  UC
# CD  95  47  84   5
# UC   0   0   0 104



Mode <- function(x) {
  # """
  # Computing the mode
  # """
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


# Read data ---------------------------------------------------------

# delim = "\t"
# header = F
data_importing <- function(path, delim = "\t", header = T){
  # """
  # Read the data in the given modality
  # """
  w0_n335_OTU_table =  fread(file = path, sep = delim)
}


# https://gitlab.com/arcgl/rmagma

# covariate = covariates_selected
# dati = otu_table_selected
importing_base_data_and_select_cropping = function(dati,covariate, crop, eliminate = "SampleList_Metadata_Prediction",
                                                   prevalance_par = 0.25, sequencing_depth_par = 500) # as in https://gitlab.com/arcgl/rmagma
{
  # """
  # The dati function must be in the format each row is an idnviidual and each column is a taxa
  # 
  # """
  # rownames(dati) = dati[,1]
  # dati = dati %>% select(- eliminate)
  
  prevalence       <- colMeans(dati>0)
  sequencing_depth <- rowSums(dati)
  
  icol <- prevalence >prevalance_par
  irow <- sequencing_depth >sequencing_depth_par
  sum(icol)
  
  otu_table <- dati[irow,icol]
  
  dim(otu_table) #[1] 335 122 ## individuals * taxas retained
  
  filt_cov = covariate[rownames(covariate) %in% rownames(otu_table),]
  
  
  return(list(otu_table, filt_cov))
}

# table(filtering_disease$Therapy_2)
# # TNF UST VDZ 
# # 140  70 125 
# table(filtering_disease$Diagnosis)
# # CD  UC 
# # 231 104 

# starting table
#    TNF UST VDZ
# CD  85  70  76
# UC  55   0  49

## after filtering with individuals with few taxas shared
#    TNF UST VDZ
# CD  79  63  75
# UC  53   0  49



# default parameters ------------------------------------------------------

# 
# diagnosi = c("CD","UC")
# terapia = c("TNF", "UST", "VDZ")
# diagnosi = c("CD")
# 
# table(full[[3]]$Diagnosis, full[[3]]$Therapy_2)
# # Unit test
# full = selecting_only_group(otu_table_eliminated, covariates_eliminated, select = filtering_disease, diagnosi =c("CD","UC") , terapia = c("TNF", "UST", "VDZ"))
# only_cd = selecting_only_group(otu_table_eliminated, covariates_eliminated, select = filtering_disease, diagnosi =c("CD") , terapia = c("TNF", "UST", "VDZ"))[[1]]
# only_cd_tnf = selecting_only_group(otu_table_eliminated, covariates_eliminated, select = filtering_disease, diagnosi =c("CD") , terapia = c("TNF"))[[1]]
# 
# dati = otu_table_eliminated
# covariate = covariates_eliminated
# select = filtering_disease


selecting_only_group = function(dati, covariate, select, diagnosi =c("CD","UC") , terapia = c("TNF", "UST", "VDZ"))
{
  # """
  # Select only the group to which we are interested into 
  # e.g. 
  # only_cd_tnf = selecting_only_group(otu_table_eliminated, 
  #             covariates_eliminated, select = filtering_disease, diagnosi =c("CD") , terapia = c("TNF"))[[1]]
  # 
  # """
  select = select[select$FC.nummer %in% rownames(dati),]
  filtering_disease_g_d = select[select$Therapy_2 %in% terapia & 
                                   select$Diagnosis %in% diagnosi  ,]
  dati_g = dati[rownames(dati) %in%  filtering_disease_g_d$FC.nummer,]
  covariate_g = covariate[rownames(covariate) %in% filtering_disease_g_d$FC.nummer,]
  
  
  final = "group"
  for (el in diagnosi){
    print(el)
    final = paste(final,el,sep = "_")
  }
  for(el in terapia){
    final = paste(final, el,sep = "_")
  }
  
  
  
  return ( list(dati_g, covariate_g, filtering_disease_g_d) )
  
}

build_igraph_from_magma = function(magma_Stool_AllCov)
{
  
  index = match(magma_Stool_AllCov$opt.lambda,magma_Stool_AllCov$lambda ) 
  binaryresult = magma_Stool_AllCov$path[[index]]
  
  gg = graph.adjacency(binaryresult)
  plot(gg)
}



build_LOO_net = function(dati, covariate, mapping_nameIndividual_SPARCC, weeks = "0")
{
  # """
  # Build LOO network, 
  # 
  # needed the dati, our otu table, 
  # covariate with same number of individual and 
  # mapping_nameIndividual_SPARCC = the mapping from individualID to real ID to write the correct ID for the individual 
  # Otu_table --> rows INDIVIDUAL, column TAXAS
  # # 
  # Covariates_selected 
  #   > covariates_selected[1:5,1:5]
  #         [,1] [,2] [,3] [,4] [,5]
  # FC-105  19.9    0 0.26    0    0
  # FC-108  23.8    0 1.68    0    0
  # FC-1147 31.9    0 0.20    0    0
  # FC-1263 61.7    1 4.26    1    0
  # FC-1287 26.4    0 0.13    1    0
  # --> mapping_nameIndividual_SPARCC
  # Individual_ID    NAME
  # 1          FC-105   Ind_1
  # 2          FC-108   Ind_2
  # 3         FC-1147   Ind_3
  # 4         FC-1247   Ind_4
  # 5         FC-1263   Ind_5
  # 6         FC-1287   Ind_6
  # """
  curr_wd = getwd()
  dir.create(file.path(curr_wd, "LooNet"))
  if (weeks == "0"){
  for (i in 1: nrow(dati)){
    ## does not work, interrupt on the 7th
    (Individual = mapping_nameIndividual_SPARCC[mapping_nameIndividual_SPARCC$Individual_ID == rownames(dati)[i], 2] )
    magma_Stool_AllCov <- magma(data = dati[-i,],X = covariate[-i,] )#, distrib = "ZIP")
    index = match(magma_Stool_AllCov$opt.lambda,magma_Stool_AllCov$lambda )
    binaryresult = magma_Stool_AllCov$path[[index]]
    colnames(binaryresult) = rownames(binaryresult) = colnames(dati)
    continuousresult = magma_Stool_AllCov$opt.icov
    colnames(continuousresult) = rownames(continuousresult) = colnames(dati)
    print(rownames(dati)[i])
    write.table(continuousresult, file = paste0("LooNet/MAGMA_continuous_Data",Individual, ".tsv" ), quote = FALSE, row.names = TRUE, col.names = TRUE,sep = "\t")
    write.table(binaryresult, file = paste0("LooNet/MAGMA_binary_Data",Individual, ".tsv" ), quote = FALSE, row.names = TRUE, col.names = TRUE,sep = "\t")
    
  }
  } else {
    for (i in 1: nrow(dati)){
      ## does not work, interrupt on the 7th
      (Individual = rownames(dati)[i])
      magma_Stool_AllCov <- magma(data = dati[-i,],X = covariate[-i,] )#, distrib = "ZIP")
      index = match(magma_Stool_AllCov$opt.lambda,magma_Stool_AllCov$lambda )
      binaryresult = magma_Stool_AllCov$path[[index]]
      colnames(binaryresult) = rownames(binaryresult) = colnames(dati)
      continuousresult = magma_Stool_AllCov$opt.icov
      colnames(continuousresult) = rownames(continuousresult) = colnames(dati)
      print(rownames(dati)[i])
      write.table(continuousresult, file = paste0("LooNet/MAGMA_continuous_Data",Individual, ".tsv" ), quote = FALSE, row.names = TRUE, col.names = TRUE,sep = "\t")
      write.table(binaryresult, file = paste0("LooNet/MAGMA_binary_Data",Individual, ".tsv" ), quote = FALSE, row.names = TRUE, col.names = TRUE,sep = "\t")
      
    }
  }
}

build_GLOB_net = function(dati, covariate)
{
  magma_Stool_AllCov <- magma(data = dati,X = covariate)
  index = match(magma_Stool_AllCov$opt.lambda,magma_Stool_AllCov$lambda ) 
  binaryresult = magma_Stool_AllCov$path[[index]]
  colnames(binaryresult) = rownames(binaryresult) = colnames(dati)
  continuousresult = magma_Stool_AllCov$opt.icov
  colnames(continuousresult) = rownames(continuousresult) = colnames(dati)
  # print(rownames(otu_table_eliminated)[i])
  write.table(continuousresult, file = paste0("MAGMA_continuous_Data_GLOBAL", ".tsv" ), quote = FALSE, row.names = TRUE, col.names = TRUE,sep = "\t")
  write.table(binaryresult, file = paste0("MAGMA_binary_Data_GLOBAL", ".tsv" ), quote = FALSE, row.names = TRUE, col.names = TRUE,sep = "\t")
  
  return(list(magma_Stool_AllCov, binaryresult, continuousresult))
  
  
}

graphical_printing = function(magma_data, layout, level = 2, dataset_match_taxa)
{
  matching_name = data.frame(Name =c("graphs_Covariates_MAGMA_2ndlevel_PHYLUM.png",
                                     "graphs_Covariates_MAGMA_3rdlevel_CLASS.png",
                                     "graphs_Covariates_MAGMA_4thlevel_ORDER.png",
                                     "graphs_Covariates_MAGMA_5thlevel_FAMILY.png",
                                     "graphs_Covariates_MAGMA_6thlevel_GENUS.png"),
                             Level = c(2,3,4,5,6) )
  name = matching_name$Name[matching_name$Level == level]
  png(name, width = 465, height = 225, units='mm', res = 300)
  plot(magma_data,layout = layout, V.color.factor = dataset_match_taxa[level,])
  dev.off()
  plot(magma_data,layout = layout, V.color.factor = dataset_match_taxa[level,])
  
}

matching_list_creator = function(matching_string= "MAGMA_continuous_Data")
{
  # """
  # We are looking if continuous. e.g. MAGMA_continuous_Data or MAGMA_binary_Data
  # """
  
  files <- list.files(path=".", pattern=paste0(matching_string,"*"), ## to have all the subsequent
                      full.names=TRUE, recursive=FALSE)
  names_file = gsub("./", "",files)
  names_file = gsub(".tsv","",names_file)
  names_file = gsub(matching_string,"", names_file)
  
  
  files_name_matching = cbind(files, names_file)
  return(list(files,files_name_matching,names_file))
}

preprocessing_global_net = function(global_net)
{
  
  Node_list = paste0("Node", seq(1,length(rownames(global_net))))
  Sequences_nodes = cbind(Node_list, rownames(global_net))
  
  try(if(sum(rownames(global_net) != Sequences_nodes[,2]) > 0 ) stop("NOT MATCHINGs"))
  
  colnames(global_net) = rownames(global_net) = Node_list
  
  global_net_vect <- c(as.matrix(global_net))
  lionessOutput <- matrix(NA, nrow(global_net) * ncol(global_net), length(files) + 
                            2)
  samples = names_file
  colnames(lionessOutput) <- c("reg", "tar", samples)
  lionessOutput[, 1] <- rep(row.names(global_net), ncol(global_net))
  lionessOutput[, 2] <- rep(colnames(global_net), each = nrow(global_net))
  lionessOutput <- as.data.frame(lionessOutput, stringsAsFactors = FALSE)
  
  
  return(list(global_net, global_net_vect,lionessOutput,Sequences_nodes))
}

ISN_computation = function(files, global_net, global_net_vect, lionessOutput, Sequences_nodes ,matching_string= "MAGMA_continuous_Data")
{
  
  Node_list = rownames(global_net)
  
  for (i in 1: length(files)){
    # for (i in 1: 10){
    recons_net = read.table(file= list.files(path=".", pattern=paste0(matching_string,"*"), 
                                             full.names=TRUE, recursive=FALSE)[i], sep = '\t', header = TRUE, row.names = 1)
    colnames(recons_net) = rownames(recons_net)
    dim(recons_net)
    try(if(sum(rownames(recons_net) != Sequences_nodes[,2]) > 0 ) stop("NOT MATCHINGs"))
    colnames(recons_net) = rownames(recons_net) = Node_list 
    corr_recons_net = c(as.matrix(recons_net)) #  scales a covariance matrix into the corresponding correlation matrix efficiently.
    head(recons_net[1:5,1:5])
    
    
    lionessOutput[, i + 2] <- length(files) * (global_net_vect - corr_recons_net) + corr_recons_net
    print(paste0(i, " = " , str(mean(abs( lionessOutput[, i + 2] )))))
  }
  
  
  edges <- paste(lionessOutput[, 1], lionessOutput[, 2], sep = "_")
  rownames(lionessOutput) = edges
  
  return(lionessOutput)
}



## THe parameter for the trial

# OTU_RESULT = continuous_result
# label_name_final = "_graphs_Covariates_MAGMA_Phylum_shrinking"
# otu_table = otu_table_eliminated
# Taxonomy_ASVs = ASV_final
# continuous_or_binary = "CONTINUOUS"
# month = ""
# import_layout = LAYOUT

plot_graphical_phylum = function(OTU_RESULT, otu_table, Taxonomy_ASVs,continuous_or_binary, label_name_final = "_graphs_Covariates_MAGMA_Phylum_shrinking", month = "6M", import_layout = F  )
{
  #   """
  #   OTU_RESULTS: or the binary or the continuous results for the OTU FINAL TABLE
  # Originated from 
  # magma_Stool_AllCov <- magma(data = otu_table9M)
  # index = match(magma_Stool_AllCov$opt.lambda,magma_Stool_AllCov$lambda ) 
  # binaryresult = magma_Stool_AllCov$path[[index]]
  # colnames(binaryresult) = rownames(binaryresult) = colnames(otu_table9M)
  # continuousresult = magma_Stool_AllCov$opt.icov
  # colnames(continuousresult) = rownames(continuousresult) = colnames(otu_table9M)
  # similar code
  #   
  #   otu_table: the otu table with the taxas abundance tused only to filter the taxonomy: format individuals on the rows, taxas on the columns
  #   Taxonomy_ASVs: the taxonomy of every individual in the form 
  #          Kingdom             Phylum                  Class                 Order                 Family               Genus           Species
  # 1   k__Bacteria  p__Proteobacteria c__Gammaproteobacteria  o__Enterobacteriales  f__Enterobacteriaceae      g__Escherichia           s__coli
  # 2   k__Bacteria  p__Actinobacteria      c__Actinobacteria  o__Bifidobacteriales  f__Bifidobacteriaceae  g__Bifidobacterium         s__longum
  # 3   k__Bacteria  p__Actinobacteria      c__Actinobacteria  o__Bifidobacteriales  f__Bifidobacteriaceae  g__Bifidobacterium               s__
  # 4   k__Bacteria  p__Actinobacteria      c__Actinobacteria  o__Bifidobacteriales  f__Bifidobacteriaceae  g__Bifidobacterium         s__longum
  # 
  # continuous_or_binary = a string "continuous" or "binary" depending on which data you put in OTU_RESULT that is used to save 
  
  
  # label_name_final = "" the name that we want to give to the graph 
  # month = "6M" : the month that we need for the graph
  # 
  
  # import_layout = F --> if we want to import layout e.g. "../6M/LAYOUT_PHYLUM_6M_graph.txt"
  # N.b we will use the LOO of the month 6 for the phylum graph
  # 
  # Return: nothing, 
  #  
  # plot the graphs
  #   """
  
  Taxonomy_ASVs[2,]
  ASV_final = Taxonomy_ASVs[Taxonomy_ASVs$ASV %in% colnames(otu_table),]
  df = (t(ASV_final))
  dim(df)
  colnames(df) = df[1,]
  df = df[-1,]
  
  Phylums = unique(ASV_final$Phylum)
  # REMOVAL PHYLUM IF UNKWNOW
  Phylums = Phylums[!is.na(Phylums)] 
  Taxas_name = colnames(OTU_RESULT)
  
  phy_sum = matrix(data = 0, nrow = length(Phylums), ncol = length(Phylums))
  colnames(phy_sum) = rownames(phy_sum) = Phylums
  phy_cont = phy_sum
  
  
  for (i in 1:nrow(OTU_RESULT)){
    for (j in 1:ncol(OTU_RESULT)){
      # IF ONE OF THE PHYLYM IS NA --> REMOVE IT 
      if (is.na(df["Phylum",Taxas_name[i]]) || is.na(df["Phylum",Taxas_name[j]]) ) {next;}
      # NA REMOVAL, REMOVE IF THE PHYLUM IS UNKNOWN 
      
      phy_sum[df["Phylum",Taxas_name[i]],df["Phylum",Taxas_name[j]] ] = phy_sum[df["Phylum",Taxas_name[i]],df["Phylum",Taxas_name[j]] ] + abs(OTU_RESULT[i,j])
      phy_cont[df["Phylum",Taxas_name[i]],df["Phylum",Taxas_name[j]] ] = phy_cont[df["Phylum",Taxas_name[i]],df["Phylum",Taxas_name[j]] ] +1 
    }
  }
  
  taxas_per_phylum = table(df["Phylum",])
  ordered_taxa_per_phylum = taxas_per_phylum[colnames(phy_cont)]
  phy_cont_less_diag = phy_cont  
  diag(phy_cont_less_diag) = diag(phy_cont) - ordered_taxa_per_phylum 
  
  
  phy_sum_norm = round((phy_sum / phy_cont_less_diag ) * 100, 1)
  
  
  # NA CREATED because 0 / 0 --> there is no signal, superimpose 0 ----------
  phy_sum_norm[is.na(phy_sum_norm) | is.infinite(phy_sum_norm)] = 0
  
  gg = graph_from_adjacency_matrix(
    phy_sum_norm,
    mode = "lower", ## lower traingolar 
    weighted = TRUE,
    diag = TRUE,
    add.colnames = NULL,
    add.rownames = NA
  )
  
  
  # FROM THE rMAGMA::plot_magma() function ----------------------------------
  
  
  # ONLY NON NAN ones
  V(gg)$color = as.numeric(factor(na.omit(unique(df[2,]))) )
  V.color.factor = (factor(na.omit(unique(df[2,]))))
  # Only NON NAN ones
  
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
  
  set.seed(123)
  ll <- igraph::layout_with_fr(gg)
  plot(gg,
       vertex.size = (diag(phy_sum_norm) + 5),
       edge.width = (E(gg)$weight / 2),
       V.color.factor = df[2,],
       layout = ll, 
       edge.label = E(gg)$weight)
  graphics::legend(x = +1.2, y = +1.3, levels(V.color.factor), 
                   pch = 21, col = V.color, pt.bg = V.color, pt.cex = 2, 
                   cex = 0.8, bty = "n", ncol = 1)
  
  
  # with edge weight
  
  
  if (is.character(import_layout))
  {
    LO = as.matrix(read.table( import_layout ))
  } else
  {
    set.seed(1234)
    LO = layout_nicely(gg)
  }
  
  png(paste0(month,label_name_final,"_",continuous_or_binary,"_","withEDGE_LABEL.png" ), width = 465, height = 225, units='mm', res = 300)
  plot(gg,
       vertex.size = (ordered_taxa_per_phylum^(4/5)+10),
       edge.width = (E(gg)$weight / 2),
       V.color.factor = df[2,],
       layout = LO, 
       edge.label = paste0("  ",round(E(gg)$weight,1)),
       #edge.color="black", 
       vertex.label = ordered_taxa_per_phylum,
       vertex.label.font=2, edge.label.font = 2, edge.label.cex = 1.2, 
       vertex.label.cex = 1 , edge.label.color = "black",vertex.label.color = "black")
  graphics::legend(x = +1.2, y = +1.3, levels(V.color.factor), 
                   pch = 21, col = V.color, pt.bg = V.color, pt.cex = 2, 
                   cex = 0.8, bty = "n", ncol = 1)
  dev.off()
  
  # without edge weight label  
  png(paste0(month,label_name_final,"_",continuous_or_binary,"_","NO_EDGE_LABEL.png" ), width = 465, height = 225, units='mm', res = 300)
  plot(gg,
       vertex.size = (ordered_taxa_per_phylum^(4/5)+10),
       edge.width = (E(gg)$weight / 2),
       V.color.factor = df[2,],
       layout = LO, 
       edge.label = NA ,# round(E(gg)$weight,1),
       # edge.color="black", 
       vertex.label = ordered_taxa_per_phylum,
       vertex.label.font=2, edge.label.font =4, edge.label.cex = 1.2, 
       vertex.label.cex = 1 , edge.label.color = "black", vertex.label.color = "black")
  graphics::legend(x = +1.2, y = +1.3, levels(V.color.factor), 
                   pch = 21, col = V.color, pt.bg = V.color, pt.cex = 2, 
                   cex = 0.8, bty = "n", ncol = 1)  
  dev.off()
  
  if (!(is.character(import_layout)))
  {
    write.table(LO,file = paste0("LAYOUT_PHYLUM_",month,"_graph.txt") ) 
  }
}

checkStrict <- function(f, silent=FALSE) {
  vars <- codetools::findGlobals(f)
  found <- !vapply(vars, exists, logical(1), envir=as.environment(2))
  if (!silent && any(found)) {
    warning("global variables used: ", paste(names(found)[found], collapse=', '))
    return(invisible(FALSE))
  }
  
  !any(found)
}
checkStrict(plot_graphical_phylum)
