
library(dplyr)
library(RankAggreg)
direttoria = args$direttoria
outcome_file = args$outcome

data(OSdata)
rowData <- DataFrame(row.names = rownames(exp), gene = rownames(exp))
colData <- DataFrame(row.names = targets$sample, sample = as.character(targets$sample), mets = targets$mets)

se <- SummarizedExperiment(assays = list(counts = as.matrix(exp)), 
                           colData = colData, rowData = rowData)
setwd(path_base)
getwd()
Resulting_net_6m= data.frame(fread(file=paste0(direttoria,"/ISNs/Resulting_net_notNULL_NONODE_MAGMACONF.txt"),sep = " "  )) # keeps the rownames

setwd(direttoria)

metadata <- read.delim("---/Data/w0_n335_metadata.txt")
mapping_nameIndividual_SPARCC <- read.csv("---/Data/mapping_nameIndividual_SPARCC.txt", sep="")

colnames(Resulting_net_6m) %in% mapping_nameIndividual_SPARCC$NAME

mapping_nameIndividual_SPARCC_filt = mapping_nameIndividual_SPARCC[mapping_nameIndividual_SPARCC$NAME %in% c( colnames(Resulting_net_6m)),]

## MERGE OUTCOME

merged_outcome =  merge(mapping_nameIndividual_SPARCC_filt,metadata, by.y = "FC.nummer", by.x = "Individual_ID"    )

metadata[1:5,1:5]
Resulting_net_6m[1:5,1:5]


rownames(Resulting_net_6m) = Resulting_net_6m$V1
name_vect = Resulting_net_6m$V1



cc       <- strsplit(name_vect,'_')
part1    <- unlist(cc)[2*(1:length(name_vect))-1]
part2    <- unlist(cc)[2*(1:length(name_vect))  ]
head(part1) ; head(part2) ; head(name_vect)


dat_with_taxa = Resulting_net_6m %>% select (-c(merged_outcome$NAME[is.na(merged_outcome %>% select(c(outcome_file))   )]) )
dat = dat_with_taxa %>% select (-c(V1))

outcome = colnames(dat) %in% merged_outcome$NAME[merged_outcome %>% select(c(outcome_file)) == 1]



# ELIMINATE data column with no variation
variance_row_after_cutting = round(apply(dat,1, var ), 4) 
dat = dat[variance_row_after_cutting != 0,]
dat_with_taxa = dat_with_taxa[variance_row_after_cutting != 0,]

temp_res = rbind(dat,outcome)
f_outcome = data.frame(t ( rbind(temp_res[nrow(temp_res),],colnames(temp_res)) ))
colnames(f_outcome) = c("Y", "Ind")

dat_with_taxa[1:5,1:5]
colnames(dat_with_taxa)[1] = "Taxas"

colnames(dat_with_taxa)



# Random Forest -----------------------------------------------------------
# https://stackoverflow.com/questions/19760169/how-to-perform-random-forest-cross-validation-in-r
# https://cran.r-project.org/web/packages/AUCRF/AUCRF.pdf
# https://stackoverflow.com/questions/31637259/random-forest-crossvalidation-in-r
# https://stats.stackexchange.com/questions/168415/random-forest-in-r-using-unbalanced-data/168424



dir.create(file.path(paste0("RandomForest_", outcome_file)), showWarnings = FALSE)
setwd(paste0("RandomForest_", outcome_file))

trans_ord = data.frame(t(dat_with_taxa))
dim(trans_ord)
colnames(trans_ord) = trans_ord[1,]
trans_ord_right = trans_ord[-1,]
try(rownames(trans_ord_right)  == f_outcome$Ind)
if( sum(rownames(trans_ord_right)  == f_outcome$Ind)!= nrow(trans_ord_right) ) stop("NON matching column")

net_with_out = cbind(trans_ord_right, f_outcome)

sample_data = net_with_out %>% dplyr::select(-c(Y,Ind))

sample_y = factor(net_with_out$Y)

sample_data_num = data.frame(apply(sample_data, 2, function(x){ as.numeric(x)}))
rownames(sample_data_num) = rownames(sample_data)


set.seed(123)






dataset_auc = data.frame("NOVEL_stratified_RFCV_AUC" = NA)




# RFCV
source(paste0(curr_wd,"/../Code/Prediction/Edge_based_onlyRFCV_strat/script_computing_Random_forest_rfcv_with_stratification.R" ) )
write.table("../../FINAL_AUC_NOVEL_stratified_RFCV_AUC.txt",x = dataset_auc,quote = F,sep = "\t",row.names = T, col.names = T)


# Here is a snippet from Breiman's official documentation:
# 
# In random forests, there is no need for cross-validation or a separate test set to get an unbiased estimate of the test set error.
# It is estimated internally, during the run, as follows:
# 
# Each tree is constructed using a different bootstrap sample from the original data. 
# About one-third of the cases are left out of the bootstrap sample and not used in the construction of the kth tree.
# 
# Put each case left out in the construction of the kth tree down the kth tree to get a classification.
# In this way, a test set classification is obtained for each case in about one-third of the trees. 
# At the end of the run, take j to be the class that got most of the votes every time case n was oob. The proportion of times that j is not equal to the true class of n averaged over all cases is the oob error estimate. This has proven to be unbiased in many tests.

