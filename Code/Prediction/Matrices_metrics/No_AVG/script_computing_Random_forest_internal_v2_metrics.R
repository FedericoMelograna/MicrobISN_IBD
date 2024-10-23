
library(data.table)
library(dplyr)
library(RankAggreg)
library(rlang)
direttoria = args$direttoria
pred_dir = gsub("A_", "", direttoria)
outcome_file = args$outcome
feat_sel = args$feat_sel 
na_treatment = args$na_treatment 

final_direttoria = paste(direttoria,feat_sel,na_treatment, sep = "_")

metrics_matr <- readRDS(paste0(direttoria, "/metricmatrices/list_result_3d.rds"))
class(metrics_matr)
aa = metrics_matr[[1]]

matr_prediction = data.frame(t(metrics_matr[feat_sel][[1]]))

rownames(matr_prediction) = gsub("_",".", rownames(matr_prediction))
rownames(matr_prediction) = gsub("-",".", rownames(matr_prediction))



setwd(paste0(path_base, "Analysis_NOAVG"))
getwd()

if (grepl("14", direttoria)){
  metadata = read.delim("---/Data/w14w24_Data/w14_n188_metadata.txt")
  metadata$FC.nummer = gsub("-",".",metadata$FC.nummer)
  merged_outcome =  metadata[metadata$FC.nummer %in% rownames(matr_prediction), ]
  

} else if (grepl("24", direttoria)){
  metadata = read.delim("---/Data/w14w24_Data/w24_n160_metadata.txt")
  metadata$FC.nummer = gsub("-",".",metadata$FC.nummer)
  merged_outcome =  metadata[metadata$FC.nummer %in% rownames(matr_prediction), ]
  
  
  
  } else { 
    metadata <- read.delim("---/Data/w0_n335_metadata.txt")
    mapping_nameIndividual_SPARCC <- read.csv("---/Data/Follow_up_data/mapping_nameIndividual_SPARCC.txt", sep="")
    mapping_nameIndividual_SPARCC$NAME = gsub("_",".", mapping_nameIndividual_SPARCC$NAME)
    mapping_nameIndividual_SPARCC_filt = mapping_nameIndividual_SPARCC[mapping_nameIndividual_SPARCC$NAME %in% c( rownames(matr_prediction)),]
    merged_outcome =  merge(mapping_nameIndividual_SPARCC_filt,metadata, by.y = "FC.nummer", by.x = "Individual_ID"    )
    pos = match( rownames(matr_prediction), merged_outcome$NAME)
    
    rownames(matr_prediction) = gsub("-",".",merged_outcome$Individual_ID[pos])
    merged_outcome$FC.nummer = merged_outcome$Individual_ID
    metadata$FC.nummer = gsub("-",".",metadata$FC.nummer)
    merged_outcome =  metadata[metadata$FC.nummer %in% rownames(matr_prediction), ]
    #abort("using a wrong directory, not 14 or 24")
}


## Check

rownames(matr_prediction) %in% merged_outcome$FC.nummer
metadata[1:5,1:5]



if (sum(rownames(matr_prediction) %in% merged_outcome$FC.nummer) != nrow(matr_prediction)){
  abort("MISMATCH --> the metadata does not have all the individual ")
}


  
dat = matr_prediction[!is.na(merged_outcome %>% select(c(outcome_file))   ),]
outcome = rownames(dat) %in% merged_outcome$FC.nummer[merged_outcome %>% select(c(outcome_file)) == 1]

if (na_treatment == "elim"){
  dat = dat[,apply(dat, 2, function(x) sum(is.na(x)) ==0)]
  
} else if (na_treatment == "zero"){
  dat[is.na(dat)] = 0
}

dir.create(final_direttoria)
setwd(final_direttoria)

kend = apply(dat,2, function(x) cor.test(x, ifelse(outcome,1,0),  method="kendall")$p.value)
spear= apply(dat,2, function(x) cor.test(x, ifelse(outcome,1,0),  method="spearman")$p.value)

sign_both = data.frame(Kend = kend, Spear = spear)
write.table(sign_both, file = "Significance_correlation.txt", quote = F, row.names = T, col.names = T, sep = "\t")

# ELIMINATE data column with no variation
variance_row_after_cutting = round(apply(dat,2, var ), 5) 
dat = dat[,variance_row_after_cutting != 0]
dat_with_taxa = data.frame(cbind(rownames(dat), dat))

temp_res = cbind(dat,outcome)
f_outcome = data.frame(Y = temp_res[,"outcome"], Ind = rownames(temp_res))



# Random Forest -----------------------------------------------------------

dir.create(file.path(paste0("RandomForest_", outcome_file)), showWarnings = FALSE)
setwd(paste0("RandomForest_", outcome_file))

fwrite(temp_res, file = "Net_with_outcome.txt" , quote = FALSE, sep = "\t", row.names = T, col.names = T)

sample_data = as.data.frame(temp_res) %>% dplyr::select(-c(outcome))

sample_y = factor(f_outcome$Y)

# GLOBAL  -----------------------------------------------------------------

# TOP feature set ---------------------------------------------------------

FeatureCorrelation = apply(sample_data, 2, function(x){abs(cor(x,as.numeric(sample_y)))})

featureSet = which(FeatureCorrelation > .2)

reduced_samples = sample_data[,colnames(sample_data) %in% names(featureSet)]

nRareSamples = min(table(sample_y)[1] , table(sample_y)[2])

rf.strata = randomForest(x = sample_data , y = sample_y , importance = TRUE,
                         sampsize=c(nRareSamples,nRareSamples))
print(rf.strata)

# # Mean decrease gini &  MeanDecreaseAccuracy------------------------------------------------------
v_strat = varImpPlot(rf.strata,type=2)
v_strat_a =varImpPlot(rf.strata,type=1)
df = data.frame(cbind(v_strat,v_strat_a, "Name" = rownames(v_strat)))
v_strat2 = df[order(df[,1], decreasing = T),]


write.table(v_strat2, file = "Variable_importance_cutoff.txt",quote = FALSE, sep = " ", row.names = T)
# # PLOTS -------------------------------------------------------------------

dataset_auc = data.frame("strata" = NA, "cutoff" = NA, "threhsold" = NA)

png("VarImp_plot_global.png")
varImpPlot(rf.strata)
dev.off()
# 
# 
imp_score = importance(rf.strata)
scaled = data.frame(cbind(imp_score[,c(1,2)],scale(imp_score[,c(3,4)])))
scaled$tot = scaled$MeanDecreaseAccuracy + scaled$MeanDecreaseGini
sorted_scale = scaled[order(scaled$tot,decreasing = T),]
top_20 = sorted_scale[1:20 ,]


top_feature = sample_data[,colnames(sample_data) %in% rownames(top_20), drop = F]
rownames(top_feature) = rownames(sample_data)
write.table("RF_top_global_feature.txt",x = top_feature,quote = F,sep = "\t",row.names = T, col.names = T)
write.table("RF_top_global_feature_importance.txt",x = top_20,quote = F,sep = "\t",row.names = T, col.names = T)


# INTERNAL-KFOLDED --------------------------------------------------------
source(paste0(path_outcome,"/Code/Prediction/Matrices_metrics/No_AVG/script_computing_Random_forest_internal_kfolds.R" ) ) 



# SVM  --------------------------------------------------------------------
rankfeature = "internal"
# Rep = 10
Rep = 5 # FASTER
MeanAUC = c()
MeanSD = c()
start = Sys.time()
for (nFeat in seq(2,20,2)){

  featureSet = order(FeatureCorrelation, decreasing = TRUE)[1:nFeat]

  X = sample_data[,featureSet]
  X = scale(X)

  ## LOOCV
  N = nrow(sample_data)
  kernelParam = .1
  Result = list()
  Result[["AUC"]] = c()

  for (rep in 1:Rep){
    y_hat_cont = rep(NA, N)
    y_hat = rep(NA, N)

    for (i in 1:N){
      if (rankfeature == "internal"){

        FeatureCorrelation = apply(sample_data[-i,], 2, function(x){abs(cor(x,as.numeric(sample_y[-i])))})
        # featureSet = which(FeatureCorrelation > .2)
        featureSet = order(FeatureCorrelation, decreasing = TRUE)[1:nFeat]
        # print(featureSet)
        X = sample_data[,featureSet]
        X = scale(X)        ## DO we really want to scale for CLR?
      }
      X_test = X[i,]
      y_test = sample_y[i]
      X_train = X[-i,]
      y_train = sample_y[-i]

      # Down sampling the majority class
      n1 = table(y_train)[1]
      n2 = table(y_train)[2]
      majorityClass = ifelse(n1>n2,1,2)
      n_resample = abs(n1-n2)

      indeces_Cm = which(y_train == levels(y_train)[majorityClass])
      indeces_Cm = sample(indeces_Cm, n_resample, replace = FALSE)

      if (n_resample > 0 ){

      X_train_resampled = X_train[-indeces_Cm,]
      y_train_resampled = y_train[-indeces_Cm]

      } else {
        X_train_resampled = X_train
        y_train_resampled = y_train
      }
      model = svm(y = as.factor(y_train_resampled), x = data.frame(X_train_resampled),
                  kernel = "radial", gamma = kernelParam, scale = F)

      y_hat_cont[i] = predict(model, newdata = data.frame(rbind(X_test)))
      y_hat[i] = ifelse(y_hat_cont[i]>=1.5,2,1)

    }

    pred = prediction(y_hat, sample_y)
    auc = performance(pred, measure = "auc")

    Result[["AUC"]] = c(Result[["AUC"]], as.numeric(auc@y.values))
  }

  MeanAUC = c(MeanAUC, mean(Result[["AUC"]]))
  MeanSD = c(MeanSD, sd(Result[["AUC"]]))

}
stops = Sys.time()
stops - start
dir.create(file.path(paste0("SVM", outcome_file)), showWarnings = FALSE)
setwd(paste0("SVM", outcome_file))


# plot(MeanAUC)
# max(MeanAUC)

FeatureCorrelation = apply(sample_data, 2, function(x){abs(cor(x,as.numeric(sample_y)))})
nFeat = which.max(MeanAUC)+1
featureSet = order(FeatureCorrelation, decreasing = TRUE)[1:nFeat]
topFeatures = data.frame(feature_top = featureSet, Taxas = colnames(sample_data[,featureSet] ))

write.table(topFeatures,  paste0("TopFeatures_internal.txt") )


## two value external or internal
write(MeanAUC, paste0("MeanAUC_int.txt"),)
write(max(MeanAUC), paste0("MaxAUC_int.txt"))
write(MeanSD, paste0("MeanSD_int.txt"),)

png(filename =  paste0("PlotAUC_int.png"),)
plot(MeanAUC)
dev.off()


png(filename =  paste0("histAUC_int.png"),)
hist(FeatureCorrelation, 10)
dev.off()


dataset_auc = rbind(dataset_auc, "SVM" = 
                      c(rep(max(MeanAUC),3)))

rownames(dataset_auc)[nrow(dataset_auc)] = "SVM"

write.table("../FINAL_AUC_kfold.txt",x = dataset_auc,quote = F,sep = "\t",row.names = T, col.names = T)

