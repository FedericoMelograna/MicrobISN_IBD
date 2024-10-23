AUC_calculator_c = function(predics, test_data_y){
  pred = prediction(predics, test_data_y)
  auc = ROCR::performance(pred, "auc")#, sample_y_FOLD_test)
  
  return(as.numeric(auc@y.values))
}


# TOP feature set ---------------------------------------------------------

library(caret)
library(plyr)
folds <- createFolds(factor(sample_y), k = 10, list = FALSE)

sample_y[folds]
prediction_y = c()
factor(prediction_y, levels = levels(sample_y))

y_pred_global_strata = c()
y_pred_global_cutoff = c()
y_pred_global_classwt = c()
y_pred_global_normal = c()

y_pred_global_strata_onl = c()
y_pred_global_cutoff_onl = c()
y_pred_global_classwt_onl = c()
y_pred_global_normal_onl = c()

y_pred_global_strata_20 = c()
y_pred_global_cutoff_20 = c()
y_pred_global_classwt_20 = c()
y_pred_global_normal_20 = c()


y_pred_global_strata_glm = c()
y_pred_global_cutoff_glm = c()
y_pred_global_classwt_glm = c()
y_pred_global_normal_glm = c()

start = Sys.time()
y_pred_glm = c()
for(fold_i in 1:10){
  #Segment your data by fold using the which() function 
  testIndexes <- which(folds==fold_i,arr.ind=TRUE)
  sample_data_FOLD_test <- sample_data[testIndexes, ]
  sample_data_FOLD <- sample_data[-testIndexes, ]
  #Use the test and train data partitions however you desire...
  sample_y_FOLD_test <- sample_y[testIndexes]
  sample_y_FOLD <- sample_y[-testIndexes]
  prediction_y = c(prediction_y, as.character(sample_y_FOLD_test))

  FeatureCorrelation = apply(sample_data_FOLD, 2, function(x){abs(cor(x,as.numeric(sample_y_FOLD)))})
  featureSet = which(FeatureCorrelation > .2)

  
  set.seed(123)
  
  # rare.class.prevalence = 0.1
  rf.normal = randomForest(x = sample_data_FOLD , y = sample_y_FOLD ,
                           importance = TRUE,)

  
  
  set.seed(123)
  
  rare.class.prevalence = as.numeric(table(sample_y_FOLD)[1]/(table(sample_y_FOLD)[2]+table(sample_y_FOLD)[1]))
  if(table(sample_y_FOLD)[1] > table(sample_y_FOLD)[2]){
    cutf = c(1-rare.class.prevalence, rare.class.prevalence)
    
  } else {
    cutf = c(rare.class.prevalence,1-rare.class.prevalence)
  }
  
  rf.cutoff = randomForest(x = sample_data_FOLD , y = sample_y_FOLD ,
                           importance = TRUE,cutoff=c(rare.class.prevalence,1-rare.class.prevalence))
  
  summary(rf.cutoff)
  
  nRareSamples = min(table(sample_y_FOLD)[1] , table(sample_y_FOLD)[2])
  
  rf.strata = randomForest(x = sample_data_FOLD , y = sample_y_FOLD , importance = TRUE,
                           sampsize=c(nRareSamples,nRareSamples))
  
  auc(roc(rf.strata$votes[,2], sample_y_FOLD)) 
  auc(roc(rf.cutoff$votes[,2], sample_y_FOLD))
  auc(roc(rf.normal$votes[,2], sample_y_FOLD)) 
  
  pred_strata = predict(rf.strata, newdata = sample_data_FOLD_test, type = "prob")
  pred_cutoff = predict(rf.cutoff, newdata = sample_data_FOLD_test, type = "prob")
  pred_normal = predict(rf.normal, newdata = sample_data_FOLD_test, type = "prob")
  
  
  y_pred_global_strata = c(y_pred_global_strata, pred_strata[,2]) 
  y_pred_global_cutoff = c(y_pred_global_cutoff, pred_cutoff[,2])
  y_pred_global_normal = c(y_pred_global_normal, pred_normal[,2])
  
  
  
  imp_score = importance(rf.strata)
  scaled = data.frame(cbind(imp_score[,c(1,2)],scale(imp_score[,c(3,4)])))
  scaled$tot = scaled$MeanDecreaseAccuracy + scaled$MeanDecreaseGini
  sorted_scale = scaled[order(scaled$tot,decreasing = T),]
  
  
  top_20 = sorted_scale[1:20,]
  
  top_feature = sample_data_FOLD[,colnames(sample_data_FOLD) %in% rownames(top_20), drop = F]
  rownames(top_feature) = rownames(sample_data_FOLD)
  

  set.seed(123)
  
  rf.normal_20 = randomForest(x = top_feature , y = sample_y_FOLD ,
                              importance = TRUE,)

  set.seed(123)
  
  rare.class.prevalence = as.numeric(table(sample_y_FOLD)[1]/(table(sample_y_FOLD)[2]+table(sample_y_FOLD)[1]))
  
  
  if(table(sample_y_FOLD)[1] > table(sample_y_FOLD)[2]){
    cutf = c(1-rare.class.prevalence, rare.class.prevalence)
    
  } else {
    cutf = c(rare.class.prevalence,1-rare.class.prevalence)
  }
  
  rf.cutoff_20 = randomForest(x = top_feature , y = sample_y_FOLD ,
                              importance = TRUE,cutoff=cutf)

  nRareSamples = min(table(sample_y_FOLD)[1] , table(sample_y_FOLD)[2])
  

  rf.strata_20 = randomForest(x = top_feature , y = sample_y_FOLD , importance = TRUE,
                              sampsize=c(nRareSamples,nRareSamples))

  set.seed(123)
  
  pred_strata_20 = predict(rf.strata_20, newdata = sample_data_FOLD_test, type = "prob")
  pred_cutoff_20 = predict(rf.cutoff_20, newdata = sample_data_FOLD_test, type = "prob")
  pred_normal_20 = predict(rf.normal_20, newdata = sample_data_FOLD_test, type = "prob")
  
  
  y_pred_global_strata_20 = c(y_pred_global_strata_20, pred_strata_20[,2]) 
  y_pred_global_cutoff_20 = c(y_pred_global_cutoff_20, pred_cutoff_20[,2])
  y_pred_global_normal_20 = c(y_pred_global_normal_20, pred_normal_20[,2])
  
  rf.strata_20
  class(rf.strata_20)
  
  

}  
stop = Sys.time()
stop-start
prediction_y = factor(prediction_y, levels = levels(sample_y))

AUC_calculator_c(y_pred_global_strata, prediction_y)



dataset_auc = rbind(dataset_auc, "GLOBAL_10fold" = 
                      c(AUC_calculator_c(y_pred_global_strata, prediction_y),
                        AUC_calculator_c(y_pred_global_cutoff, prediction_y),
                        # AUC_calculator_c(y_pred_global_classwt, prediction_y),
                        AUC_calculator_c(y_pred_global_normal, prediction_y) ) )#,main="black strata, red normal, green classwt")

rownames(dataset_auc)[nrow(dataset_auc)] = "GLOBAL_RF_10fold"


dataset_auc = rbind(dataset_auc, "top20VarImp_10fold" = 
                      c(AUC_calculator_c(y_pred_global_strata_20, prediction_y),
                        AUC_calculator_c(y_pred_global_cutoff_20, prediction_y),
                        # AUC_calculator_c(y_pred_global_classwt_20, prediction_y),
                        AUC_calculator_c(y_pred_global_normal_20, prediction_y) ))

rownames(dataset_auc)[nrow(dataset_auc)] = "topRF_RF_10fold"



write.table("FINAL_AUC_kfold.txt",x = dataset_auc,quote = F,sep = "\t",row.names = T, col.names = T)
  
  
