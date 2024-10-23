AUC_calculator_c = function(predics, test_data_y){
  pred = prediction(predics, test_data_y)
  auc = ROCR::performance(pred, "auc")#, sample_y_FOLD_test)
  
  return(as.numeric(auc@y.values))
}



rfcv_custom_strat = function (trainx, trainy, cv.fold = 5, scale = "log", step = 0.5, 
                        mtry = function(p) max(1, floor(sqrt(p))), recursive = FALSE,sampsize = sampsize, 
                        ...) 
{
  classRF <- is.factor(trainy)
  n <- nrow(trainx)
  p <- ncol(trainx)
  if (scale == "log") {
    k <- floor(log(p, base = 1/step))
    n.var <- round(p * step^(0:(k - 1)))
    same <- diff(n.var) == 0
    if (any(same)) 
      n.var <- n.var[-which(same)]
    if (!1 %in% n.var) 
      n.var <- c(n.var, 1)
  }
  else {
    n.var <- seq(from = p, to = 1, by = step)
  }
  k <- length(n.var)
  cv.pred <- vector(k, mode = "list")
  cv.AUC <- vector(k, mode = "list")
  for (i in 1:k) cv.pred[[i]] <- trainy
  if (classRF) {
    f <- trainy
  }   else {
    f <- factor(rep(1:5, length = length(trainy))[order(order(trainy))])
  }
  nlvl <- table(f)
  idx <- numeric(n)
  for (i in 1:length(nlvl)) {
    idx[which(f == levels(f)[i])] <- sample(rep(1:cv.fold, 
                                                length = nlvl[i]))
  }
  for (i in 1:cv.fold) {
    all.rf <- randomForest(trainx[idx != i, , drop = FALSE], 
                           trainy[idx != i], trainx[idx == i, , drop = FALSE], 
                           trainy[idx == i], mtry = mtry(p), importance = TRUE,sampsize =  rep(min(table(trainy[idx != i])),2), ...)#  )#, 
    # ...)
    cv.pred[[1]][idx == i] <- all.rf$test$predicted
    cv.AUC[[1]][idx == i] <- all.rf$test$votes[,2]
    
    impvar <- (1:p)[order(importance(all.rf, type = 1), 
                          decreasing = TRUE)]
    for (j in 2:k) {
      imp.idx <- impvar[1:n.var[j]]
      sub.rf <- randomForest(trainx[idx != i, imp.idx, 
                                    drop = FALSE], trainy[idx != i], trainx[idx == 
                                                                              i, imp.idx, drop = FALSE], trainy[idx == i], 
                             mtry = mtry(n.var[j]), importance = recursive, sampsize = rep(min(table(trainy[idx != i])),2), ...) # sampsize = sampsize ) 
      cv.pred[[j]][idx == i] <- sub.rf$test$predicted
      cv.AUC[[j]][idx == i] <- sub.rf$test$votes[,2]
      if (recursive) {
        impvar <- (1:length(imp.idx))[order(importance(sub.rf, 
                                                       type = 1), decreasing = TRUE)]
      }
      NULL
    }
    NULL
  }
  if (classRF) {
    error.cv <- sapply(cv.pred, function(x) mean(trainy != 
                                                   x))
    AUC_performance.cv = sapply(cv.AUC, function(x) AUC_calculator_c(x,trainy))
  }
  else {
    error.cv <- sapply(cv.pred, function(x) mean((trainy - 
                                                    x)^2))
    
  }
  names(error.cv) <- names(cv.pred) <- n.var
  names(AUC_performance.cv) <- names(cv.AUC) <- n.var
  
  
  list(n.var = n.var, error.cv = error.cv, predicted = cv.pred, prob.prediction = cv.AUC, AUC.performance = AUC_performance.cv, importance_var = impvar)
}


set.seed(123)
# CV - custom validated

dir.create(file.path(paste0("RFCV", outcome_file)), showWarnings = FALSE)
setwd(paste0("RFCV", outcome_file))

dir.create(file.path("RFCV_novel_strata"), showWarnings = FALSE)
setwd("RFCV_novel_strata")


rf.strat.cv <- replicate(10, rfcv_custom_strat(sample_data_num, sample_y, cv.fold=5, scale="log", step=0.7,
                                          mtry=function(p) max(1, floor(sqrt(p))), recursive=FALSE )#,
, simplify =  F)  

rf.strat.cv.importance_var <- sapply(rf.strat.cv, "[[", "importance_var")
rf.strat.cv.AUC_perf <- sapply(rf.strat.cv, "[[", "AUC.performance")


png("novel_STRAT_graph_AUC_performances.png")
par(mfrow = c(1,1))
matplot(rf.strat.cv[[1]]$n.var, cbind(rowMeans(rf.strat.cv.AUC_perf), rf.strat.cv.AUC_perf), type="l",
        lwd=c(2, rep(1, ncol(rf.strat.cv.AUC_perf))), col=1, lty=1, log="x",
        xlab="Number of variables", ylab="AUC perf")
dev.off()

rf.strat.cv.averaged_AUC = (apply(rf.strat.cv.AUC_perf,1, mean))
t_rf.strat.cv.averaged_AUC = t(rf.strat.cv.averaged_AUC)

dff.novel_STRAT.cv.averaged_AUC = data.frame("Features" = colnames(t_rf.strat.cv.averaged_AUC), "AUC" = rf.strat.cv.averaged_AUC )
write.table(dff.novel_STRAT.cv.averaged_AUC, "FEATURES_results_novel_STRAT.txt",quote = F,sep = "\t",row.names = F, col.names = T)
write.table(dff.novel_STRAT.cv.averaged_AUC[dff.novel_STRAT.cv.averaged_AUC$AUC == max(dff.novel_STRAT.cv.averaged_AUC$AUC),], "MAX_FEATURES_results_novel_STRAT.txt",quote = F,sep = "\t",row.names = F, col.names = T)
max_auc_novel_STRAT = dff.novel_STRAT.cv.averaged_AUC[dff.novel_STRAT.cv.averaged_AUC$AUC == max(dff.novel_STRAT.cv.averaged_AUC$AUC),"AUC"]


top25CE <- RankAggreg(t(rf.strat.cv.importance_var), 25, seed=100, rho=.01)

topFeatures_names = colnames(sample_data_num)[top25CE$top.list]
topFeature = sample_data_num[, colnames(sample_data_num)[top25CE$top.list]]

write.table("top25_global_feature_novel_STRAT.txt",x = topFeature,quote = F,sep = "\t",row.names = T, col.names = T)

top_s = as.integer(dff.novel_STRAT.cv.averaged_AUC[dff.novel_STRAT.cv.averaged_AUC$AUC == max(dff.novel_STRAT.cv.averaged_AUC$AUC),"Features"])

top_s = ifelse(top_s == 1, 2, top_s)
top_s = ifelse(top_s > 30, 30, top_s)

tokCE = tryCatch(
  expr = {

    tokCE <- RankAggreg(t(rf.strat.cv.importance_var), top_s , seed=100, rho=.01)
    return(tokCE)
  },
  error = function(e){ 
    print("error01")
    skipping = F
    tryCatch(
      expr = {

        tokCE <- RankAggreg(t(rf.strat.cv.importance_var), top_s , seed=100, rho=.08)
        return(tokCE)
      },
      error = function(e){ 
        print("error08")
        skipping = F
        tokCE <- RankAggreg(t(rf.strat.cv.importance_var), top_s , seed=100, rho=.12)
        return(tokCE)
        
      } 
    )
    
    
    
  } 
)


topFeatures_names = colnames(sample_data_num)[tokCE$top.list]
topFeature = sample_data_num[, colnames(sample_data_num)[tokCE$top.list]]



write.table("top_best_sel_global_feature_novel_STRAT.txt",x = topFeature,quote = F,sep = "\t",row.names = T, col.names = T)
dataset_auc[1,1] = max(dff.novel_STRAT.cv.averaged_AUC$AUC)

