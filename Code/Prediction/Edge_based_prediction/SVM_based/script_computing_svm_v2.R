library(data.table)
library(dplyr)
library(RankAggreg)
library(rlang)

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

if (grepl("14", direttoria)){
  metadata = read.delim("---/Data/w14w24_Data/w14_n188_metadata.txt")
} else if (grepl("24", direttoria)){
  metadata = read.delim("---/Data/w14w24_Data/w24_n160_metadata.txt")
  
  } else {
    abort("using a wrong directory, not 14 or 24")
}

mapping_nameIndividual_SPARCC <- read.csv("---/Data/mapping_nameIndividual_SPARCC.txt", sep="")



## MERGE OUTCOME

merged_outcome =  metadata 
merged_outcome$FC.nummer = gsub("-", ".", merged_outcome$FC.nummer)
colnames(Resulting_net_6m) %in% merged_outcome$FC.nummer

metadata[1:5,1:5]
Resulting_net_6m[1:5,1:5]


rownames(Resulting_net_6m) = Resulting_net_6m$V1
name_vect = Resulting_net_6m$V1



cc       <- strsplit(name_vect,'_')
part1    <- unlist(cc)[2*(1:length(name_vect))-1]
part2    <- unlist(cc)[2*(1:length(name_vect))  ]
head(part1) ; head(part2) ; head(name_vect)


if (sum(colnames(Resulting_net_6m) %in% merged_outcome$FC.nummer) != ncol(Resulting_net_6m) -1){
  abort("MISMATCH --> the metadata does not have all the individual ")
}

merged_outcome = merged_outcome[ merged_outcome$FC.nummer %in% colnames(Resulting_net_6m),]

dat_with_taxa = Resulting_net_6m %>% select (-c(merged_outcome$FC.nummer[is.na(merged_outcome %>% select(c(outcome_file))   )]) )
dat = dat_with_taxa %>% select (-c(V1))

# we have few data --> outcome T or F if there is the disease --> 1 disease 
outcome = colnames(dat) %in% merged_outcome$FC.nummer[merged_outcome %>% select(c(outcome_file)) == 1]

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


# 
set.seed(123)

# SVM  --------------------------------------------------------------------
rankfeature = "internal"
Rep = 10
MeanAUC = c()
sdAUC = c()

for (nFeat in seq(2,20,1)){

  featureSet = order(FeatureCorrelation, decreasing = TRUE)[1:nFeat]

  X = sample_data_num[,featureSet]
  X = scale(X)

  ## LOOCV
  N = nrow(sample_data_num)
  kernelParam = .1
  Result = list()
  Result[["AUC"]] = c()

  for (rep in 1:Rep){
    y_hat_cont = rep(NA, N)
    y_hat = rep(NA, N)

    for (i in 1:N){
      if (rankfeature == "internal"){

        FeatureCorrelation = apply(sample_data_num[-i,], 2, function(x){abs(cor(x,as.numeric(sample_y[-i])))})
        featureSet = which(FeatureCorrelation > .2)
        featureSet = order(FeatureCorrelation, decreasing = TRUE)[1:nFeat]
        X = sample_data_num[,featureSet]
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
  sdAUC = c(sdAUC, sd(Result[["AUC"]]))
  
}


dir.create(file.path(paste0("SVM", outcome_file)), showWarnings = FALSE)
setwd(paste0("SVM", outcome_file))


plot(MeanAUC)
max(MeanAUC)

FeatureCorrelation = apply(sample_data_num, 2, function(x){abs(cor(x,as.numeric(sample_y)))})
nFeat = which.max(MeanAUC)+1
featureSet = order(FeatureCorrelation, decreasing = TRUE)[1:nFeat]
topFeatures = data.frame(feature_top = featureSet, Taxas = colnames(sample_data_num[,featureSet] ))

write.table(topFeatures,  paste0("TopFeatures_internal.txt") )


## two value external or internal
write(MeanAUC, paste0("MeanAUC_int.txt"))
write(max(MeanAUC), paste0("MaxAUC_int.txt"))
write(sdAUC,paste0("SdAUC_int.txt"))
write(sdAUC[which.max(MeanAUC)],paste0("MaxSdAUC_int.txt"))


png(filename =  paste0("PlotAUC_int.png"),)
plot(MeanAUC)
dev.off()


png(filename =  paste0("histAUC_int.png"),)
hist(FeatureCorrelation, 20)
dev.off()
