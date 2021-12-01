##---------------------------------------------------------
## Control vs case analysis
## gsea code adapted from 
## http://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
##---------------------------------------------------------

## libraries 
library(fgsea)
library(data.table)
library(ggplot2)
library(caret)

load("ctrl_vs_case/ctrl_data")

##---------------------------------------------------------
## Split into train/test
##---------------------------------------------------------
# target variable
target_var <- "CtrlVsCase_Classifier"

# case = 1, ctrl = 0
case_rows <- which(df_case[, target_var] == 1)
ctrl_rows <- which(df_case[, target_var] == 0)

# stratified test/train sampling from both classes
set.seed(99)
percentageTest <- 30

nTestCase <- floor(percentageTest * length(case_rows) / 100)
case_test <- sort(sample(case_rows, nTestCase, replace=FALSE))

nTestCtrl <- floor(percentageTest * length(ctrl_rows) / 100)
ctrl_test <- sort(sample(ctrl_rows, nTestCtrl, replace=FALSE))

test_rows <- sort(c(case_test, ctrl_test))

## discard columns with near zero variance 
nzv <- nearZeroVar(df_case)
df_case <- df_case[, -nzv]

## standardize (z-transform) the columns 
df_case_scaled <- preProcess(df_case[, c(-1, -2)], 
                               method=c("center", "scale"))

df_test <- predict(df_case_scaled, df_case[test_rows, ])
df_train <- predict(df_case_scaled, df_case[-test_rows, ])

##---------------------------------------------------------
## Running GSEA
##---------------------------------------------------------

## Note: may want to upsample the controls (0's in phenotype vector)
## Using SMOTE algorithm: https://www.rdocumentation.org/packages/performanceEstimation/versions/1.1.0/topics/smote
library(performanceEstimation)

df_case_no_ID <- df_train[, 2:ncol(df_train)]
df_case_no_ID[, target_var] <- as.factor(df_case_no_ID[, target_var])
df_upsample <- smote(CtrlVsCase_Classifier~., df_case_no_ID, 
                     perc.over=6, perc.under=1, k=5)
df_upsample[,target_var] <- as.integer(df_upsample[,target_var])-1

## Split dataset into [Feature, vector]
X <- df_upsample[, 2:ncol(df_upsample)]
y <- df_upsample[, 1]

## function to compute correlation between gene expression and phenotype
gene_pheno_corr <- function(x, y) {
  if (sd(x) == 0) {
    return(0)
  }
  else {
    return(cor(x, y, method="pearson"))
  }
}

## Compute correlations 
ranks <- sapply(X, function(x) gene_pheno_corr(x, y))

## Compute t statistics (not using for now)
#tranks <- sapply(X, function(x) t.test(x ~ y)$statistic)

## Read in gene set
gmt.file <- "genesets.gmt"
pathways <- gmtPathways(gmt.file)

set.seed(99)
fgseaRes <- fgsea(pathways=pathways, stats=ranks)
# plotEnrichment(pathways[["GOBP_NEURON_INTRINSIC_APOPTOTIC
#       _SIGNALING_PATHWAY_IN_RESPONSE_TO
#       _OXIDATIVE_STRESS"]], 
#                ranks) + labs(title="ALS")
# plotEnrichment(pathways[["GOBP_MOTOR_NEURON_APOPTOTIC_PROCESS"]], 
#                ranks) + labs(title="ALS")

#Visualize the pathways and their performance (based on GSEA example code)
nPathways <- 3
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=nPathways), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=nPathways), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways[topPathways], ranks, fgseaRes, 
              gseaParam=0.5)

#Grab a new dataframe containing only the GSEA results for the top performing pathways
topPathwaysRes <- fgseaRes[fgseaRes$pathway %in% topPathways, ]
#print(topPathwaysRes[order(fgseaRes$pval),]$leadingEdge[[17]])

#Grab the leading edge of each of the top pathway
features <- topPathwaysRes[order(topPathwaysRes$pval),]$leadingEdge[[1]]

#Alternatively, grab the *best* nFeatureSelect features from each top pathway
#This is similar to that paper we did in class
top_features_temp <- topPathwaysRes$leadingEdge
nFeatureSelect <- 10
i = 1
for (f in top_features_temp) {
  top_features_temp[[i]] <- f[1:min(length(f), nFeatureSelect)]
  i <- i+1
}

#Turn this into a vector because R datatypes R fun :)
top_features_w_vec <- unlist(top_features_temp)

#Grab unique values from the above vector
top_features <- unique(top_features_w_vec)

#Keep only our chosen features to train the model on
df_upsample_features <- df_upsample[, c(target_var, top_features)] #or features

#Factorize the predictor column (because 1 and 0 aren't good enough for R)
df_upsample_features[,target_var] <- factor(df_upsample_features[,target_var])

#Create a trainControl object to carry out 5fold cross-validation
train_control <- trainControl(method="cv",number=5)

#Split out the training features and predictor column
upsample_noclass <- subset(df_upsample_features,select=-CtrlVsCase_Classifier)

#Train the model *****Specify different classifiers using the method parameter*****
model <- train(upsample_noclass, df_upsample_features[,target_var], method = "adaboost", trControl = train_control)

#Prepare the test dataframe with the correct features
test_df <- df_test[, c(target_var, top_features)]  # or features
test_df[,target_var] <- factor(test_df[,target_var])

#Predict the values of the target variable and store in predictions
predictions <- predict(model, test_df)

print(predictions)

#Compute and print the confusion matrix
summ_df <- confusionMatrix(predictions, test_df[,target_var])
print(summ_df)

## Plot ROC curve
library(pROC)

rocCurve <- pROC::roc(response=test_df[,1],
                      predictor=as.numeric(predictions),
                      level=levels(test_df[,1]))
plot(rocCurve, legacy.axes=TRUE)

##---------------------------------------------------------
## Use a random forest on the leading edge gene set to see if 
## we can determine some variable importance
##---------------------------------------------------------
## Code taken from https://compgenomr.github.io/book/trees-and-forests-random-forests-in-action.html#trees-to-forests
set.seed(99)
rfFit <- train(CtrlVsCase_Classifier~.,
               data=df_upsample_features,  ## or use df_upsample_features for GSEA gene set
               method="ranger", 
               trControl=train_control,
               importance="permutation", 
               tuneGrid=data.frame(mtry=7,
                                   min.node.size=1,
                                   splitrule="gini")
)

predictions <- predict(rfFit, test_df)

print(predictions)

#Compute and print the confusion matrix
summ_df <- confusionMatrix(predictions, test_df[,target_var])
print(summ_df)

## variable importance 
plot(varImp(rfFit),top=20)

rocCurve <- pROC::roc(response=test_df[,1],
                      predictor=as.numeric(predictions),
                      level=levels(test_df[,1]))
plot(rocCurve, legacy.axes=TRUE)

##---------------------------------------------------------
## Try using all of the genes from all of the leading edge 
## gene sets
##---------------------------------------------------------
## Code taken from https://compgenomr.github.io/book/trees-and-forests-random-forests-in-action.html#trees-to-forests
set.seed(99)
df_all_leading_genes <- df_upsample[, c(target_var, unique(unlist(topPathwaysRes$leadingEdge)))]
df_all_leading_genes[,target_var] <- factor(df_all_leading_genes[,target_var])
rfFit <- train(CtrlVsCase_Classifier~.,
               data=df_all_leading_genes,
               method="ranger", 
               trControl=train_control,
               importance="permutation", 
               tuneGrid=data.frame(mtry=17,
                                   min.node.size=1,
                                   splitrule="gini")
)

## OOB Error
rfFit$finalModel$prediction.error


## variable importance 
plot(varImp(rfFit),top=20)

## Fit a second model on the top nTopFeatures random variables 
nTopFeatures <- 10
rfFitImportance <- varImp(rfFit)$importance 
rfFitImportance$Gene <- rownames(rfFitImportance)
nTopFeatures <- min(nTopFeatures, nrow(rfFitImportance))
topFeatures <- rfFitImportance[order(rfFitImportance$Overall, decreasing=TRUE),][1:nTopFeatures, "Gene"]

df_topFeatures <- df_upsample[, c(target_var, topFeatures)]
df_topFeatures[,target_var] <- factor(df_topFeatures[,target_var])

## Using Adaboost for second classifier
model <- train(CtrlVsCase_Classifier~., data=df_topFeatures, method = "adaboost", trControl = train_control)

## Test Error
test_df <- df_test[, c(target_var, colnames(df_all_leading_genes))]
predictions <- predict(rfFit, test_df)

print(predictions)

#Compute and print the confusion matrix
summ_df <- confusionMatrix(predictions, factor(test_df[,target_var]))
print(summ_df)

rocCurve <- pROC::roc(response=test_df[,1],
                      predictor=as.numeric(predictions),
                      level=levels(factor(test_df[,1])))
plot(rocCurve, legacy.axes=TRUE)

