##---------------------------------------------------------
## Bulbar vs limb analysis
## code adapted from 
## http://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
##---------------------------------------------------------

## set seed
set.seed(99)   # Wayne Gretzky
#set.seed(14)  # Austin Matthews

## libraries 
library(fgsea)
library(data.table)
library(ggplot2)
library(caret)

## Note - we assume preprocessing steps from bulbar_vs_limb_code.R
load("bulbar_vs_limb/bulbar_data_preprocessed")

##---------------------------------------------------------
## Running GSEA
##---------------------------------------------------------

## Note: may want to upsample the bulbars (0's in phenotype vector)
## Using SMOTE algorithm: https://www.rdocumentation.org/packages/performanceEstimation/versions/1.1.0/topics/smote
library(performanceEstimation)

df_bulbar_no_ID <- df_train[, 2:ncol(df_train)]
df_bulbar_no_ID$SiteOnset_Class <- as.factor(df_bulbar_no_ID$SiteOnset_Class)
df_upsample <- smote(SiteOnset_Class~., df_bulbar_no_ID, 
                     perc.over=6, perc.under=1, k=5)
df_upsample$SiteOnset_Class <- as.integer(df_upsample$SiteOnset_Class)-1

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

##---------------------------------------------------------
## Try using all of the genes from all of the leading edge 
## gene sets
##---------------------------------------------------------

#Create a trainControl object to carry out 5fold cross-validation
train_control <- trainControl(method="cv",number=5)

## Code taken from https://compgenomr.github.io/book/trees-and-forests-random-forests-in-action.html#trees-to-forests
df_all_leading_genes <- df_upsample[, c("SiteOnset_Class", unique(unlist(topPathwaysRes$leadingEdge)))]
df_all_leading_genes$SiteOnset_Class <- factor(df_all_leading_genes$SiteOnset_Class)
rfFit <- train(SiteOnset_Class~.,
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

df_topFeatures <- df_upsample[, c("SiteOnset_Class", topFeatures)]
df_topFeatures$SiteOnset_Class <- factor(df_topFeatures$SiteOnset_Class)

## Using Adaboost for second classifier
model <- train(SiteOnset_Class~., data=df_topFeatures, method = "adaboost", trControl = train_control)

## Test Error
test_df <- df_test[, c("SiteOnset_Class", colnames(df_all_leading_genes))]
predictions <- predict(model, test_df)

print(predictions)

#Compute and print the confusion matrix
summ_df <- confusionMatrix(predictions, factor(test_df$SiteOnset_Class))
print(summ_df)

rocCurve <- pROC::roc(response=test_df[,1],
                      predictor=as.numeric(predictions),
                      level=levels(factor(test_df[,1])))
plot(rocCurve, legacy.axes=TRUE)

