##---------------------------------------------------------
## Bulbar vs limb analysis
## code adapted from 
## http://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
##---------------------------------------------------------

## libraries 
library(fgsea)
library(data.table)
library(ggplot2)

df_bulbar <- read.csv("bulbar_vs_limb.csv")
#load(gene_ex_data)

##---------------------------------------------------------
## Split into train/test
##---------------------------------------------------------
# limb = 1, bulbar = 0
limb_rows <- which(df_bulbar[, c("SiteOnset_Class")] == 1)
bulbar_rows <- which(df_bulbar[, c("SiteOnset_Class")] == 0)

# stratified test/train sampling from both classes
set.seed(99)
percentageTest <- 30

nTestLimb <- floor(percentageTest * length(limb_rows) / 100)
limb_test <- sort(sample(limb_rows, nTestLimb, replace=FALSE))

nTestBulbar <- floor(percentageTest * length(bulbar_rows) / 100)
bulbar_test <- sort(sample(bulbar_rows, nTestBulbar, replace=FALSE))

test_rows <- sort(c(limb_test, bulbar_test))

df_test <- df_bulbar[test_rows,]
df_train <- df_bulbar[-test_rows,]

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

## Compute t statistics
tranks <- sapply(X, function(x) t.test(x ~ y)$statistic)

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
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways[topPathways], ranks, fgseaRes, 
              gseaParam=0.5)

#Grab a new dataframe containing only the GSEA results for the top performing pathways
topPathwaysRes <- fgseaRes[fgseaRes$pathway %in% topPathways, ]
#print(topPathwaysRes[order(fgseaRes$pval),]$leadingEdge[[17]])

##---------------------------------------------------------
## Train a classifier on the leading edge gene set
##---------------------------------------------------------
library(caret)

#Grab the leading edge of each of the top pathway
features <- topPathwaysRes[order(topPathwaysRes$pval),]$leadingEdge[[1]]

#Alternatively, grab the *best* singular feature from each top pathway
#This is similar to that paper we did in class
top_features_temp <- topPathwaysRes$leadingEdge
i = 1
for (f in top_features_temp) {
  top_features_temp[[i]] <- f[1]
  i <- i+1
}

#Turn this into a vector because R datatypes R fun :)
top_features_w_vec <- unlist(top_features_temp)

#Grab unique values from the above vector
top_features <- unique(top_features_w_vec)

#Keep only our chosen features to train the model on
df_upsample_features <- df_upsample[, c("SiteOnset_Class", top_features)] #or features

#Factorize the predictor column (because 1 and 0 aren't good enough for R)
df_upsample_features$SiteOnset_Class <- factor(df_upsample_features$SiteOnset_Class)

#Create a trainControl object to carry out 5fold cross-validation
train_control <- trainControl(method="cv",number=5)

#Split out the training features and predictor column
upsample_noclass <- subset(df_upsample_features,select=-SiteOnset_Class)

#Train the model *****Specify different classifiers using the method parameter*****
model <- train(upsample_noclass, df_upsample_features$SiteOnset_Class, method = "adaboost", trControl = train_control)

#Prepare the test dataframe with the correct features
test_df <- df_test[, c("SiteOnset_Class", features)]
test_df$SiteOnset_Class <- factor(test_df$SiteOnset_Class)

#Predict the values of SiteOnset_Class and store in predictions
predictions <- predict(model, test_df)

print(predictions)

#Compute and print the confusion matrix
summ_df <- confusionMatrix(predictions, test_df$SiteOnset_Class)
print(summ_df)