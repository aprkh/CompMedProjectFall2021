##---------------------------------------------------------
## Bulbar vs limb analysis
## code adapted from 
## http://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
##---------------------------------------------------------

## libraries 
library(fgsea)
library(data.table)
library(ggplot2)

gene_ex_data <- "bulbar_vs_limb/bulbar_data"
load(gene_ex_data)

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
plotEnrichment(pathways[["GOBP_REGULATION_OF_NEURONAL_SYNAPTIC_PLASTICITY"]], 
               ranks) + labs(title="ALS")
plotEnrichment(pathways[["GOBP_MOTOR_NEURON_APOPTOTIC_PROCESS"]], 
               ranks) + labs(title="ALS")

##---------------------------------------------------------
## Train a classifier on the leading edge gene set
##---------------------------------------------------------

# feature selection
features <- fgseaRes[order(fgseaRes$pval),]$leadingEdge[[1]]
df_upsample_features <- df_upsample[, c("SiteOnset_Class", features)]

logit_model <- glm(SiteOnset_Class~., family=binomial, data=df_upsample_features)
summary(logit_model)
