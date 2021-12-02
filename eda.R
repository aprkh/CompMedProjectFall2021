###########################################################
## Exploratory Data Analysis of bulbar vs limb data
###########################################################

library(DESeq2)
library(pheatmap)
library(ggfortify)

load("bulbar_vs_limb/bulbar_data")  # no preprocessing

## cluster based on gene expressions, top 100 variance
df_bulbar <- df_bulbar[order(df_bulbar$SiteOnset_Class),]
X <- t(as.matrix(subset(df_bulbar, select=c(-Participant_ID, -SiteOnset_Class))))
V <- apply(X, 1, var)

# annotate the observations by columns (each column is a person)
colnames(X) <- df_bulbar$Participant_ID
source_name <- ifelse(df_bulbar$SiteOnset_Class==1, "limb", "bulbar")
group <- ifelse(df_bulbar$SiteOnset_Class==1, 1, 0)
annot <- data.frame(source_name=source_name, group=group)
rownames(annot) <- df_bulbar$Participant_ID

# select the genes with the top 100 variance accross samples
selectedGenes <- names(V[order(V, decreasing=TRUE)][1:100])
pheatmap(X[selectedGenes,], scale="row", show_rownames=FALSE,
         annotation_col=annot)

## PCA plot 
library(stats)
library(ggplot2)

#transpose the matrix 
M <- t(X[selectedGenes,])

#compute PCA 
pcaResults <- prcomp(M)

#plot PCA results making use of ggplot2's autoplot function 
#ggfortify is needed to let ggplot2 know about PCA data structure
autoplot(pcaResults, data=annot, colour="group")
