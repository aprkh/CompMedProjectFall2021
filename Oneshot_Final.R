##---------------------------------------------------------
## ctrl vs limb analysis
## code adapted from 
## http://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
##---------------------------------------------------------

## libraries 
library(fgsea)
library(data.table)
library(ggplot2)
library(caret)

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

df_ctrl <- read.csv("ctrl_vs_case.csv")

df_oneshot <- df_train[FALSE,]

#Set up parallel backend
library(parallel)
library(doParallel)
library(foreach)

core_cluster <- makeCluster(
  7,
  type = "FORK"
)

doParallel::registerDoParallel(cl = core_cluster)

#Save off indices for case & ctrl training datapoints
df_train_p <- which(df_train[, c("CtrlVsCase_Classifier")] == 1)
df_train_n <- which(df_train[, c("CtrlVsCase_Classifier")] == 0)

#Set the fraction of rows used to generate differential samples from
pct_comparison <- 0.5

#Set the proportion of case:control differential samples
p_frac <- 0.5
set.seed(99)

#If reproducing, load df_oneshot_50pct.RData instead to reduce wait
for (row_1 in 1:nrow(df_train)){
  cat('\r',row_1)
  random_positive <- sample(df_train_p, size=(floor(length(df_train_p)*p_frac*pct_comparison)))
  random_negative <- sample(df_train_n, size=(floor(length(df_train_n)*(1-p_frac)*pct_comparison)))
  row_temp <- foreach(row_2 = random_positive, .combine="rbind") %dopar% {
    if (row_1 != row_2){
      df_train[row_1,]-df_train[row_2,]
    }
  }
  df_oneshot <- rbind(df_oneshot, row_temp)
  
  row_temp <- foreach(row_2 = random_negative, .combine="rbind") %dopar% {
    if (row_1 != row_2){
      df_train[row_1,]-df_train[row_2,]
    }
  }
  df_oneshot <- rbind(df_oneshot, row_temp)
}

#Stop the cluster & unregister the backend
stopCluster(cl=core_cluster)
unregister_dopar()

#Treat 0-1 same as 1-0
df_oneshot[,1] <- abs(df_oneshot[,1])

test_processed <- df_test_le[FALSE,]

#Set up the parallel backend
core_cluster <- makeCluster(
  7,
  type = "FORK"
)

doParallel::registerDoParallel(cl = core_cluster)

#Generate n_predictors*2 one-shot samples for each test point
cls_true <- rep(0,0)
n_predictors <- 3

for (row_1 in 1:nrow(df_test_le)){
  cat('\r',row_1)
  random_positive <- sample(df_train_p, size=n_predictors)
  random_negative <- sample(df_train_n, size=n_predictors)
  row_temp <- foreach(row_2 = random_positive, .combine="rbind") %dopar% {
    df_test_le[row_1,]-df_train_le[row_2,]
  }
  test_processed <- rbind(test_processed, row_temp)
  
  row_temp <- foreach(row_2 = random_negative, .combine="rbind") %dopar% {
    df_test_le[row_1,]-df_train_le[row_2,]
  }
  test_processed <- rbind(test_processed, row_temp)
}

#Stop the cluster and unregister dopar
stopCluster(cl=core_cluster)
unregister_dopar()

#Treat 1-0 same as 0-1
test_processed[,1] <- abs(test_processed[,1])
test_processed[,1] <- as.factor(test_processed[,1])
