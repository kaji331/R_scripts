library(mlr)
library(parallelMap)
load("~/Projects/R/R_scripts/数据/classify_3_4.task.rda")
generateFilterValuesData(classify_3_4.task,"oneR")$data %>% 
{ggplot(.,aes(name,oneR)) + geom_col()}

ctrl <- makeFeatSelControlExhaustive(maxit=100)
rdesc <- makeResampleDesc("Subsample",predict="both",iters=100,stratify=T)
parallelStartMulticore(4)
l <- c("C50","ctree","featureless","fnn","lda","multinom","naiveBayes","nnet",
       "oneR","randomForest","rda","rpart","svm")
L <- vector()
for (i in l) {
  sfeats <- selectFeatures(learner=paste("classif",i,sep="."),
                           task=classify_3_4.task,resampling=rdesc,control=ctrl,
                           measures=list(mmce,ber))
  L <- c(L,i,sfeats$x,sfeats$y)
  print(L)
}
print(L)
parallelStop()
