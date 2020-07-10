setwd("Simulation")

library(tidyverse)    # dplyr, ggplot2, etc
library(doParallel)   # parallel computing
library(fdapace)      # functional PCA
library(e1071)        # SVM, Naive bayes
library(MASS)         # LDA, QDA
library(data.table)   # list rbind
library(gridExtra)    # subplot in ggplot2
library(reshape2)     # melt function
source("R/bagFPCA.R")

# parallel computing setting
ncores <- detectCores() - 2
registerDoParallel(ncores)

packages <- c("fdapace","e1071","MASS","tidyverse")   # foreach에서 사용할 package 정의

### load data and preprocessing
data <- read.csv("data/spnbmd.csv", header=T)

ind <- data %>%
  # filter(ethnic == "Hispanic") %>%
  group_by(idnum) %>%
  summarise(n=n()) %>%
  filter(n >= 2) %>%
  dplyr::select(idnum)
data <- data[which(data$idnum %in% ind$idnum), ]


# gender classification
data <- data %>% 
  mutate(y = factor(ifelse(sex == "fem", 1, 0), levels=c(0, 1))) %>% 
  dplyr::select(y, idnum, age, spnbmd)
length(unique(data$idnum))   # 280

colnames(data) <- c("y","id","time","val")


### construct classification models
result <- list()
set.seed(1000)
seed <- sample(1:10000, 100)
for (simm in 1:100) {
  # simm <- 10
  print( paste(simm, ":", seed[simm]) )
  set.seed(seed[simm])
  
  start.time <- Sys.time()
  
  ### train, test split
  train_test <- train_test_split(data, train.prop = 2/3)
  X.train <- train_test$X.train
  X.test <- train_test$X.test
  y.train <- train_test$y.train
  y.test <- train_test$y.test
  
  ### Single classifier
  err.single <- get_single_err(X.train, X.test, y.train, y.test)
  
  ## Bagging
  # Bootstrap aggregating
  err.bag <- get_bag_err(X.train, X.test, y.train, y.test, B = 50, packages = packages, hyper.para = err.single$hyper.para)
  
  end.time <- Sys.time()
  print(end.time - start.time)
  
  # save result
  res <- as.data.frame(rbind(err.single = err.single$err.single,
                             err.majority = err.bag$err.majority,
                             err.oob = err.bag$err.oob))
  colnames(res) <- c("Logit","SVM(Linear)","SVM(Gaussian)","LDA","QDA","NaiveBayes")
  
  result[[simm]] <- res
}

## 결과 정리
result <- result[!sapply(result, is.null)]

res <- sapply(1:3, function(i){
  paste(lapply(result, function(x){ x[i, ]*100 }) %>% 
          rbindlist %>% 
          colMeans %>% 
          round(2),
        " (",
        apply(lapply(result[!sapply(result, is.null)], 
                     function(x){ x[i, ]*100 }) %>% 
                rbindlist, 2, sd) %>% 
          round(2),
        ")",
        sep="")
}) %>% 
  t() %>% 
  as.data.frame
colnames(res) <- c("Logit","SVM(Linear)","SVM(Gaussian)","LDA","QDA","NaiveBayes")
rownames(res) <- c("Single","Majority vote","OOB weight")
res


plot(FPC1, FPC2, data=train.fpc)

ggplot(train.fpc, aes(FPC1, FPC2, col=y)) + 
  geom_point()




  
result[[1]][1:3, ]

library(ranger)
rf.fit <- ranger(y~., train.fpc, num.trees=100)
pred <- predict(rf.fit, test.fpc)
mean(test.fpc$y != pred$predictions)

library(randomForest)
set.seed(100)
rf.fit <- randomForest(y~., train.fpc, ntree=100)
rf.fit
plot(rf.fit)
pred <- predict(rf.fit, test.fpc)
mean(test.fpc$y != pred)

library(gbm)
gbm.fit <- gbm(y~., data=train.fpc, distribution="bernoulli", n.trees = 1000, shrinkage = 0.1,             
               interaction.depth = 2, verbose=T)
best.iter <- gbm.perf(gbm.fit, method = "OOB")
pred <- predict(gbm.fit, test.fpc, n.trees = best.iter)
mean(test.fpc$y != pred)
