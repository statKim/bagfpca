#########################################
### Simulation
###   Functional Robust Support Vector Machines for Sparse and Irregular Longitudinal Data
###   Link: https://www.tandfonline.com/doi/pdf/10.1080/10618600.2012.680823?needAccess=true
#########################################
setwd("C:/Users/user/Desktop/KHS/bagfpca")

library(tidyverse)    # dplyr, ggplot2, etc
library(doParallel)   # parallel computing
library(fdapace)      # functional PCA
library(e1071)        # SVM, Naive bayes
library(MASS)         # LDA, QDA
library(data.table)   # list rbind
library(gridExtra)    # subplot in ggplot2
library(reshape2)     # melt function
library(xtable)
source("R/bagFPCA.R")

# parallel computing setting
ncores <- detectCores() - 3
registerDoParallel(ncores)

packages <- c("fdapace","e1071","MASS","tidyverse")   # foreach에서 사용할 package 정의


### construct classification models
sim_model <- "A"   # different mean and variance
# sim_model <- "B"   # different mean
# sim_model <- "C"   # different variance
p <- 0.5   # class proportion
result <- list()   # error rate result
simm <- 0      # loop index
num.sim <- 0   # number of simulations
while (num.sim < 500) {
  start.time <- Sys.time()
  
  simm <- simm + 1
  print(simm)
  set.seed(simm)
  
  ### generate simulated data
  data <- sim.curve(200, sparsity=5:10, model=sim_model, prop=p)
  
  ### train, test split
  train_test <- train_test_split(data, train.prop = 1/2)
  X.train <- train_test$X.train
  X.test <- train_test$X.test
  y.train <- train_test$y.train
  y.test <- train_test$y.test
  
  ### Single classifier
  err.single <- tryCatch({ 
    get_single_err(X.train, X.test, y.train, y.test) 
  }, error = function(e) {
    print(e)
    return(FALSE) 
  })
  # check that error occurs
  if ( isFALSE(err.single) ) {
    next
  }
  
  ### Bagging
  err.bag <- tryCatch({ 
    get_bag_err(X.train, X.test, y.train, y.test, B = 100, packages = packages, hyper.para = err.single$hyper.para) 
  }, error = function(e) {
    print(e)
    return(FALSE) 
  })
  # check that error occurs
  if ( isFALSE(err.bag) ) {
    next
  } else {
    num.sim <- num.sim + 1
  }
  
  end.time <- Sys.time()
  print(end.time - start.time)
  
  # save result
  res <- as.data.frame(rbind(err.single = err.single$err.single,
                             err.majority = err.bag$err.majority,
                             err.oob = err.bag$err.oob))
  colnames(res) <- c("Logit","SVM(Linear)","SVM(Gaussian)","LDA","QDA","NaiveBayes")
  
  result[[simm]] <- list(error = res,
                         single.K = err.single$K,   # number of PCs
                         bag.K = sapply(err.bag$y.pred, function(x){ x$K }))   # number of PCs
  
  # print number of repetitions
  print(paste("# of sim =", num.sim))
  
  # 10개마다 RData 저장
  if (num.sim %% 10 == 0) {
    save(result, file=paste("RData/sim_", sim_model, ".RData", sep=""))
  }
}

save(result, file=paste("RData/sim_", sim_model, ".RData", sep=""))



# avg error rate and standard error
res <- sapply(1:3, function(i){
  paste(lapply(result[!sapply(result, is.null)], function(x){ x$error[i, ]*100 }) %>% 
          rbindlist %>% 
          colMeans %>% 
          round(2),
        " (",
        apply(lapply(result[!sapply(result, is.null)], 
                     function(x){ x$error[i, ]*100 }) %>% 
                rbindlist, 2, sd) %>% 
          round(2),
        ")",
        sep="")
}) %>% 
  t() %>% 
  as.data.frame
colnames(res) <- c("Logit","SVM(Linear)","SVM(Gaussian)","LDA","QDA","NaiveBayes")
rownames(res) <- c("Single","Majority vote","OOB weight")
print(res)


# number of selected PCs
single.K <- sapply(result[!sapply(result, is.null)],
                   function(x){ x$single.K })
bag.K <- lapply(result[!sapply(result, is.null)],
                function(x){ x$bag.K })



