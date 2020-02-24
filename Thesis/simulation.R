#########################################
### Simulation
###   Functional Robust Support Vector Machines for Sparse and Irregular Longitudinal Data
###   Link: https://www.tandfonline.com/doi/pdf/10.1080/10618600.2012.680823?needAccess=true
#########################################
setwd("C:\\Users\\user\\Desktop\\KHS\\FDA-Lab\\Thesis")

library(tidyverse)    # dplyr, ggplot2, etc
library(doParallel)   # parallel computing
library(fdapace)      # functional PCA
library(e1071)        # SVM, Naive bayes
library(MASS)         # LDA, QDA
library(data.table)   # list rbind
library(gridExtra)    # subplot in ggplot2
library(reshape2)     # melt function
library(xtable)
source("bagFPCA.R")

# parallel computing setting
ncores <- detectCores() - 2
registerDoParallel(ncores)

packages <- c("fdapace","e1071","MASS","tidyverse")   # foreach에서 사용할 package 정의


### different sample size
err <- list()
# for (n in c(100, 200, 400, 600, 800, 1000)) {
for (n in c(400, 600, 800, 1000)) {
  print(paste("n =", n))
  ### construct classification models
  result <- list()   # error rate result
  simm <- 0      # loop index
  num.sim <- 0   # number of simulations
  while (num.sim < 100) {
    simm <- simm + 1
    print(simm)
    set.seed(simm)
    
    ### generate simulated data
    data <- sim.curve(n, sparsity=5:10, model="A")   # different mean and variance
    # data <- sim.curve(200, sparsity=5:10, model="B")   # different mean
    # data <- sim.curve(200, sparsity=5:10, model="C")   # different variance
    
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
    
    # save result
    res <- as.data.frame(rbind(err.single = err.single$err.single,
                               err.majority = err.bag$err.majority,
                               err.oob = err.bag$err.oob))
    colnames(res) <- c("Logit","SVM(Linear)","SVM(Gaussian)","LDA","QDA","NaiveBayes")
    
    result[[simm]] <- res
    
    # print number of repetitions
    print(paste("# of sim =", num.sim))
    
    # # 10개마다 RData 저장
    # if (num.sim %% 10 == 0) {
    #   # save(result, file="RData/sim_A.RData")
    #   # save(result, file="RData/sim_B.RData")
    #   # save(result, file="RData/sim_C.RData")
    # }
  }
  
  err[[paste("n=", n, sep="")]] <- result

  # save(result, file="RData/sim_A.RData")
  # save(result, file="RData/sim_B.RData")
  # save(result, file="RData/sim_C.RData")
  save(err, file="RData/sim_A_n.RData")
}





## 결과 정리
result <- result[!sapply(result, is.null)]

res <- sapply(1:3, function(i){
  paste(lapply(result, function(x){ x[i, ]*100 }) %>% 
          rbindlist %>% 
          colMeans %>% 
          round(1),
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

xtable(res)





