#########################################
### Simulation
###   Functional Robust Support Vector Machines for Sparse and Irregular Longitudinal Data
###   Link: https://www.tandfonline.com/doi/pdf/10.1080/10618600.2012.680823?needAccess=true
#########################################
setwd("C:\\Users\\user\\Desktop\\KHS\\bagfpca\\Thesis")

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


### different class proportion
err <- list()
# load("RData/sim_C.RData")
# err[["P(A) = 0.5"]] <- result
for (p in c(0.5, 0.4, 0.3, 0.2)) {
  print(paste("P(A) =", p))
  ### construct classification models
  result <- list()   # error rate result
  simm <- 0      # loop index
  num.sim <- 0   # number of simulations
  while (num.sim < 100) {
    start.time <- Sys.time()
    
    simm <- simm + 1
    print(simm)
    set.seed(simm)
    
    ### generate simulated data
    # data <- sim.curve(200, sparsity=5:10, model="A", prop=p)   # different mean and variance
    # data <- sim.curve(200, sparsity=5:10, model="B", prop=p)   # different mean
    data <- sim.curve(200, sparsity=5:10, model="C", prop=p)   # different variance
    
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
                               err.oob = err.bag$err.oob,
                               # sensitivity and specificity
                               sens.single = err.single$sensitivitiy,
                               spec.single = err.single$specificity,
                               sens.major = err.bag$sens.major,
                               spec.major = err.bag$spec.major,
                               sens.oob = err.bag$sens.oob,
                               spec.oob = err.bag$spec.oob))
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
  
  err[[paste("P(A) =", p)]] <- result

  # save(result, file="RData/sim_A.RData")
  # save(result, file="RData/sim_B.RData")
  # save(result, file="RData/sim_C.RData")
  
  # save(err, file="RData/sim_A_p.RData")
  # save(err, file="RData/sim_B_p.RData")
  # save(err, file="RData/sim_C_p.RData")
  
  # save(err, file="RData/sim_A_modify.RData")
  # save(err, file="RData/sim_B_modify.RData")
  save(err, file="RData/sim_C_modify.RData")
  
  # error, se check
  res <- sapply(1:3, function(i){
    paste(lapply(result[!sapply(result, is.null)], function(x){ x[i, ]*100 }) %>% 
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
  print(res)
}





## 결과 정리
for (i in 1:4) {
  result <- err[[i]]
  result <- result[!sapply(result, is.null)]
  
  res <- sapply(1:9, function(i){
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
  rownames(res) <- c("Single","Majority vote","OOB weight",
                     "Sens(Single)","Spec(Single)","Sens(Majority)","Spec(Majority)","Sens(OOB)","Spec(OOB)")
  
  print(names(err)[i])
  print(xtable(res))
}






