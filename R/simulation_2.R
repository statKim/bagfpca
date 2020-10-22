#########################################
### Simulation
###   Probability-enhanced effective dimension reduction
###   for classifying sparse functional data
###   Yao et al.
###   Link: https://link.springer.com/content/pdf/10.1007%2Fs11749-015-0470-2.pdf
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
source("R/bagFPCA.R")

# parallel computing setting
ncores <- detectCores() - 3
registerDoParallel(ncores)

packages <- c("fdapace","e1071","MASS","tidyverse")   # foreach에서 사용할 package 정의


# res.sim <- list()  # classification result of each bootstrap resamples
sim_model <- "A"
# sim_model <- "B"
# sim_model <- "C"
result <- list()   # error rate result
simm <- 0          # loop index
num.sim <- 0       # number of simulations
while (num.sim < 500) {
  simm <- simm + 1
  print(simm)
  
  start.time <- Sys.time()
  
  ### generate the simulated dataset
  set.seed(simm)
  n <- 200   # number of observations
  j <- 1:50  # number of basis
  
  # parameters when generate class label
  b <- matrix(ifelse(j <= 2, 1, (j-2)^(-3)), ncol=1)
  
  ## generate curves
  data <- list()
  for (i in 1:n) {
    # random sparsify
    num.obs <- sample(10:20, 1)
    t <- sort(runif(num.obs, 0, 10))
    
    # 101~200 are test set
    if (i > 100) {
      range.train <- range(unlist(data$Lt[1:100]))
      while( (max(t) > range.train[2]) | (min(t) < range.train[1]) ) {
        t <- sort(runif(num.obs, 0, 10))
      }
    }
    
    # eigenfunctions
    phi <- sapply(j, function(x){
      if (x %% 2 == 0) {
        sqrt(1/5)*sin(pi*t*x/5)
      } else {
        sqrt(1/5)*cos(pi*t*x/5)
      }
    })
    
    # generate PC scores
    xi <- sapply(j, function(x){ rnorm(1, 0, sqrt( x^(-1.5) )) })
    
    # parameters when generate class label
    beta.1 <- phi %*% b
    beta.2 <- matrix(sqrt(3/10)*(t/5-1), ncol=1)
    
    # measurement error
    eps <- rnorm(num.obs, 0, sqrt(0.1))
    
    # generate the curve
    X <- xi %*% t(phi) + eps
    
    # generate class label
    eps <- rnorm(1, 0, sqrt(0.1))   # model error
    if (sim_model == "A") {
      fx <- exp((X %*% beta.1)/2)-1   # model 2
    } else if (sim_model == "B") {
      fx <- atan(pi*(X %*% beta.1)) + exp((X %*% beta.2)/3) - 1    # model 4
    } else {
      fx <- atan(pi*(X %*% beta.1)/4)
    }
    # fx <- sin(pi*(X %*% beta.1)/4)  # model 1
    # fx <- sin(pi*(X %*% beta.1)/3) + exp((X %*% beta.2)/3) - 1   # model 3
    # fx <- atan(pi*(X %*% beta.1)) + exp((X %*% beta.2)/3) - 1    # model 4
    y <- factor(ifelse(fx + eps < 0, 0, 1), levels=c(0, 1))
    
    data$id[[i]] <- rep(i, num.obs)
    data$y[[i]] <- rep(y, num.obs)
    data$Lt[[i]] <- t
    data$Ly[[i]] <- X
  }
  
  # # plot the generated curves
  # sapply(data, unlist) %>% 
  #   data.frame() %>% 
  #   mutate(y = ifelse(unlist(data$y) == 0, "G1", "G2")) %>% 
  #   ggplot(aes(x=Lt, y=Ly, group=id, color=y)) +
  #   geom_line() +
  #   xlab("Time") +
  #   ylab("") +
  #   theme_bw() +
  #   theme(legend.title = element_blank())
  
  
  ### train, test split => train: 1~100, test: 101~200
  id <- 1:n
  ind.train <- 1:100
  ind.test <- 101:200
  
  # split to train, test set
  y <- sapply(data$y, unique)
  X.train <- list(Lid = lapply(data$id[ind.train], unique),
                  Lt = data$Lt[ind.train],
                  Ly = data$Ly[ind.train])
  X.test <- list(Lid = lapply(data$id[ind.test], unique),
                 Lt = data$Lt[ind.test],
                 Ly = data$Ly[ind.test])
  y.train <- y[ind.train]
  y.test <- y[ind.test]
  
  
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
    save(result, file=paste("RData/simulation_2_", sim_model, ".RData", sep=""))
  }
}

save(result, file=paste("RData/simulation_2_", sim_model, ".RData", sep=""))



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

