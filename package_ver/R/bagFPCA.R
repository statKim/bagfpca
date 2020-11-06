#---------------------------------------------------------------
# Bootstrap aggregated classification models using sparse FPCA
#---------------------------------------------------------------
require(fdapace)
require(tidyverse)
require(doParallel)
require(data.table)
require(e1071)
require(class)
require(MASS)


# The packages to use "foreach" statement
packages <- c("fdapace","e1071","class","MASS")


### Fit the bagged classifier using sparse FPCA
# - Provided classification models
#   - Logistic regression("glm" function)
#   - SVM("svm" function in "e1071" package)
#   - KNN("knn" function in "class" package)
#   - LDA("lda" function in "MASS" package)
#   - qDA("qda" function in "MASS" package)
#   - Naive bayes("naiveBayes" function in "e1071" package)
# - Input parameters
#   - X.train : input data for "FPCA" function("fdapace" package)
#   - y.train : factor type vector of train set
#   - selectK : the method of select the number of FPCs("FVE","AIC","BIC"), "FVE" is default.
#   - method : classification method("logit", "svm", "knn)
#   - kernel : (only method="svm")option of "svm" function("radial","linear","polynomial","sigmoid"), "radial" is default.
#   - B : the number of bootstrap replication
#   - ncores : the number of cores to use parallel computing
#   - seed : seed number for reproducing
bagFPCA <- function(X.train, y.train, selectK="FVE", method="logit", kernel=NULL, B=10, ncores=NULL, seed=777) {
  # check input data type
  if (check_data(X.train) == FALSE) {
    stop("Input data of X.train is not appropriate.(Not input type of FPCA in fdapace package)")
  }
  if (!is.factor(y.train)) {
    stop("Input data of y.train is not appropriate.(Not factor)")
  }
  
  # FPCA option
  dtype <- check_dtype(X.train)
  if (dtype == "Sparse") {
    FPC.method <- "CE"
  } else {
    FPC.method <- "IN"
  }
  
  # parallel computing setting
  if (is.null(ncores)) {
    ncores <- ceiling(detectCores() / 3)
  }
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  # for reproducing
  set.seed(seed)
  
  # Bootstrap aggregating
  N <- length(y.train)
  boot.mat <- matrix(sample(1:N, N*B, replace=T), N, B)
  y.pred <- foreach(b=1:B, .packages=packages) %dopar% {
    # fit FPCA for bootstrapped data
    boot.ind <- boot.mat[, b]   # bootstraped sample index
    fpca.fit <- FPCA(X.train$Ly[boot.ind], 
                     X.train$Lt[boot.ind], 
                     optns=list(dataType = dtype,
                                methodXi = FPC.method,
                                methodSelectK = selectK,
                                FVEthreshold = 0.99,
                                verbose = F))
    k <- fpca.fit$selectK   # optimal number of PCs
    fpc.score <- fpca.fit$xiEst
    
    # train FPC scores
    train.fpc <- data.frame(y = y.train[boot.ind], 
                            x = fpc.score)
    colnames(train.fpc) <- c("y", paste("FPC", 1:k, sep=""))
    
    # fit the classification model
    if (method == "svm") {
      if (kernel == "linear") {
        tune.linear <- tune.svm(y ~ ., 
                                data = train.fpc, 
                                kernel = "linear",
                                cost = c(10^(-3:1), 2^(5:10)) )
        fit.class <- svm(y ~ ., 
                         data = train.fpc, 
                         kernel = "linear", 
                         cost = tune.linear$best.model$cost)
      } else if (kernel == "sigmoid") {
        tune.sigmoid <- tune.svm(y ~ ., 
                                 data = train.fpc, 
                                 kernel = "sigmoid",
                                 cost = c(10^(-3:1), 2^(5:10)), 
                                 gamma = c(10^(-3:1), 2^(5:10)))
        fit.class <- svm(y ~ ., 
                         data = train.fpc, 
                         kernel = "sigmoid", 
                         cost = tune.sigmoid$best.model$cost,
                         gamma = tune.sigmoid$best.model$gamma)
      } else if (kernel == "polynomial") {
        tune.poly <- tune.svm(y ~ ., 
                              data = train.fpc, 
                              kernel = "polynomial",
                              cost = c(10^(-3:1), 2^(5:10)), 
                              gamma = c(10^(-3:1), 2^(5:10)))
        fit.class <- svm(y ~ ., 
                         data = train.fpc,
                         kernel = "polynomial",
                         cost = tune.poly$best.model$cost,
                         gamma = tune.poly$best.model$gamma,
                         degree = tune.poly$best.model$degree)
      } else {
        tune.radial <- tune.svm(y ~ ., 
                                data = train.fpc, 
                                kernel = "radial",
                                cost = c(10^(-3:1), 2^(5:10)), 
                                gamma = c(10^(-3:1), 2^(5:10)))
        fit.class <- svm(y ~ ., 
                         data = train.fpc, 
                         kernel = "radial", 
                         cost = tune.radial$best.model$cost,
                         gamma = tune.radial$best.model$gamma)
      }
    } else if (method == "logit") {
      fit.class <- glm(y~., train.fpc, family=binomial)
    } else if (method == "knn") {
      fit.class <- tune.knn(x = train.fpc[, -1],
                            y = train.fpc[, 1],
                            k = 1:ceiling(N/2))$best.model
    } else if (method == "lda") {
      fit.class <- lda(y~., train.fpc)
    } else if (method == "qda") {
      fit.class <- qda(y~., train.fpc)
    } else if (method == "naivebayes") {
      fit.class <- naiveBayes(y~., train.fpc)
    }
    
    # OOB FPC scores
    id.oob <- get_oob_index(X.train, boot.ind)   # predictable out-of-bag curve ID
    oob.ind <- which(X.train$Lid %in% id.oob)    # predictable out-of-bag curve index
    fpc.score.oob <- predict(fpca.fit, 
                             X.train$Ly[oob.ind], 
                             X.train$Lt[oob.ind], 
                             K = k, 
                             xiMethod = FPC.method)
    oob.fpc <- data.frame(y = y.train[oob.ind], 
                          x = fpc.score.oob)
    colnames(oob.fpc) <- c("y", paste("FPC", 1:k, sep=""))
    
    # OOB error rate
    if (method == "svm") {
      pred.oob <- predict(fit.class, oob.fpc)
    } else if (method == "logit") {
      pred.oob <- predict(fit.class, oob.fpc, type="response")
      pred.oob <- factor(ifelse(pred.oob > 0.5, 1, 0), levels=c(0, 1))
    } else if (method == "knn") {
      pred.oob <- knn(train = train.fpc, test = oob.fpc, cl = train.fpc$y, k = fit.class$k)
    } else if (method %in% c("lda","qda")) {
      pred.oob <- predict(fit.class, oob.fpc)$class
    } else if (method == "naivebayes") {
      pred.oob <- predict(fit.class, oob.fpc, type='class')
    }
    oob.error <- mean(pred.oob != oob.fpc$y)
    
    # training error rate
    if (method == "svm") {
      pred.train <- predict(fit.class, train.fpc)
    } else if (method == "logit") {
      pred.train <- predict(fit.class, train.fpc, type="response")
      pred.train <- factor(ifelse(pred.train > 0.5, 1, 0), levels=c(0, 1))
    } else if (method == "knn") {
      pred.train <- knn(train = train.fpc, test = train.fpc, cl = train.fpc$y, k = fit.class$k)
    } else if (method %in% c("lda","qda")) {
      pred.train <- predict(fit.class, train.fpc)$class
    } else if (method == "naivebayes") {
      pred.train <- predict(fit.class, train.fpc, type='class')
    }
    train.error <- mean(pred.train != train.fpc$y)
    
    return( list(fpca.fit = fpca.fit,
                 K = k,
                 boot.index = boot.ind,
                 model = fit.class,
                 oob.error = oob.error,
                 train.error = train.error) )
  }
  
  # # define method of SVM for kernel
  # if (method == "svm") {
  #   if (is.null(kernel)) {
  #     method <- "svm.radial"
  #   } else {
  #     method <- paste(method, kernel, sep="")
  #   }
  # }
  
  
  K <- sapply(y.pred, function(x) { x$K }) 
  fpca.fit <- lapply(y.pred, function(x) { x$fpca.fit })
  boot.index <- lapply(y.pred, function(x) { x$boot.index })
  model <- lapply(y.pred, function(x) { x$model })
  OOB.error <- sapply(y.pred, function(x) { x$oob.error }) 
  train.error <- sapply(y.pred, function(x) { x$train.error }) 
  
  result <- list(dtype = dtype,
                 FPC.method = FPC.method,
                 method = method,
                 kernel = kernel,
                 fpca.fit = fpca.fit,
                 K = K,
                 boot.index = boot.index,
                 model = model,
                 B = B,
                 OOB.error = OOB.error,
                 train.error = train.error)
  class(result) <- "bagFPCA"   # define class of object
  
  stopCluster(cl)    # End the parallel computing
  
  return( result )
}
