##################################
### bagging classifier with FPCA
##################################

require(tidyverse)    # dplyr, ggplot2, etc
require(doParallel)   # parallel computing
require(fdapace)      # functional PCA
require(e1071)        # SVM, Naive bayes
require(MASS)         # LDA, QDA
require(data.table)   # list rbind

### majority voting functon
majority_vote <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


### sensistivity and specificity
get_sens_spec <- function(pred, y.test) {
  tb <- table(pred, y.test)
  
  # print(tb)
  if (length(unique(pred)) == 1) {
    if (unique(pred) == 0) {
      sens <- tb[1, 1] / sum(tb[, 1]) 
      spec <- 0
    } else {
      sens <- 0
      spec <- tb[1, 2] / sum(tb[, 2])
    }
  } else {
    sens <- tb[1, 1] / sum(tb[, 1]) 
    spec <- tb[2, 2] / sum(tb[, 2])
  }
  
  return( c(sens, spec) )
}


### train, test split
##  - input data form : cbind(id, time, val, y)
train_test_split <- function(data, train.prop=2/3) {
  # # train, test split => range(test set) > range(train set) 인 경우의 id 제거
  # min.max.grid <- data %>% 
  #   filter(time %in% range(time)) %>% 
  #   dplyr::select(id) %>% 
  #   unique
  # id <- setdiff(unique(data$id), min.max.grid$id)
  # 
  # N <- length(unique(data$id))
  # N.train <- ceiling(N * train.prop)
  # id.train <- c(sample(id, N.train-length(unique(min.max.grid$id))),
  #               min.max.grid$id)
  # id.test <- id[-which(id %in% id.train)]
  
  ### train, test split => train set에 없는 time point 포함시 에러
  # range(test set) > range(train set) => 넘는 부분에 해당하는 id 제거
  id <- unique(data$id)
  N <- length(unique(data$id))
  N.train <- ceiling(N * train.prop)
  id.train <- sample(id, N.train)
  id.test <- NULL
  range.train <- range(data$time[which(data$id %in% id.train)])
  range.test <- range(data$time[-which(data$id %in% id.train)])
  if (range.test[1] < range.train[1]) {
    over.ind <- which(data$time[-which(data$id %in% id.train)] < range.train[1] )
    over.ind <- data$id[-which(data$id %in% id.train)][over.ind]
    id.test <- id[-which(id %in% c(id.train, over.ind))]
  }
  if (range.test[2] > range.train[2]) {
    over.ind <- which(data$time[-which(data$id %in% id.train)] > range.train[2] )
    over.ind <- data$id[-which(data$id %in% id.train)][over.ind]
    if (is.numeric(id.test)) {
      id.test <- intersect(id.test,
                           id[-which(unique(data$id) %in% c(id.train, over.ind))])
    } else {
      id.test <- id[-which(id %in% c(id.train, over.ind))]
    }
  }
  if (!is.numeric(id.test)) {
    id.test <- id[-which(id %in% id.train)]
  }
  
  # transform to FPCA input
  train <- data[which(data$id %in% id.train), ]
  test <- data[which(data$id %in% id.test), ]
  
  X.train <- MakeFPCAInputs(IDs = train$id,
                            tVec = train$time,
                            yVec = train$val)
  X.test <- MakeFPCAInputs(IDs = test$id,
                           tVec = test$time,
                           yVec = test$val)
  y.train <- unique(train[, c("id", "y")])[, "y"]
  y.test <- unique(test[, c("id", "y")])[, "y"]
  
  return(list(X.train = X.train,
              X.test = X.test,
              y.train = y.train,
              y.test = y.test))
}


### Get ID of test and OOB samples
get_id_test_oob <- function(X.train, X.test, boot.ind) {
  # bootstrap sample의 time range를 벗어나는 test set sample 제거
  id.test <- NULL
  max.boot <- max(unlist(X.train$Lt[boot.ind]))
  min.boot <- min(unlist(X.train$Lt[boot.ind]))
  if (min(unlist(X.test$Lt)) < min.boot) {
    # train set에서의 index(Not id)
    over.ind <- which(sapply(X.test$Lt, 
                             function(x){ sum(x < min.boot) }) != 0)
    id.test <- unlist(X.test$Lid)[-over.ind]
  }
  if (max(unlist(X.test$Lt)) > max.boot) {
    # train set에서의 index(Not id)
    over.ind <- which(sapply(X.test$Lt, 
                             function(x){ sum(x > max.boot) }) != 0)
    if (is.numeric(id.test)) {
      id.test <- intersect(id.test,
                           unlist(X.test$Lid)[-over.ind])
    } else {
      id.test <- unlist(X.test$Lid)[-over.ind]
    }
  }
  if (!is.numeric(id.test)) {
    id.test <- unlist(X.test$Lid)
  }
  
  # bootstrap sample의 time range를 벗어나는 OOB sample 제거
  id.oob <- NULL
  if (min(unlist(X.train$Lt[-unique(boot.ind)])) < min(unlist(X.train$Lt[boot.ind]))) {
    # train set에서의 index(Not id)
    over.ind <- which(sapply(X.train$Lt[-unique(boot.ind)], 
                             function(x){ sum( x < min(unlist(X.train$Lt[boot.ind])) ) }) != 0)
    id.oob <- unlist(X.train$Lid)[-unique(boot.ind)][-over.ind]
  }
  if (max(unlist(X.train$Lt[-unique(boot.ind)])) > max(unlist(X.train$Lt[boot.ind]))) {
    # train set에서의 index(Not id)
    over.ind <- which(sapply(X.train$Lt[-unique(boot.ind)], 
                             function(x){ sum( x > max(unlist(X.train$Lt[boot.ind])) ) }) != 0)
    if (is.numeric(id.oob)) {
      id.oob <- intersect(id.oob,
                          unlist(X.train$Lid)[-unique(boot.ind)][-over.ind])
    } else {
      id.oob <- unlist(X.train$Lid)[-unique(boot.ind)][-over.ind]
    }
  }
  if (!is.numeric(id.oob)) {
    id.oob <- unlist(X.train$Lid)[-unique(boot.ind)]
  }
  
  return(list(id.test = id.test,
              id.oob = id.oob))
}


### Fit FPCA and make data using FPC scores
fit_FPCA <- function(X.train, X.test, y.train, y.test, boot.ind=NULL) {
  if (is.null(boot.ind)) {
    ## Not bootstrap (Don't have boot.ind)
    # fit FPCA for training set
    fpca.fit <- FPCA(X.train$Ly, 
                     X.train$Lt, 
                     optns=list(dataType="Sparse",
                                methodXi="CE",
                                # methodSelectK="BIC",
                                FVEthreshold=0.99,
                                verbose=F))
    k <- fpca.fit$selectK   # optimal number of PCs
    fpc.score <- fpca.fit$xiEst
    
    # train FPC scores
    train.fpc <- data.frame(y = y.train, 
                            x = fpc.score)
    colnames(train.fpc) <- c("y", paste("FPC", 1:k, sep=""))
    
    # test FPC scores
    fpc.score.test <- predict(fpca.fit, X.test$Ly, X.test$Lt, K=k, xiMethod="CE")
    test.fpc <- data.frame(y = y.test, 
                           x = fpc.score.test)
    colnames(test.fpc) <- c("y", paste("FPC", 1:k, sep=""))
    
    return(list(train.fpc = train.fpc,
                test.fpc = test.fpc,
                fpca.fit = fpca.fit,
                k = k))
  } else {
    ## Bootstrap (Have boot.ind)
    # get index of test and OOB samples without out-of-range time points
    id_test_oob <- get_id_test_oob(X.train, X.test, boot.ind)
    id.test <- id_test_oob$id.test
    id.oob <- id_test_oob$id.oob
    
    # fit FPCA for bootstrapped data
    fpca.fit <- FPCA(X.train$Ly[boot.ind], 
                     X.train$Lt[boot.ind], 
                     optns=list(dataType="Sparse",
                                methodXi="CE",
                                # methodSelectK="BIC",
                                FVEthreshold=0.99,
                                verbose=F))
    k <- fpca.fit$selectK   # optimal number of PCs
    fpc.score <- fpca.fit$xiEst
    
    # train FPC scores
    train.fpc <- data.frame(y = y.train[boot.ind], 
                            x = fpc.score)
    colnames(train.fpc) <- c("y", paste("FPC", 1:k, sep=""))
    
    # test FPC scores
    test.ind <- which(unlist(X.test$Lid) %in% id.test)
    fpc.score.test <- predict(fpca.fit, 
                              X.test$Ly[test.ind], 
                              X.test$Lt[test.ind], 
                              K=k, 
                              xiMethod="CE")
    test.fpc <- data.frame(y = y.test[test.ind],
                           x = fpc.score.test)
    colnames(test.fpc) <- c("y", paste("FPC", 1:k, sep=""))
    
    # OOB FPC scores
    oob.ind <- which(X.train$Lid %in% id.oob)
    fpc.score.oob <- predict(fpca.fit,
                             X.train$Ly[oob.ind],
                             X.train$Lt[oob.ind],
                             K=k,
                             xiMethod="CE")
    oob.fpc <- data.frame(y = y.train[oob.ind],
                          x = fpc.score.oob)
    colnames(oob.fpc) <- c("y", paste("FPC", 1:k, sep=""))
    
    return(list(train.fpc = train.fpc,
                test.fpc = test.fpc,
                oob.fpc = oob.fpc,
                id.test = id.test,
                test.ind = test.ind,
                fpca.fit = fpca.fit,
                k = k))
  }
}


### Fit classification models
fit_classif <- function(formula, train.fpc, hyper.para=NULL) {
  # hyperparameters
  tune.linear <- hyper.para$tune.linear
  tune.radial <- hyper.para$tune.radial
  
  # fit classifiers
  fit.logit <- glm(formula, train.fpc, family=binomial)
  fit.svm.linear <- svm(formula,
                        data=train.fpc, 
                        kernel="linear", 
                        cost=tune.linear$best.model$cost)
  fit.svm.radial <- svm(formula,
                        data=train.fpc, 
                        kernel="radial", 
                        cost=tune.radial$best.model$cost,
                        gamma=tune.radial$best.model$gamma)
  fit.lda <- lda(formula, train.fpc)
  fit.qda <- qda(formula, train.fpc)
  fit.nb <- naiveBayes(formula, train.fpc)
  
  return(list(fit.logit = fit.logit,
              fit.svm.linear = fit.svm.linear,
              fit.svm.radial = fit.svm.radial,
              fit.lda = fit.lda,
              fit.qda = fit.qda,
              fit.nb = fit.nb))
}


### Predict obtained by fit_classif
predict_classif <- function(fit_classif, test.fpc) {
  pred <- data.frame(logit = factor(ifelse(predict(fit_classif$fit.logit, test.fpc, type="response") > 0.5, 1, 0), levels=c(0, 1)),
                     svm.linear = predict(fit_classif$fit.svm.linear, test.fpc),
                     svm.radial = predict(fit_classif$fit.svm.radial, test.fpc),
                     lda = predict(fit_classif$fit.lda, test.fpc)$class,
                     qda = predict(fit_classif$fit.qda, test.fpc)$class,
                     nb = predict(fit_classif$fit.nb, test.fpc, type="class"))
  
  return(pred)
}


### prediction of aggregated classifier and its classification error rate
agg_classif <- function(pred, method) {
  if (method == "majority") {
    # majority vote
    pred_agg <- pred %>% 
      group_by(id, y) %>% 
      summarise(logit = majority_vote(logit),
                svm.linear = majority_vote(svm.linear),
                svm.radial = majority_vote(svm.radial),
                lda = majority_vote(lda),
                qda = majority_vote(qda),
                nb = majority_vote(nb),
                .groups="drop")   # remove grouping information
  } else if (method == "oob") {
    # OOB error weighted vote
    pred_agg <- pred %>% 
      # if OOB error = 0, w = min(other OOB error)
      mutate(w.logit = 1/ifelse(oob.logit == 0, min(oob.logit[oob.logit > 0]), oob.logit),
             w.svm.linear = 1/ifelse(oob.svm.linear == 0, min(oob.svm.linear[oob.svm.linear > 0]), oob.svm.linear),
             w.svm.radial = 1/ifelse(oob.svm.radial == 0, min(oob.svm.radial[oob.svm.radial > 0]), oob.svm.radial),
             w.lda = 1/ifelse(oob.lda == 0, min(oob.lda[oob.lda > 0]), oob.lda),
             w.qda = 1/ifelse(oob.qda == 0, min(oob.qda[oob.qda > 0]), oob.qda),
             w.nb = 1/ifelse(oob.nb == 0, min(oob.nb[oob.nb > 0]), oob.nb)) %>%
      # # if OOB error = 0, w = runif(1, 0, min(other OOB error))
      # mutate(w.logit = 1/ifelse(oob.logit == 0, runif(1, 0, min(oob.logit[oob.logit > 0])), oob.logit),
      #        w.svm.linear = 1/ifelse(oob.svm.linear == 0, runif(1, 0, min(oob.svm.linear[oob.svm.linear > 0])), oob.svm.linear),
      #        w.svm.radial = 1/ifelse(oob.svm.radial == 0, runif(1, 0, min(oob.svm.radial[oob.svm.radial > 0])), oob.svm.radial),
      #        w.lda = 1/ifelse(oob.lda == 0, runif(1, 0, min(oob.lda[oob.lda > 0])), oob.lda),
      #        w.qda = 1/ifelse(oob.qda == 0, runif(1, 0, min(oob.qda[oob.qda > 0])), oob.qda),
      #        w.nb = 1/ifelse(oob.nb == 0, runif(1, 0, min(oob.nb[oob.nb > 0])), oob.nb)) %>% 
      group_by(id, y) %>% 
      summarise(logit = factor(ifelse(sum(w.logit*as.numeric(logit)) / sum(w.logit) > 1.5, 1, 0), 
                               levels=c(0, 1)),
                svm.linear = factor(ifelse(sum(w.svm.linear*as.numeric(svm.linear)) / sum(w.svm.linear) > 1.5, 1, 0), 
                                    levels=c(0, 1)),
                svm.radial = factor(ifelse(sum(w.svm.radial*as.numeric(svm.radial)) / sum(w.svm.radial) > 1.5, 1, 0), 
                                    levels=c(0, 1)),
                lda = factor(ifelse(sum(w.lda*as.numeric(lda)) / sum(w.lda) > 1.5, 1, 0), 
                             levels=c(0, 1)),
                qda = factor(ifelse(sum(w.qda*as.numeric(qda)) / sum(w.qda) > 1.5, 1, 0), 
                             levels=c(0, 1)),
                nb = factor(ifelse(sum(w.nb*as.numeric(nb)) / sum(w.nb) > 1.5, 1, 0), 
                            levels=c(0, 1)),
                .groups="drop")   # remove grouping information
  }
  
  # classification error rate of aggregated classsifier
  err_agg <- apply(pred_agg[, 3:8], 2, function(x){ mean(x != pred_agg$y) })
  
  return(list(pred = pred_agg,
              error = err_agg))
}



### Single classifier
get_single_err <- function(X.train, X.test, y.train, y.test) {
  # obtain FPC scores of train and test data
  fit.fpc <- fit_FPCA(X.train, X.test, y.train, y.test)
  train.fpc <- fit.fpc$train.fpc
  test.fpc <- fit.fpc$test.fpc
  k <- fit.fpc$k
  
  # tune the hyperparmeters
  tune.linear <- tune.svm(y ~ .,
                          data = train.fpc,
                          kernel = "linear",
                          cost = c(10^(-3:1), 2^(5:10)))
  tune.radial <- tune.svm(y ~ .,
                          data = train.fpc,
                          kernel = "radial",
                          cost = c(10^(-3:1), 2^(5:10)),
                          gamma = c(10^(-3:1), 2^(5:10)))
  
  # fit classifiers
  fit <- fit_classif(y~., train.fpc, hyper.para=list(tune.linear = tune.linear,
                                                     tune.radial = tune.radial))
  
  # predict
  pred <- predict_classif(fit, test.fpc)
  
  # save the error rate
  err.single <- apply(pred, 2, function(x){ mean(x != y.test) })
  
  # # sensitivity and specificity
  # sens_spec <- apply(pred, 2, get_sens_spec, y.test)
  
  return(list(err.single = err.single,
              model = fit,   # fitted model objects
              pred = pred,   # prediction
              hyper.para = list(tune.linear = tune.linear,
                                tune.radial = tune.radial)))
}



### Bagging classifier
get_bag_err <- function(X.train, X.test, y.train, y.test, B = 100, packages = c("fdapace","e1071","MASS"), hyper.para) {
  # list of manual functions
  ftns <- c("fit_FPCA","get_id_test_oob","fit_classif","predict_classif")
  
  N <- length(X.train$Ly)
  boot.mat <- matrix(sample(1:N, N*B, replace=T), N, B)   # index of bootstrap samples
  
  # make predictions from bootstrap samples
  y.pred <- foreach(b=1:B, .packages=packages, .export=ftns) %dopar% {
    boot.ind <- boot.mat[, b]
    
    # obtain FPC scores of train, test, OOB data
    fit.fpc <- fit_FPCA(X.train, X.test, y.train, y.test, boot.ind=boot.ind)
    train.fpc <- fit.fpc$train.fpc
    test.fpc <- fit.fpc$test.fpc
    oob.fpc <- fit.fpc$oob.fpc
    id.test <- fit.fpc$id.test
    test.ind <- fit.fpc$test.ind
    k <- fit.fpc$k
    
    # classifiers with bootstrapped FPC scores
    fit <- fit_classif(y~., train.fpc, hyper.para=hyper.para)
    
    # predict
    pred <- predict_classif(fit, test.fpc)
    
    # OOB error rate
    oob.error <- apply(predict_classif(fit, oob.fpc),
                       2,
                       function(x){ mean(x != oob.fpc$y) })
    
    # train error rate
    train.error <- apply(predict_classif(fit, train.fpc),
                         2,
                         function(x){ mean(x != train.fpc$y) })
    
    return( list(oob.error = oob.error,
                 train.error = train.error,
                 model = fit,   # fitted model objects
                 pred = pred,   # prediction
                 boot = data.frame(id = id.test,
                                   y = y.test[test.ind],
                                   pred,   # prediction of test set (data.frame)
                                   # OOB error
                                   oob.logit = oob.error[1],
                                   oob.svm.linear = oob.error[2],
                                   oob.svm.radial = oob.error[3],
                                   oob.lda = oob.error[4],
                                   oob.qda = oob.error[5],
                                   oob.nb = oob.error[6])) )
  }
  
  # calculate error rate for different number of aggregated models
  for (b in 1:B) {
    # make predictions of classifiers into row-wise dataframe
    pred <- lapply(y.pred[1:b], function(x){ x$boot }) %>% 
      rbindlist() %>% 
      as.data.frame()
    
    # majority vote
    pred.majoriry <- agg_classif(pred, method="majority")
    
    # oob error weighted vote
    pred.oob <- agg_classif(pred, method="oob")
    
    if (b == 1) {
      err.majority <- pred.majoriry$error
      err.oob <- pred.oob$error
    } else {
      err.majority <- rbind(err.majority, 
                            pred.majoriry$error)
      err.oob <- rbind(err.oob, 
                       pred.oob$error)
    }
  }
  rownames(err.majority) <- paste("B=", 1:B, sep="")
  rownames(err.oob) <- paste("B=", 1:B, sep="")
  
  # # make predictions of classifiers into row-wise dataframe
  # pred <- lapply(y.pred, function(x){ x$boot }) %>% 
  #   rbindlist() %>% 
  #   as.data.frame()
  # 
  # # majority vote
  # pred.majoriry <- agg_classif(pred, method="majority")
  # err.majority <- pred.majoriry$error
  # # sens.spec.major <- apply(pred.majoriry$pred[, 3:8], 2, get_sens_spec, y.test)   # sensitivity and specificity
  # 
  # # oob error weighted vote
  # pred.oob <- agg_classif(pred, method="oob")
  # err.oob <- pred.oob$error
  # # sens.spec.oob <- apply(pred.oob$pred[, 3:8], 2, get_sens_spec, y.test)   # sensitivity and specificity
  
  return(list(err.majority = err.majority[B, ],
              err.oob = err.oob[B, ],
              err = list(majority = err.majority,
                         oob = err.oob),
              y.pred = y.pred))
}


### generate sparse functional data
# model A : different mean and variance 
# model B : different mean
# model C : different variance
sim.curve <- function(n, sparsity=5:10, model="A", prop=0.5) {
  n_A <- ceiling(n*prop)
  n_B <- n - n_A
  ## generate curves
  y <- factor(c(rep(0, n_A), rep(1, n_B)), levels=c(0, 1))
  data <- list()
  # random sparsify => generate curve
  for (i in 1:n) {
    # random sparsify
    num.obs <- sample(sparsity, 1)
    t <- sort(runif(num.obs, 0, 10))
    
    # mean function and variance
    if (model == "A") {
      if (i <= n_A) {
        mu <- t + sin(t)
        lambda <- c(4, 2, 1)
      } else {
        mu <- t + cos(t)
        lambda <- c(16, 8, 4)
      }
    } else if (model == "B") {
      if (i <= n_A) {
        mu <- t + sin(t)
        lambda <- c(4, 2, 1)
      } else {
        mu <- t + cos(t)
        lambda <- c(4, 2, 1)
      }
    } else if (model == "C") {
      if (i <= n_A) {
        mu <- t + sin(t)
        lambda <- c(4, 2, 1)
      } else {
        mu <- t + sin(t)
        lambda <- c(16, 8, 4)
      }
    } else {
      stop("Model is not existed. 3 models(A, B, C) can be applied.")
    }
    
    # eigenfunctions
    phi <- cbind(cos(pi*t/5) / sqrt(5),
                 sin(2*pi*t/5) / sqrt(5),
                 cos(3*pi*t/5) / sqrt(5))
    
    # FPC scores
    xi <- cbind(rnorm(1, 0, sqrt(lambda[1])),
                rnorm(1, 0, sqrt(lambda[2])),
                rnorm(1, 0, sqrt(lambda[3])))
    
    # measurement error
    eps <- rnorm(num.obs, 0, 0.5)
    
    # generate the curve with measurement error
    X <- mu + xi %*% t(phi) + eps
    
    data$id[[i]] <- rep(i, num.obs)
    data$y[[i]] <- rep(y[i], num.obs)
    data$Lt[[i]] <- t
    data$Ly[[i]] <- X
  }
  
  # return type of data
  id <- unlist( mapply(function(x, y){ rep(x, y) },
                       1:n,
                       sapply(data$Lt, length)) )
  time <- unlist(data$Lt)
  val <- unlist(data$Ly)
  y <- unlist( mapply(function(x, z){ rep(x, z) },
                      y,
                      sapply(data$Lt, length)) )
  data <- data.frame(id, time, val, y)

  return(data)
}
