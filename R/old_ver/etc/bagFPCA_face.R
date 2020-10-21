##################################
### bagging classifier with FPCA
##################################

require(tidyverse)    # dplyr, ggplot2, etc
require(doParallel)   # parallel computing
require(face)         # functional PCA
require(e1071)        # SVM, Naive bayes
require(MASS)         # LDA, QDA
require(data.table)   # list rbind

### majority voting functon
majority_vote <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


### train, test split
##  - input data form : cbind(id, time, val, y)
train_test_split <- function(data, train.prop=2/3, num.grid=51) {
  id <- unique(data$id)
  N <- length(unique(data$id))
  N.train <- ceiling(N * train.prop)
  id.train <- sample(id, N.train)
  id.test <- id[-which(id %in% id.train)]
  
  # transform data for FACE
  train <- data[which(data$id %in% id.train), ]
  test <- data[which(data$id %in% id.test), ]
  
  X.train <- data.frame(y = train$val,
                        argvals = train$time,
                        subj = train$id)
  X.test <- data.frame(y = test$val,
                       argvals = test$time,
                       subj = test$id)
  y.train <- unique(train[, c("id", "y")])[, "y"]
  y.test <- unique(test[, c("id", "y")])[, "y"]
  
  # grids of time points
  t.grid <- seq(min(data$time), max(data$time), length.out=num.grid)
  
  return(list(X.train = X.train,
              X.test = X.test,
              y.train = y.train,
              y.test = y.test,
              t.grid = t.grid))
}


### Fit FPCA
fit_FPCA <- function(X, t.grid) {
  ### estimate covariance using face algorithm
  fpca.fit <- face.sparse(X,
                          newdata=X,
                          argvals.new=t.grid,
                          calculate.scores=T)
  
  return(fpca.fit)
}


### make data using FPC scores
get_FPCscore <- function(X.train, X.test, y.train, y.test, t.grid, boot.ind=NULL) {
  if (is.null(boot.ind)) {
    ## Not bootstrap (Don't have boot.ind)
    # fit FPCA for training set
    fpca.fit <- fit_FPCA(X.train,
                         t.grid)
    fpc.score <- fpca.fit$scores$scores
    k <- ncol(fpc.score)   # optimal number of PCs
    
    # train FPC scores
    train.fpc <- data.frame(y = y.train, 
                            x = fpc.score)
    colnames(train.fpc) <- c("y", paste("FPC", 1:k, sep=""))
    
    # test FPC scores
    fpc.score.test <- predict(fpca.fit, X.test)
    test.fpc <- data.frame(y = y.test, 
                           x = fpc.score.test$scores$scores)
    colnames(test.fpc) <- c("y", paste("FPC", 1:k, sep=""))
    
    return(list(train.fpc = train.fpc,
                test.fpc = test.fpc,
                fpca.fit = fpca.fit,
                k = k))
  } else {
    ## Bootstrap (Have boot.ind)
    id.train <- unique(X.train$subj)
    id.test <- unique(X.test$subj)
    id.boot <- id.train[boot.ind]
    id.oob <- id.train[-unique(boot.ind)]
    
    # get bootstrap sample and OOB sample
    X.boot <- data.frame(subj=1:length(id.boot),
                         id=id.boot) %>% 
      left_join(X.train, by=c("id"="subj")) %>% 
      dplyr::select(y, argvals, subj)
    X.oob <- data.frame(subj=id.oob) %>% 
      left_join(X.train, by="subj") %>% 
      dplyr::select(y, argvals, subj)
    
    y.boot <- y.train[boot.ind]
    y.oob <- y.train[-boot.ind]
    
    # fit FPCA for bootstrapped data
    fpca.fit <- fit_FPCA(X.boot,
                         t.grid)
    fpc.score <- fpca.fit$scores$scores
    k <- ncol(fpc.score)   # optimal number of PCs
    
    # train FPC scores
    train.fpc <- data.frame(y = y.boot, 
                            x = fpc.score)
    colnames(train.fpc) <- c("y", paste("FPC", 1:k, sep=""))
    
    # test FPC scores
    fpc.score.test <- predict(fpca.fit, X.test)
    test.fpc <- data.frame(y = y.test, 
                           x = fpc.score.test$scores$scores)
    colnames(test.fpc) <- c("y", paste("FPC", 1:k, sep=""))
    
    # OOB FPC scores
    fpc.score.oob <- predict(fpca.fit, X.oob)
    oob.fpc <- data.frame(y = y.oob, 
                          x = fpc.score.oob$scores$scores)
    colnames(oob.fpc) <- c("y", paste("FPC", 1:k, sep=""))
    
    return(list(train.fpc = train.fpc,
                test.fpc = test.fpc,
                oob.fpc = oob.fpc,
                id.test = id.test,
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
get_single_err <- function(X.train, X.test, y.train, y.test, t.grid) {
  # obtain FPC scores of train and test data
  fit.fpc <- get_FPCscore(X.train, X.test, y.train, y.test, t.grid)
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
get_bag_err <- function(X.train, X.test, y.train, y.test, t.grid, B = 100, packages = c("face","e1071","MASS"), hyper.para) {
  # list of manual functions
  ftns <- c("fit_FPCA","get_FPCscore","fit_classif","predict_classif")
  
  N <- length(unique(X.train$subj))
  boot.mat <- matrix(sample(1:N, N*B, replace=T), N, B)   # index of bootstrap samples
  
  # make predictions from bootstrap samples
  y.pred <- foreach(b=1:B, .packages=packages, .export=ftns) %dopar% {
    boot.ind <- boot.mat[, b]
    
    # obtain FPC scores of train, test, OOB data
    fit.fpc <- get_FPCscore(X.train, X.test, y.train, y.test, t.grid, boot.ind=boot.ind)
    train.fpc <- fit.fpc$train.fpc
    test.fpc <- fit.fpc$test.fpc
    oob.fpc <- fit.fpc$oob.fpc
    id.test <- fit.fpc$id.test
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
                                   y = y.test,
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
