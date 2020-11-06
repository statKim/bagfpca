### predict function of "fpca.bag.classifier" object(the bagged classifier using FPCA)
# - bagFPCA.obj : output object of "bagFPCA" function
# - X.test : input data for "FPCA" function("fdapace" package)
# - y.test : vector contained data labels
# - ncores : the number of cores to use parallel computing
predict.bagFPCA <- function(bagFPCA.obj, X.test, y.test, ncores=NULL) {
  # the number of bootstrap resamples
  B <- bagFPCA.obj$B
  
  # parallel computing setting
  if (is.null(ncores)) {
    ncores <- ceiling(detectCores() / 3)
  }
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  y.pred <- foreach(b=1:B, .packages=packages) %dopar% {
    fpca.fit <- bagFPCA.obj$fpca.fit[[b]]
    K <- bagFPCA.obj$K[b]
    model <- bagFPCA.obj$model[[b]]
    boot.ind <- bagFPCA.obj$boot.index[[b]]
    oob.error <-  bagFPCA.obj$OOB.error[b]
    
    # get predictable test set ID
    id.test <- get_test_index(fpca.fit, X.test)      # predictable test set ID
    ind <- which(unlist(X.test$Lid) %in% id.test)    # predictable test set index
    
    # test FPC scores
    fpc.score.test <- predict(fpca.fit, 
                              X.test$Ly[ind], 
                              X.test$Lt[ind], 
                              K = K, 
                              xiMethod = FPC.method)
    test.fpc <- data.frame(y = y.test[ind], 
                           x = fpc.score.test)
    colnames(test.fpc) <- c("y", paste("FPC", 1:K, sep=""))
   
    # predict for each model
    if (bagFPCA.obj$method == "svm") {
      pred <- predict(model, test.fpc, probability=probability)
    } else if (bagFPCA.obj$method == "logit") {
      pred <- predict(model, test.fpc, type="response")
      pred <- factor(ifelse(pred > 0.5, 1, 0), levels=c(0, 1))
    } else if (bagFPCA.obj$method == "knn") {
      pred <- knn(train = cbind(y=model$cl,
                                    model$train), test = test.fpc, cl = model$cl, k = model$k)
    } else if (bagFPCA.obj$method %in% c("lda","qda")) {
      pred <- predict(model, test.fpc)$class
    } else if (bagFPCA.obj$method == "naivebayes") {
      pred <- predict(model, test.fpc, type="class")
    }
    
    return( data.frame(id = id.test,
                       y = y.test[ind],
                       pred = pred,   # predicted value
                       pred.oob = as.numeric(pred) * (1 - oob.error)  # predicted value * OOB accuracy
                       ) )
  }
  
  # Ensemble predicted values
  res <- as.data.frame( rbindlist(y.pred) )     # predicted values of all bootstrap samples
  test <- data.frame(id = unlist(X.test$Lid),   # test set
                     y = y.test) 
  # majority voting
  pred.majority <- res %>% 
    group_by(id, y) %>% 
    summarise(pred = majority_vote(pred))
  pred.majority <- left_join(test, pred.majority, by = id)
  # err.majority <- 1 - mean(pred.majority$pred == pred.majority$y)
  
  # oob error
  oob.acc <- 1 - sum(bagFPCA.obj$OOB.error)
  pred.oob <- res %>% 
    group_by(id, y) %>% 
    summarise(pred = factor(ifelse(sum(pred.oob)/oob.acc > 1.5, 1, 0), 
                            levels = c(0, 1)))
  pred.oob <- left_join(test, pred.oob, by = id)
  # err.oob <- 1 - mean(pred.oob$pred == pred.oob$y)
  
  stopCluster(cl)   # End the parallel computing
  
  return( data.frame(majority = pred.majority$pred,
                     oob = pred.oob$pred) )
}
