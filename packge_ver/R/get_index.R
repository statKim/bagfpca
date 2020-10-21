### The functions obtained the curve's ID for bagged classifier with sparse FPCA
### - exclude curves out-of-range time points 

# return ID of test curves
get_test_index <- function(fpca.fit, X.test) {
  time.train <- fpca.fit$obsGrid   # time points of FPCA (train set)
  time.test <- unlist(X.test$Lt)   # time points of test set
  
  ### exclude test curves out-of-range of time points 
  id.test <- NULL
  
  # min(time of test) < min(time of train)
  if (min(time.test) < min(time.train)) {
    # index of train set(Not id)
    over.ind <- which(sapply(X.test$Lt, 
                             function(x){ sum( x < min(time.train) ) }) != 0)
    id.test <- unlist(X.test$Lid)[-over.ind]
  }
  
  # max(time of test) > max(time of train)
  if (max(time.test) > max(time.train)) {
    # index of train set(Not id)
    over.ind <- which(sapply(X.test$Lt, 
                             function(x){ sum( x > max(time.train) ) }) != 0)
    if (is.numeric(id.test)) {
      id.test <- intersect(id.test,
                           unlist(X.test$Lid)[-over.ind])
    } else {
      id.test <- unlist(X.test$Lid)[-over.ind]
    }
  }
  
  # Not out-of-range
  if (!is.numeric(id.test)) {
    id.test <- unlist(X.test$Lid)
  }

  return(id.test)
}


# return ID of Out-of-Bag curves 
get_oob_index <- function(X.train, boot.ind) {
  ### exclude OOB curves out-of-range of time points 
  id.oob <- NULL
  
  # min(time of test) < min(time of train)
  if (min(unlist(X.train$Lt[-unique(boot.ind)])) < min(unlist(X.train$Lt[boot.ind]))) {
    # index of train set(Not id)
    over.ind <- which(sapply(X.train$Lt[-unique(boot.ind)], 
                             function(x){ sum( x < min(unlist(X.train$Lt[boot.ind])) ) }) != 0)
    id.oob <- unlist(X.train$Lid)[-unique(boot.ind)][-over.ind]
  }
  
  # max(time of test) > max(time of train)
  if (max(unlist(X.train$Lt[-unique(boot.ind)])) > max(unlist(X.train$Lt[boot.ind]))) {
    # index of train set(Not id)
    over.ind <- which(sapply(X.train$Lt[-unique(boot.ind)], 
                             function(x){ sum( x > max(unlist(X.train$Lt[boot.ind])) ) }) != 0)
    if (is.numeric(id.oob)) {
      id.oob <- intersect(id.oob,
                          unlist(X.train$Lid)[-unique(boot.ind)][-over.ind])
    } else {
      id.oob <- unlist(X.train$Lid)[-unique(boot.ind)][-over.ind]
    }
  }
  
  # Not out-of-range
  if (!is.numeric(id.oob)) {
    id.oob <- unlist(X.train$Lid)[-unique(boot.ind)]
  }
  
  return(id.oob)
}

