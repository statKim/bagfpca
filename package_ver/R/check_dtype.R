# check input data type is "sparse" or "dense"
check_dtype <- function(X.train) {
  time.len <- sapply(X.train$Lt, length)   # obtain length of time points in each curves
  uniq.len <- unique(time.length)          # unique number of time points length
  if (length(uniq.len) == 1) {
    dtype <- "Dense"
  } else {
    dtype <- "Sparse"
  }
  
  return(dtype)
}

# check input data is appropriate
check_data <- function(X.train) {
  id.len <- sapply(X.train$Lid, length)   # obtain length of ID in each curves
  time.len <- sapply(X.train$Lt, length)   # obtain length of time points in each curves
  y.len <- sapply(X.train$Ly, length)   # obtain length of observations in each curves
  
  # check length of curves
  len.TF <- (length(id.len) == length(time.len)) & (length(time.len) == length(y.len))
  
  # check length of observations
  y.time.len.TF <- identical(time.len, y.len)
  
  # check all length is 1
  id.len.TF <- identical(id.len, rep(1, length(id.len)))
  
  return( len.TF & y.time.len.TF & id.len.TF )
}