### print function of "bagFPCA" object
print.bagFPCA <- function(bagFPCA.obj) {
  dtype <- bagFPCA.obj$dtype
  B <- bagFPCA.obj$B
  FPC.method <- bagFPCA.obj$FPC.method
  avg.K <- round(mean(bagFPCA.obj$K), 2)
  method <- bagFPCA.obj$method
  oob.error <- round(mean(bagFPCA.obj$OOB.error), 2)
  train.error <- round(mean(bagFPCA.obj$train.error), 2)
  
  cat( "The bootstrap aggregating classifier using", dtype, "FPC scores with", B, "iterations \n",
       "FPCA method :", FPC.method, "\n" ,
       "Avg selected K:", avg.K, "\n",
       "Classification method :", method, "\n",
       "OOB error rate :", oob.error, "\n",
       "Training error rate :", train.error )
}
