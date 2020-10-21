library(tidyverse)
library(face)

### load data and preprocessing
data <- read.csv("Simulation/data/spnbmd.csv", header=T)

ind <- data %>%
  # filter(ethnic == "Hispanic") %>%
  group_by(idnum) %>%
  summarise(n=n()) %>%
  filter(n >= 2) %>%
  dplyr::select(idnum)
data <- data[which(data$idnum %in% ind$idnum), ]


# gender classification
data <- data %>% 
  mutate(y = factor(ifelse(sex == "fem", 1, 0), levels=c(0, 1))) %>% 
  dplyr::select(y, idnum, age, spnbmd)
length(unique(data$idnum))   # 280

colnames(data) <- c("y","id","time","val")

# transform data for face
X <- data.frame(y = data$val,
                argvals = data$time,
                subj = data$id)

# grids of time points
t_grid <- seq(min(X$argvals), max(X$argvals), length.out=51)

### estimate covariance using face algorithm
fit_face <- face.sparse(X, argvals.new=t_grid)

# fitted grid of time point (same with t_grid)
tnew <- fit_face$argvals.new

## scatter plots
Xlab <- "Age (years)"
Ylab <- "Spinal bone mineral density"
par(mfrow=c(1, 1))
id <- X$subj
uid <- unique(id)
plot(X$argvals, X$y,
     type="n", xlab=Xlab, ylab=Ylab)

for (i in 1:length(uid)) {
  seq <- which(id == uid[i])
  lines(X$argvals[seq], X$y[seq], lty=1, col="gray")
  #points(X$argvals[seq],X$y[seq],col=1,lty=1,pch=1)
}

Sample <- sample(uid, 5)
for (i in Sample) {
  seq <- which(id == i)
  lines(X$argvals[seq], X$y[seq], lty=1, col="black", lwd=1.5)
}
lines(tnew, fit_face$mu.new, lwd=2, lty=2, col="red")

## plots of variance/correlation functions
Cov <- fit_face$Chat.new
Cov_diag <- diag(Cov)
Cor <- fit_face$Cor.new

par(mfrow=c(1, 2))
plot(tnew, Cov_diag, 
     type="l", xlab=Xlab, ylab="", main="Variance function", lwd=2)

library(fields)
image.plot(tnew, tnew,Cor,
           xlab=Xlab, ylab=Xlab,
           main="Correlation function",
           # cex.axis=1.25, cex.lab=1.25,cex.main=1.25,
           axis.args=list(at = c(0,0.2,0.4,0.6,0.8,1.0)),
           legend.shrink=0.75, legend.line=-1.5)

## prediction of several subjects
par(mfrow=c(2, 2))
Sample <- sample(uid, 4)
for (i in Sample) {
  sel <- which(id == i)
  dati <- X[sel, ]
  
  seq <- t_grid
  k <- length(seq)
  dati_pred <- data.frame(y = rep(NA, nrow(dati) + k ),
                          argvals = c(rep(NA, nrow(dati)), seq),
                          subj = rep(dati$subj[1], nrow(dati) + k))
  
  dati_pred[1:nrow(dati), ] <- dati
  yhat2 <- predict(fit_face, dati_pred)
  
  data3 <- dati
  Ylim <- range(c(data3$y, yhat2$y.pred))
  
  plot(data3$argvals, data3$y,
       xlim=range(X$argvals), ylim=range(X$y),
       xlab=Xlab, ylab=Ylab, main=paste("Trajectory ", i, sep=""), pch=1)
  
  Ord <- nrow(dati) + 1:k
  lines(dati_pred$argvals[Ord], yhat2$y.pred[Ord], col="red", lwd=2)
  lines(dati_pred$argvals[Ord], yhat2$y.pred[Ord] - 1.96*yhat2$se.pred[Ord], col="red", lwd=1, lty=2)
  lines(dati_pred$argvals[Ord], yhat2$y.pred[Ord] + 1.96*yhat2$se.pred[Ord], col="red", lwd=1, lty=2)
  
  lines(tnew,fit_face$mu.new, lty=3, col="black", lwd=2)
  legend("topleft", c("mean","prediction"), lty=c(3,1), col=1:2, lwd=2, bty="n")
}


eig <- eigen(Cov)
cumsum(eig$values[1:5]) / sum(eig$values)
fit_face$eigenvalues
fit_face$eigenfunctions
