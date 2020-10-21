####################
### PACE Example
####################
setwd("Simulation")
library(fdapace)
library(tidyverse)

### load data and preprocessing
data <- read.csv("data/spnbmd.csv", header=T)
ind <- data %>%
  # filter(ethnic == "Hispanic") %>%
  group_by(idnum) %>%
  summarise(n=n()) %>%
  filter(n >= 2) %>%
  dplyr::select(idnum)
data <- data[which(data$idnum %in% ind$idnum), ]
data <- MakeFPCAInputs(IDs = data$idnum,
                       tVec = data$age,
                       yVec = data$spnbmd)

### fit FPCA
fpca.fit <- FPCA(data$Ly, 
                 data$Lt, 
                 optns=list(dataType="Sparse",
                            methodXi="CE",
                            # kernel="epan",
                            # methodSelectK="BIC",
                            FVEthreshold=0.99,
                            verbose=T))
fpca.fit$timings

### PACE
k <- fpca.fit$selectK
obs_grid <- fpca.fit$obsGrid
work_grid <- fpca.fit$workGrid
fitted_cov <- fpca.fit$fittedCov
sigma_2 <- fpca.fit$sigma2
mu <- fpca.fit$mu
phi <- fpca.fit$phi
lambda <- fpca.fit$lambda


j <- 2
Y <- data$Ly[[j]]
t <- data$Lt[[j]]


# source file load
for (fname in list.files("R/source")) {
  source( paste("R/source/", fname, sep="") )
}

## mean function
bw_mu <- unlist( GCVLwls1D1(yy=data$Ly, tt=data$Lt, kernel="gauss", npoly=1, nder=0, dataType="Sparse") )[1]
xin <- unlist(data$Lt)
yin <- unlist(data$Ly)[order(xin)]
xin <- sort(xin)
win <- rep(1, length(xin))
mu_obs <- Lwls1D(bw_mu, kernel_type="gauss", npoly=1, nder=0, xin=xin, yin=yin, xout=obs_grid, win=win)
# mu_obs <- ConvertSupport(obs_grid, obs_grid, mu=mu_obs)

## covariance function
rcov <- GetRawCov(data$Ly, data$Lt, obs_grid, mu_obs, dataType="Sparse", error=T)
bwCov <- GCVLwls2DV2(obs_grid, work_grid, kern="gauss", rcov=rcov, verbose=T, t=data$Lt)$h
smoothCov <- Lwls2D(bwCov, kern="gauss", 
                    xin=rcov$tPairs, yin=rcov$cxxn,
                    xout1=work_grid, xout2=work_grid)



## PACE test 1
mu_obs <- ConvertSupport(work_grid, obs_grid, mu=mu)
phi_obs <- ConvertSupport(work_grid, obs_grid, phi=phi)
Sigma_Y <- ConvertSupport(work_grid, obs_grid, Cov=fitted_cov) + diag(sigma_2, nrow=length(obs_grid))

t_ind <- which(obs_grid %in% t)
xi_1 <- diag(lambda) %*% t(phi_obs[t_ind, ]) %*% solve(Sigma_Y[t_ind, t_ind], 
                                                       matrix(Y - mu_obs[t_ind], ncol=1)) 

## PACE test 2
muVec <- approx(obs_grid, mu_obs, t)$y
phiMat <- matrix(apply(phi_obs, 2,
                       function(phivec){ approx(obs_grid, phivec, t)$y }),
                 nrow=length(t))
Sigma_Yi <- matrix(pracma::interp2(obs_grid, obs_grid, Sigma_Y,
                                   sapply(t, function(x){ rep(x, length(t)) }),
                                   rep(t, length(t))),
                   length(t), length(t))
xi_2 <- diag(lambda) %*% t(phiMat) %*% solve(Sigma_Yi, 
                                             matrix(Y - muVec, ncol=1)) 

data.frame(manual_1 = xi_1,
           manual_2 = xi_2,
           FPCA = fpca.fit$xiEst[j, ])

# xi variance
diag(lambda) - diag(lambda) %*% t(phi_obs[t_ind, ]) %*% t(diag(lambda) %*% t(phi_obs[t_ind, ]) %*% solve(Sigma_Y[t_ind, t_ind]))
diag(lambda) - diag(lambda) %*% t(phiMat) %*% t(diag(lambda) %*% t(phiMat) %*% solve(Sigma_Yi))
fpca.fit$xiVar[[j]]



### fitted curve
obs_grid <- fpca.fit$obsGrid
work_grid <- fpca.fit$workGrid
mu_est <- fpca.fit$mu
phi_est <- fpca.fit$phi
xi_est <- fpca.fit$xiEst
lambda_est <- fpca.fit$lambda
k <- fpca.fit$selectK
cov_est <- fpca.fit$fittedCov
sigma_2_est <- fpca.fit$sigma2

xi_manual <- pace(data, cov_est, mu_est, lambda_est, phi_est, sigma_2_est, k, work_grid, obs_grid)

par(mfrow=c(2, 2))
for (i in 1:4) {
  plot(data$Lt[[i]], data$Ly[[i]], xlab="", ylab="", xlim=range(obs_grid), ylim=range(unlist(data$Ly)))
  X_est <- mu_est + matrix(xi_est[i, ], nrow=1) %*% t(phi_est)
  lines(work_grid, X_est, col="red")
  
  # X_est <- mu_est + matrix(xi_manual[i, ], nrow=1) %*% t(phi_est)
  # lines(work_grid, X_est, col="green")

  X_est <- ret$mu + matrix(ret$xiEst[i, ], nrow=1) %*% t(ret$phi)
  lines(work_grid, X_est, col="green")
}













### eigen analysis
fit_eig <- function(cov, work_grid) {
  ### eigen analysis
  eig <- eigen(cov)
  positive_ind <- which(eig$values > 0)
  lambda <- eig$values[positive_ind]
  phi <- eig$vectors[, positive_ind]
  
  ### normalize eigenvalues and eigenvectors
  grid_interval <- diff(work_grid)[1]
  muWork <- 1:nrow(phi)
  lambda <- lambda * grid_interval
  phi <- apply(phi, 2, function(x) {
    x <- x / sqrt(trapzRcpp(work_grid, x^2))   # trapezoidal rule for numerical integration
    if (sum(x * muWork) >= 0) {
      return(x)
    } else {
      return(-x)
    }
  })
  
  return(list(lambda = lambda,
              phi = phi))
}

eig <- fit_eig(fpca.fit$smoothedCov, work_grid)

eig$lambda[1:k]
fpca.fit$lambda

eig$phi[, 2]
fpca.fit$phi[, 2]


### PACE
pace <- function(data, cov, mu, lambda, phi, sigma_2, k, work_grid, obs_grid) {
  N <- length(data$Ly)
  phi <- phi[, 1:k]
  lambda <- lambda[1:k]
  
  mu_obs <- ConvertSupport(work_grid, obs_grid, mu=mu)
  phi_obs <- ConvertSupport(work_grid, obs_grid, phi=phi)
  Sigma_Y <- ConvertSupport(work_grid, obs_grid, Cov=cov) + diag(sigma_2, nrow=length(obs_grid))
  
  xi <- matrix(0, N, k)
  for (n in 1:N) {
    Y <- data$Ly[[n]]
    t <- data$Lt[[n]]
    
    # t_ind <- which(obs_grid %in% t)
    # # numerical problem (https://stackoverflow.com/questions/9508518/why-are-these-numbers-not-equal)
    # if (length(t) != length(t_ind)) {
    #   not_in <- t[which(!(t %in% obs_grid[t_ind]))]
    #   for (i in 1:length(not_in)) {
    #     eq_ind <- sapply(obs_grid, function(x) {
    #       all.equal(x, not_in[i])
    #     })
    #     t_ind <- c(t_ind, which(eq_ind == TRUE))
    #   }
    #   t_ind <- sort(t_ind)
    # }
    # xi[n, ] <- diag(lambda) %*% t(phi_obs[t_ind, ]) %*% solve(Sigma_Y[t_ind, t_ind], 
    #                                                           matrix(Y - mu_obs[t_ind], ncol=1))
    
    muVec <- approx(obs_grid, mu_obs, t)$y
    phiMat <- matrix(apply(phi_obs, 2,
                           function(phivec){ approx(obs_grid, phivec, t)$y }),
                     nrow=length(t))
    Sigma_Yi <- matrix(pracma::interp2(obs_grid, obs_grid, Sigma_Y,
                                       sapply(t, function(x){ rep(x, length(t)) }),
                                       rep(t, length(t))),
                       length(t), length(t))
    xi[n, ] <- diag(lambda) %*% t(phiMat) %*% solve(Sigma_Yi, 
                                                    matrix(Y - muVec, ncol=1))
  }
  
  return(xi)
}










##### Simulation
sim.curve <- function(n) {
  t <- seq(0, 1, length.out=51)
  num.obs <- length(t)
  
  # mean function and variance
  mu <- t + cos(t)
  lambda <- c(16, 8, 4)
  
  # eigenfunctions
  phi <- cbind(cos(pi*t/5) / sqrt(5),
               sin(2*pi*t/5) / sqrt(5),
               cos(3*pi*t/5) / sqrt(5))
  
  X <- matrix(0, n, num.obs)
  scores <- matrix(0, n, 3)
  for (i in 1:n) {
    # FPC scores
    xi <- cbind(rnorm(1, 0, sqrt(lambda[1])),
                rnorm(1, 0, sqrt(lambda[2])),
                rnorm(1, 0, sqrt(lambda[3])))
    
    # measurement error
    eps <- rnorm(num.obs, 0, 0.5)
    
    # generate the curve with measurement error
    X[i, ] <- mu + xi %*% t(phi) + eps
    
    scores[i, ] <- xi
  }
  
  return(list(X = X,
              time = t,
              scores = scores,
              mu = mu,
              eigfunc = phi,
              eigval = lambda))  
}



set.seed(100)
data <- sim.curve(100)
t <- data$time
xi <- data$scores
phi <- data$eigfunc

# dense
X_dense <- MakeFPCAInputs(tVec = t,
                          yVec = data$X)
# fpca.fit <- FPCA(X_dense$Ly, 
#                  X_dense$Lt, 
#                  optns=list(dataType="Dense",
#                             # methodXi="CE",
#                             # kernel="epan",
#                             # methodSelectK="BIC",
#                             FVEthreshold=0.99,
#                             verbose=T))

# sparse
X_sparse <- Sparsify(data$X, t, sparsity=c(5, 10))
fpca.fit <- FPCA(X_sparse$Ly, 
                 X_sparse$Lt, 
                 optns=list(dataType="Sparse",
                            methodXi="CE",
                            # kernel="epan",
                            # methodSelectK="BIC",
                            FVEthreshold=0.99,
                            verbose=T))
fpca.fit$timings

obs_grid <- fpca.fit$obsGrid
work_grid <- fpca.fit$workGrid
mu_est <- fpca.fit$mu
phi_est <- fpca.fit$phi
xi_est <- fpca.fit$xiEst
lambda_est <- fpca.fit$lambda
k <- fpca.fit$selectK
cov_est <- fpca.fit$fittedCov
sigma_2_est <- fpca.fit$sigma2

# manual PC scores
xi_manual <- pace(X_sparse, cov_est, mu_est, lambda_est[1:k], phi_est[, 1:k], sigma_2_est, k, work_grid, obs_grid)
  
# true eigenfunctions and estimated eigenfunctions
par(mfrow=c(2, 2))
plot(t, phi[, 1], type="l", xlab="", ylab="", main="1st eigenfunction")
lines(work_grid, phi_est[, 1], col="red")
plot(t, phi[, 2], type="l", xlab="", ylab="", main="2nd eigenfunction")
lines(work_grid, phi_est[, 2], col="red")
plot(t, phi[, 3], type="l", xlab="", ylab="", main="3rd eigenfunction")
lines(work_grid, phi_est[, 3], col="red")


sum((unlist(xi) - unlist(xi_est))^2)


par(mfrow=c(2, 2))
for (i in 1:4) {
  plot(X_sparse$Lt[[i]], X_sparse$Ly[[i]], xlab="", ylab="", xlim=c(0, 1))
  X_est <- mu_est + matrix(xi_est[i, ], nrow=1) %*% t(phi_est)
  lines(work_grid, X_est, col="red")
  
  X_est <- mu_est + matrix(xi_manual[i, ], nrow=1) %*% t(phi_est)
  lines(work_grid, X_est, col="green")
  
  # lines(X_dense$Lt[[i]], X_dense$Ly[[i]], col="blue")
}


