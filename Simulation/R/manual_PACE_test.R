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



### eigen analysis
eig <- eigen(fpca.fit$smoothedCov)
positive_ind <- which(eig$values > 0)
lambda <- eig$values[positive_ind]
phi <- eig$vectors[, positive_ind]

### normalize eigenvalues and eigenvectors
grid_interval <- diff(work_grid)[1]
muWork <- 1:nrow(phi)
lambda <- lambda * grid_interval
phi <- apply(phi, 2, function(x) {
  x <- x / sqrt(sum(x^2 * grid_interval)) 
  if (sum(x * muWork) >= 0) {
    return(x)
  } else {
    return(-x)
  }
})

lambda[1:k]
fpca.fit$lambda

phi[, 2]
fpca.fit$phi[, 2]

