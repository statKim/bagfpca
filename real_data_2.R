#########################################
### Real data analysis
###   Berkely growth data
#########################################

library(tidyverse)
library(doParallel)   # parallel computing
library(fdapace)
library(e1071)
library(MASS)
library(data.table)
library(reshape2)
library(fda)
library(xtable)
source("bagFPCA.R")

# parallel computing setting
ncores <- detectCores() - 3
registerDoParallel(ncores)

packages <- c("fdapace","e1071","MASS")   # foreach에서 사용할 package 정의


### load data and preprocessing
data <- rbind(t(growth$hgtm),
              t(growth$hgtf))
y <- factor(c(rep(0, ncol(growth$hgtm)), 
              rep(1, ncol(growth$hgtf))), levels = c(0, 1))
time <- growth$age
N <- length(y)

# graph of data
cbind(as.data.frame(data),
      y = ifelse(y == 1, "Female", "Male"),
      id = 1:nrow(data)) %>% 
  reshape2::melt(id.vars = c("id", "y"), 
                 variable.name = "time",
                 value.name = "height") %>% 
  mutate(time = as.numeric(as.character(time))) %>% 
  ggplot(aes(x=time, y=height, group=id, color=y)) +
  geom_point(size = 1.5) +
  geom_line(size = 1) +
  xlab("Time") +
  ylab("Height") +
  theme_bw() +
  theme(legend.position = c(0.15, 0.85),
        legend.title = element_blank(),
        legend.background = element_rect(color = "black", 
                                         # size = 1,
                                         linetype = "solid"),
        panel.border = element_rect(colour = "black", 
                                    # size = 1.5
                                    fill = NA)
  ) +
  scale_x_continuous(breaks = seq(2, 18, 2))

# random sparsify
set.seed(1000)
data <- Sparsify(data, time, sparsity = c(12, 15))
id <- unlist( mapply(function(x, y){ rep(x, y) },
                     1:N,
                     sapply(data$Lt, length)) )
time <- unlist(data$Lt)
val <- unlist(data$Ly)
y <- unlist( mapply(function(x, z){ rep(x, z) },
                    y,
                    sapply(data$Lt, length)) )
data <- data.frame(id, time, val, y)


### construct classification models
result <- list()
num.sim <- 0
set.seed(1000)
seed <- sample(1:10000, 1000)
for (simm in 1:500) {
  print( paste(simm, ":", seed[simm]) )
  set.seed(seed[simm])
  
  start.time <- Sys.time()
  
  ### train, test split
  train_test <- train_test_split(data, train.prop = 2/3)
  X.train <- train_test$X.train
  X.test <- train_test$X.test
  y.train <- train_test$y.train
  y.test <- train_test$y.test
  
  ### Single classifier
  err.single <- get_single_err(X.train, X.test, y.train, y.test)
  
  ## Bagging
  # Bootstrap aggregating
  err.bag <- get_bag_err(X.train, X.test, y.train, y.test, B = 100, packages = packages, hyper.para = err.single$hyper.para)
  
  end.time <- Sys.time()
  print(end.time - start.time)
  
  # save result
  res <- as.data.frame(rbind(err.single = err.single$err.single,
                             err.majority = err.bag$err.majority,
                             err.oob = err.bag$err.oob))
  colnames(res) <- c("Logit","SVM(Linear)","SVM(Gaussian)","LDA","QDA","NaiveBayes")
  
  result[[simm]] <- list(error = res,
                         single.K = err.single$K,   # number of PCs
                         bag.K = sapply(err.bag$y.pred, function(x){ x$K }))   # number of PCs
  
  # save RData per 10 times
  num.sim <- num.sim + 1
  if (num.sim %% 10 == 0) {
    save(result, file="RData/result_growth.RData")
  }
}

save(result, file="RData/result_growth.RData")


# avg error rate and standard error
res <- sapply(1:3, function(i){
  paste(lapply(result[!sapply(result, is.null)], function(x){ x$error[i, ]*100 }) %>% 
          rbindlist %>% 
          colMeans %>% 
          round(2),
        " (",
        apply(lapply(result[!sapply(result, is.null)], 
                     function(x){ x$error[i, ]*100 }) %>% 
                rbindlist, 2, sd) %>% 
          round(2),
        ")",
        sep="")
}) %>% 
  t() %>% 
  as.data.frame
colnames(res) <- c("Logit","SVM(Linear)","SVM(Gaussian)","LDA","QDA","NaiveBayes")
rownames(res) <- c("Single","Majority vote","OOB weight")
print(res)


# number of selected PCs
single.K <- sapply(result[!sapply(result, is.null)],
                   function(x){ x$single.K })
bag.K <- lapply(result[!sapply(result, is.null)],
                function(x){ x$bag.K })
range(single.K)
mean(single.K)
range(unlist(bag.K))
mean(unlist(bag.K))

