#########################################
### Real data analysis
###   Spinal bone mineral density data
#########################################
setwd("C:\\Users\\user\\Desktop\\KHS\\bagfpca\\Thesis")

library(tidyverse)
library(doParallel)   # parallel computing
library(fdapace)
library(e1071)
library(MASS)
library(data.table)
library(xtable)
source("R/bagFPCA.R")

# parallel computing setting
ncores <- detectCores() - 3
registerDoParallel(ncores)

packages <- c("fdapace","e1071","MASS")   # foreach에서 사용할 package 정의

### load data and preprocessing
data <- read.csv("data/spnbmd.csv", header=T)

ind <- data %>%
  # filter(ethnic == "Hispanic") %>%
  group_by(idnum) %>%
  summarise(n=n()) %>%
  filter(n >= 2) %>%
  dplyr::select(idnum)
data <- data[which(data$idnum %in% ind$idnum), ]

data %>% 
  mutate(sex = ifelse(sex == 'fem', "Female", "Male")) %>% 
  ggplot(aes(x=age, y=spnbmd, group=idnum, color=sex)) +
  geom_point(size = 1.5) +
  geom_line(size = 1) +
  theme_bw() +
  ylab("Spinal bone mineral density") +
  xlab("Age (years)") +
  theme(legend.position = c(0.15, 0.85),
        legend.title = element_blank(),
        legend.background = element_rect(color = "black", 
                                         # size = 1,
                                         linetype = "solid"),
        panel.border = element_rect(colour = "black", 
                                    # size = 1.5
                                    fill = NA)
        )
  # theme(legend.title = element_blank())
ggsave(file="img/spnbmd.eps") 


# gender classification
data <- data %>% 
  mutate(y = factor(ifelse(sex == "fem", 1, 0), levels=c(0, 1))) %>% 
  dplyr::select(y, idnum, age, spnbmd)
length(unique(data$idnum))   # 280

colnames(data) <- c("y","id","time","val")


### construct classification models
result <- list()
set.seed(1000)
seed <- sample(1:10000, 100)
for (simm in 1:100) {
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
                             err.oob = err.bag$err.oob,
                             # sensitivity and specificity
                             sens.single = err.single$sensitivitiy,
                             spec.single = err.single$specificity,
                             sens.major = err.bag$sens.major,
                             spec.major = err.bag$spec.major,
                             sens.oob = err.bag$sens.oob,
                             spec.oob = err.bag$spec.oob))
  colnames(res) <- c("Logit","SVM(Linear)","SVM(Gaussian)","LDA","QDA","NaiveBayes")
  
  result[[simm]] <- res
}

save(result, file="RData/real_data_1_modify.RData")

## 결과 정리
result <- result[!sapply(result, is.null)]

res <- sapply(1:3, function(i){
  paste(lapply(result, function(x){ x[i, ]*100 }) %>% 
          rbindlist %>% 
          colMeans %>% 
          round(2),
        " (",
        apply(lapply(result[!sapply(result, is.null)], 
                     function(x){ x[i, ]*100 }) %>% 
                rbindlist, 2, sd) %>% 
          round(2),
        ")",
        sep="")
}) %>% 
  t() %>% 
  as.data.frame
colnames(res) <- c("Logit","SVM(Linear)","SVM(Gaussian)","LDA","QDA","NaiveBayes")
rownames(res) <- c("Single","Majority vote","OOB weight")
res

xtable(res)


