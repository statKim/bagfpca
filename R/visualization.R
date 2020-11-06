library(tidyverse)
library(gridExtra)
library(fdapace)
######################################
### Trajactories of simulated data
######################################

### Simulation 1
p <- 0.5
simm <- 1

# all trajactories
for (sim_model in c("A","B","C")) {
  set.seed(simm)
  data <- sim.curve(200, sparsity=5:10, model=sim_model, prop=p)
  
  # mean curve for each groups
  for (i in unique(data$y)) {
    data2 <- data[which(data$y == i), ]
    data2 <- MakeFPCAInputs(IDs=data2$id,
                            tVec=data2$time,
                            yVec=data2$val)
    mean_y <- GetMeanCurve(data2$Ly, data2$Lt, 
                           optns=list(dataType="Sparse"))
    if (i == unique(data$y)[1]) {
      mean_curve <- data.frame(time=mean_y$workGrid,
                               val=mean_y$mu,
                               y=i)
    } else {
      mean_curve <- rbind(mean_curve,
                          data.frame(time=mean_y$workGrid,
                                     val=mean_y$mu,
                                     y=i))
    }
  }
  
  # each trajectories and mean curve for each groups  
  assign(paste("fig_", sim_model, sep=""),
         ggplot() +
           geom_line(data=data, aes(x=time, y=val, group=id, color=y), linetype="dashed") +
           geom_point(data=data, aes(x=time, y=val, group=id, color=y), size=1.3) +
           geom_line(data=mean_curve, aes(x=time, y=val, color=y), size=1.3) +
           xlab("") +
           ylab("") +
           ggtitle(paste("Model", sim_model)) +
           theme_bw() +
           theme(legend.position="none")
  )
}
grid.arrange(fig_A,
             fig_B,
             fig_C,
             nrow=1)
g <- arrangeGrob(fig_A,
                 fig_B,
                 fig_C,
                 nrow=1)
ggsave(file="img/sim_1.eps", g, width=15, height=5)


# randomly selected trajectories
num_curve <- 10
for (sim_model in c("A","B","C")) {
  set.seed(simm)
  data <- sim.curve(200, sparsity=5:10, model=sim_model, prop=p)
  
  # select curves for each groups
  ind <- c(sample(data %>% 
                    filter(y == 0) %>% 
                    dplyr::select(id) %>% 
                    unique() %>% 
                    unlist(),
                  num_curve),
           sample(data %>% 
                    filter(y == 1) %>% 
                    dplyr::select(id) %>% 
                    unique() %>% 
                    unlist(),
                  num_curve))
  
  data <- data %>% 
    filter(id %in% ind)
  
  # mean curve for each groups
  for (i in unique(data$y)) {
    data2 <- data[which(data$y == i), ]
    data2 <- MakeFPCAInputs(IDs=data2$id,
                            tVec=data2$time,
                            yVec=data2$val)
    mean_y <- GetMeanCurve(data2$Ly, data2$Lt, 
                           optns=list(dataType="Sparse"))
    if (i == unique(data$y)[1]) {
      mean_curve <- data.frame(time=mean_y$workGrid,
                               val=mean_y$mu,
                               y=i)
    } else {
      mean_curve <- rbind(mean_curve,
                          data.frame(time=mean_y$workGrid,
                                     val=mean_y$mu,
                                     y=i))
    }
  }
  
  # each trajectories and mean curve for each groups  
  assign(paste("fig_", sim_model, sep=""),
         ggplot() +
           geom_line(data=data, aes(x=time, y=val, group=id, color=y), linetype="dashed") +
           geom_point(data=data, aes(x=time, y=val, group=id, color=y), size=1.3) +
           geom_line(data=mean_curve, aes(x=time, y=val, color=y), size=1.3) +
           xlab("") +
           ylab("") +
           ggtitle(paste("Model", sim_model)) +
           theme_bw() +
           theme(legend.position="none")
  )
}
grid.arrange(fig_A,
             fig_B,
             fig_C,
             nrow=1)
g <- arrangeGrob(fig_A,
                 fig_B,
                 fig_C,
                 nrow=1)
ggsave(file="img/sim_1-1.eps", g, width=15, height=5)


# # randomly selected trajectories for each groups
# num_curve <- 10
# for (sim_model in c("A","B","C")) {
#   set.seed(simm)
#   data <- sim.curve(200, sparsity=5:10, model=sim_model, prop=p)
#   
#   # select curves for each groups
#   ind <- c(sample(data %>% 
#                     filter(y == 0) %>% 
#                     dplyr::select(id) %>% 
#                     unique() %>% 
#                     unlist(),
#                   num_curve),
#            sample(data %>% 
#                     filter(y == 1) %>% 
#                     dplyr::select(id) %>% 
#                     unique() %>% 
#                     unlist(),
#                   num_curve))
#   
#   data <- data %>% 
#     filter(id %in% ind)
#   
#   for (g in 0:1) {
#     assign(paste("fig_", sim_model, "_", g, sep=""),
#            ggplot(data[data$y == g, ], aes(x=time, y=val, group=id)) +
#              geom_line(linetype=g+1) +
#              xlab("") +
#              ylab("") +
#              ggtitle(paste("Model", sim_model)) +
#              theme_bw()
#     )
#   }
# }
# grid.arrange(fig_A_0, fig_A_1,
#              fig_B_0, fig_B_1,
#              fig_C_0, fig_C_1,
#              nrow=3)
# g <- arrangeGrob(fig_A_0, fig_A_1,
#                  fig_B_0, fig_B_1,
#                  fig_C_0, fig_C_1,
#                  nrow=3)
# ggsave(file="img/sim_1-2.eps", g, width=10, height=10)



### Simulation 2
simm <- 1
# all trajactories
for (sim_model in c("A","B","C")) {
  set.seed(simm)
  data <- sim.curve.2(200, sparsity=10:20, model=sim_model, split.prop=0.5)
  data <- data$data
  data <- data.frame(id=unlist(data$id),
                     time=unlist(data$Lt),
                     val=unlist(data$Ly),
                     y=unlist(data$y))
  
  # mean curve for each groups
  for (i in unique(data$y)) {
    data2 <- data[which(data$y == i), ]
    data2 <- MakeFPCAInputs(IDs=data2$id,
                            tVec=data2$time,
                            yVec=data2$val)
    mean_y <- GetMeanCurve(data2$Ly, data2$Lt, 
                           optns=list(dataType="Sparse"))
    if (i == unique(data$y)[1]) {
      mean_curve <- data.frame(time=mean_y$workGrid,
                               val=mean_y$mu,
                               y=i)
    } else {
      mean_curve <- rbind(mean_curve,
                          data.frame(time=mean_y$workGrid,
                                     val=mean_y$mu,
                                     y=i))
    }
  }
  
  # each trajectories and mean curve for each groups  
  assign(paste("fig_", sim_model, sep=""),
         ggplot() +
           geom_line(data=data, aes(x=time, y=val, group=id, color=y), linetype="dashed") +
           geom_point(data=data, aes(x=time, y=val, group=id, color=y), size=1.3) +
           geom_line(data=mean_curve, aes(x=time, y=val, color=y), size=1.3) +
           xlab("") +
           ylab("") +
           ggtitle(paste("Model", sim_model)) +
           theme_bw() +
           theme(legend.position="none")
  )
}
grid.arrange(fig_A,
             fig_B,
             fig_C,
             nrow=1)
g <- arrangeGrob(fig_A,
                 fig_B,
                 fig_C,
                 nrow=1)
ggsave(file="img/sim_2.eps", g, width=15, height=5)


# randomly selected trajectories
num_curve <- 10
for (sim_model in c("A","B","C")) {
  set.seed(simm)
  data <- sim.curve.2(200, sparsity=10:20, model=sim_model, split.prop=0.5)
  data <- data$data
  data <- data.frame(id=unlist(data$id),
                     time=unlist(data$Lt),
                     val=unlist(data$Ly),
                     y=unlist(data$y))
  
  # select curves for each groups
  ind <- c(sample(data %>% 
                    filter(y == 0) %>% 
                    dplyr::select(id) %>% 
                    unique() %>% 
                    unlist(),
                  num_curve),
           sample(data %>% 
                    filter(y == 1) %>% 
                    dplyr::select(id) %>% 
                    unique() %>% 
                    unlist(),
                  num_curve))
  
  data <- data %>% 
    filter(id %in% ind)
  
  # mean curve for each groups
  for (i in unique(data$y)) {
    data2 <- data[which(data$y == i), ]
    data2 <- MakeFPCAInputs(IDs=data2$id,
                            tVec=data2$time,
                            yVec=data2$val)
    mean_y <- GetMeanCurve(data2$Ly, data2$Lt, 
                           optns=list(dataType="Sparse"))
    if (i == unique(data$y)[1]) {
      mean_curve <- data.frame(time=mean_y$workGrid,
                               val=mean_y$mu,
                               y=i)
    } else {
      mean_curve <- rbind(mean_curve,
                          data.frame(time=mean_y$workGrid,
                                     val=mean_y$mu,
                                     y=i))
    }
  }
  
  # each trajectories and mean curve for each groups  
  assign(paste("fig_", sim_model, sep=""),
         ggplot() +
           geom_line(data=data, aes(x=time, y=val, group=id, color=y), linetype="dashed") +
           geom_point(data=data, aes(x=time, y=val, group=id, color=y), size=1.3) +
           geom_line(data=mean_curve, aes(x=time, y=val, color=y), size=1.3) +
           xlab("") +
           ylab("") +
           ggtitle(paste("Model", sim_model)) +
           theme_bw() +
           theme(legend.position="none")
           # theme(legend.position = c(0.15, 0.85),
           #       legend.title = element_blank(),
           #       legend.background = element_rect(color = "black", 
           #                                        # size = 1,
           #                                        linetype = "solid"),
           #       panel.border = element_rect(colour = "black", 
           #                                   # size = 1.5
           #                                   fill = NA))
  )
}
grid.arrange(fig_A,
             fig_B,
             fig_C,
             nrow=1)
g <- arrangeGrob(fig_A,
                 fig_B,
                 fig_C,
                 nrow=1)
ggsave(file="img/sim_2-1.eps", g, width=15, height=5)


# # randomly selected trajectories for each groups
# num_curve <- 10
# for (sim_model in c("A","B","C")) {
#   set.seed(simm)
#   n <- 200   # number of observations
#   j <- 1:50  # number of basis
#   
#   # parameters when generate class label
#   b <- matrix(ifelse(j <= 2, 1, (j-2)^(-3)), ncol=1)
#   
#   ## generate curves
#   data <- list()
#   for (i in 1:n) {
#     # random sparsify
#     num.obs <- sample(10:20, 1)
#     t <- sort(runif(num.obs, 0, 10))
#     
#     # 101~200 are test set
#     if (i > 100) {
#       range.train <- range(unlist(data$Lt[1:100]))
#       while( (max(t) > range.train[2]) | (min(t) < range.train[1]) ) {
#         t <- sort(runif(num.obs, 0, 10))
#       }
#     }
#     
#     # eigenfunctions
#     phi <- sapply(j, function(x){
#       if (x %% 2 == 0) {
#         sqrt(1/5)*sin(pi*t*x/5)
#       } else {
#         sqrt(1/5)*cos(pi*t*x/5)
#       }
#     })
#     
#     # generate PC scores
#     xi <- sapply(j, function(x){ rnorm(1, 0, sqrt( x^(-1.5) )) })
#     
#     # parameters when generate class label
#     beta.1 <- phi %*% b
#     beta.2 <- matrix(sqrt(3/10)*(t/5-1), ncol=1)
#     
#     # measurement error
#     eps <- rnorm(num.obs, 0, sqrt(0.1))
#     
#     # generate the curve
#     X <- xi %*% t(phi) + eps
#     
#     # generate class label
#     eps <- rnorm(1, 0, sqrt(0.1))   # model error
#     if (sim_model == "A") {
#       fx <- exp((X %*% beta.1)/2)-1   # model 2
#     } else if (sim_model == "B") {
#       fx <- atan(pi*(X %*% beta.1)) + exp((X %*% beta.2)/3) - 1    # model 4
#     } else {
#       fx <- atan(pi*(X %*% beta.1)/4)
#     }
#     # fx <- sin(pi*(X %*% beta.1)/4)  # model 1
#     # fx <- sin(pi*(X %*% beta.1)/3) + exp((X %*% beta.2)/3) - 1   # model 3
#     # fx <- atan(pi*(X %*% beta.1)) + exp((X %*% beta.2)/3) - 1    # model 4
#     y <- factor(ifelse(fx + eps < 0, 0, 1), levels=c(0, 1))
#     
#     data$id[[i]] <- rep(i, num.obs)
#     data$y[[i]] <- rep(y, num.obs)
#     data$Lt[[i]] <- t
#     data$Ly[[i]] <- X
#   }
#   
#   
#   # select curves for each groups
#   ind <- sapply(data$y, function(x){ x[1] })
#   ind <- c(sample(which(ind == 0), num_curve),
#            sample(which(ind == 1), num_curve))
#   
#   data <- lapply(data, function(x){ x[ind] })
#   
#   # plot the generated curves
#   for (g in 0:1) {
#     assign(paste("fig_", sim_model, "_", g, sep=""),
#            sapply(data, unlist) %>%
#              data.frame() %>%
#              mutate(y = ifelse(unlist(data$y) == 0, "0", "1")) %>%
#              filter(y == g) %>% 
#              ggplot(aes(x=Lt, y=Ly, group=id)) +
#              geom_line(linetype=g+1) +
#              xlab("") +
#              ylab("") +
#              ggtitle(paste("Model", sim_model)) +
#              theme_bw() +
#              theme(legend.title = element_blank(),
#                    legend.position = c(0.8, 0.1))
#     )
#   }
# }
# grid.arrange(fig_A_0, fig_A_1,
#              fig_B_0, fig_B_1,
#              fig_C_0, fig_C_1,
#              nrow=3)
# g <- arrangeGrob(fig_A_0, fig_A_1,
#                  fig_B_0, fig_B_1,
#                  fig_C_0, fig_C_1,
#                  nrow=3)
# ggsave(file="img/sim_2-2.eps", g, width=10, height=10)




### box plot
library(tidyverse)
library(gridExtra)

model_name <- c("Logit","SVM(Linear)","SVM(Gaussian)","LDA","QDA","NaiveBayes")

load("Simulation/RData/real_data_2_modify.RData")

for (i in 1:6) {
  res <- lapply(result, function(x){ t(x)[i, 1:3]*100 %>% matrix(nrow=1, ncol=3) %>% as.data.frame }) %>% 
    rbindlist %>% 
    rename(Single=V1,
           Majority=V2,
           OOB_weight=V3) %>% 
    gather(key="Model", value="Error")
  
  assign(paste("p_", i, sep=""),
         ggplot(res, aes(x=Model, y=Error)) +
           geom_boxplot() +
           geom_hline(yintercept = median(res$Error[res$Model == "OOB_weight"]), color="red", size=0.7) +
           scale_x_discrete(limits=c("Single", "Majority", "OOB_weight")) +
           xlab("") +
           # ylab("Classification error rate") +
           ylab("") +
           ggtitle(model_name[i]) +
           theme_bw())
}

grid.arrange(p_1,
             p_2,
             p_3,
             p_4,
             p_5,
             p_6,
             nrow=2)

