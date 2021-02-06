library(tidyverse)
library(gridExtra)
library(fdapace)
source("bagFPCA.R")


### Simulation 1
p <- 0.5
simm <- 1

# randomly selected trajectories
num_curve <- 5
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
  
  # each trajectories for each groups  
  assign(paste("traj_", sim_model, sep=""),
         ggplot() +
           geom_line(data=data, aes(x=time, y=val, group=id, color=y)) +
           geom_point(data=data, aes(x=time, y=val, group=id, color=y), size=1.3) +
           xlab("") +
           ylab("") +
           ggtitle(paste("Model", sim_model)) +
           theme_bw() +
           theme(legend.position="none")
  )
  
  # mean curve for each groups  
  assign(paste("mean_", sim_model, sep=""),
         ggplot() +
           geom_line(data=mean_curve, aes(x=time, y=val, color=y), size=1) +
           xlab("") +
           ylab("") +
           ggtitle(paste("Model", sim_model)) +
           theme_bw() +
           theme(legend.position="none")
  )
}
grid.arrange(traj_A,
             traj_B,
             traj_C,
             nrow=1)
grid.arrange(mean_A,
             mean_B,
             mean_C,
             nrow=1)
p_curve <- arrangeGrob(traj_A,
                       traj_B,
                       traj_C,
                       nrow=1)
p_mean <- arrangeGrob(mean_A,
                      mean_B,
                      mean_C,
                      nrow=1)
ggsave(file="img/sim_1_curve.eps", p_curve, width=15, height=5)
ggsave(file="img/sim_1_mean.eps", p_mean, width=15, height=5)




### Simulation 2
simm <- 1
# randomly selected trajectories
num_curve <- 5
for (sim_model in c("A","B","C")) {
  set.seed(simm)
  data <- sim.curve.2(200, sparsity=10:20, model=sim_model, split.prop=0.5)
  data <- data$data
  data <- data.frame(id=unlist(data$id),
                     time=unlist(data$Lt),
                     val=unlist(data$Ly),
                     y=unlist(data$y))
  print(table(unique(data[, c("id","y")])$y))
  
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
  
  # each trajectories for each groups  
  assign(paste("traj_", sim_model, sep=""),
         ggplot() +
           geom_line(data=data, aes(x=time, y=val, group=id, color=y)) +
           geom_point(data=data, aes(x=time, y=val, group=id, color=y), size=1.3) +
           xlab("") +
           ylab("") +
           ggtitle(paste("Model", sim_model)) +
           theme_bw() +
           theme(legend.position="none")
  )
  
  # mean curve for each groups  
  assign(paste("mean_", sim_model, sep=""),
         ggplot() +
           geom_line(data=mean_curve, aes(x=time, y=val, color=y), size=1) +
           xlab("") +
           ylab("") +
           ggtitle(paste("Model", sim_model)) +
           theme_bw() +
           theme(legend.position="none")
  )
}
grid.arrange(traj_A,
             traj_B,
             traj_C,
             nrow=1)
grid.arrange(mean_A,
             mean_B,
             mean_C,
             nrow=1)
p_curve <- arrangeGrob(traj_A,
                       traj_B,
                       traj_C,
                       nrow=1)
p_mean <- arrangeGrob(mean_A,
                      mean_B,
                      mean_C,
                      nrow=1)
ggsave(file="img/sim_2_curve.eps", p_curve, width=15, height=5)
ggsave(file="img/sim_2_mean.eps", p_mean, width=15, height=5)
