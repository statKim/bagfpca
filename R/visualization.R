# box plot
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

