library(UpSetR)
library(tidyverse)

######
# plot upset function
plotUpset2 <- function(data, prefix){ 
  data <- read.csv(data, sep = "\t")
  data <- data[, -1]
  data.melt <- reshape2::melt(data, id.vars="None")
  as_tibble(data.melt) %>% group_by(None, variable) %>% summarise(total=sum(value)) %>% 
    pivot_wider(names_from = variable, values_from = total)
  
  as_tibble(data.melt) %>% group_by(None, variable) %>% summarise(total=sum(value)) %>% 
    mutate(t2 = if_else(total>0,1,0)) %>% select(!total)  %>% 
    pivot_wider(names_from = variable, values_from = t2) -> data.melt2
  
  data.melt2 <- as.data.frame(data.melt2)
  
  rownames(data.melt2) <- data.melt2$None
  pdf(paste0(prefix,"_piRNA_intersetion_m2.pdf"), width = 7, height = 6)
  print(upset(data.melt2[, 2:13], nsets = 12, nintersects = 20, mb.ratio = c(0.6, 0.4), 
              order.by = c("freq", "degree"), decreasing = c(TRUE, FALSE),
              matrix.color = "gray23", main.bar.color = "gray23",  
              mainbar.y.label = "Intersection of piRNA species",  sets.x.label = "Total piRNA species",
              scale.intersections = "identity", scale.sets = "identity"
  ))
  dev.off()
}

# seed and length
plotUpset2("cluster1_piRNA_seq_2.txt", "cluster1") # piC-TBX5-RBM19 
plotUpset2("cluster2_piRNA_seq_2.txt", "cluster2") # piC-ZNF518B-WDR1


