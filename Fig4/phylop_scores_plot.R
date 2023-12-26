library(wig)
library(tidyverse)
library(cowplot)

# phyloP scores cut interval ----------------------------------------------
plot.cluster.interval <- function(data1, data2, data3, t1, t2, t3, output){ 
  
  wig_data <- import_wig(data1)
  wig_data$position <-as.numeric(as.factor(cut_interval(wig_data$pos, 100)))
  wig_data <- wig_data %>% group_by(position) %>% summarise(value=mean(val),Group=t1)
  
  wig_data2 <- import_wig(data2)
  wig_data2$position <-as.numeric(as.factor(cut_interval(wig_data2$pos, 100)))
  wig_data2 <- wig_data2 %>% group_by(position) %>% summarise(value=mean(val),Group=t2)
  
  wig_data3 <- import_wig(data3)
  wig_data3$position <-as.numeric(as.factor(cut_interval(wig_data3$pos, 100)))
  wig_data3 <- wig_data3 %>% group_by(position) %>% summarise(value=mean(val),Group=t3)
  
  wig_data_merge1 <- rbind(wig_data, wig_data2, wig_data3)
  
  p <- ggplot(wig_data_merge1, aes(x=position, y=value, group=Group, color=Group)) + geom_line(linewidth=1)+
    theme_cowplot() +
    coord_cartesian(ylim = c(-3, 3)) +
    labs(x="Position", y="PhyloP scores")
  
  pdf(output, width = 6, height = 3)
  print(p)
  dev.off()
}


# Cluster piC-ZNF518B-WDR1
plot.cluster.interval(
  "Cluster1.phyloPscores.alnFa.wig", "ZNF518B_WDR1_flank10k.phyloPscores.2.wig", "ZNF518B.phyloPscores.wig",
  "piC-ZNF518B-WDR1", "Flank intergenic region", "Flank gene ZNF518B",
  "Plot_piC-ZNF518B-WDR1_phyloP_scores.pdf"
)

# Cluster piC-TBX5-RBM19
plot.cluster.interval(
  "Cluster2.phyloPscores.alnFa.wig", "TBX5_RBM19_flank10k.phyloPscores.wig", "TBX5.phyloPscores.wig",
  "piC-TBX5-RBM19", "Flank intergenic region", "Flank gene TBX5",
  "Plot_piC-TBX5-RBM191_phyloP_scores.pdf"
)


