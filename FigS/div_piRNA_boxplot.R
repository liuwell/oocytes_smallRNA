library(ggpubr)
library(ggplot2)
library(cowplot)

Div.plot2 <- function(data, species){ 
  
  data1 <- read.table(data)
  data1$Div_Group <- factor(data1$Div_Group, levels = c("Low", "Medium-low","Medium-high", "High"))
  my_comparisons <- list( c("Low", "Medium-low"), c("Low", "Medium-high"), c("Low", "High") )
  p <- ggplot(data1, aes(x=Div_Group, y=piRNA, color=Div_Group, fill=Div_Group)) + 
    geom_boxplot(outlier.size = 0.5) + cowplot::theme_cowplot() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(x="Divergence", y="piRNA abundance (log10)", title = species) +
    scale_color_brewer(palette = "Set1", direction = 1) + 
    scale_fill_brewer(palette = "Set1", direction = 1)+
    stat_compare_means(aes(label = ..p.signif..), comparisons = my_comparisons, size=3, method = "t.test")
  
  return(p)

}


p1 <- Div.plot2("Div.Human.txt", "Human")
p2 <- Div.plot2("Div.Monkey.txt", "Monkey")
p3 <- Div.plot2("Div.Rabbit.txt", "Rabbit")
p4 <- Div.plot2("Div.Guinea pig.txt", "Guinea pig")

p5 <- Div.plot2("Div.Mouse.txt", "Mouse")
p6 <- Div.plot2("Div.Rat.txt", "Rat")
p7 <- Div.plot2("Div.Golden hamster.txt", "Golden hamster")
p8 <- Div.plot2("Div.Chinese hamster.txt", "Chinese hamster")

p9 <- Div.plot2("Div.Dog.txt", "Dog")
p10 <- Div.plot2("Div.Pig.txt", "Pig")
p11 <- Div.plot2("Div.Goat.txt", "Goat")
p12 <- Div.plot2("Div.Zebrafish.txt", "Zebrafish")

rcombined <- plot_grid(p1,p2,p3,p4, 
                       p5,p6,p7,p8,
                       p9,p10,p11,p12,
                       nrow = 3)

pdf("Div.piRNA.combined.pdf", width = 13, height = 11)
print(rcombined)
dev.off()
