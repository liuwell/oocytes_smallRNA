library(ggbreak)
library(ggplot2)
library(cowplot)
library(tidyverse)


##########
barplot <- function(data){ 
  data <- read.csv(data, sep="\t", header = F)
  colnames(data) <- c("length",'miRNA', 'piRNA', 'tsRNA', 'rsRNA', 'snoRNA', 'lncRNA', 'mRNA', 'others')
  
  data.melt <- reshape2::melt(data, id.vars="length")
  data.melt$ratio <- data.melt$value/sum(data.melt$value)*100
  #data.melt$variable <- factor(data.melt$variable, levels = c("piRNA", "miRNA"))
  data.melt$variable <- factor(data.melt$variable, levels = c('miRNA', 'piRNA', 'tsRNA', 'rsRNA', 'snoRNA', 'lncRNA', 'mRNA', 'others'))
  
  ggplot(data.melt, aes(x=length, y=ratio, fill=variable)) + geom_bar(position = position_stack(reverse = T), stat="identity") +
    coord_cartesian(xlim = c(17, 40), ylim = c(0,30)) + scale_fill_brewer(palette = "Set1", direction = 1) +
    scale_y_continuous(expand = expansion(mult=c(0,0.05))) + labs(x="Length", y="% of total small RNA") + 
    theme_half_open() + theme(legend.position = "none")
}

# normalized by total spikein
barplot.spikein <- function(data, spikein, ylim=100){ 
  data <- read.csv(data, sep="\t", header = F)
  # 'miRNA', 'piRNA', 'tsRNA', 'rsRNA', 'snoRNA', 'lncRNA', 'mRNA', 'others'
  colnames(data) <- c("length",'miRNA', 'piRNA', 'tsRNA', 'rsRNA', 'snoRNA', 'lncRNA', 'mRNA', 'others')
  
  data.melt <- reshape2::melt(data, id.vars="length")
  data.melt$ratio <- data.melt$value/spikein
  data.melt$variable <- factor(data.melt$variable, levels = c('miRNA', 'piRNA', 'tsRNA', 'rsRNA', 'snoRNA', 'lncRNA', 'mRNA', 'others'))
  
  p <- ggplot(data.melt, aes(x=length, y=ratio, fill=variable)) + geom_bar(stat="identity", position=position_stack(reverse = T)) +
    coord_cartesian(xlim = c(17, 40), ylim = c(0, ylim)) + scale_fill_brewer(palette = "Set1", direction = 1) +
    scale_y_continuous(expand = expansion(mult=c(0,0.05))) + labs(x="Length", y="Normalized by spiekin") + 
    theme_half_open() + theme(legend.position = "none") #+ scale_y_break(c(10,90))
  return(p)
}

barplot.spikein2 <- function(data, spikein, c1, c2, ylim=100){ 
  data <- read.csv(data, sep="\t", header = F)
  colnames(data) <- c("length",'miRNA', 'piRNA', 'tsRNA', 'rsRNA', 'snoRNA', 'lncRNA', 'mRNA', 'others')
  data.melt <- reshape2::melt(data, id.vars="length")
  data.melt$ratio <- data.melt$value/spikein
  #data.melt$variable <- factor(data.melt$variable, levels = c("piRNA", "miRNA"))
  data.melt$variable <- factor(data.melt$variable, levels = c('miRNA', 'piRNA', 'tsRNA', 'rsRNA', 'snoRNA', 'lncRNA', 'mRNA', 'others'))
  data.melt <- data.melt %>% as_tibble() %>% dplyr::filter(length<=40)
  
  p<-  ggplot(data.melt, aes(x=length, y=ratio, fill=variable)) + geom_bar(position = position_stack(reverse = T), stat="identity") +
    #coord_cartesian(xlim = c(17, 40), ylim = c(0,30)) + 
    scale_fill_brewer(palette = "Set1", direction = 1) +
    scale_y_continuous(expand = expansion(mult=c(0,0.05))) + labs(x="Length", y="Normalized by spiekin") + 
    scale_y_break(c(c1, c2), scales = 0.5) + coord_cartesian(xlim = c(17, 40), ylim = c(0,ylim)) +
    theme_half_open() + theme(legend.position = "none")
  return(p)
}

# guinea pig
gp.input <- barplot.spikein("GuineaPig-input.length_RNA_counts.txt", 251741, 150)
gp.lgG <- barplot.spikein("GuineaPig-IP-IgG.length_RNA_counts.txt", 802257, 150)
gp.P1 <- barplot.spikein("GuineaPig-IP-P1.length_RNA_counts.txt", 541880, 150)
gp.P3 <- barplot.spikein("GuineaPig-IP-P3.length_RNA_counts.txt", 47798, 150)
prow <- plot_grid(
  gp.input, gp.lgG, gp.P1, gp.P3, 
  align = 'vh',
  labels = c("Input", "IgG IP", "Piwil1 IP", "Piwil3 IP"),
  hjust = -0.1,label_x = 0.25,
  nrow = 1, ncol=4
)
pdf("Lenght_GuineaPig_normalized.pdf", width = 11, height = 2.5)
print(prow)
dev.off()

### pig
pig.input <- barplot.spikein("Pig-input.length_RNA_counts.txt", 391299, 2)
pig.lgG <- barplot.spikein("Pig-IP-IgG.length_RNA_counts.txt", 1445305, 2)
pig.P1 <- barplot.spikein("Pig-IP-P1.length_RNA_counts.txt", 3250391, 2)
pig.P3 <- barplot.spikein("Pig-IP-P3.length_RNA_counts.txt", 986984, 2)
prow <- plot_grid(
  pig.input, pig.lgG, pig.P1, pig.P3, 
  align = 'vh',
  labels = c("Input", "IgG IP", "Piwil1 IP", "Piwil3 IP"),
  hjust = -0.1,label_x = 0.25,
  nrow = 1, ncol=4
)
pdf("Lenght_Pig_normalized.pdf", width = 11, height = 2.5)
prow
dev.off()

### Goat
goat.input <- barplot.spikein("Goat-input.length_RNA_counts.txt", 138264, 40)
goat.lgG <- barplot.spikein("Goat-IP-IgG.length_RNA_counts.txt", 1885821, 40)
goat.P1 <- barplot.spikein("Goat-IP-P1.length_RNA_counts.txt", 1662899, 40)
goat.P3 <- barplot.spikein("Goat-IP-P3.length_RNA_counts.txt", 63171, 40)
prow <- plot_grid(
  goat.input, goat.lgG, goat.P1, goat.P3, 
  align = 'vh',
  labels = c("Input", "IgG IP", "Piwil1 IP", "Piwil3 IP"),
  hjust = -0.1,label_x = 0.25,
  nrow = 1, ncol=4
)
pdf("Lenght_Goat_normalized.pdf", width = 11, height = 2.5)
prow
dev.off()


### Rat
rat.input <- barplot.spikein("Rat-input.length_RNA_counts.txt", 1105883, 3)
rat.lgG <- barplot.spikein("Rat-IP-IgG.length_RNA_counts.txt", 7427410, 3)
rat.P1 <- barplot.spikein("Rat-IP-P1.length_RNA_counts.txt", 1888394, 3)
prow <- plot_grid(
  rat.input, rat.lgG, rat.P1,
  align = 'vh',
  labels = c("Input", "IgG IP", "Piwil1 IP"),
  hjust = -0.15, label_x = 0.25,
  nrow = 1, ncol=4
)
pdf("Lenght_Rat_normalized.pdf", width = 11, height = 2.5)
prow
dev.off()

