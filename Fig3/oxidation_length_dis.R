library(ggplot2)
library(reshape2)
library(cowplot)


##########
# normalized by total spikein
barplot.spikein <- function(data, spikein, ylim=100){ 
  data <- read.csv(data, sep="\t", header = F)
  # 'miRNA', 'piRNA', 'tsRNA', 'rsRNA', 'snoRNA', 'lncRNA', 'mRNA', 'others'
  data <- data[, 1:3]
  colnames(data) <- c("length", "miRNA", "piRNA")
  
  data.melt <- reshape2::melt(data, id.vars="length")
  data.melt$ratio <- data.melt$value/spikein
  data.melt$variable <- factor(data.melt$variable, levels = c("piRNA", "miRNA"))
  
  ggplot(data.melt, aes(x=length, y=ratio, fill=variable)) + geom_bar(position = "stack", stat="identity") +
    coord_cartesian(xlim = c(17, 40), ylim = c(0, ylim)) + scale_fill_brewer(palette = "Set1", direction = -1) +
    scale_y_continuous(expand = expansion(mult=c(0,0.05))) + labs(x="Length", y="Normalized by spiekin") + 
    theme_half_open() + theme(legend.position = "none")
  
}

barplot.mouse.spikein <- function(data, spikein, ylim=100){ 
  # miRNA', 'piRNA', 'siRNA','tsRNA', 'rsRNA', 'snoRNA', 'lncRNA', 'mRNA'
  data <- read.csv(data, sep="\t", header = F)
  data <- data[, 1:4]
  colnames(data) <- c("length", "miRNA", "piRNA", "endo-siRNA")
  
  data.melt <- reshape2::melt(data, id.vars="length")
  data.melt$ratio <- data.melt$value/spikein
  data.melt$variable <- factor(data.melt$variable, levels = c( "miRNA","piRNA", "endo-siRNA"))
  
  ggplot(data.melt, aes(x=length, y=ratio, fill=variable)) + geom_bar(position = "stack", stat="identity") +
    coord_cartesian(xlim = c(17, 40), ylim = c(0, ylim)) + scale_fill_brewer(palette = "Set1", direction = 1) +
    scale_y_continuous(expand = expansion(mult=c(0,0.05))) + labs(x="Length", y="Normalized by spiekin") + 
    theme_half_open() + theme(legend.position = "none")
  
}

# B: no oxidation
# O: oxidation
h1 <- barplot.spikein("HumanB.length_RNA_counts.txt", 2008, 100)
h2 <- barplot.spikein("HumanO.length_RNA_counts.txt", 471, 100)
mk1 <- barplot.spikein("MonkeyB.length_RNA_counts.txt", 5523, 30)
mk2 <- barplot.spikein("MonkeyO.length_RNA_counts.txt", 18537, 30)
rb1 <- barplot.spikein("RabbitB.length_RNA_counts.txt", 37911, 15)
rb2 <- barplot.spikein("RabbitO.length_RNA_counts.txt", 44315, 15)
gp1 <- barplot.spikein("GuineaPigB.length_RNA_counts.txt", 59412, 15)
gp2 <- barplot.spikein("GuineaPigO.length_RNA_counts.txt", 54930, 15)

gh1 <- barplot.spikein("GoldenHamsterB.length_RNA_counts.txt", 38427, 30)
gh2 <- barplot.spikein("GoldenHamsterO.length_RNA_counts.txt", 14925, 30)
m1 <- barplot.mouse.spikein("MouseB.length_RNA_counts.txt", 626, 300)
m2 <- barplot.mouse.spikein("MouseO.length_RNA_counts.txt", 2703, 300)
rat1 <- barplot.spikein("RatB.length_RNA_counts.txt", 6002, 60)
rat2 <- barplot.spikein("RatO.length_RNA_counts.txt", 2890, 60)
ch1 <- barplot.spikein("ChineseHamsterB.length_RNA_counts.txt", 10380, 200)
ch2 <- barplot.spikein("ChineseHamsterO.length_RNA_counts.txt", 6368, 200)

dog1 <- barplot.spikein("DogB.length_RNA_counts.txt", 40474,  200)
dog2 <- barplot.spikein("DogO.length_RNA_counts.txt", 21945, 200)
pig1 <- barplot.spikein("PigB.length_RNA_counts.txt", 2518, 120)
pig2 <- barplot.spikein("PigO.length_RNA_counts.txt", 5219, 120)
goat1 <- barplot.spikein("GoatB.length_RNA_counts.txt", 7539, 160)
goat2 <- barplot.spikein("GoatO.length_RNA_counts.txt", 12619, 160)
zf1 <- barplot.spikein("ZebrafishB.length_RNA_counts.txt", 147749.9, 10)
zf2 <- barplot.spikein("ZebrafishO.length_RNA_counts.txt", 72015, 10)

prow <- plot_grid(
  h1,h2, mk1, mk2, rb1,rb2, gp1,gp2, gh1,gh2, m1,m2, 
  align = 'vh',
  labels = c("Human ctrl", "Human treat", "Monkey ctrl", "Monkey treat", "Rabbit ctrl", "Rabbit treat",
             "GP ctrl","GP treat", "GH ctrl", "GH treat", "Mouse ctrl", "Mouse treat"),
  hjust = -0.1,label_x = 0.25,
  nrow = 3, ncol=4
)
prow2 <- plot_grid(
  rat1,rat2, ch1, ch2, dog1,dog2, pig1,pig2, goat1,goat2, zf1,zf2, 
  align = 'vh',
  labels = c("Rat ctrl", "Rat treat", "CH ctrl", "CH treat", "Dog ctrl", "Dog treat", 
             "Pig ctrl", "Pig treat", "Goat ctrl", "Goat treat", "ZF ctrl", "ZF treat"),
  hjust = -0.1,label_x = 0.25,
  nrow = 3, ncol=4
)

pdf("Lenght_normalized.oxidation.pdf", width = 11, height = 7.5)
prow
prow2
dev.off()

