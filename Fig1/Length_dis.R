library(ggplot2)
library(reshape2)
library(cowplot)
library(RColorBrewer)


# Length distribution -----------------------------------------------------
colors <- brewer.pal(9, "Set1")
colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#FDBF6F", "#999999")
barplot <- function(data){ 
  data <- read.csv(data, sep="\t", header = F)
  # 'miRNA', 'piRNA', 'tsRNA', 'rsRNA', 'snoRNA', 'lncRNA', 'mRNA', 'others'
  data$endo.siRNA <- 0
  colnames(data) <- c("length", "miRNA", "piRNA", "tsRNA", "rsRNA", "snoRNA", "lncRNA", "mRNA", "others", "endo-siRNA")
  
  data.melt <- reshape2::melt(data, id.vars="length")
  data.melt$ratio <- data.melt$value/sum(data.melt$value)*100
  data.melt$variable <- factor(data.melt$variable, levels = c("miRNA", "piRNA", "endo-siRNA","tsRNA", "rsRNA", "snoRNA", "lncRNA", "mRNA", "others"))
  
  ggplot(data.melt, aes(x=length, y=ratio, fill=variable)) + geom_bar(position = position_stack(reverse = TRUE), stat="identity") +
    coord_cartesian(xlim = c(16.5, 40), ylim = c(0,30)) + #scale_fill_brewer(palette = "Set1", direction = 1) +
    scale_fill_manual(values = colors) +
    #scale_y_continuous(expand = expansion(mult=c(0,0.05))) + labs(x="Length", y="Percent(%)") + 
    scale_y_continuous(expand = expansion(mult=c(0,0.05))) + labs(x="Length", y="% of total small RNA") + 
    theme_half_open() + 
    theme(legend.position = "none") 
  #theme(legend.position=c(0.6, 0.8), legend.title = element_blank(), legend.text = element_text(size = 18, face = "bold"))
}
barplot.mouse <- function(data){ 
  data <- read.csv(data, sep="\t", header = F)
  # miRNA', 'piRNA', 'siRNA','tsRNA', 'rsRNA', 'snoRNA', 'lncRNA', 'mRNA'
  colnames(data) <- c("length", "miRNA", "piRNA", "endo-siRNA", "tsRNA", "rsRNA", "snoRNA", "lncRNA", "mRNA", "others")
  
  data.melt <- reshape2::melt(data, id.vars="length")
  data.melt$ratio <- data.melt$value/sum(data.melt$value)*100
  data.melt$variable <- factor(data.melt$variable, 
                               levels = c("miRNA", "piRNA", "endo-siRNA", "tsRNA", "rsRNA", "snoRNA", "lncRNA", "mRNA", "others"))
  
  ggplot(data.melt, aes(x=length, y=ratio, fill=variable)) + geom_bar(position = position_stack(reverse = TRUE), stat="identity") +
    coord_cartesian(xlim = c(16.5, 40), ylim = c(0,30)) + scale_fill_manual(values = colors) +
    #    scale_y_continuous(expand = expansion(mult=c(0,0.05))) + labs(x="Length", y="Percent(%)") + 
    scale_y_continuous(expand = expansion(mult=c(0,0.05))) + labs(x="Length", y="% of total small RNA") + 
    theme_half_open() + 
    theme(legend.position = "none")
  #theme(legend.position=c(0.6, 0.8), legend.title = element_blank(), legend.text = element_text(size = 15, face = "bold"))
}


human  <- barplot("Human.length.txt")
monkey <- barplot("Monkey.length.txt")
rabbit <- barplot("Rabbit.length.txt")
gp <- barplot("GuineaPig.length.txt")

mouse  <- barplot.mouse("Mouse.length.txt")
rat <- barplot("Rat.length.txt")
gh  <- barplot("GoldenHamster.length.txt")
ch  <- barplot("ChineseHamster.length.txt") 

dog <- barplot("Dog.length.txt")
pig <- barplot("Pig.length.txt")
goat <- barplot("Goat.lengths.txt")
zf  <- barplot("Zebrafish.length.txt")


prow <- plot_grid(
  human, monkey, rabbit, gp, mouse, rat, gh,ch, dog,pig,goat,zf, 
  align = 'vh',
  labels = c("Human", "Monkey", "Rabbit", "Guinea pig", 
             "Mouse", "Rat","Golden hamster","Chinese hamster", 
             "Dog", "Pig", "Goat", "Zebrafish"),
  label_x = 0.25,
  hjust = -0.1, #vjust = -0.1,
  nrow = 3, ncol=4
)

pdf("Length.distribution.pdf", width = 10, height = 7.5)
prow
dev.off()


# phylogenetic tree ---------------------------------------------------------------

library(ggtree)
library(ggimage)

data <- read.tree("species.nwk")
ggtree(data) + geom_tiplab()
ggtree(data, size=2) + geom_tiplab(aes(label=str_replace(label, "_", " ")), offset = 0.05, font="italic")
ggtree(log10(data), size=2) + geom_tiplab(aes(label=str_replace(label, "_", " ")), offset = 0.05, font="italic")

pdf("species_tree.pdf", width = 10, height = 8)
ggtree(data, size=1) + geom_tiplab(aes(label=str_replace(label, "_", " ")), offset = 0.1) +
  theme_tree() + #geom_treescale(x=0, y=10, width=1) +
  geom_tiplab(aes(image=paste0("img/", label, '.png')), geom="image", offset=150, align=0) +
  xlim(NA, 800)
dev.off()



