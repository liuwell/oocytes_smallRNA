library(ggplot2)
library(cowplot)
library(tidyverse)


# ping-pong, divergence ---------------------------------------------------

plot.pp.div <- function(data1, data2, prefix){ 
  data <- read.csv(data1, sep = "\t", header = T, 
                   col.names = c("Name", "Sense", "Antisense", "piRNA", "TE", "Div", "Copies", "Lenght"))
  data <- na.omit(data)
  
  pingpong <- read.csv(data2, sep = "\t", header = F,col.names = c("Name", "pingpong"))
  data <- merge(data, pingpong, by="Name", all.x = T)
  data[is.na(data)] <- 0
  
  data$RPM <- data$piRNA/sum(data$piRNA)*1000000
  data <- data[order(data$RPM, decreasing = T), ]
  #print(data$Name[duplicated(data$Name)])
  print(dim(data))
  data$Name <- factor(data$Name, levels = rev(data$Name))
  p <- ggplot(data[1:30,], aes(x=Name, y=RPM, color=Div, 
                               size=pingpong)) + 
    geom_point() + theme_cowplot(font_size = 10) +
    coord_flip() +
    scale_fill_gradientn(colors = c("#A50026", "orange", "#313695")) +
    scale_color_gradientn(colors = c("#A50026", "orange", "#313695")) +
    scale_size_continuous(breaks = seq(0,30,5)) +
    labs(x=NULL, y="RPM", title = prefix) 
    #theme(legend.position = "none")
  
  return(p)
}

####
plot.pp.div2 <- function(data1, data2, prefix){ 
  data <- read.csv(data1, sep = "\t", header = T, 
                   col.names = c("TE", "Name", "Sense", "Antisense", "piRNA", "TE", "Div", "Copies", "Lenght"))
  data <- na.omit(data)
  
  pingpong <- read.csv(data2, sep = "\t", header = F,col.names = c("TE", "pingpong"))
  data <- merge(data, pingpong, by="TE")
  
  data$RPM <- data$piRNA/sum(data$piRNA)*1000000
  data <- data[order(data$RPM, decreasing = T), ]
  data$Name <- factor(data$Name, levels = rev(unique(data$Name)))
  p <- ggplot(data[1:30,], aes(x=Name, y=RPM, color=Div, fill=Div, size=pingpong)) + 
    geom_point() + theme_cowplot(font_size = 10) +
    coord_flip() +
    scale_fill_gradientn(colors = c("#A50026", "orange", "#313695")) +
    scale_color_gradientn(colors = c("#A50026", "orange", "#313695")) +
    labs(x=NULL, y="RPM", title = prefix) +
    scale_size_continuous(breaks = seq(0,30,5)) +theme(legend.position = "none")
  
  return(p)
}


p1 <- plot.pp.div("Human.div.txt", "Human.pingpong.txt", "Human")
p2 <- plot.pp.div("Monkey.div.txt", "Monkey.pingpong.txt", "Monkey")
p3 <- plot.pp.div("Rabbit.div.txt", "Rabbit.pingpong.txt", "Rabbit")
p4 <- plot.pp.div("GuineaPig.div.txt", "GuineaPig.pingpong.txt", "Guinea pig")

p5 <- plot.pp.div("Mouse.div.txt", "Mouse.pingpong.txt", "Mouse")
p6 <- plot.pp.div("Rat.div.txt", "Rat.pingpong.txt", "Rat")
p7 <- plot.pp.div2("GoldenHamster.div.txt", "GoldenHamster.pingpong.txt", "Golden hamster")
p8 <- plot.pp.div("ChineseHamster.div.txt", "ChineseHamster.pingpong.txt", "Chinese hamster")

p9 <- plot.pp.div("Dog.div.txt", "Dog.pingpong.txt", "Dog")
p10 <- plot.pp.div("Pig.div.txt", "Pig.pingpong.txt", "Pig")
p11 <- plot.pp.div("Goat.div.txt", "Goat.pingpong.txt", "Goat")
p12 <- plot.pp.div("Zebrafish.div.txt", "Zebrafish.pingpong.txt", "Zebrafish")

p.combeind <- plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, nrow = 2)
pdf("pingpong_div.pdf", width = 17, height = 8)
print(p.combeind)
dev.off()


# heatmap -----------------------------------------------------------------
library(pheatmap)

piRNA.corr <- read.csv("heatmap_piRNA_corr.csv", row.names = 1)
TE.corr <- read.csv("heatmap_TE_corr.csv", row.names = 1)

# piRNA
pheatmap(piRNA.corr[1:3,], filename = "Heatmap.piRNA.corr.pdf", width = 10, height = 3.5, 
         scale ="none" , cluster_cols = F, cluster_rows = F,
         color = rev(colorRampPalette(brewer.pal(n = 9, name = "RdBu"))(9)), 
         breaks=c(seq(-0.8, 0.8,length=9)), 
         display_numbers = T, number_format = "%.2f", 
         angle_col = 45, fontsize_row = 12,
         border_color = "black", na_col = "#DDDDDD")
# TE
pheatmap(TE.corr[1:3,], filename = "Heatmap.TE.corr.pdf", width = 10, height = 3.5, 
         scale ="none" , cluster_cols = F, cluster_rows = F,
         color = rev(colorRampPalette(brewer.pal(n = 9, name = "RdBu"))(9)), 
         breaks=c(seq(-0.8, 0.8,length=9)), 
         display_numbers = T, number_format = "%.2f",
         angle_col = 45, fontsize_row = 12,
         border_color = "black", na_col = "#DDDDDD")



