library(ggplot2)
library(reshape2)
library(cowplot)
library(tidyverse)


### normalize by size factor
calc_sf <- function (expr_mat){
  geomeans <- exp(rowMeans(log(expr_mat)))
  SF <- function(cnts){
    median((cnts/geomeans)[(is.finite(geomeans) & geomeans >0)])
  }
  norm_factor <- apply(expr_mat,2,SF)
  return(t(t(expr_mat)/norm_factor))
}

data.sf <- data.frame(calc_sf(as.matrix(data)))

###
sample.freq <- function(data){ 
  a <- NULL
  for(i in 1:nrow(data)){
    x <- 0
    for(j in 1:ncol(data)){
      if(data[i,j]>0){ 
        x <- x+1
      }
    }
    a <- append(a,x)
  }
  return(a)
}

data <- read.csv("miRNA_families_expression_counts.csv", row.names = 1)
a <- sample.freq(data)
data.freq <- data[a>2,]
data.sf <- data.frame(calc_sf(as.matrix(data.freq)))
data.melt <- reshape2::melt(data.sf)

### DEseq2
library(DESeq2)
species <- factor(rep(c("ChineseHamster","Dog", "Goat", "GoldenHamster", "GuineaPig", "Human", "Monkey", "Mouse", "Pig","Rabbit","Rat","Zebrafish"), 
                      c(5, 3, 3, 4, 4, 4, 5, 5, 5, 4, 3, 3)), 
                  levels = c("Human", "Monkey", "Rabbit", "GuineaPig", "Mouse", "Rat",  
                             "GoldenHamster", "ChineseHamster", "Dog", "Pig", "Goat", "Zebrafish"))
meta <- data.frame(sampletype=species)
rownames(meta) <- colnames(data)
dds <- DESeqDataSetFromMatrix(countData = data.freq, colData = meta, design = ~ sampletype)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="miRNA_normalized_counts.txt", sep="\t", quote=F, col.names=NA)

c12 <- c("#E31A1C", "#FB9A99", "#1F78B4", "#A6CEE3", "#8C510A", "#BF812D",  
         "#FF7F00", "#FDBF6F", "#6A3D9A", "#CAB2D6", "#252525", "#666666")

p <- data.melt %>% as_tibble()%>% filter(value >0)%>%dplyr::count(variable) %>% 
  mutate(Species=species) %>%
  ggplot(aes(x=Species, y=n, fill=Species, color=Species)) + 
  geom_boxplot(outlier.shape = 18,outlier.size = 1)+
  #scale_fill_brewer(palette = "Paired", direction = 1) + scale_color_brewer(palette = "Paired", direction = 1) +
  scale_fill_manual(values = c12) + scale_color_manual(values = c12)+ 
  coord_cartesian(ylim = c(0,200)) + labs(x=NULL,y="miRNA family numbers") + 
  scale_y_continuous(expand = expansion(mult=c(0,0.05))) +
  theme_half_open() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=10), legend.position = "right") 

pdf("miRNA_family_numbers.pdf", width = 6,height = 4)
p
dev.off()

### miRNA number
miRNA.family.averge <- data.melt %>% as_tibble()%>% filter(value >0)%>%dplyr::count(variable) %>% 
  mutate(Species=species) %>% group_by(Species) %>% summarise(m=mean(n))
write.table(miRNA.family.averge, file = 'miRNA_family_number_average.txt', sep = "\t", quote=F)

data.melt %>% as_tibble()%>% filter(value >0)%>%dplyr::count(variable) %>% mutate(Species=species) -> miRNA.number
write.table(miRNA.number, file = 'miRNA_family_number.txt', sep = "\t", quote=F)

### heatmap
library(pheatmap)
library(RColorBrewer)

data.sf <- data.sf[order(rowSums(data.sf), decreasing = T),]
data.log = log10(data.sf+1)
anno_col = data.frame(row.names = colnames(data.log), 
                      Species=species)
# creat colours for each group
#color = brewer.pal(12, "Paired")
#names(color) <- unique(anno_col$Species)
#color <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", 
#           "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")
color <- c12
names(color) <- c("Human", "Monkey", "Rabbit", "GuineaPig", "Mouse", "Rat", 
                  "GoldenHamster", "ChineseHamster", "Dog", "Pig", "Goat", "Zebrafish")
annoCol <- list(Species = color)
#anno_colors = list(Species = c(Human=color[1], Monkey=color[2], Rabbit=color[3], GuineaPig=color[4],
#                               Mouse=color[5], Rat=color[6], ChineseHamster=color[7], GoldenHamster=color[8],
#                               Dog=color[9], Pig=color[10], Sheep=color[11], Zebrafish=color[12]))

# heatmap
data2 <- data.log[, c(20:28, 39:42, 16:19, 29:33, 43:45, 12:15, 1:8, 34:38, 9:11, 46:48)]
pheatmap(data2[1:100,], filename = "Plot_miRNA_family_heatmap.pdf",  #border_color = NA,
         color = colorRampPalette(brewer.pal(n = 9, name = "Reds"))(100), cluster_cols = F,
         height = 8, width = 10, fontsize_row = 5, fontsize_col = 8,
         annotation_col = anno_col, annotation_colors = annoCol)

### PCA
pca<-prcomp(t(data.log), center = T, scale. =F)
pca2 <- as.data.frame(pca$x)
pca2$Species <- species

labx <- paste("PC1", sprintf('(%0.1f%%)', 100 * pca$sdev[1]^2/sum(pca$sdev^2)))
laby <- paste("PC2", sprintf('(%0.1f%%)', 100 * pca$sdev[2]^2/sum(pca$sdev^2)))

c12 <- c("#E31A1C", "#FB9A99", "#1F78B4", "#A6CEE3", "#8C510A", "#BF812D",  
         "#FF7F00", "#FDBF6F", "#6A3D9A", "#CAB2D6", "#252525", "#666666")

p <- ggplot(pca2, aes(x=PC1, y=PC2, color=Species)) +
  geom_point(size=6, alpha=0.8) +labs(x=labx, y=laby) + 
  scale_fill_manual(values = c12) + scale_color_manual(values = c12)+ 
  theme_cowplot(font_size = 23, font_family = "", line_size = 1)

output.2d <- "miRNA_family_PCA.pdf"
pdf(output.2d, width = 10, height = 6)
print(p)
dev.off()
