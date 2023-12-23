library(pheatmap)
library(RColorBrewer)

### piRNA pathway gene
data <- read.csv("piRNA_pathway_gene.txt", sep="\t", row.names = 1)
data = log10(data+1)
# change the order of ChineseHamster and GoldenHamster
data <- data[, c(1:26, 32:35, 27:31, 36:49)]
anno_col = data.frame(row.names = colnames(data), 
                      Species=rep(c("Human", "Monkey", "Rabbit", "GuineaPig", "Mouse", "Rat",  "GoldenHamster", "ChineseHamster",
                                    "Dog", "Pig", "Goat", "Zebrafish"), c(5, 4, 4, 3, 5, 5, 4, 5, 3, 5, 3, 3)))
# creat colours for each group
color = brewer.pal(12, "Paired")
names(color) <- unique(anno_col$Species)
annoCol <- list(Species = color)
color = brewer.pal(12, "Paired")
color <- c("#E31A1C", "#FB9A99", "#1F78B4", "#A6CEE3", "#33A02C", "#B2DF8A",  
           "#FF7F00", "#FDBF6F", "#6A3D9A", "#CAB2D6", "#B15928", "#666666")
anno_colors = list(Species = c(Human=color[1], Monkey=color[2], Rabbit=color[3], GuineaPig=color[4],
                               Mouse=color[5], Rat=color[6], GoldenHamster=color[7], ChineseHamster=color[8], 
                               Dog=color[9], Pig=color[10], Goat=color[11], Zebrafish=color[12]))

pheatmap(data, filename = "piRNA_pathway_gene.pdf",  na_col = "black",  #border_color = F,
         color = colorRampPalette(brewer.pal(n = 9, name = "Blues"))(100),
         height = 6.8, width = 7.5, fontsize_row = 10, fontsize_col = 8, cluster_cols = F,
         annotation_col = anno_col, annotation_colors = anno_colors)



### miRNA pathway gene
data <- read.csv("miRNA_pathway_gene.txt", sep="\t", row.names = 1)
data = log10(data+1)

anno_col = data.frame(row.names = colnames(data), 
                      Species=rep(c("Human", "Monkey", "Rabbit", "GuineaPig", "Mouse", "Rat",  "GoldenHamster", "ChineseHamster",
                                    "Dog", "Pig", "Goat", "Zebrafish"), c(5, 4, 4, 3, 5, 5, 4, 5, 3, 5, 3, 3)))

# creat colours for each group
color = brewer.pal(12, "Paired")
names(color) <- unique(anno_col$Species)
annoCol <- list(Species = color)

color = brewer.pal(12, "Paired")
color <- c("#E31A1C", "#FB9A99", "#1F78B4", "#A6CEE3", "#33A02C", "#B2DF8A",  
           "#FF7F00", "#FDBF6F", "#6A3D9A", "#CAB2D6", "#B15928", "#666666")
anno_colors = list(Species = c(Human=color[1], Monkey=color[2], Rabbit=color[3], GuineaPig=color[4],
                               Mouse=color[5], Rat=color[6],  GoldenHamster=color[7], ChineseHamster=color[8],
                               Dog=color[9], Pig=color[10], Goat=color[11], Zebrafish=color[12]))
#x <- zoo::na.fill(data, 0)
#rownames(x) <- rownames(data)

pheatmap(data, filename = "miRNA_pathway_gene.pdf",  na_col = "black",  #border_color = "white", 
         color = colorRampPalette(brewer.pal(n = 9, name = "Blues"))(100),
         height = 6.8, width = 7.5, fontsize_row = 10, fontsize_col = 8, cluster_cols = F,
         annotation_col = anno_col, annotation_colors = anno_colors)

