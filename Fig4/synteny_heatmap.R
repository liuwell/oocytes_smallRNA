library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(zoo)

#########################
### RPM
data <- read.csv("synteny_cluster.txt", sep = "\t", header = T, row.names = 1)
mapped <- c(21931027,57134403,33657884,47773155,
            53457193,22578199,13034058,27655905,
            48151599,36655335,12605704,8891766)
data <- t(t(data)/mapped * 1000000)

data <- t(log10(data+1))
data <- as.data.frame(data)

# reorder columns
data2 <- data[, c(1:12, 15,25,29,30,36:41,47,14,43,35,13,45,48, 49:53, 46,42,44,17,19,21,24,26,27,28,16,18,20,22,23,31:34, 54:56)]
pheatmap(data2, filename = "Cluster_synteny_RPM.pdf", width = 16, height = 6, 
         scale ="none" , cluster_cols = F, cluster_rows = F,
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")[1:4]))(11), 
         angle_col = 90, fontsize_row = 12,
         border_color = "black", na_col = "#DDDDDD")


#########################
### percent
total_piRNA <- c(14939334,34279045,31173757,45703454, 
                 10676160,16461440,11142708,23859831, 
                 36392271,24620450,9516378,2589642)
data <- read.csv("synteny_cluster.txt", sep = "\t", header = T, row.names = 1)
data <- t(data)/total_piRNA * 100
data <- as.data.frame(data)

# reorder columns
data2 <- data[, c(1:12, 15,25,29,30,36:41,47,14,43,35,13,45,48, 49:53, 46,42,44,17,19,21,24,26,27,28,16,18,20,22,23,31:34, 54:56)]
pheatmap(data2, filename = "Cluster_synteny_percent.pdf", width = 16, height = 6, 
         scale ="none" , cluster_cols = F, cluster_rows = F,
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")[1:4]))(11), 
         display_numbers = T, number_format = "%.1f",
         angle_col = 90, fontsize_row = 12,
         border_color = "black", na_col = "#DDDDDD")

### total percent
data2 <- na.fill(data2, 0)
x <- rowSums(data2)

data3 <- data.frame(species=rownames(data), sum=x)
write.table(data3, file = "total_percent_sum.txt", sep = "\t", quote = F)
