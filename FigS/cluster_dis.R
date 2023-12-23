library(ggplot2)
library(cowplot)


# cluster dis -------------------------------------------------------------

combine1 <- read.csv("cluster_dis_data1.csv", row.names = 1)
combine2 <- read.csv("cluster_dis_data2.csv", row.names = 1)
combine3 <- read.csv("cluster_dis_data3.csv", row.names = 1)

c12 <- c("#E31A1C", "#FB9A99", "#1F78B4", "#A6CEE3", "#33A02C", "#B2DF8A",  
         "#FF7F00", "#FDBF6F", "#6A3D9A", "#CAB2D6", "#B15928", "#666666")
p1 <- ggplot(combine1, aes(x=num, y=cum, color=species, fill=species)) + 
  geom_line(size=1) + geom_point(size=1.2) +
  coord_cartesian(ylim = c(0,100)) + scale_fill_manual(values = c12[1:4]) + 
  scale_color_manual(values = c12[1:4])+
  labs(x="Top expressed piRNA clusters", y="% of cluster \nderived piRNA reads") + theme_cowplot()+ 
  geom_vline(xintercept = 10, linetype="dashed")

p1
p2 <- ggplot(combine2, aes(x=num, y=cum, color=species, fill=species)) + 
  geom_line(size=1) + geom_point(size=1.2) +
  coord_cartesian(ylim = c(0,100))+ scale_fill_manual(values = c12[5:8]) + 
  scale_color_manual(values = c12[5:8])+
  labs(x="Top expressed piRNA clusters", y="% of cluster \nderived piRNA reads") + theme_cowplot()+ 
  geom_vline(xintercept = 10, linetype="dashed")
p2
p3 <- ggplot(combine3, aes(x=num, y=cum, color=species, fill=species)) + 
  geom_line(size=1) + geom_point(size=1.2) +
  coord_cartesian(ylim = c(0,100))+ scale_fill_manual(values = c12[9:12]) + 
  scale_color_manual(values = c12[9:12])+
  labs(x="Top expressed piRNA clusters", y="% of cluster \nderived piRNA reads") + theme_cowplot()+ 
  geom_vline(xintercept = 10, linetype="dashed")
p3

pdf("piCluster_dis.pdf", width = 10, height = 6)
plot_grid(p1,p2,p3, ncol = 2)
dev.off()


# cluster number ----------------------------------------------------------

data <- read.csv("cluster_number.txt", header = F, sep="\t")
data$V1 <- factor(data$V1, levels=c("Human", "Monkey", "Rabbit", "GuineaPig", 
                                    "Mouse","Rat","GoldenHamster", "ChineseHamster",
                                    "Dog", "Pig", "Goat", "Zebrafish"))
pdf("Cluster_numbers.pdf", width = 4, height = 4)
ggplot(data, aes(x=V1,y=V2, fill=V1, color=V1)) + #geom_bar(stat = "identity", position = "stack", fill="#E31A1C") +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c12) + scale_color_manual(values = c12)+
  labs(x="", y="Cluster numbers") + theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "none") 
dev.off()
