library(tidyverse)
library(cowplot)
library(aplot)
library(RColorBrewer)


# TE derived piRNA of cluster1 and cluster2 -------------------------------

plot1<- function(data, outfile, topN=20, width=12){  
  data <- read.csv(data, sep = "\t", header = F)
  data$V7 <- gsub("LINE", "1_LINE", data$V7)
  data$V7 <- gsub("SINE", "2_SINE", data$V7)
  data$V7 <- gsub("LTR", "3_LTR", data$V7)
  data$V7 <- gsub("DNA", "4_DNA", data$V7)
  
  data2 <- data[1:100, 5:8] %>% as_tibble() %>% filter(V6>0) %>% 
    dplyr::filter(V7 %in% c("1_LINE", "2_SINE", "3_LTR", "4_DNA")) %>%
    group_by(V7) %>% top_n(topN, V6) %>% arrange(V7)
  data2$V8 <- factor(data2$V8, levels = data2$V8)
  
  p1 <- ggplot(data2, aes(x=V8, y=V6, fill=V7, color=V7)) + geom_bar(stat = "identity") +
    cowplot::theme_cowplot() + labs(x=NULL, y="RPM") + 
    scale_y_continuous(expand = expansion(mult=c(0,0.05))) +
    scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
    theme(axis.text.x = element_blank())
  
  p2 <- ggplot(data2, aes(x=V8, y=0)) +geom_tile(aes(fill=V5)) +
  #  scale_fill_gradientn(name="Strand bias", low = "#A50026", high="#313695") +
    scale_fill_gradientn(colors = brewer.pal(11, "RdBu"), name="Strand bias") +
    cowplot::theme_cowplot() + labs(x=NULL, y=NULL) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.9, hjust = 1), 
          axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  pdf(outfile, width = width, height = 5)
  print(p1 %>% insert_bottom(p2, height = 0.2))
  dev.off()

}

###
plot2<- function(data, outfile, topN=20, height=12, width=7){  
  data <- read.csv(data, sep = "\t", header = F)
  data$V7 <- gsub("LINE", "1_LINE", data$V7)
  data$V7 <- gsub("SINE", "2_SINE", data$V7)
  data$V7 <- gsub("LTR", "3_LTR", data$V7)
  data$V7 <- gsub("DNA", "4_DNA", data$V7)
  
  data2 <- data[1:100, 5:8] %>% as_tibble() %>% filter(V6>0)  %>% 
    dplyr::filter(V7 %in% c("1_LINE", "2_SINE", "3_LTR", "4_DNA")) %>%
    group_by(V7) %>% top_n(topN, V6)  %>% arrange(V6)  %>% arrange(desc(V7))
  data2$V8 <- factor(data2$V8, levels = data2$V8)
  
  p1 <- ggplot(data2, aes(x=V8, y=V6, fill=V7, color=V7)) + geom_bar(stat = "identity") +
    cowplot::theme_cowplot() + labs(x=NULL, y="RPM") +
    scale_y_continuous(expand = expansion(mult=c(0,0.05))) +
    scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
    theme(legend.position = "none", axis.text.y = element_blank()) +coord_flip() 

  p2 <- ggplot(data2, aes(x=V8, y=0)) +geom_tile(aes(fill=V5)) +
    scale_fill_gradientn(colors = brewer.pal(11, "RdBu"), name="Strand bias") +
    cowplot::theme_cowplot() + labs(x=NULL, y=NULL) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") + 
    coord_flip() 

  pcol <- plot_grid(
    p2,p1, align = "h", axis = "l",
    ncol = 2, rel_widths = c(0.8,1)
  )
  
  pdf(outfile, width = width, height = height)
  print(pcol)
  dev.off()
  
}

# Human
plot1("Human_cluster1_TE.txt", "Human_cluster1.pdf")
plot1("Human_cluster2_TE.txt", "Human_cluster2.pdf")

# Monkey
plot1("Monkey_cluster1_TE.txt", "Monkey_cluster1.pdf")
plot1("Monkey_cluster2_TE.txt", "Monkey_cluster2.pdf")

# Rabbit
plot1("Rabbit_cluster1_TE.txt", "Rabbit_cluster1.pdf", topN=10, width = 6)
plot1("Rabbit_cluster2_TE.txt", "Rabbit_cluster2.pdf", topN=10, width = 6)

# GuineaPig
plot1("GuineaPig_cluster1_TE.txt", "GuineaPig_cluster1.pdf", topN=10, width = 8)
plot1("GuineaPig_cluster2_TE.txt", "GuineaPig_cluster2.pdf", topN=10, width = 8)


# Pig
plot1("Pig_cluster1_TE.txt", "Pig_cluster1.pdf", topN=10, width = 8)
plot1("Pig_cluster2_TE.txt", "Pig_cluster2.pdf", topN=10, width = 8)


# dog
plot1("Dog2_cluster1_TE.txt", "Dog_cluster1.pdf", topN=10, width = 8)
plot1("Dog2_cluster2_TE.txt", "Dog_cluster2.pdf", topN=10, width = 8)


# Goat
plot1("Goat_cluster1_TE.txt", "Goat_cluster1.pdf", topN=10, width = 10)
plot1("Goat_cluster2_TE.txt", "Goat_cluster2.pdf", topN=10, width = 10)


# CH
plot1("CH_cluster1_TE.txt", "CH_cluster1.pdf", topN=10, width = 10)
plot1("CH_cluster2_TE.txt", "CH_cluster2.pdf", topN=10, width = 10)

# GH
plot1("GH_cluster1_TE.txt", "GH_cluster1.pdf", topN=10, width = 10)
plot1("GH_cluster2_TE.txt", "GH_cluster2.pdf", topN=10, width = 10)





