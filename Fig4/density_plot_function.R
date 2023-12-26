library(aplot)
library(tidyverse)
library(cowplot)
library(stringr)
library(ggplot2)
library(gggenes)
library(RColorBrewer)

gene_cluster_plot1 <- function(data, TE.bed, chr, start, end, outname, strand.rev=FALSE, width=100, ylim1 = 0, ylim2 = 100){ 
  ### add ylim1, ylim2
  ### TE.bed, 7 columns
  if(strand.rev){ p1 <- data %>% as_tibble() %>% filter(V1==chr) %>% filter(V2 > start) %>% filter(V3 < end) %>%
    ggplot(aes(x=V2, y=abs(V4))) + geom_bar(stat="identity", width = width, color="black", fill="black") + labs(x=NULL, y=NULL) + 
    theme_cowplot(font_size = 10) + #background_grid(size.major = 0.1) + 
    coord_cartesian(xlim = c(start-1000, end+1000), ylim = c(ylim1, ylim2))+
    #scale_y_continuous(expand = expansion(mult=c(0,0.05))) +
    theme(axis.text.x = element_blank(), legend.position = "none")
  }else{p1 <- data %>% as_tibble() %>% filter(V1==chr) %>% filter(V2 > start) %>% filter(V3 < end) %>%
    ggplot(aes(x=V2, y=V4)) + geom_bar(stat="identity", width = width,  color="black", fill="black") + labs(x=NULL, y=NULL) + 
    theme_cowplot(font_size = 10) + #background_grid(size.major = 0.1) + 
    coord_cartesian(xlim = c(start-1000, end+1000), ylim = c(ylim1, ylim2))+
    #scale_y_continuous(expand = expansion(mult=c(0,0.05))) +
    theme(axis.text.x = element_blank(), legend.position = "none")}
  
  TE.bed2 <- TE.bed %>% as_tibble() %>% filter(V1==chr) %>% filter(V2 > start) %>% filter(V3 < end)
  p3 <- ggplot(TE.bed2, aes(xmin = V2, xmax = V3, y = V1, fill=V7, color=V7)) +
    geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(0.1, "mm")) + #geom_gene_label(aes(label = gene)) +
    labs(x=NULL, y=NULL) +
    #scale_fill_brewer(palette = "Set1") + scale_color_brewer(palette = "Set1") +
    scale_fill_manual(values = c("#E41A1C","#377EB8", "#FF7F00")) + scale_color_manual(values = c("#E41A1C","#377EB8", "#FF7F00")) +
    theme_genes() #+ theme(legend.position = "none")
  
  pdf(outname, height = 2, width = 8)
  print(p1 %>% insert_bottom(p3, height = 0.3))
  dev.off()
}

#####
plot_TE <- function(data, TE.bed, chr, start, end, outname, strand.rev=FALSE, width=100, p.height=2, ylim1 = 0, ylim2 = 100){ 
  ### TE.bed, 7 columns
  if(strand.rev){ p1 <- data %>% as_tibble() %>% filter(V1==chr) %>% filter(V2 > start) %>% filter(V3 < end) %>%
    ggplot(aes(x=V2, y=abs(V4))) + geom_bar(stat="identity", width = width, color="black", fill="black") + labs(x=NULL, y=NULL) + 
    theme_cowplot(font_size = 10) + #background_grid(size.major = 0.1) + 
    coord_cartesian(xlim = c(start-1000, end+1000), ylim = c(ylim1, ylim2))+
    #scale_y_continuous(expand = expansion(mult=c(0,0.05))) +
    theme(axis.text.x = element_blank(), legend.position = "none")
  }else{p1 <- data %>% as_tibble() %>% filter(V1==chr) %>% filter(V2 > start) %>% filter(V3 < end) %>%
    ggplot(aes(x=V2, y=V4)) + geom_bar(stat="identity", width = width,  color="black", fill="black") + labs(x=NULL, y=NULL) + 
    theme_cowplot(font_size = 10) + #background_grid(size.major = 0.1) + 
    coord_cartesian(xlim = c(start-1000, end+1000), ylim = c(ylim1, ylim2))+
    #scale_y_continuous(expand = expansion(mult=c(0,0.05))) +
    theme(axis.text.x = element_blank(), legend.position = "none")}
  
  TE.bed2 <- TE.bed %>% as_tibble() %>% filter(V1==chr) %>% filter(V2 > start) %>% filter(V3 < end)
  TE.bed2$V6 <- as.numeric(TE.bed2$V6)
  TE.bed2$V6 <- base::gsub(1, FALSE, TE.bed2$V6)
  TE.bed2$V6 <- base::gsub(2, TRUE, TE.bed2$V6)
  p3 <- ggplot(TE.bed2, aes(xmin = V2, xmax = V3, y = V1, fill=V4, color=V4, forward=V6)) +
    geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(0.5, "mm")) + 
    geom_gene_label(aes(label = V4), size=2) +
    labs(x=NULL, y=NULL) +
    scale_fill_brewer(palette = "Paired") + scale_color_brewer(palette = "Paired") +
    theme_genes() #+ theme(legend.position = "none")
  
  pdf(outname, height = p.height, width = 8)
  print(p1 %>% insert_bottom(p3, height = 0.3))
  dev.off()
}

gene_cluster_plot2 <- function(data, TE.bed, TE.bed3, chr, start, end, outname, strand.rev=FALSE, width=100, ylim1 = 0, ylim2 = 100){ 
  ### add ylim1, ylim2
  ### TE.bed, 7 columns
  if(strand.rev){ p1 <- data %>% as_tibble() %>% filter(V1==chr) %>% filter(V2 > start) %>% filter(V3 < end) %>%
    ggplot(aes(x=V2, y=abs(V4))) + geom_bar(stat="identity", width = width, color="black", fill="black") + labs(x=NULL, y=NULL) + 
    theme_cowplot(font_size = 10) + #background_grid(size.major = 0.1) + 
    coord_cartesian(xlim = c(start-1000, end+1000), ylim = c(ylim1, ylim2))+
    #scale_y_continuous(expand = expansion(mult=c(0,0.05))) +
    theme(axis.text.x = element_blank(), legend.position = "none")
  }else{p1 <- data %>% as_tibble() %>% filter(V1==chr) %>% filter(V2 > start) %>% filter(V3 < end) %>%
    ggplot(aes(x=V2, y=V4)) + geom_bar(stat="identity", width = width,  color="black", fill="black") + labs(x=NULL, y=NULL) + 
    theme_cowplot(font_size = 10) + #background_grid(size.major = 0.1) + 
    coord_cartesian(xlim = c(start-1000, end+1000), ylim = c(ylim1, ylim2))+
    #scale_y_continuous(expand = expansion(mult=c(0,0.05))) +
    theme(axis.text.x = element_blank(), legend.position = "none")}
  
  TE.bed2 <- TE.bed %>% as_tibble() %>% filter(V1==chr) %>% filter(V2 > start) %>% filter(V3 < end)
  p3 <- ggplot(TE.bed2, aes(xmin = V2, xmax = V3, y = V1, fill=V7, color=V7)) +
    geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(0.1, "mm")) + #geom_gene_label(aes(label = gene)) +
    labs(x=NULL, y=NULL) +
    #scale_fill_brewer(palette = "Set1") + scale_color_brewer(palette = "Set1") +
    scale_fill_manual(values = c("#E41A1C","#377EB8", "#FF7F00")) + scale_color_manual(values = c("#E41A1C","#377EB8", "#FF7F00")) +
    theme_genes() #+ theme(legend.position = "none")
  
  TE.bed3 <- TE.bed3 %>% as_tibble() %>% filter(V1==chr) %>% filter(V2 > start) %>% filter(V3 < end)
  TE.bed3$V6 <- as.numeric(TE.bed3$V6)
  TE.bed3$V6 <- base::gsub(1, FALSE, TE.bed3$V6)
  TE.bed3$V6 <- base::gsub(2, TRUE, TE.bed3$V6)
  p4 <- ggplot(TE.bed3, aes(xmin = V2, xmax = V3, y = V1, fill=V7, color=V7, forward=V6)) + 
    geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(3, "mm")) + #geom_gene_label(aes(label = gene)) +
    labs(x=NULL, y=NULL) +
    #scale_fill_brewer(palette = "Set1") + scale_color_brewer(palette = "Set1") +
    scale_fill_manual(values = c("#E41A1C","#377EB8", "#FF7F00")) + scale_color_manual(values = c("#E41A1C","#377EB8", "#FF7F00")) +
    theme_genes() #+ theme(legend.position = "none")
  pdf(outname, height = 2, width = 8)
  print(p1 %>% insert_bottom(p3, height = 0.3) %>% insert_bottom(p4, height = 0.3))
  dev.off()
}

