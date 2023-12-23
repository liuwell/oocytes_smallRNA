library(Hmisc)
library(reshape2)
library(cowplot)
library(ggrepel)
library(corrplot)


# corrlation plot of TE features ------------------------------------------

col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
plot1 <- function(fi, prefix){ 
  
  data <- read.csv(fi, sep = "\t", header = T, 
                   col.names = c("Name", "Sense", "Antisense", "piRNA", "TE", "Div", "Copies", "Lenght"))
  data <- na.omit(data)
  
  ### correlation
  data1 <- data[,4:8]
  data1$piRNA <- log10(data1$piRNA+1)
  data1$TE <- log10(data1$TE+1)
  data1$Copies <- log10(data1$Copies+1)
  data1$Lenght <- log10(data1$Lenght+1)
  #data.corr <- as.data.frame(cor(data, method = "pearson"))
  #data.corr.p <- cor.mtest(data, conf.level=0.95)
  data.corr <-rcorr(as.matrix(data1[, 1:4]), type = "pearson")
  pdf(paste0(prefix, ".corr.pdf"), width = 5.5, height = 5)
  corrplot(data.corr$r, type="upper", tl.pos = "d", p.mat = data.corr$P, insig = "label_sig", #method = "color",
           #col = colorRampPalette("blue", "white", "red")(200),
           col=rev(col2(200)),
           sig.level = c(0.001, 0.01, 0.05), pch.cex = 2, pch.col = "black", mar = c(1,1,1,1))
  corrplot(data.corr$r, add=TRUE, method = "number", col="black", type="lower", diag=FALSE,tl.pos="n", cl.pos="n")
  dev.off()
  
  ### All ### Scatterplot, piRNA, RNAseq, consensus
  data$RPM <- log10(data$piRNA/sum(data$piRNA)*1000000 + 1)
  data$FPKM <- log10(data$TE+1)
  data$total <- rowSums(data[,9:10])
  subd <- data[order(data$total, decreasing = T), ][1:10,]
  x <- paste0("cor=", round(cor(data$RPM, data$FPKM),2))
  p <- ggplot(data, aes(x=RPM, y=FPKM, color=Div, fill=Div)) + geom_point(color="black", size=3, shape=21) + theme_cowplot(font_size = 12) + 
    coord_cartesian(ylim = c(0, 6), xlim = c(0, 6)) + 
    labs(x="piRNA abundance (log10, RPM)", y="TE abundance (log10, FPKM)", title = prefix) +
    #scale_fill_gradientn(colours =rainbow(3)) + 
    ggplot2::annotate("text", x=1, y=5.5, label=x, size=5) +
    scale_fill_gradientn(colors = c("#A50026", "#FFFFBF", "#313695")) +
    geom_abline(intercept = 0, slope = 1,size=0.5) +
    geom_abline(intercept = log10(2), slope = 1,size=0.3,linetype="dashed")+
    geom_abline(intercept = -log10(2), slope = 1,size=0.3,linetype="dashed") +
    geom_text_repel(data=subd, aes(label=Name), color="blue", size=3)
  pdf(paste0(prefix, ".piRNA.TE.div.pdf"), width = 5, height = 4)
  print(p)
  dev.off()
  
  # All ### Scatterplot, Sense, Antisense, consensus
  data$x1 <- log10(data$Sense/sum(data$piRNA)*1000000+1)
  data$x2 <- log10(data$Antisense/sum(data$piRNA)*1000000+1)
  subd <- data[order(data$piRNA, decreasing = T), ][1:10,]
  x <- paste0("cor=", round(cor(data$x1, data$x2),2))
  p <- ggplot(data, aes(x=x1, y=x2, color=Div, fill=Div)) + geom_point(color="black", size=3, shape=21) + theme_cowplot(font_size = 12) + 
    coord_cartesian(ylim = c(0, 6), xlim = c(0, 6)) + 
    labs(x="Sense piRNA abundance(log10, RPM)", y="Antisense piRNA abundance(log10, RPM)", title = prefix) +
    #scale_fill_gradientn(colours =rainbow(3)) + 
    ggplot2::annotate("text", x=1, y=5.5, label=x, size=5) +
    scale_fill_gradientn(colors = c("#A50026", "#FFFFBF", "#313695")) +
    geom_abline(intercept = 0, slope = 1,size=0.5) +
    geom_abline(intercept = log10(2), slope = 1,size=0.3,linetype="dashed")+
    geom_abline(intercept = -log10(2), slope = 1,size=0.3,linetype="dashed") +
    geom_text_repel(data=subd, aes(label=Name), color="blue", size=3)
  pdf(paste0(prefix, ".piRNA.strand.div.pdf"), width = 5, height = 4)
  print(p)
  dev.off()
}

# plot correlation
plot1("Human.piRNA.TE.div.txt", "Human")
plot1("Monkey.piRNA.TE.div.txt", "Monkey")
plot1("Rabbit.piRNA.TE.div.txt", "Rabbit")
plot1("GuineaPig.piRNA.TE.div.txt", "GuineaPig")

plot1("Mouse.piRNA.TE.div.txt", "Mouse") 
plot1("Rat.piRNA.TE.div.txt", "Rat")
plot1("GH.piRNA.TE.div.txt", "GH")
plot1("CH.piRNA.TE.div.txt", "CH")

plot1("Dog.piRNA.TE.div.txt", "Dog")
plot1("Pig.piRNA.TE.div.txt", "Pig")
plot1("Goat.piRNA.TE.div.txt", "Goat")
plot1("Zebrafish.piRNA.TE.div.txt", "Zebrafish")


