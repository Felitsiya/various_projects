# Barplots with modified y-axis
# https://stackoverflow.com/a/44697832
library(readxl)
library(reshape2)
library(dplyr)
library(ggplot2)

setwd("F:/barplot")

tissue <- read_excel("RQ values interrupted histograms.xlsx", col_names = TRUE)

trans <- function(x){pmin(x, 2) + 0.05 * pmax(x - 2, 0)}
tissue_mod <- as.data.frame(lapply(tissue[2:14], trans))
tissue_mod <- cbind(tissue$Samples, tissue_mod)
colnames(tissue_mod) <- colnames(tissue)

tissue$Samples <- factor(tissue$Samples, levels = tissue$Samples)
yticks <- c(0, 1, 2, 10, 20, 30, 40)

len <- ncol(tissue_mod)
for (i in 2:len) {
  nam <- noquote(colnames(tissue_mod)[i])
  x <- ggplot(data = tissue_mod[,c(1,i)], aes(x = Samples, y = get(nam))) +
    geom_col(position = "dodge", fill = "#FF9966") +
    ylab(paste("RQ", nam, sep = " ")) +
    theme(text = element_text(size = 16),
          axis.line = element_line(colour = "black"), 
          panel.background = element_blank(),
          panel.grid.major = element_line(colour = "#FFFFFF"), 
          legend.title=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y = element_text(size = 10, face = "bold"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7, face = "bold", color = "black"),
          axis.text.y = element_text(size = 7, color = "black", face = "bold")) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), breaks = trans(yticks), labels=yticks) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", size = 1) +
    geom_hline(yintercept = 2, linetype = "dashed", color = "black", size = 1)
  print(x)
  ggsave(paste('barplot',i, '.png'))
}
