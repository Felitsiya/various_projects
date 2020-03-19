library(readxl)
library(reshape2)
library(dplyr)
library(ggplot2)

setwd("I:/boxplot")

df <- as.data.frame(read_excel("cancer_data.xlsx", col_names = TRUE))
df <- df[,2:14]
df_t <- melt(df[1:82,])
df_t$type <- "tumour tissue"
df_n <- melt(df[83:164,])
df_n$type <- "normal tissue"
df_melted <- as.data.frame(rbind(df_t, df_n))
df_melted <- na.omit(df_melted)

df_melted <- df_melted %>% mutate(variable =  factor(variable, levels = colnames(df)),
                                  type = factor(type, levels = c("tumour tissue", "normal tissue"),
                                                ordered = TRUE)) %>%
                           arrange(variable) %>% as.data.frame()

ggplot(df_melted, aes(x = variable, y = value)) + 
  geom_boxplot(aes(fill = type), width = 0.5, outlier.size = 0.3) + 
  ylab("Normalized dCt Expression Levels") + 
  theme(text = element_text(size = 16),
        axis.line = element_line(colour = "black"), 
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5)) +
  scale_x_discrete(breaks = colnames(df)) +
  scale_fill_manual(values = c("#FF6666", "#3399CC")) +
  annotate("text", x = 1:13,y = 15, label = c("***","***","*","***","**","*","*","***",
                                            "*","*","***", "**","**"), color = "black")
