library(dplyr)
library(tidyr)
library(ggplot2)


PCs   <-read.csv("C:\\Users\\Kim\\Documents\\Romberg_eQTL\\ImputationII_withX\\all_chr_LDpruned.eigenvec", header=T, sep="")

ggplot(data=PCs, aes(x=PC1, y=PC2)) + geom_point(size=3, alpha = 0.5) + xlab("PC1") + ylab("PC2")
ggsave("C:\\Users\\Kim\\Documents\\Romberg_eQTL\\ImputationII_withX\\final_LDpruned_PC12.png", device = "png", units = "in", dpi = 600)


ggplot(data=PCs, aes(x=PC3, y=PC4)) + geom_point(size=3, alpha = 0.5) + xlab("PC3") + ylab("PC4")
ggsave("C:\\Users\\Kim\\Documents\\Romberg_eQTL\\ImputationII_withX\\final_LDpruned_PC34.png", device = "png", units = "in", dpi = 600)
