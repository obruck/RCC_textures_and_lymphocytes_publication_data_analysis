rm(list=ls())

print("Start Ly/Scripts/Ly_PBRM1_Nature_Medicine_study.R")

# Load libraries
library(tidyverse)
library(readxl)
library(ggrepel)
library(ggpubr)
library(ggplot2)

# Load data
a <- read_xlsx("./../data/pbrm1_nature_medicine_study.xlsx") %>%
  filter(!is.na(PBRM1)) %>%
  filter(!PBRM1 == "NA") %>%
  mutate(TC_CD8_Density = as.numeric(TC_CD8_Density),
         PBRM1 = factor(PBRM1) )

# Plot
g <- ggplot(a, aes(x = PBRM1, y = TC_CD8_Density)) +
  geom_jitter(size=5, width = 0.2, aes(fill=PBRM1), shape = 21, color = "black") +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
  labs(x="PBRM1", y="Intratumoral CD8+ Density") +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, colour = "black"),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title.y = element_text(size=14, face="bold", colour = "black"),
        axis.title.x = element_text(size=14, face="bold", colour = "black"),
        legend.position = "none") +
  scale_fill_brewer(palette = c("Set1")) +
  stat_compare_means(method = "wilcox.test",
                     # label = "p.signif",
                     label.x = 1.25,
                     label.y = 1.05*max(a$TC_CD8_Density, na.rm = TRUE),
                     size = 5)
ggsave(plot = g, filename = "./Ly/Images/Mut/Ly_pbrm1_nature_medicine.png", width = 5, height = 5, units = 'in', dpi = 300)


