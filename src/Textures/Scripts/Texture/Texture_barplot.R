rm(list=ls())

print("Start Textures/Scripts/Texture/Texture_barplot.R")

# Load libraries
library(jsonlite)
library(readxl)
library(tidyverse)
library(data.table)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggpubr)


############################# LOAD DATA ##########################################################################################################


# Image data
tcga_kirc <- read_xlsx("../data/image_analysis_results_final.xlsx")


# Normalize by removing empty
tcga_kirc$`texture_blood_%` <- 100*tcga_kirc$texture_blood / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`texture_cancer_%` <- 100*tcga_kirc$texture_cancer / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`texture_normal_%` <- 100*tcga_kirc$texture_normal / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`texture_stroma_%` <- 100*tcga_kirc$texture_stroma / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`texture_other_%` <- 100*tcga_kirc$texture_other / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)

tcga_kirc <- tcga_kirc %>%
  dplyr::filter(`texture_cancer_%` > 5) %>%
  dplyr::mutate(tissue_source_site = gsub("-[[:print:]]{4}", "", gsub("TCGA-", "", tcga_id)))

  
############################# PLOT ##########################################################################################################


# Rename
colnames(tcga_kirc)[grep(pattern = "^texture_[[:print:]]*%", colnames(tcga_kirc))] <- c("Blood", "Cancer", "Normal", "Stroma", "Other")
textures <- c("Blood", "Cancer", "Normal", "Stroma", "Other")
# sapply(tcga_kirc[textures], quantile)


# Melt longer
tcga_kirc_long <- tcga_kirc %>%
  dplyr::select("tcga_id", "Blood", "Cancer", "Normal", "Stroma", "Other") %>%
  dplyr::rename("ID" = "tcga_id") %>%
  reshape2::melt() %>%
  dplyr::rename("Texture" = "variable") %>%
  arrange(desc(value))

# Order
tcga_kirc_long <- tcga_kirc_long %>%
  dplyr::mutate(orderid = ifelse(Texture=="Cancer", value, NA)) %>%
  dplyr::filter(!is.na(orderid)) %>%
  dplyr::select(ID, orderid) %>%
  dplyr::right_join(tcga_kirc_long)


# Plot
g <- ggplot(tcga_kirc_long, aes(x = reorder(ID, -orderid), y = value, fill = Texture)) +
  geom_bar(stat = "identity") +
  labs(y="Proportion (%)") +
  scale_y_continuous(limits = c(0,100.1), expand = c(0, 0)) +
  scale_fill_brewer(palette="Set1") +
  # scale_fill_manual(values=c("#4daf4a", "#e41a1c", "#a6cee3", "#984ea3", "#ffff33")) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title.y=element_text(size=14, face="bold", colour = "black"),
        axis.title.x=element_blank(),
        legend.title = element_text(size=14, face="bold", colour = "black"),
        legend.text = element_text(size=12, colour = "black"),
        legend.key = element_rect(color="black"),
        legend.position = "bottom")
ggsave(plot = g, filename = "Textures/Images/TCGA_texture_barplot.png", width = 18, height = 5, units = 'in', dpi = 300)


# Plot
## Plot parameters
a <- max(tcga_kirc_long$value, na.rm=TRUE)
## Adjust p values
pairwise.test = tcga_kirc_long %>%
  wilcox_test(value~Texture) %>%
  adjust_pvalue(method = 'BH')

g <- ggplot(tcga_kirc_long, aes(x = Texture, y = value)) +
  geom_jitter(size=5, width = 0.2, aes(fill=Texture), shape = 21, color = "black") +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
  labs(y="Texture proportion (%)", x="") +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, colour = "black"),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title.y = element_text(size=14, face="bold", colour = "black"),
        axis.title.x = element_text(size=14, face="bold", colour = "black"),
        legend.position = "none") +
  scale_fill_brewer(palette = c("Set1")) +
  scale_y_continuous(breaks = c(0, 50, 100), labels = c(0, 50, 100)) +
  stat_compare_means(method = "kruskal.test",
                     label.y = 1.67*a,
                     size = 5) +
  stat_pvalue_manual(
    pairwise.test,
    label = "p.adj.signif",
    bracket.size = 1.0,
    size = 5,
    y.position = c(1.21*a, 1.42*a, 1.56*a,
                   1.63*a, 1.14*a, 1.35*a,
                   1.49*a, 1.07*a, 1.28*a,
                   1.00*a)
  )
ggsave(plot = g, filename = "Textures/Images/TCGA_texture_scatterplot.png", width = 8, height = 5.5, units = 'in', dpi = 300)
