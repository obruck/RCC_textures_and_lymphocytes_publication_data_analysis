rm(list=ls())

print("Start Textures/Scripts/Margin/Margin_Ly_Margin_Cancer_correlation.R")

# Load libraries
library(readxl)
library(tidyverse)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(corrplot)
library(psych)
library(ggpubr)


############################# LOAD DATA ##########################################################################################################

# Image data
tcga_kirc <- read_xlsx("../data/image_analysis_results_final.xlsx")

# Normalize by removing empty
tcga_kirc$`texture_cancer_%` <- 100*tcga_kirc$texture_cancer / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`texture_normal_%` <- 100*tcga_kirc$texture_normal / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`margin_texture_blood_%` <- 100*tcga_kirc$margin_texture_blood / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)
tcga_kirc$`margin_texture_normal_%` <- 100*tcga_kirc$margin_texture_normal / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)
tcga_kirc$`margin_texture_stroma_%` <- 100*tcga_kirc$margin_texture_stroma / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)
tcga_kirc$`margin_texture_other_%` <- 100*tcga_kirc$margin_texture_other / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)

tcga_kirc <- tcga_kirc %>%
  dplyr::filter(`texture_cancer_%` > 5) %>%
  dplyr::filter(`texture_normal_%` >= 1) %>%
  dplyr::mutate(tissue_source_site = gsub("-[[:print:]]{4}", "", gsub("TCGA-", "", tcga_id)))


############################# PLOT ##########################################################################################################


## TEXTURES
## Rename
colnames(tcga_kirc)[grep(pattern = "^inf_bin_lymphocytes_", colnames(tcga_kirc))] <- c("Ly/Blood", "Ly/Cancer", "Ly/Normal", "Ly/Stroma", "Ly/Other")
colnames(tcga_kirc)[grep(pattern = "^inf_margin_bin_lymphocytes_", colnames(tcga_kirc))] <- c("Ly/Blood_Margin", "Ly/Normal_Margin", "Ly/Stroma_Margin", "Ly/Other_Margin")
colnames(tcga_kirc)[grep(pattern = "^margin_texture_[[:print:]]*%", colnames(tcga_kirc))] <- c("Blood_Margin", "Cancer_Margin", "Normal_Margin", "Stroma_Margin", "Other_Margin")
tcga_kirc$`Ly/Total` = (tcga_kirc$bin_lymphocytes_blood + tcga_kirc$bin_lymphocytes_normal + tcga_kirc$bin_lymphocytes_other + tcga_kirc$bin_lymphocytes_stroma + tcga_kirc$bin_lymphocytes_cancer) / (tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_blood + tcga_kirc$texture_other)
tcga_kirc$`Ly/Total_Margin` = (tcga_kirc$margin_bin_lymphocytes_blood + tcga_kirc$margin_bin_lymphocytes_normal + tcga_kirc$margin_bin_lymphocytes_stroma + tcga_kirc$margin_bin_lymphocytes_other) / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)
textures <- c("Blood_Margin", "Normal_Margin", "Stroma_Margin", "Other_Margin", "Ly/Total", "Ly/Cancer", "Ly/Blood", "Ly/Normal", "Ly/Stroma", "Ly/Other", "Ly/Total_Margin", "Ly/Blood_Margin", "Ly/Normal_Margin", "Ly/Stroma_Margin", "Ly/Other_Margin")


## Quantity
## Normalize in all textures (= values to equal 100%)
tcga_kirc[textures] <- sapply(tcga_kirc[textures], function(x) x*100)  # 100 = 100%, 5 = 5 textures so max is 500% which looks weird in the plot

# Group by infiltration
tcga_kirc$Infiltration = ifelse(tcga_kirc$`Ly/Cancer` > quantile(tcga_kirc$`Ly/Cancer`)[4], "High", "Low")


# Plot
g <- ggplot(tcga_kirc %>% dplyr::rename(Ly_Cancer = "Ly/Cancer", Ly_Total_Margin = "Ly/Total_Margin"), aes(x = Ly_Cancer, y = Ly_Total_Margin)) +
  geom_point(size=5, color = "black") +
  geom_smooth(formula = y ~ x + 0, method='lm', fullrange = T, size=2, color="#984ea3") +
  stat_cor(size = 5, color="#984ea3", p.accuracy = 0.001, r.accuracy = 0.01, method = "spearman", label.x = 10, label.y = 60, show.legend = FALSE) +
  # geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x="Lymphocytes in Cancer (%)", y="Lymphocytes in Margin (%)") +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0), limits = c(-1, 80)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-1, 80)) +
  theme(axis.text.x = element_text(size=12, colour = "black"),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title.y = element_text(size=12, face="bold", colour = "black"),
        axis.title.x = element_text(size=12, face="bold", colour = "black"),
        legend.title = element_text(size=12, face="bold", colour = "black"),
        legend.text = element_text(size=12, colour = "black"),
        legend.position = "bottom",
        legend.box = "vertical",
        # legend.key.size = unit(0.1, "cm"),
        legend.margin = margin(),
        legend.spacing.y = unit(0.1, "cm"))
ggsave(plot = g, filename = paste0("Textures/Images/Margin/TCGA_ly_correlation_margin_cancer1.png"), width = 6, height = 6, units = 'in', dpi = 300)

g <- ggplot(tcga_kirc %>% dplyr::rename(Ly_Cancer = "Ly/Cancer", Ly_Total_Margin = "Ly/Total_Margin"), aes(x = Ly_Cancer, y = Ly_Total_Margin)) +
  geom_point(size=5, color = "black") +
  geom_smooth(formula = y ~ x + 0, method='lm', aes(color=Infiltration), fullrange = T, size=2) +
  stat_cor(aes(color=Infiltration), size = 5, p.accuracy = 0.001, r.accuracy = 0.01, method = "spearman", label.x = 10, label.y = c(60, 56), show.legend = FALSE) +
  # geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x="Lymphocytes in Cancer (%)", y="Lymphocytes in Margin (%)") +
  theme_bw() +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(expand = c(0, 0), limits = c(-1, 80)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-1, 80)) +
  theme(axis.text.x = element_text(size=12, colour = "black"),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title.y = element_text(size=12, face="bold", colour = "black"),
        axis.title.x = element_text(size=12, face="bold", colour = "black"),
        legend.title = element_text(size=12, face="bold", colour = "black"),
        legend.text = element_text(size=12, colour = "black"),
        legend.position = "bottom",
        legend.box = "vertical",
        # legend.key.size = unit(0.1, "cm"),
        legend.margin = margin(),
        legend.spacing.y = unit(0.1, "cm"))
ggsave(plot = g, filename = paste0("Textures/Images/Margin/TCGA_ly_correlation_margin_cancer2.png"), width = 6, height = 6, units = 'in', dpi = 300)

