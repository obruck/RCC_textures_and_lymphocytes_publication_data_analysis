rm(list=ls())

print("Start Textures/Scripts/Texture/Texture_by_sample_normal_proportion.R")

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
library(rstatix)


############################# LOAD DATA ##########################################################################################################


# Image data
tcga_kirc <- read_xlsx("../data/image_analysis_results_final.xlsx")


# Normalize by removing empty
tcga_kirc$`texture_blood_%` <- 100*tcga_kirc$texture_blood / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`texture_cancer_%` <- 100*tcga_kirc$texture_cancer / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`texture_normal_%` <- 100*tcga_kirc$texture_normal / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`texture_stroma_%` <- 100*tcga_kirc$texture_stroma / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`texture_other_%` <- 100*tcga_kirc$texture_other / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)

# Rename
colnames(tcga_kirc)[grep(pattern = "^inf_bin_lymphocytes_[[:print:]]*", colnames(tcga_kirc))] <- c("Ly_Blood", "Ly_Cancer", "Ly_Normal", "Ly_Stroma", "Ly_Other")
textures <- c("Ly_Blood", "Ly_Cancer", "Ly_Stroma", "Ly_Other")

# Proportions
# Normalize in all textures (= values to equal 100%)
tcga_kirc$Ly_Blood <- 100*tcga_kirc$Ly_Blood / (tcga_kirc$Ly_Blood + tcga_kirc$Ly_Cancer + tcga_kirc$Ly_Normal + tcga_kirc$Ly_Stroma + tcga_kirc$Ly_Other)
tcga_kirc$Ly_Cancer <- 100*tcga_kirc$Ly_Cancer / (tcga_kirc$Ly_Blood + tcga_kirc$Ly_Cancer + tcga_kirc$Ly_Normal + tcga_kirc$Ly_Stroma + tcga_kirc$Ly_Other)
tcga_kirc$Ly_Normal <- 100*tcga_kirc$Ly_Normal / (tcga_kirc$Ly_Blood + tcga_kirc$Ly_Cancer + tcga_kirc$Ly_Normal + tcga_kirc$Ly_Stroma + tcga_kirc$Ly_Other)
tcga_kirc$Ly_Stroma <- 100*tcga_kirc$Ly_Stroma / (tcga_kirc$Ly_Blood + tcga_kirc$Ly_Cancer + tcga_kirc$Ly_Normal + tcga_kirc$Ly_Stroma + tcga_kirc$Ly_Other)
tcga_kirc$Ly_Other <- 100*tcga_kirc$Ly_Other / (tcga_kirc$Ly_Blood + tcga_kirc$Ly_Cancer + tcga_kirc$Ly_Normal + tcga_kirc$Ly_Stroma + tcga_kirc$Ly_Other)

a <- 100/(tcga_kirc$Ly_Blood + tcga_kirc$Ly_Cancer + tcga_kirc$Ly_Normal + tcga_kirc$Ly_Stroma + tcga_kirc$Ly_Other)
tcga_kirc[textures] <- sapply(tcga_kirc[textures], function(x) x*a)
 
tcga_kirc <- tcga_kirc %>%
  dplyr::filter(`texture_cancer_%` > 5) %>%
  dplyr::mutate(tissue_source_site = gsub("-[[:print:]]{4}", "", gsub("TCGA-", "", tcga_id)))
tcga_kirc1 <- read_xlsx("../data/TCGA_source_sites.xlsx")
tcga_kirc <- tcga_kirc %>% dplyr::left_join(tcga_kirc1)
a <- tcga_kirc %>% group_by(ClinicalCenter) %>% summarise(me = median(Ly_Cancer))


############################# PLOT ##########################################################################################################


# Rename
colnames(tcga_kirc)[grep(pattern = "^texture_[[:print:]]*%", colnames(tcga_kirc))] <- c("Blood", "Cancer", "Normal", "Stroma", "Other")
textures <- c("Blood", "Cancer", "Stroma", "Other")

# Divide by normal tissue proportion
tcga_kirc <- tcga_kirc %>%
  dplyr::mutate(
    Normal_Present = ifelse(Normal < "1", "No", "Yes")
  )

quantile(tcga_kirc[tcga_kirc$Normal_Present=="No",]$Cancer)
quantile(tcga_kirc[tcga_kirc$Normal_Present=="Yes",]$Cancer)

# Plot
for (texture1 in textures) {
  a <- max(tcga_kirc[[texture1]], na.rm=TRUE)
  g <- ggplot(tcga_kirc, aes_string(x = "Normal_Present", y = texture1)) +
    geom_jitter(size=5, width = 0.2, aes(fill=Normal_Present), shape = 21, color = "black") +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
    labs(y=paste0(str_to_sentence(texture1), " proportion (%)"), x="Normal tissue in sample") +
    theme_bw() +
    theme(axis.text.x = element_text(size=12, colour = "black"),
          axis.text.y = element_text(size=12, colour = "black"),
          axis.title=element_text(size=14, face="bold", colour = "black"),
          legend.position = "none") +
    scale_fill_brewer(palette = c("Set1")) +
    stat_compare_means(method = "wilcox.test",
                       comparisons = list(c("Yes", "No")),
                       label.y = 1.0*a,
                       label.x = 1.5,
                       label = "p.signif",
                       bracket.size = 1.5,
                       size = 5)
  ggsave(plot = g, filename = paste0("Textures/Images/Texture_normal/Scatter_", texture1, "_normal_status.png"), width = 5, height = 5, units = 'in', dpi = 300, pointsize = 12) #original pointsize = 12
  
}
