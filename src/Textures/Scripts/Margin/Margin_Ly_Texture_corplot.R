rm(list=ls())

print("Start Textures/Scripts/Margin/Margin_Ly_Texture_corplot.R")

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

tcga_kirc1 <- tcga_kirc


############################# PLOT ##########################################################################################################



## TEXTURES
## Rename
colnames(tcga_kirc)[grep(pattern = "^inf_margin_bin_lymphocytes_", colnames(tcga_kirc))] <- c("Ly/Blood", "Ly/Normal", "Ly/Stroma", "Ly/Other")
colnames(tcga_kirc)[grep(pattern = "^margin_texture_[[:print:]]*%", colnames(tcga_kirc))] <- c("Blood", "Cancer", "Normal", "Stroma", "Other")
tcga_kirc$`Ly/Total` = tcga_kirc$`Ly/Blood` + tcga_kirc$`Ly/Normal` + tcga_kirc$`Ly/Stroma` + tcga_kirc$`Ly/Other`
textures <- c("Blood", "Normal", "Stroma", "Other", "Ly/Total", "Ly/Blood", "Ly/Normal", "Ly/Stroma", "Ly/Other")


## Quantity
## Normalize in all textures (= values to equal 100%)
tcga_kirc[textures] <- sapply(tcga_kirc[textures], function(x) x*100)  # 100 = 100%, 5 = 5 textures so max is 500% which looks weird in the plot



## Correlation matrix
cor_table <- corr.test(tcga_kirc[textures], adjust = "none", ci=F, method = "spearman")
r <- cor_table$r
p <- cor_table$p



# Plot
png("Textures/Images/Margin/TCGA_margin_texture_ly_corrplot.png", width = 6, height = 6, units = 'in', res = 300, pointsize = 12) #original pointsize = 12
corrplot(r %>% round(2),
         cl.ratio = 0.2,
         method="color",
         col= colorRampPalette(c("blue","white", "red"))(100),
         type="upper", 
         addrect=4,
         outline = TRUE,
         addgrid.col = "white",
         order="hclust",
         addCoef.col = "black", # Add coefficient of correlation
         hclust.method="ward.D2",
         cl.pos="r",  # b, l, r, t
         tl.cex=0.8,
         tl.srt=45, #Text label color and rotation
         p.mat = p,
         sig.level = 0.05,
         insig = "blank",
         diag=FALSE,   # hide correlation coefficient on the principal diagonal
         tl.col = "black")
dev.off()


