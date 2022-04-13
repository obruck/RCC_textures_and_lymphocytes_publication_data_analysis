rm(list=ls())

print("Start Clinical_center/Scripts/Clinical_center_lymphocyte.R")

# Load libraries
library(jsonlite)
library(readxl)
library(tidyverse)
library(data.table)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)


############################# LOAD DATA ##########################################################################################################


# Image data
tcga_kirc <- read_xlsx("../data/image_analysis_results_final.xlsx")


# Normalize by removing empty
tcga_kirc$`texture_cancer_%` <- 100*tcga_kirc$texture_cancer / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)


# Add TCGA id
tcga_kirc <- tcga_kirc %>%
  dplyr::mutate(tissue_source_site = gsub("-[[:print:]]{4}", "", gsub("TCGA-", "", tcga_id)))


# TCGA clinical sources
clinical_sources <- read_xlsx("../data/TCGA_source_sites.xlsx")


############################# JOIN ##########################################################################################################


# Join Otso's data and TCGA centers
tcga_kirc <- tcga_kirc %>%
  dplyr::left_join(clinical_sources %>% dplyr::select(-StudyName, -BCR))


############################# PLOT ##########################################################################################################


# Ly/Total %
tcga_kirc$bin_total <- tcga_kirc$bin_lymphocytes_blood + tcga_kirc$bin_lymphocytes_cancer + tcga_kirc$bin_lymphocytes_normal +
  tcga_kirc$bin_lymphocytes_stroma + tcga_kirc$bin_lymphocytes_other
tcga_kirc$Ly_Total <- tcga_kirc$bin_total / tcga_kirc$texture_total


# Rename
colnames(tcga_kirc)[grep(pattern = "^inf_bin_lymphocytes_[[:print:]]*", colnames(tcga_kirc))] <- c("Ly_Blood", "Ly_Cancer", "Ly_Normal", "Ly_Stroma", "Ly_Other")
textures <- c("Ly_Total", "Ly_Blood", "Ly_Cancer", "Ly_Normal", "Ly_Stroma", "Ly_Other")
tcga_kirc[textures] <- tcga_kirc[textures]*100


# For loop
for (texture in textures) {
  # texture = textures[1]
  g <- ggplot(tcga_kirc, aes_string(x="ClinicalCenter", y=texture)) +
    geom_jitter(size=3, width = 0.2, aes(fill=ClinicalCenter), shape = 21, color = "black") +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    labs(x="Clinical Center", y=paste0("Proportion of ", gsub("_", "/", texture), " (%)")) +
    theme_bw() +
    theme(axis.text.x = element_text(size=12, colour = "black", angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(size=12, colour = "black"),
          axis.title=element_text(size=14, face="bold", colour = "black"),
          legend.position = "none",
          plot.margin = margin(0.1, 0.1, 0.1, 0.3, "in")) +
    # scale_fill_brewer(palette = c("Set1")) +
    stat_compare_means(method = "kruskal.test",
                       # label.y = max(tcga_kirc[texture], na.rm = TRUE),
                       label.x = 4.75,
                       # label.y = a5,
                       size = 5)
  ggsave(plot = g,
         filename = paste0("Clinical_center/Images/TCGA_ly_", texture, ".png"),
         width = 7, height = 7,
         units = "in",
         dpi = 300)
}








# Heatmap
tcga_kirc_hm <- tcga_kirc %>%
  dplyr::select("ClinicalCenter", "Ly_Blood", "Ly_Cancer", "Ly_Normal", "Ly_Stroma", "Ly_Other") %>%
  dplyr::group_by(ClinicalCenter) %>%
  summarise_all(
    funs(median(.))
  )


# Annotation data
anno <- data.frame(Center=tcga_kirc_hm$ClinicalCenter)

## Annotation colors 
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(length(unique(anno$Center)))
col_anno1 = list(Center = structure(mycolors))
col_anno1 <- setNames(object = col_anno1$Center, nm = unique(anno$Center))
col_anno <- list()
col_anno$Center <- col_anno1

## Heatmap annotation
hm_anno <- HeatmapAnnotation(df = anno,
                             name = "Top Annotation",
                             col = col_anno,
                             show_annotation_name = T,
                             gp = gpar(col = "black"),
                             height = unit(2, "cm"),
                             # annotation_name_gp = gpar(fontsize=8, fontface="bold"),
                             annotation_name_gp = gpar(fontface="bold"),
                             annotation_legend_param = list(
                               title_gp = gpar(fontsize = 13),
                               labels_gp = gpar(fontsize = 13),
                               grid_height = unit(5, "mm"),
                               border=TRUE),
                             na_col = "white")


## Remove unnecessary metadata
tcga_kirc_hm <- tcga_kirc_hm %>%
  dplyr::select(-ClinicalCenter) %>%
  as.matrix()


################################# Colors and scaling #############################################################################


# Color ColorBrewerilla
cols <- colorRampPalette(brewer.pal(11, "RdBu"))(256)

# Scale data using median as the center and the maximum value as the scale = (x-median)/max value
tcga_kirc_hm1 <- scale(tcga_kirc_hm, center = apply(tcga_kirc_hm, 2, median, na.rm=TRUE), scale=apply(tcga_kirc_hm, 2, max, na.rm=TRUE))
tcga_kirc_hm1[is.nan(tcga_kirc_hm1)] <- NA


# Plot
png("Clinical_center/Images/TCGA_texture_heatmap.png", width = 12, height = 5, units = 'in', res = 300, pointsize = 12) #original pointsize = 12
g <- Heatmap(t(tcga_kirc_hm1),
             # Heatmap(t(tcga_kirc_hm),
             top_annotation = hm_anno,
             col=rev(cols),
             na_col = "lightgrey",
             column_dend_height = unit(1,"cm"),
             column_dend_gp = grid::gpar(lwd = 2),
             row_dend_width = unit(1,"cm"),
             row_dend_gp = grid::gpar(lwd = 2),
             name = "Scale",
             # row_names_gp = gpar(fontsize=5, fontface="bold"),
             row_names_gp = gpar(fontface = "bold"),  # fontsize=7,
             cluster_columns = FALSE,
             clustering_method_rows = "ward.D2",
             clustering_distance_rows = "spearman")

print(g)
dev.off()

