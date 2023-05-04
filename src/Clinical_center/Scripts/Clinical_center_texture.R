rm(list=ls())

print("Start Clinical_center/Scripts/Clinical_center_texture.R")

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

# Rename
colnames(tcga_kirc)[grep(pattern = "^texture_[[:print:]]*%", colnames(tcga_kirc))] <- c("Blood", "Cancer", "Normal", "Stroma", "Other")
textures <- c("Blood", "Cancer", "Normal", "Stroma", "Other")


tcga_kirc %>%
  group_by(ClinicalCenter) %>%
  summarise(m = median(Cancer)) %>%
  arrange(m)
tcga_kirc %>%
  group_by(ClinicalCenter) %>%
  summarise(m = median(Normal)) %>%
  arrange(m)
sapply(tcga_kirc[textures], quantile)

tcga_kirc %>%
  filter(Cancer>50) %>%
  nrow()


# For loop
for (texture in textures) {
  g <- ggplot(tcga_kirc, aes_string(x="ClinicalCenter", y=texture)) +
    geom_jitter(size=3, width = 0.2, aes(fill=ClinicalCenter), shape = 21, color = "black") +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    labs(x="Clinical Center", y=paste0("Proportion of ", texture, " texture (%)")) +
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
  # stat_compare_means(method = "wilcox.test",
  #                    # label = "p.signif",
  #                    comparisons = my_comparisons,
  #                    # label.y = c(1.13*max(panels_prop_merged_modified2$hlaabc, na.rm = TRUE), 1.10*max(panels_prop_merged_modified2$hlaabc, na.rm = TRUE),
  #                    #             1.07*max(panels_prop_merged_modified2$hlaabc, na.rm = TRUE), 1.04*max(panels_prop_merged_modified2$hlaabc, na.rm = TRUE)),
  #                    label.y = c(a1*max(panels_prop_merged_modified2$hlaabc, na.rm = TRUE), a2*max(panels_prop_merged_modified2$hlaabc, na.rm = TRUE),
  #                                a3*max(panels_prop_merged_modified2$hlaabc, na.rm = TRUE), a4*max(panels_prop_merged_modified2$hlaabc, na.rm = TRUE)),
  #                    size = 6)
  ggsave(plot = g,
         filename = paste0("Clinical_center/Images/TCGA_texture_", texture, ".png"),
         width = 7, height = 7,
         units = "in",
         dpi = 300)
}



# Heatmap
# a <- tcga_kirc %>% dplyr::group_by(ClinicalCenter) %>% summarise(n = n()) %>% filter(n>=5)
tcga_kirc_hm <- tcga_kirc %>%
  dplyr::select("ClinicalCenter", "Blood", "Cancer", "Normal", "Stroma", "Other") %>%
  dplyr::group_by(ClinicalCenter) %>%
  summarise_all(
    funs(median(.))
  )


# Annotation data
anno <- data.frame(Center=tcga_kirc_hm$ClinicalCenter)

## Dark2 colors
### #1b9e77 green #d95f02 orange #7570b3 violet #e7298a pink #66a61e light green #e6ab02 yellow #a6761d brown #666666 grey
## Set1 colors
### #e41a1c #377eb8 #4daf4a #984ea3 #ff7f00 #ffff33 #a65628 #f781bf #999999
## Paired colors
#a6cee3, #1f78b4, #b2df8a, #33a02c, #fb9a99, #e31a1c, #fdbf6f, #ff7f00, #cab2d6, #6a3d9a

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
                             annotation_name_gp = gpar(fontface="bold", fontsize = 16),
                             annotation_legend_param = list(
                               title_gp = gpar(fontsize = 10),
                               ncol = 2,
                               labels_gp = gpar(fontsize = 10),
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

## Remove columns with more than 50% NA
# scaled_hm <- scaled_hm[-which(rowMeans(is.na(scaled_hm)) > 0.50),]


# Plot
png("Clinical_center/Images/TCGA_texture_heatmap.png", width = 7, height = 4, units = 'in', res = 300, pointsize = 12) #original pointsize = 12
# Heatmap(t(tcga_kirc_hm1),
g <- Heatmap(t(tcga_kirc_hm1),
             top_annotation = hm_anno,
             col=rev(cols),
             na_col = "lightgrey",
             column_dend_height = unit(1,"cm"),
             column_dend_gp = grid::gpar(lwd = 2),
             row_dend_width = unit(1,"cm"),
             row_dend_gp = grid::gpar(lwd = 2),
             name = "Scale",
             # row_names_gp = gpar(fontsize=5, fontface="bold"),
             row_names_gp = gpar(fontsize = 16, fontface = "bold"),  # fontsize=7,
             cluster_columns = FALSE,
             clustering_method_rows = "ward.D2",
             clustering_distance_rows = "spearman")

# print(g)
draw(g, heatmap_legend_side="left", annotation_legend_side="bottom")
dev.off()

