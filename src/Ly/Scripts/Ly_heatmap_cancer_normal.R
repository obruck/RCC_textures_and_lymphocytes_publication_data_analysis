rm(list=ls())

print("Start Ly/Scripts/Ly_heatmap_cancer_normal.R")

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


# Load data
cancer_gexp_pathways0 <- readxl::read_xlsx("Ly/Images/Gexp/Ly_Cancer_gexp_pathways.xlsx")
normal_gexp_pathways0 <- readxl::read_xlsx("Ly/Images/Gexp/Ly_Normal_gexp_pathways.xlsx")

# Filter by p
cancer_gexp_pathways <- cancer_gexp_pathways0 %>% 
  dplyr::filter(padj<0.10) %>%
  dplyr::mutate(Cancer = ifelse(is.na(NES), 0, NES))
normal_gexp_pathways <- normal_gexp_pathways0 %>%
  dplyr::filter(padj<0.10) %>%
  dplyr::mutate(Normal = ifelse(is.na(NES), 0, NES))

# Rename variables
colnames(normal_gexp_pathways)[2:(ncol(normal_gexp_pathways))] <- paste0(colnames(normal_gexp_pathways)[2:(ncol(normal_gexp_pathways))], "_normal")
normal_gexp_pathways = normal_gexp_pathways %>% dplyr::rename(Normal = Normal_normal)
colnames(cancer_gexp_pathways)[2:(ncol(cancer_gexp_pathways))] <- paste0(colnames(cancer_gexp_pathways)[2:(ncol(cancer_gexp_pathways))], "_cancer")
cancer_gexp_pathways = cancer_gexp_pathways %>% dplyr::rename(Cancer = Cancer_cancer)

# Select top 5 and last 5 rows by NES value
cancer_gexp_pathways1 <- cancer_gexp_pathways %>%
  dplyr::filter(!is.na(Cancer)) %>%
  dplyr::arrange(NES_cancer) %>%
  slice(c(1:5, (n()-4):n()))
normal_gexp_pathways1 <- normal_gexp_pathways %>%
  dplyr::filter(!is.na(Normal)) %>%
  dplyr::arrange(NES_normal) %>%
  slice(c(1:5, (n()-4):n()))
# Do not remove pathways from cancer_gexp_pathways if these are adj<0.10 in normal_gexp_pathways and vice versa
cancer_gexp_pathways2 <- cancer_gexp_pathways %>%
  dplyr::filter(!is.na(Cancer)) %>%
  dplyr::filter(pathway %in% normal_gexp_pathways1$pathway)
normal_gexp_pathways2 <- normal_gexp_pathways %>%
  dplyr::filter(!is.na(Normal)) %>%
  dplyr::filter(pathway %in% cancer_gexp_pathways1$pathway)

# Join
gexp_pathways <- cancer_gexp_pathways1 %>%
  full_join(cancer_gexp_pathways2) %>%
  full_join(normal_gexp_pathways1) %>% 
  full_join(normal_gexp_pathways2)
gexp_pathways <- gexp_pathways %>%
  group_by(pathway) %>%
  summarise_all(funs(max), na.rm = TRUE)
gexp_pathways[sapply(gexp_pathways, is.infinite)] <- 0

# Make sure the other is at least adjp<0.05
gexp_pathways <- gexp_pathways %>%
  dplyr::filter(!(padj_cancer==0 & padj_normal>0.05)) %>%
  dplyr::filter(!(padj_normal==0 & padj_cancer>0.05))


# Fill NES with 0 if NA
# gexp_pathways <- gexp_pathways %>%
#   dplyr::mutate(Cancer = ifelse(is.na(NES_cancer), 0, NES_cancer),
#                 Normal = ifelse(is.na(NES_normal), 0, NES_normal))

# Filter pathways that are padj >0.05 but padj<0.10
# normal_gexp_pathways1 <- normal_gexp_pathways0 %>%
#   dplyr::filter(padj < 0.10 & padj >= 0.05) %>%
#   dplyr::select(pathway)
# cancer_gexp_pathways1 <- cancer_gexp_pathways0 %>%
#   dplyr::filter(padj < 0.10 & padj >= 0.05) %>%
#   dplyr::select(pathway)
# exclude1 <- gexp_pathways %>% dplyr::filter(Cancer == 0) %>% dplyr::select(pathway) %>% dplyr::filter(pathway %in% cancer_gexp_pathways1$pathway)
# exclude2 <- gexp_pathways %>% dplyr::filter(Normal == 0) %>% dplyr::select(pathway) %>% dplyr::filter(pathway %in% normal_gexp_pathways1$pathway)
# gexp_pathways <- gexp_pathways %>% dplyr::filter(!pathway %in% c(exclude1$pathway, exclude2$pathway))

# Transpose
gexp_pathways <- gexp_pathways %>% 
  dplyr::select(pathway, Cancer, Normal) %>%
  pivot_longer(cols = c(Cancer, Normal), names_to = "col1", values_to = "col2") %>% 
  pivot_wider(names_from = "pathway", values_from = "col2") %>%
  dplyr::rename(Texture = col1)

# Identify FGSEA type
a <- colnames(gexp_pathways)[2:ncol(gexp_pathways)]
a = ifelse(str_detect(a, "HALLMARK"), "HALLMARK",
           ifelse(str_detect(a, "KEGG"), "KEGG",
                  ifelse(str_detect(a, "PID"), "PID",
                         ifelse(str_detect(a, "REACTOME"), "REACTOME", 
                                ifelse(str_detect(a, "chr"), "CHROMOSOME", NA)))))


# Rename colnames
colnames(gexp_pathways)[2:ncol(gexp_pathways)] <- gsub("HALLMARK_", "", colnames(gexp_pathways)[2:ncol(gexp_pathways)])
colnames(gexp_pathways)[2:ncol(gexp_pathways)] <- gsub("KEGG_", "", colnames(gexp_pathways)[2:ncol(gexp_pathways)])
colnames(gexp_pathways)[2:ncol(gexp_pathways)] <- gsub("PID_", "", colnames(gexp_pathways)[2:ncol(gexp_pathways)])
colnames(gexp_pathways)[2:ncol(gexp_pathways)] <- gsub("CYTOKINE_CYTOKINE_RECEPTOR", "CYTOKINE_RECEPTOR", colnames(gexp_pathways)[2:ncol(gexp_pathways)])
colnames(gexp_pathways)[2:ncol(gexp_pathways)] <- gsub("REACTOME_", "", colnames(gexp_pathways)[2:ncol(gexp_pathways)])
colnames(gexp_pathways)[2:ncol(gexp_pathways)] <- gsub("CHROMOSOME_", "", colnames(gexp_pathways)[2:ncol(gexp_pathways)])
colnames(gexp_pathways)[2:ncol(gexp_pathways)] <- gsub("TRIACYLGLYCEROL_AND_KETONE_BODY", "TAG & KETONE", colnames(gexp_pathways)[2:ncol(gexp_pathways)])
colnames(gexp_pathways)[2:ncol(gexp_pathways)] <- gsub("VALINE", "VAL", colnames(gexp_pathways)[2:ncol(gexp_pathways)])
colnames(gexp_pathways)[2:ncol(gexp_pathways)] <- gsub("ISOLEUCINE", "ILE", colnames(gexp_pathways)[2:ncol(gexp_pathways)])
colnames(gexp_pathways)[2:ncol(gexp_pathways)] <- gsub("LEUCINE", "LEU", colnames(gexp_pathways)[2:ncol(gexp_pathways)])
colnames(gexp_pathways)[2:ncol(gexp_pathways)] <- gsub("EPITHELIAL_MESENCHYMAL", "E-M", colnames(gexp_pathways)[2:ncol(gexp_pathways)])
colnames(gexp_pathways)[2:ncol(gexp_pathways)] <- gsub("INTERFERON GAMMA", "IFNG", colnames(gexp_pathways)[2:ncol(gexp_pathways)])
colnames(gexp_pathways)[2:ncol(gexp_pathways)] <- gsub("ATP_SYNTHESIS_BY[[:PRINT:]]*", "ATP_SYNTHESIS", colnames(gexp_pathways)[2:ncol(gexp_pathways)])
colnames(gexp_pathways)[2:ncol(gexp_pathways)] <- gsub("RESPIRATORY", "RESP", colnames(gexp_pathways)[2:ncol(gexp_pathways)])
colnames(gexp_pathways)[2:ncol(gexp_pathways)] <- gsub("ELECTRON", "E-", colnames(gexp_pathways)[2:ncol(gexp_pathways)])
colnames(gexp_pathways)[2:ncol(gexp_pathways)] <- gsub("TRANSPORT", "TRANSP", colnames(gexp_pathways)[2:ncol(gexp_pathways)])
colnames(gexp_pathways)[2:ncol(gexp_pathways)] <- gsub("_", " ", colnames(gexp_pathways)[2:ncol(gexp_pathways)])


# Annotation data
anno <- data.frame(Texture=gexp_pathways$Texture)
anno_row <- data.frame(GSEA = a)
col_anno = list(Texture = c(
  "Cancer" = "#377eb8",
  "Normal" = "#4daf4a"))
unique(a)
col_anno_row = list(GSEA = c(
  "HALLMARK"="#e41a1c",
  "KEGG" = "#984ea3",
  "PID" = "#ff7f00",
  "REACTOME" = "#ffff33",
  "CHROMOSOME" = "#a65628"))

## Heatmap annotation
hm_anno <- HeatmapAnnotation(df = anno,
                             name = "Top Annotation",
                             col = col_anno,
                             show_annotation_name = T,
                             gp = gpar(col = "black"),
                             height = unit(2, "cm"),
                             show_legend=FALSE,
                             # annotation_name_gp = gpar(fontsize=8, fontface="bold"),
                             annotation_name_gp = gpar(fontface="bold"),
                             annotation_legend_param = list(
                               title_gp = gpar(fontsize = 13),
                               labels_gp = gpar(fontsize = 13),
                               grid_height = unit(5, "mm"),
                               border=TRUE),
                             na_col = "white")
hm_anno_row <- rowAnnotation(df = anno_row,
                             name = "Row Annotation",
                             col = col_anno_row,
                             show_annotation_name = T,
                             gp = gpar(col = "black"),
                             width = unit(2, "cm"),
                             # annotation_name_gp = gpar(fontsize=8, fontface="bold"),
                             annotation_name_gp = gpar(fontsize = 0, fontface="bold"),
                             annotation_legend_param = list(
                               title_gp = gpar(fontsize = 13, fontface="bold"),
                               labels_gp = gpar(fontsize = 13),
                               grid_height = unit(5, "mm"),
                               border=TRUE),
                             na_col = "white")

## Remove unnecessary metadata
gexp_pathways_hm <- gexp_pathways %>%
  dplyr::select(-Texture) %>%
  as.matrix()


################################# Colors and scaling #############################################################################


# Color ColorBrewerilla
cols <- colorRampPalette(brewer.pal(11, "RdBu"))(256)


# Plot
png("Ly/Images/Gexp/Cancer_vs_normal_heatmap.png", width = 7, height = 4.5, units = 'in', res = 300, pointsize = 12) #original pointsize = 12
g <- Heatmap(t(gexp_pathways_hm),
             width = ncol(gexp_pathways_hm)*unit(0.12, "lines"), 
             top_annotation = hm_anno,
             right_annotation = hm_anno_row,
             col=rev(cols),
             na_col = "lightgrey",
             column_dend_height = unit(1,"cm"),
             column_dend_gp = grid::gpar(lwd = 2),
             row_dend_width = unit(1,"cm"),
             row_dend_gp = grid::gpar(lwd = 2),
             rect_gp = gpar(col = "white", lwd = 2),
             name = "Scale",
             # row_names_gp = gpar(fontsize=5, fontface="bold"),
             row_names_gp = gpar(fontface = "bold"),  # fontsize=7,
             cluster_columns = FALSE,
             clustering_method_rows = "ward.D2",
             clustering_distance_rows = "euclidean")

# print(g)

draw(g, annotation_legend_side="left", heatmap_legend_side="left")
dev.off()


