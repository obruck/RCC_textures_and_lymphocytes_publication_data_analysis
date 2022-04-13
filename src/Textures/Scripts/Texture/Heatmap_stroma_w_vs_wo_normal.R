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


# Load data
stroma_w_normal_gexp_pathways0 <- readxl::read_xlsx("/Users/oscarbruck/OneDrive - University of Helsinki/RCC/Otso/Analysis/Textures/Images/Stroma_WithNormal/Gexp/gexp_gene_pathways.xlsx")
stroma_wo_normal_gexp_pathways0 <- readxl::read_xlsx("/Users/oscarbruck/OneDrive - University of Helsinki/RCC/Otso/Analysis/Textures/Images/Stroma/Gexp/gexp_gene_pathways.xlsx")

# Filter by p
stroma_w_normal_gexp_pathways <- stroma_w_normal_gexp_pathways0 %>% 
  dplyr::filter(padj<0.05)
  # dplyr::mutate(source = "stroma_w_normal")
stroma_wo_normal_gexp_pathways <- stroma_wo_normal_gexp_pathways0 %>% 
  dplyr::filter(padj<0.05)
  # dplyr::mutate(source = "stroma_wo_normal")

# Rename variables
colnames(stroma_wo_normal_gexp_pathways)[2:(ncol(stroma_wo_normal_gexp_pathways))] <- paste0(colnames(stroma_wo_normal_gexp_pathways)[2:(ncol(stroma_wo_normal_gexp_pathways))], "_wo_normal")
colnames(stroma_w_normal_gexp_pathways)[2:(ncol(stroma_w_normal_gexp_pathways))] <- paste0(colnames(stroma_w_normal_gexp_pathways)[2:(ncol(stroma_w_normal_gexp_pathways))], "_w_normal")
 
# Join
stroma_gexp_pathways <- stroma_w_normal_gexp_pathways %>% full_join(stroma_wo_normal_gexp_pathways)

# Fill NES with 0 if NA
stroma_gexp_pathways <- stroma_gexp_pathways %>%
  dplyr::mutate(W_NORMAL = ifelse(is.na(NES_w_normal), 0, NES_w_normal),
                WO_NORMAL = ifelse(is.na(NES_wo_normal), 0, NES_wo_normal))

# Filter pathways that are padj >0.05 but padj<0.10
stroma_wo_normal_gexp_pathways1 <- stroma_wo_normal_gexp_pathways0 %>%
  dplyr::filter(padj < 0.10 & padj >= 0.05) %>%
  dplyr::select(pathway)
stroma_w_normal_gexp_pathways1 <- stroma_w_normal_gexp_pathways0 %>%
  dplyr::filter(padj < 0.10 & padj >= 0.05) %>%
  dplyr::select(pathway)
exclude1 <- stroma_gexp_pathways %>% dplyr::filter(W_NORMAL == 0) %>% dplyr::select(pathway) %>% dplyr::filter(pathway %in% stroma_w_normal_gexp_pathways1$pathway)
exclude2 <- stroma_gexp_pathways %>% dplyr::filter(WO_NORMAL == 0) %>% dplyr::select(pathway) %>% dplyr::filter(pathway %in% stroma_wo_normal_gexp_pathways1$pathway)
stroma_gexp_pathways <- stroma_gexp_pathways %>% dplyr::filter(!pathway %in% c(exclude1$pathway, exclude2$pathway))


# Transpose
stroma_gexp_pathways <- stroma_gexp_pathways %>% 
  dplyr::select(pathway, W_NORMAL, WO_NORMAL) %>%
  pivot_longer(cols = c(W_NORMAL, WO_NORMAL), names_to = "col1", values_to = "col2") %>% 
  pivot_wider(names_from = "pathway", values_from = "col2") %>%
  dplyr::rename(Stroma = col1)


# Identify FGSEA type
a <- colnames(stroma_gexp_pathways)[2:ncol(stroma_gexp_pathways)]
a = ifelse(str_detect(a, "HALLMARK"), "HALLMARK",
           ifelse(str_detect(a, "KEGG"), "KEGG",
                  ifelse(str_detect(a, "PID"), "PID",
                         ifelse(str_detect(a, "REACTOME"), "REACTOME", 
                                ifelse(str_detect(a, "chr"), "CHROMOSOME", NA)))))

# Rename colnames
colnames(stroma_gexp_pathways)[2:ncol(stroma_gexp_pathways)] <- gsub("HALLMARK_", "", colnames(stroma_gexp_pathways)[2:ncol(stroma_gexp_pathways)])
colnames(stroma_gexp_pathways)[2:ncol(stroma_gexp_pathways)] <- gsub("KEGG_", "", colnames(stroma_gexp_pathways)[2:ncol(stroma_gexp_pathways)])
colnames(stroma_gexp_pathways)[2:ncol(stroma_gexp_pathways)] <- gsub("PID_", "", colnames(stroma_gexp_pathways)[2:ncol(stroma_gexp_pathways)])
colnames(stroma_gexp_pathways)[2:ncol(stroma_gexp_pathways)] <- gsub("REACTOME_", "", colnames(stroma_gexp_pathways)[2:ncol(stroma_gexp_pathways)])
colnames(stroma_gexp_pathways)[2:ncol(stroma_gexp_pathways)] <- gsub("CHROMOSOME_", "", colnames(stroma_gexp_pathways)[2:ncol(stroma_gexp_pathways)])
colnames(stroma_gexp_pathways)[2:ncol(stroma_gexp_pathways)] <- gsub("TRIACYLGLYCEROL_AND_KETONE_BODY", "TAG & KETONE", colnames(stroma_gexp_pathways)[2:ncol(stroma_gexp_pathways)])
colnames(stroma_gexp_pathways)[2:ncol(stroma_gexp_pathways)] <- gsub("_", " ", colnames(stroma_gexp_pathways)[2:ncol(stroma_gexp_pathways)])


# Annotation data
anno <- data.frame(Stroma=stroma_gexp_pathways$Stroma)
anno_row <- data.frame(GSEA = a)

## Set1 colors
### #e41a1c #377eb8 #4daf4a #984ea3 #ff7f00 #ffff33 #a65628 #f781bf #999999

## Annotation colors 
# col_anno = list(Center = c(
#   "Christiana Healthcare"="#66a61e",
#   "Cureline" = "#d95f02",
#   "Harvard" = "#e7298a",
#   "International Genomics Consortium"="#1b9e77",
#   "MD Anderson Cancer Center" = "#666666",
#   "MSKCC"="#a6761d",
#   "NCI Urologic Oncology Branch"="#000000",
#   "UNC" = "#e6ab02",
#   "Fox Chase" = "#d95f02",
#   "University of Pittsburgh"="#7570b3"))
# mycolors <- colorRampPalette(brewer.pal(2, "Set1"))(length(unique(anno$Center)))
col_anno = list(Stroma = c(
  "W_NORMAL"="#e41a1c",
  "WO_NORMAL" = "#377eb8"))
col_anno_row = list(GSEA = c(
  "HALLMARK"="#4daf4a",
  "KEGG" = "#984ea3",
  "PID" = "#ff7f00",
  "REACTOME" = "#ffff33",
  "CHROMOSOME" = "#a65628"))


# col_anno1 = list(Center = structure(mycolors))
# col_anno1 <- setNames(object = col_anno1$Center, nm = unique(anno$Center))
# col_anno$Center <- col_anno1

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
hm_anno_row <- rowAnnotation(df = anno_row,
                             name = "Row Annotation",
                             col = col_anno_row,
                             show_annotation_name = T,
                             gp = gpar(col = "black"),
                             width = unit(2, "cm"),
                             # annotation_name_gp = gpar(fontsize=8, fontface="bold"),
                             annotation_name_gp = gpar(fontface="bold"),
                             annotation_legend_param = list(
                               title_gp = gpar(fontsize = 13, fontface="bold"),
                               labels_gp = gpar(fontsize = 13),
                               grid_height = unit(5, "mm"),
                               border=TRUE),
                             na_col = "white")

## Remove unnecessary metadata
stroma_gexp_pathways_hm <- stroma_gexp_pathways %>%
  dplyr::select(-Stroma) %>%
  as.matrix()


################################# Colors and scaling #############################################################################


# Color ColorBrewerilla
cols <- colorRampPalette(brewer.pal(11, "RdBu"))(256)

# Scale data using median as the center and the maximum value as the scale = (x-median)/max value
# stroma_gexp_pathways_hm1 <- scale(stroma_gexp_pathways_hm, center = apply(stroma_gexp_pathways_hm, 2, median, na.rm=TRUE), scale=apply(stroma_gexp_pathways_hm, 2, max, na.rm=TRUE))
# stroma_gexp_pathways_hm1[is.nan(stroma_gexp_pathways_hm1)] <- NA

## Remove columns with more than 50% NA
# scaled_hm <- scaled_hm[-which(rowMeans(is.na(scaled_hm)) > 0.50),]


# Plot
png("/Users/oscarbruck/OneDrive - University of Helsinki/RCC/Otso/Analysis/Textures/Images/Stroma/Stroma_w_vs_wo_normal_heatmap.png", width = 15, height = 6, units = 'in', res = 300, pointsize = 12) #original pointsize = 12
g <- Heatmap(t(stroma_gexp_pathways_hm),
             width = ncol(stroma_gexp_pathways_hm)*unit(0.35, "mm"), 
             top_annotation = hm_anno,
             right_annotation = hm_anno_row,
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
             clustering_distance_rows = "euclidean")

# print(g)

draw(g, annotation_legend_side="left", heatmap_legend_side="left")
dev.off()


