rm(list=ls())

print("Start Textures/Scripts/Margin/Margin_lymphocyte_gexp.R")

# Load libraries
library(readxl)
library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(RColorBrewer)
library(GSA)
library(fgsea)
library(ggpubr)
library(ggrepel)


############################# LOAD DATA ##########################################################################################################


# Image data
tcga_kirc <- read_xlsx("../data/image_analysis_results_final.xlsx")


# Normalize by removing empty
tcga_kirc$`texture_cancer_%` <- 100*tcga_kirc$texture_cancer / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`texture_normal_%` <- 100*tcga_kirc$texture_normal / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`margin_texture_normal_%` <- 100*tcga_kirc$margin_texture_normal / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_cancer + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)
tcga_kirc$`margin_texture_blood_%` <- 100*tcga_kirc$margin_texture_blood / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_cancer + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)
tcga_kirc$`margin_texture_cancer_%` <- 100*tcga_kirc$margin_texture_cancer / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_cancer + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)
tcga_kirc$`margin_texture_normal_%` <- 100*tcga_kirc$margin_texture_normal / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_cancer + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)
tcga_kirc$`margin_texture_stroma_%` <- 100*tcga_kirc$margin_texture_stroma / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_cancer + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)
tcga_kirc$`margin_texture_other_%` <- 100*tcga_kirc$margin_texture_other / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_cancer + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)


# Filter samples with >5% cancer normal tissue
tcga_kirc <- tcga_kirc %>%
  dplyr::filter(`texture_cancer_%` > 5) %>%
  dplyr::mutate(tissue_source_site = gsub("-[[:print:]]{4}", "", gsub("TCGA-", "", tcga_id)))

# Keep only centers with >20 samples
n <- tcga_kirc %>%
  group_by(tissue_source_site) %>%
  summarise(n=n()) %>%
  dplyr::filter(n>=20) %>%
  dplyr::filter(!tissue_source_site %in% c("AK")) 

tcga_kirc <- tcga_kirc %>%
  dplyr::filter(tissue_source_site %in% n$tissue_source_site)

# Filter centers with too much/little lymphocytes than expected
# 3Z	Mary Bird Perkins Cancer Center - Our Lady of the Lake
# 6D	University of Oklahoma HSC
# A3	International Genomics Consortium
# AK	Fox Chase
# AS	St. Joseph's Medical Center-(MD)
# B0	University of Pittsburgh
# B2	Christiana Healthcare
# B4	Cureline
# B8	UNC
# BP	MSKCC
# CB	ILSBio
# CJ	MD Anderson Cancer Center
# CW	Mayo Clinic - Rochester
# CZ	Harvard
# DV	NCI Urologic Oncology Branch
# EU	CHI-Penrose Colorado
# G6	Roswell Park
# GK	ABS - IUPUI
# MM	BLN - Baylor
# MW	University of Miami
# T7	Molecular Response
# V8	Medical College of Georgia
# WM	University of Kansas

# Save
tcga_kirc0 <- tcga_kirc


############################# PLOT ##########################################################################################################


# Rename
colnames(tcga_kirc)[grep(pattern = "^inf_margin_bin_lymphocytes", colnames(tcga_kirc))] <- c("Margin_Blood", "Margin_Normal", "Margin_Stroma", "Margin_Other")
colnames(tcga_kirc)[grep(pattern = "^inf_non_margin_bin_lymphocytes", colnames(tcga_kirc))] <- c("NonMargin_Blood", "NonMargin_Cancer", "NonMargin_Normal", "NonMargin_Stroma", "NonMargin_Other")
margin <- c("Margin_Blood", "Margin_Normal", "Margin_Stroma", "Margin_Other")
nonmargin <- c("NonMargin_Blood", "NonMargin_Cancer", "NonMargin_Normal", "NonMargin_Stroma", "NonMargin_Other")


# To percent
tcga_kirc[margin] <- sapply(tcga_kirc[margin], function(x) 100*x)
tcga_kirc[nonmargin] <- sapply(tcga_kirc[nonmargin], function(x) 100*x)

# Melt longer
tcga_kirc_long <- tcga_kirc %>%
  dplyr::select("tcga_id", "Margin_Blood", "Margin_Normal", "Margin_Stroma", "Margin_Other", "NonMargin_Blood", "NonMargin_Normal", "NonMargin_Stroma", "NonMargin_Other") %>%
  dplyr::rename("ID" = "tcga_id") %>%
  reshape2::melt() %>%
  arrange(desc(value)) %>%
  dplyr::mutate(Margin = gsub("_[[:print:]]*", "", variable),
                Texture = gsub("[[:print:]]*_", "", variable),
                Margin = factor(Margin),
                Texture = factor(Texture),
                variable = factor(variable))






# Calculate whether there are more lymphocytes in the margin or non-margin by textures

# Divide the dataframe into margin and non-margin
tcga_kirc_margin <- tcga_kirc_long %>%
  dplyr::select(-variable) %>%
  dplyr::filter(Margin == "Margin") %>%
  dplyr::select(-Margin) %>%
  dplyr::rename(Margin = value)

tcga_kirc_nonmargin <- tcga_kirc_long %>%
  dplyr::select(-variable) %>%
  dplyr::filter(Margin == "NonMargin") %>%
  dplyr::select(-Margin) %>%
  dplyr::rename(NonMargin = value)

# Join margin and non-margin dataframes
tcga_kirc_wide <- tcga_kirc_margin %>%
  full_join(tcga_kirc_nonmargin) %>%
  mutate(FC = Margin/NonMargin,
         FC = ifelse(FC == "Inf", (Margin+1)/(NonMargin+1),
                     ifelse(FC > 10, 10,
                            ifelse(FC < 0.1, 0.1, FC))))

# Read gexp data
tcga_kirc <- readRDS("../data/clinical_transcriptome.rds") %>%
  dplyr::select(one_of(tcga_kirc_wide$ID)) %>%
  rownames_to_column() %>%
  dplyr::filter(str_detect(rowname, "GEXP"))

# Clean gexp data
tcga_kirc <- tcga_kirc %>% 
  pivot_longer(!rowname, names_to = "col1", values_to = "col2") %>% 
  pivot_wider(names_from = "rowname", values_from = "col2") %>%
  dplyr::rename(tcga_id = col1)

# Keep only genes with highest median
a <- sapply(tcga_kirc[2:ncol(tcga_kirc)], function(x) median(x, na.rm = TRUE) > 8)
a <- a[a]
tcga_kirc <- tcga_kirc %>% dplyr::select(tcga_id, one_of(names(a)))

# Join
tcga_kirc <- tcga_kirc_wide %>%
  dplyr::inner_join(tcga_kirc %>% dplyr::rename(ID = tcga_id))



############################# Genesets ##########################################################################################################



# Load pathways
h <- GSA.read.gmt("../data/gsea_pathways/h.all.v6.2.symbols.gmt")
c1 <- GSA.read.gmt("../data/gsea_pathways/c1.all.v6.2.symbols.gmt")
c2 <- GSA.read.gmt("../data/gsea_pathways/c2.all.v6.2.symbols.gmt")
c5 <- GSA.read.gmt("../data/gsea_pathways/c5.all.v6.2.symbols.gmt")
c6 <- GSA.read.gmt("../data/gsea_pathways/c6.all.v6.2.symbols.gmt")

a <- grep(x = c2$geneset.names, pattern = "^KEGG|^PID|^REACTOME|^BIOCARTA", invert = FALSE)
c2$genesets <- c2$genesets[a]
c2$geneset.names <- c2$geneset.names[a]
c2$geneset.descriptions <- c2$geneset.descriptions[a]


c1$genesets <- c(c1$genesets, h$genesets, c2$genesets)
c1$geneset.names <- c(c1$geneset.names, h$geneset.names, c2$geneset.names)
c1$geneset.descriptions <- c(c1$geneset.descriptions, h$geneset.descriptions, c2$geneset.descriptions)
pathways <- c1


############################# PLOT ##########################################################################################################


# Binary
tcga_kirc$FC <- ifelse(tcga_kirc$FC > 1, "High", 
                       ifelse(tcga_kirc$FC < 1, "Low", NA))

# GEXP data
findcolnumber <- function(df, thecolumnname){
  which(colnames(df) == deparse(substitute(thecolumnname)))
}


# Join normal texture
tcga_kirc <- tcga_kirc %>% dplyr::left_join(tcga_kirc0 %>% dplyr::select(ID = tcga_id, `texture_normal_%`))

# Loop over all textures
## Save data
tcga_kirc0 <- tcga_kirc

for (texture1 in c("All", as.list(unique(tcga_kirc0$Texture))) ) {
  print(texture1)
  
  # Reset data
  if (texture1 == "All") {
    tcga_kirc <- tcga_kirc0 %>%
      dplyr::filter(Texture %in% unique(tcga_kirc0$Texture)) %>%
      dplyr::filter(!(Texture == "Normal" & `texture_normal_%` <= 1))
    
  } else if (texture1 == "Normal") {
    
    tcga_kirc <- tcga_kirc0 %>%
      dplyr::filter(`texture_normal_%` > 1) %>%
      dplyr::filter(Texture %in% texture1)
    
  } else {
    tcga_kirc <- tcga_kirc0 %>%
      dplyr::filter(Texture %in% texture1)
  }
  
  # Multiple comparison (two-sided, unpaired, Mann-Whitney U test)
  # gexp_data <- tcga_kirc[(findcolnumber(tcga_kirc, tissue_source_site)+1):ncol(tcga_kirc)]
  multiple_t_tests_p_value <- lapply(tcga_kirc[(findcolnumber(tcga_kirc, FC)+1):(ncol(tcga_kirc)-1)],
                                     function(x) wilcox.test(x ~ tcga_kirc$FC, correct=FALSE, exact=FALSE, estimate=TRUE, na.rm=TRUE))
  ## P-values can be extracted from the result object
  pvalue <- data.frame(p.value = sapply(multiple_t_tests_p_value, getElement, name = "p.value"))
  ## Create a matrix and dataframe of the p-values
  pvalue_mat <- as.matrix(pvalue)
  pvalue_df <- data.frame(pvalue)
  ## Add the p values to a new dataframe
  pvalue_df$p_adjusted <- p.adjust(p = pvalue_mat, method = "BH")
  ## Add also the t values, 95%CI to the same dataframe
  pvalue_df$median1 <- unlist(lapply(tcga_kirc[(findcolnumber(tcga_kirc, FC)+1):(ncol(tcga_kirc)-1)], function(x) median(x[tcga_kirc$FC=="High"], na.rm=TRUE)))
  pvalue_df$median2 <- unlist(lapply(tcga_kirc[(findcolnumber(tcga_kirc, FC)+1):(ncol(tcga_kirc)-1)], function(x) median(x[tcga_kirc$FC=="Low"], na.rm=TRUE)))
  ## Rownames to column
  pvalue_df <- pvalue_df %>%
    rownames_to_column() %>%
    rename(genes = rowname,
           pvalue = p.value
    ) %>%
    arrange(pvalue)
  rownames(pvalue_df) = NULL
  
  
  
  # Export data
  writexl::write_xlsx(pvalue_df, paste0("Textures/Images/Margin/Margin_vs_non_margin/Ly/Gexp/Ly_in_margin_vs_nonmargin_", texture1, "_gexp.xlsx"))
  
  
  
  
  # Rename variables
  pvalue_df <- pvalue_df %>%
    mutate(genes = gsub("N:GEXP:", "", genes)
    )
  
  
  # GSEA
  Ranks_tmp <- pvalue_df %>% dplyr::mutate(FC = log(median1/median2)) %>% dplyr::filter(!is.na(FC) & pvalue<0.05 & ((median1+median2)/2)>8)
  # Ranks_tmp <- pvalue_df %>% dplyr::mutate(FC = 1/pvalue) %>% dplyr::filter(!is.na(FC))
  Ranks <- Ranks_tmp %>% dplyr::select(FC)
  Ranks <- setNames(object = Ranks$FC, nm = Ranks_tmp$genes)
  names(pathways$genesets) <- pathways$geneset.names
  fgseaRes_blood <- fgsea(pathways = pathways$genesets,
                          stats = Ranks,
                          minSize=15,
                          maxSize=500,
                          nperm=10000)
  # top 6 enriched pathways
  head(fgseaRes_blood[order(pval), ])
  
  
  topPathwaysUp <- fgseaRes_blood[ES > 0][padj < 0.05][head(order(-NES), n=10), pathway]
  topPathwaysDown <- fgseaRes_blood[ES < 0][padj < 0.05][head(order(NES), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  
  png(paste0("Textures/Images/Margin/Margin_vs_non_margin/Ly/Gexp/", texture1, "_gexp_top20_selectedpathways.png"), width = 11.5, height = 3.5, units = 'in', res = 300)
  plotGseaTable(pathways$genesets[topPathways], Ranks, fgseaRes_blood,
                colwidths = c(5, 4, 0.7, 0.7, 0.7),
                gseaParam = 0.5)
  dev.off()
  
  
  # Barplot
  fgseaRes_res2 <- fgseaRes_blood[fgseaRes_blood$pathway %in% topPathways] %>%
    dplyr::mutate(pathway = gsub("_", " ", pathway),
                  Sig = ifelse(padj < 0.001, "***",
                               ifelse(padj < 0.01, "**",
                                      ifelse(padj < 0.05, "*", "p<0.1"))),
                  Sig = factor(Sig, levels = c("p<0.1", "*", "**", "***")))
  g <- ggplot(data = fgseaRes_res2, aes(x = reorder(pathway, NES), y = NES, fill = Sig)) +
    geom_bar(stat = "identity", width = 0.8, color="black") +
    coord_flip() +
    xlab("") +
    theme_bw() +
    # geom_hline(yintercept = 0, size = 2, linetype = 0) +
    scale_fill_brewer(palette="Reds", limits=c("*","**","***"), name = "Sig", labels=c("*","**","***")) +
    # labs(x = "Immune marker",
    #      # title = "                  Immune Markers in AML BM vs. PB",
    #      y = "Difference in Proportion (%)") +
    theme(
      # plot.title = element_text(size=20, face="bold"),
      axis.text.x = element_text(face="bold", colour = "black"),
      axis.text.y = element_text(face="bold", colour = "black"),
      axis.title.x = element_text(size=15, face="bold", colour = "black"), #,angle=90
      axis.title.y = element_text(size=15, face="bold", colour = "black"), #,angle=90
      legend.title = element_text(size=15, face="bold"),
      legend.text = element_text(size=15, face="bold", vjust = 0.0),
      # legend.justification = "top",
      # legend.position = c(0.85, 0.1) )
      legend.position = "bottom" )
  ggsave(plot = g, filename = paste0("Textures/Images/Margin/Margin_vs_non_margin/Ly/Gexp/", texture1, "_pathways_barplot.png"), width = max(5.5, 5+max(nchar(fgseaRes_res2$pathway))/23), height = max(2, (1+nrow(fgseaRes_res2)/5.5)), units = 'in', dpi = 300)
  
  
  
  
  # Prepare data for volcanoplot
  pvalue_df1 <- pvalue_df %>% mutate(
    fold_change1 = ifelse(median1 == 0 & median2 == 0, 1,
                          ifelse(median2 == 0, median1+1,
                                 ifelse(median1 == 0, 1/(median2+1), median1/median2))),
    # fold_change2 = (median2-median1)/((median2+median1)/2),
    sig1 = ifelse(pvalue<0.01, 1, 0),
    sig2 = ifelse(pvalue<0.01 & fold_change1 > 1, "p<0.01 & FC>1",
                  ifelse(pvalue<0.01 & fold_change1 < 1, "p<0.01 & FC<1", "ns"))
    # sig_adj1 = ifelse(p_adjusted<0.05, "FDR", "ns"),
    # sig_adj2 = ifelse(p_adjusted<0.05 & fold_change1 > 1, "FDR<0.05 & FC>1",
    # ifelse(p_adjusted<0.05 & fold_change1 < 1, "FDR<0.05 & FC<1", "ns"))
  ) %>%
    mutate(sig1 = as.factor(sig1))
  
  pvalue_df1 <- pvalue_df1 %>% mutate(
    label = ifelse((-log10(pvalue)>1 & log(fold_change1)>0), 2,    # 1.30103
                   ifelse((-log10(pvalue)>1 & log(fold_change1)<0), 1, 0)  # 1.30103
    )) %>% mutate(
      label = factor(label, levels = c(0,1,2))
    )
  
  
  # Colors
  # if ( unique(pvalue_df1$label %in% 2) & max(pvalue_df1$label) == 2 ) {
  if ( max(as.numeric(as.character(pvalue_df1$label))) == 2 & length(unique(as.numeric(as.character(pvalue_df1$label)))) == 2) {
    values = c("grey", "red", "blue")
  } else {
    values = c("grey", "blue", "red")
  }
  # values = c("white", "blue", "red")
  
  # Plot
  g1 <- ggplot(pvalue_df1, aes(log(fold_change1), -log10(pvalue))) +
    geom_point(aes(fill=sig2), pch = 21, color = "black", size = 4) +
    labs(x="Fold change (Log)", y="Log10 p-value") +
    theme_bw() +
    theme(axis.title.y = element_text(size = 15, face = "bold"),
          axis.title.x = element_text(size = 15, face = "bold"),
          axis.text.y = element_text(size = 15, face = "bold"),
          axis.text.x = element_text(size = 15, face = "bold"),
          legend.title = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 15, face = "bold"),
          legend.direction = 'horizontal', 
          legend.position = 'bottom',
          legend.key = element_rect(size = 5),
          legend.key.size = unit(1.5, 'lines')) +
    scale_fill_manual(name="Significance", values=values) +
    geom_vline(xintercept=0, linetype="dashed", color = "black", size=1) +
    geom_label_repel(data=pvalue_df1[pvalue_df1$pvalue<0.005,],
                     # label.padding = 0.5,
                     segment.size = 1,
                     aes(label=genes), size = 5, fontface = "bold")
  
  g2 <- ggplot(pvalue_df1, aes(log(fold_change1), -log10(pvalue))) +
    geom_point(aes(fill=sig2), pch = 21, color = "black", size = 4) +
    labs(x="Fold change (Log)", y="Log10 p-value") +
    theme_bw() +
    theme(axis.title.y = element_text(size = 15, face = "bold"),
          axis.title.x = element_text(size = 15, face = "bold"),
          axis.text.y = element_text(size = 15, face = "bold"),
          axis.text.x = element_text(size = 15, face = "bold"),
          legend.title = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 15, face = "bold"),
          legend.direction = 'horizontal', 
          legend.position = 'bottom',
          legend.key = element_rect(size = 5),
          legend.key.size = unit(1.5, 'lines')) +
    scale_fill_manual(name="Significance", values=values) +
    geom_vline(xintercept=0, linetype="dashed", color = "black", size=1) +
    geom_label_repel(data=pvalue_df1[pvalue_df1$pvalue<0.001,],
                     # label.padding = 0.5,
                     segment.size = 1,
                     aes(label=genes), size = 5, fontface = "bold")
  
  g3 <- ggplot(pvalue_df1, aes(log(fold_change1), -log10(pvalue))) +
    geom_point(aes(fill=sig2), pch = 21, color = "black", size = 4) +
    labs(x="Fold change (Log)", y="Log10 p-value") +
    theme_bw() +
    theme(axis.title.y = element_text(size = 15, face = "bold"),
          axis.title.x = element_text(size = 15, face = "bold"),
          axis.text.y = element_text(size = 15, face = "bold"),
          axis.text.x = element_text(size = 15, face = "bold"),
          legend.title = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 15, face = "bold"),
          legend.direction = 'horizontal', 
          legend.position = 'bottom',
          legend.key = element_rect(size = 5),
          legend.key.size = unit(1.5, 'lines')) +
    scale_fill_manual(name="Significance", values=values) +
    geom_vline(xintercept=0, linetype="dashed", color = "black", size=1) +
    geom_label_repel(data=pvalue_df1[pvalue_df1$pvalue<0.0005,],
                     # label.padding = 0.5,
                     segment.size = 1,
                     aes(label=genes), size = 5, fontface = "bold")
  
  g4 <- ggplot(pvalue_df1, aes(log(fold_change1), -log10(pvalue))) +
    geom_point(aes(fill=sig2), pch = 21, color = "black", size = 4) +
    labs(x="Fold change (Log)", y="Log10 p-value") +
    theme_bw() +
    theme(axis.title.y = element_text(size = 15, face = "bold"),
          axis.title.x = element_text(size = 15, face = "bold"),
          axis.text.y = element_text(size = 15, face = "bold"),
          axis.text.x = element_text(size = 15, face = "bold"),
          legend.title = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 15, face = "bold"),
          legend.direction = 'horizontal', 
          legend.position = 'bottom',
          legend.key = element_rect(size = 5),
          legend.key.size = unit(1.5, 'lines')) +
    scale_fill_manual(name="Significance", values=values) +
    geom_vline(xintercept=0, linetype="dashed", color = "black", size=1) +
    geom_label_repel(data=pvalue_df1[pvalue_df1$pvalue<0.0001,],
                     # label.padding = 0.5,
                     segment.size = 1,
                     aes(label=genes), size = 5, fontface = "bold")
  
  ggsave(plot = g1, filename = paste0("Textures/Images/Margin/Margin_vs_non_margin/Ly/Gexp/", texture1, "_gexp_volcano_p0005.png"), width = 7, height = 7, units = 'in', dpi = 300)
  ggsave(plot = g2, filename = paste0("Textures/Images/Margin/Margin_vs_non_margin/Ly/Gexp/", texture1, "_gexp_volcano_p0001.png"), width = 7, height = 7, units = 'in', dpi = 300)
  ggsave(plot = g3, filename = paste0("Textures/Images/Margin/Margin_vs_non_margin/Ly/Gexp/", texture1, "_gexp_volcano_p00005.png"), width = 7, height = 7, units = 'in', dpi = 300)
  ggsave(plot = g4, filename = paste0("Textures/Images/Margin/Margin_vs_non_margin/Ly/Gexp/", texture1, "_gexp_volcano_p00001.png"), width = 7, height = 7, units = 'in', dpi = 300)
  
}

