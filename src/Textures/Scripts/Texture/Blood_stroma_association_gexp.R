rm(list = ls())

print("Start Textures/Scripts/Texture/Blood_stroma_association_gexp.R")


# Load libraries
library(readxl)
library(tidyverse)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(rstatix)
library(GSA)
library(fgsea)
library(ggrepel)
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


# Petri's data
df <- readRDS("../data/clinical_transcriptome.rds") %>%
  # Filter patients
  dplyr::select(one_of(tcga_kirc$tcga_id))

# Clean clinical data
df1 <- df[grep(pattern = "N:GEXP", x = rownames(df)),]
df1 <- df1 %>% rownames_to_column()
df2 <- df1 %>% 
  pivot_longer(!rowname, names_to = "col1", values_to = "col2") %>% 
  pivot_wider(names_from = "rowname", values_from = "col2") %>%
  dplyr::rename(tcga_id = col1)

# Keep only genes with highest median
a <- sapply(df2[2:ncol(df2)], function(x) median(x, na.rm = TRUE) > 8)
a <- a[a]
df2 <- df2 %>% dplyr::select(tcga_id, one_of(names(a)))


# Join
tcga_kirc <- tcga_kirc %>%
  dplyr::left_join(df2)

# Rename
colnames(tcga_kirc)[grep(pattern = "^texture_[[:print:]]*%", colnames(tcga_kirc))] <- c("Blood", "Cancer", "Normal", "Stroma", "Other")
textures <- c("Blood", "Cancer", "Normal", "Stroma", "Other")

# Save
tcga_kirc0 <- tcga_kirc


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


###########################
##### BLOOD WO NORMAL #####
###########################


# Reset
tcga_kirc <- tcga_kirc0

# Filter only samples with ≥1% normal tissue
tcga_kirc <- tcga_kirc %>% dplyr::filter(Normal < 1)

# Binary
tcga_kirc$Blood <- ifelse(tcga_kirc$Blood > median(tcga_kirc$Blood, na.rm = TRUE), "High", "Low")

# GEXP data
findcolnumber <- function(df, thecolumnname){
  which(colnames(df) == deparse(substitute(thecolumnname)))
}

# Multiple comparison (two-sided, unpaired, Mann-Whitney U test)
# gexp_data <- tcga_kirc[(findcolnumber(tcga_kirc, tissue_source_site)+1):ncol(tcga_kirc)]
multiple_t_tests_p_value <- lapply(tcga_kirc[(findcolnumber(tcga_kirc, tissue_source_site)+1):ncol(tcga_kirc)],
                                   function(x) wilcox.test(x ~ tcga_kirc$Blood, correct=FALSE, exact=FALSE, estimate=TRUE, na.rm=TRUE))
## P-values can be extracted from the result object
pvalue <- data.frame(p.value = sapply(multiple_t_tests_p_value, getElement, name = "p.value"))
## Create a matrix and dataframe of the p-values
pvalue_mat <- as.matrix(pvalue)
pvalue_df <- data.frame(pvalue)
## Add the p values to a new dataframe
pvalue_df$p_adjusted <- p.adjust(p = pvalue_mat, method = "BH")
## Add also the t values, 95%CI to the same dataframe
pvalue_df$median1 <- unlist(lapply(tcga_kirc[(findcolnumber(tcga_kirc, tissue_source_site)+1):ncol(tcga_kirc)], function(x) median(x[tcga_kirc$Blood=="High"], na.rm=TRUE)))
pvalue_df$median2 <- unlist(lapply(tcga_kirc[(findcolnumber(tcga_kirc, tissue_source_site)+1):ncol(tcga_kirc)], function(x) median(x[tcga_kirc$Blood=="Low"], na.rm=TRUE)))
## Rownames to column
pvalue_df <- pvalue_df %>%
  rownames_to_column() %>%
  rename(genes = rowname,
         pvalue = p.value
  ) %>%
  arrange(pvalue)
rownames(pvalue_df) = NULL


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
                        minSize=10,
                        maxSize=500,
                        nperm=10000)
# top 6 enriched pathways
head(fgseaRes_blood[order(pval), ])

topPathwaysUp <- fgseaRes_blood[ES > 0][padj < 0.05][head(order(-NES), n=10), pathway]
topPathwaysDown <- fgseaRes_blood[ES < 0][padj < 0.05][head(order(NES), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

png("Textures/Images/Blood/Gexp/Blood_gexp_top20_selectedpathways.png", width = 10, height = 1.25, units = 'in', res = 300)
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
  theme(
    # plot.title = element_text(size=20, face="bold"),
    axis.text.x = element_text(face="bold", colour = "black"),
    axis.text.y = element_text(face="bold", colour = "black"),
    axis.title.x = element_text(size=15, face="bold", colour = "black"), #,angle=90
    axis.title.y = element_text(size=15, face="bold", colour = "black"), #,angle=90
    legend.title = element_text(size=15, face="bold"),
    legend.text = element_text(size=15, face="bold", vjust = 0.0),
    legend.position = "bottom" )
ggsave(plot = g, filename = "Textures/Images/Blood/Gexp/Blood_gexp_top20_selectedpathways_barplot.png", width = max(5.5, 5+max(nchar(fgseaRes_res2$pathway))/23), height = max(2, (1+nrow(fgseaRes_res2)/5.5)), units = 'in', dpi = 300)




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
# png(paste0("/Users/oscarbruck/OneDrive - University of Helsinki/Tutkimus/Projekteja/MDS/Analysis/HE/Results/Volcanoplot_mIHC_singlemarker_", dico, ".png"),
# width = 9, height = 9, units = 'in', res = 300)
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
  geom_label_repel(data=pvalue_df1[pvalue_df1$pvalue<0.005 & pvalue_df1$fold_change1 != 1,],
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
  geom_label_repel(data=pvalue_df1[pvalue_df1$pvalue<0.001 & pvalue_df1$fold_change1 != 1,],
                   # label.padding = 0.5,
                   segment.size = 1,
                   aes(label=genes), size = 5, fontface = "bold")
ggsave(plot = g1, filename = "Textures/Images/Blood/Blood_gexp_volcano_p0005.png", width = 7, height = 7, units = 'in', dpi = 300)
ggsave(plot = g2, filename = "Textures/Images/Blood/Blood_gexp_volcano_p0001.png", width = 7, height = 7, units = 'in', dpi = 300)



#############################
##### BLOOD WITH NORMAL #####
#############################


# Reset
tcga_kirc <- tcga_kirc0

# Filter only samples with ≥1% normal tissue
tcga_kirc <- tcga_kirc %>% dplyr::filter(Normal >= 1)

# Binary
tcga_kirc$Blood <- ifelse(tcga_kirc$Blood > median(tcga_kirc$Blood, na.rm = TRUE), "High", "Low")

# GEXP data
findcolnumber <- function(df, thecolumnname){
  which(colnames(df) == deparse(substitute(thecolumnname)))
}

# Multiple comparison (two-sided, unpaired, Mann-Whitney U test)
# gexp_data <- tcga_kirc[(findcolnumber(tcga_kirc, tissue_source_site)+1):ncol(tcga_kirc)]
multiple_t_tests_p_value <- lapply(tcga_kirc[(findcolnumber(tcga_kirc, tissue_source_site)+1):ncol(tcga_kirc)],
                                   function(x) wilcox.test(x ~ tcga_kirc$Blood, correct=FALSE, exact=FALSE, estimate=TRUE, na.rm=TRUE))
## P-values can be extracted from the result object
pvalue <- data.frame(p.value = sapply(multiple_t_tests_p_value, getElement, name = "p.value"))
## Create a matrix and dataframe of the p-values
pvalue_mat <- as.matrix(pvalue)
pvalue_df <- data.frame(pvalue)
## Add the p values to a new dataframe
pvalue_df$p_adjusted <- p.adjust(p = pvalue_mat, method = "BH")
## Add also the t values, 95%CI to the same dataframe
pvalue_df$median1 <- unlist(lapply(tcga_kirc[(findcolnumber(tcga_kirc, tissue_source_site)+1):ncol(tcga_kirc)], function(x) median(x[tcga_kirc$Blood=="High"], na.rm=TRUE)))
pvalue_df$median2 <- unlist(lapply(tcga_kirc[(findcolnumber(tcga_kirc, tissue_source_site)+1):ncol(tcga_kirc)], function(x) median(x[tcga_kirc$Blood=="Low"], na.rm=TRUE)))
## Rownames to column
pvalue_df <- pvalue_df %>%
  rownames_to_column() %>%
  rename(genes = rowname,
         pvalue = p.value
  ) %>%
  arrange(pvalue)
rownames(pvalue_df) = NULL


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
                        minSize=10,
                        maxSize=500,
                        nperm=10000)
# top 6 enriched pathways
head(fgseaRes_blood[order(pval), ])

topPathwaysUp <- fgseaRes_blood[ES > 0][padj < 0.05][head(order(-NES), n=10), pathway]
topPathwaysDown <- fgseaRes_blood[ES < 0][padj < 0.05][head(order(NES), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

png("Textures/Images/Blood_WithNormal/Gexp/Blood_gexp_top20_selectedpathways.png", width = 10, height = 1.25, units = 'in', res = 300)
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
  # geom_text(data = flow_median[flow_median$Disease=="AML",], aes(label = round(BM,1), x = name, angle = 0, y = -75), fontface="bold", size=4) + # angle=0,
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
ggsave(plot = g, filename = "Textures/Images/Blood_WithNormal/Gexp/Blood_gexp_top20_selectedpathways_barplot.png", width = max(5.5, 5+max(nchar(fgseaRes_res2$pathway))/23), height = max(2, (1+nrow(fgseaRes_res2)/5.5)), units = 'in', dpi = 300)




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
  geom_label_repel(data=pvalue_df1[pvalue_df1$pvalue<0.0005 & pvalue_df1$fold_change1 != 1,],
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
  geom_label_repel(data=pvalue_df1[pvalue_df1$pvalue<0.001 & pvalue_df1$fold_change1 != 1,],
                   # label.padding = 0.5,
                   segment.size = 1,
                   aes(label=genes), size = 5, fontface = "bold")
ggsave(plot = g1, filename = "Textures/Images/Blood_WithNormal/Blood_gexp_volcano_p0005.png", width = 7, height = 7, units = 'in', dpi = 300)
ggsave(plot = g2, filename = "Textures/Images/Blood_WithNormal/Blood_gexp_volcano_p001.png", width = 7, height = 7, units = 'in', dpi = 300)





############################
##### STROMA WO NORMAL #####
############################


# Reset
tcga_kirc <- tcga_kirc0

# Filter only samples with ≥1% normal tissue
tcga_kirc <- tcga_kirc %>% dplyr::filter(Normal < 1)

# Binary
tcga_kirc$Stroma <- ifelse(tcga_kirc$Stroma > median(tcga_kirc$Stroma, na.rm = TRUE), "High", "Low")

# GEXP data
findcolnumber <- function(df, thecolumnname){
  which(colnames(df) == deparse(substitute(thecolumnname)))
}

# Multiple comparison (two-sided, unpaired, Mann-Whitney U test)
# gexp_data <- tcga_kirc[(findcolnumber(tcga_kirc, tissue_source_site)+1):ncol(tcga_kirc)]
multiple_t_tests_p_value <- lapply(tcga_kirc[(findcolnumber(tcga_kirc, tissue_source_site)+1):ncol(tcga_kirc)],
                                   function(x) wilcox.test(x ~ tcga_kirc$Stroma, correct=FALSE, exact=FALSE, estimate=TRUE, na.rm=TRUE))
## P-values can be extracted from the result object
pvalue <- data.frame(p.value = sapply(multiple_t_tests_p_value, getElement, name = "p.value"))
## Create a matrix and dataframe of the p-values
pvalue_mat <- as.matrix(pvalue)
pvalue_df <- data.frame(pvalue)
## Add the p values to a new dataframe
pvalue_df$p_adjusted <- p.adjust(p = pvalue_mat, method = "BH")
## Add also the t values, 95%CI to the same dataframe
pvalue_df$median1 <- unlist(lapply(tcga_kirc[(findcolnumber(tcga_kirc, tissue_source_site)+1):ncol(tcga_kirc)], function(x) median(x[tcga_kirc$Stroma=="High"], na.rm=TRUE)))
pvalue_df$median2 <- unlist(lapply(tcga_kirc[(findcolnumber(tcga_kirc, tissue_source_site)+1):ncol(tcga_kirc)], function(x) median(x[tcga_kirc$Stroma=="Low"], na.rm=TRUE)))
## Rownames to column
pvalue_df <- pvalue_df %>%
  rownames_to_column() %>%
  rename(genes = rowname,
         pvalue = p.value
  ) %>%
  arrange(pvalue)
rownames(pvalue_df) = NULL


# Rename variables
pvalue_df <- pvalue_df %>%
  mutate(genes = gsub("N:GEXP:", "", genes)
  )

# Export
writexl::write_xlsx(pvalue_df, "Textures/Images/Stroma/Gexp/gexp_genes.xlsx")


# GSEA
Ranks_tmp <- pvalue_df %>% dplyr::mutate(FC = log(median1/median2)) %>% dplyr::filter(!is.na(FC) & pvalue<0.05 & ((median1+median2)/2)>8)
# Ranks_tmp <- pvalue_df %>% dplyr::mutate(FC = 1/pvalue) %>% dplyr::filter(!is.na(FC))
Ranks <- Ranks_tmp %>% dplyr::select(FC)
Ranks <- setNames(object = Ranks$FC, nm = Ranks_tmp$genes)
names(pathways$genesets) <- pathways$geneset.names
fgseaRes_Stroma <- fgsea(pathways = pathways$genesets,
                         stats = Ranks,
                         minSize=10,
                         maxSize=500,
                         nperm=10000)
# Export
writexl::write_xlsx(fgseaRes_Stroma, "Textures/Images/Stroma/Gexp/gexp_gene_pathways.xlsx")
# top 6 enriched pathways
head(fgseaRes_Stroma[order(pval), ])


topPathwaysUp <- fgseaRes_Stroma[ES > 0][padj < 0.05][head(order(-NES), n=10), pathway]
topPathwaysDown <- fgseaRes_Stroma[ES < 0][padj < 0.05][head(order(NES), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

png("Textures/Images/Stroma/Gexp/Stroma_gexp_top20_selectedpathways.png", width = 14, height = 4, units = 'in', res = 300)
plotGseaTable(pathways$genesets[topPathways], Ranks, fgseaRes_Stroma,
              colwidths = c(5, 4, 0.7, 0.7, 0.7),
              gseaParam = 0.5)
dev.off()


# Barplot
fgseaRes_res2 <- fgseaRes_Stroma[fgseaRes_Stroma$pathway %in% topPathways] %>%
  dplyr::mutate(pathway = gsub("_", " ", pathway),
                Sig = ifelse(padj < 0.001, "***",
                             ifelse(padj < 0.01, "**",
                                    ifelse(padj < 0.05, "*", "p<0.1"))),
                Sig = factor(Sig, levels = c("p<0.1", "*", "**", "***")))
g <- ggplot(data = fgseaRes_res2, aes(x = reorder(pathway, NES), y = NES, fill = Sig)) +
  geom_bar(stat = "identity", width = 0.8, color="black") +
  # geom_text(data = flow_median[flow_median$Disease=="AML",], aes(label = round(BM,1), x = name, angle = 0, y = -75), fontface="bold", size=4) + # angle=0,
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
ggsave(plot = g, filename = "Textures/Images/Stroma/Gexp/Stroma_gexp_top20_selectedpathways_barplot.png", width = max(5.5, 5+max(nchar(fgseaRes_res2$pathway))/23), height = max(2, (1+nrow(fgseaRes_res2)/5.5)), units = 'in', dpi = 300)




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
# png(paste0("/Users/oscarbruck/OneDrive - University of Helsinki/Tutkimus/Projekteja/MDS/Analysis/HE/Results/Volcanoplot_mIHC_singlemarker_", dico, ".png"),
# width = 9, height = 9, units = 'in', res = 300)
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
  geom_label_repel(data=pvalue_df1[pvalue_df1$pvalue<0.0001 & pvalue_df1$fold_change1 != 1,],
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
  geom_label_repel(data=pvalue_df1[pvalue_df1$pvalue<0.00001 & pvalue_df1$fold_change1 != 1,],
                   # label.padding = 0.5,
                   segment.size = 1,
                   aes(label=genes), size = 5, fontface = "bold")
ggsave(plot = g1, filename = "Textures/Images/Stroma/Gexp/Stroma_gexp_volcano_p00001.png", width = 7, height = 7, units = 'in', dpi = 300)
ggsave(plot = g2, filename = "Textures/Images/Stroma/Gexp/Stroma_gexp_volcano_p000001.png", width = 7, height = 7, units = 'in', dpi = 300)







##############################
##### STROMA WITH NORMAL #####
##############################


# Reset
tcga_kirc <- tcga_kirc0

# Filter only samples with ≥1% normal tissue
tcga_kirc <- tcga_kirc %>% dplyr::filter(Normal >= 1)

# Binary
tcga_kirc$Stroma <- ifelse(tcga_kirc$Stroma > median(tcga_kirc$Stroma, na.rm = TRUE), "High", "Low")

# GEXP data
findcolnumber <- function(df, thecolumnname){
  which(colnames(df) == deparse(substitute(thecolumnname)))
}

# Multiple comparison (two-sided, unpaired, Mann-Whitney U test)
# gexp_data <- tcga_kirc[(findcolnumber(tcga_kirc, tissue_source_site)+1):ncol(tcga_kirc)]
multiple_t_tests_p_value <- lapply(tcga_kirc[(findcolnumber(tcga_kirc, tissue_source_site)+1):ncol(tcga_kirc)],
                                   function(x) wilcox.test(x ~ tcga_kirc$Stroma, correct=FALSE, exact=FALSE, estimate=TRUE, na.rm=TRUE))
## P-values can be extracted from the result object
pvalue <- data.frame(p.value = sapply(multiple_t_tests_p_value, getElement, name = "p.value"))
## Create a matrix and dataframe of the p-values
pvalue_mat <- as.matrix(pvalue)
pvalue_df <- data.frame(pvalue)
## Add the p values to a new dataframe
pvalue_df$p_adjusted <- p.adjust(p = pvalue_mat, method = "BH")
## Add also the t values, 95%CI to the same dataframe
pvalue_df$median1 <- unlist(lapply(tcga_kirc[(findcolnumber(tcga_kirc, tissue_source_site)+1):ncol(tcga_kirc)], function(x) median(x[tcga_kirc$Stroma=="High"], na.rm=TRUE)))
pvalue_df$median2 <- unlist(lapply(tcga_kirc[(findcolnumber(tcga_kirc, tissue_source_site)+1):ncol(tcga_kirc)], function(x) median(x[tcga_kirc$Stroma=="Low"], na.rm=TRUE)))
## Rownames to column
pvalue_df <- pvalue_df %>%
  rownames_to_column() %>%
  rename(genes = rowname,
         pvalue = p.value
  ) %>%
  arrange(pvalue)
rownames(pvalue_df) = NULL


# Rename variables
pvalue_df <- pvalue_df %>%
  mutate(genes = gsub("N:GEXP:", "", genes)
  )

# Export
writexl::write_xlsx(pvalue_df, "Textures/Images/Stroma_WithNormal/Gexp/gexp_genes.xlsx")


# GSEA
Ranks_tmp <- pvalue_df %>% dplyr::mutate(FC = log(median1/median2)) %>% dplyr::filter(!is.na(FC) & pvalue<0.05 & ((median1+median2)/2)>8)
# Ranks_tmp <- pvalue_df %>% dplyr::mutate(FC = 1/pvalue) %>% dplyr::filter(!is.na(FC))
Ranks <- Ranks_tmp %>% dplyr::select(FC)
Ranks <- setNames(object = Ranks$FC, nm = Ranks_tmp$genes)
names(pathways$genesets) <- pathways$geneset.names
fgseaRes_Stroma <- fgsea(pathways = pathways$genesets,
                         stats = Ranks,
                         minSize=10,
                         maxSize=500,
                         nperm=10000)
# Export
writexl::write_xlsx(fgseaRes_Stroma, "Textures/Images/Stroma_WithNormal/Gexp/gexp_gene_pathways.xlsx")
# top 6 enriched pathways
head(fgseaRes_Stroma[order(pval), ])


topPathwaysUp <- fgseaRes_Stroma[ES > 0][padj < 0.05][head(order(-NES), n=10), pathway]
topPathwaysDown <- fgseaRes_Stroma[ES < 0][padj < 0.05][head(order(NES), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

png("Textures/Images/Stroma_WithNormal/Gexp/Stroma_gexp_top20_selectedpathways.png", width = 16, height = 4, units = 'in', res = 300)
plotGseaTable(pathways$genesets[topPathways], Ranks, fgseaRes_Stroma,
              colwidths = c(5, 4, 0.7, 0.7, 0.7),
              gseaParam = 0.5)
dev.off()



# Barplot
fgseaRes_res2 <- fgseaRes_Stroma[fgseaRes_Stroma$pathway %in% topPathways] %>%
  dplyr::mutate(pathway = gsub("_", " ", pathway),
                Sig = ifelse(padj < 0.001, "***",
                             ifelse(padj < 0.01, "**",
                                    ifelse(padj < 0.05, "*", "p<0.1"))),
                Sig = factor(Sig, levels = c("p<0.1", "*", "**", "***")))
g <- ggplot(data = fgseaRes_res2, aes(x = reorder(pathway, NES), y = NES, fill = Sig)) +
  geom_bar(stat = "identity", width = 0.8, color="black") +
  # geom_text(data = flow_median[flow_median$Disease=="AML",], aes(label = round(BM,1), x = name, angle = 0, y = -75), fontface="bold", size=4) + # angle=0,
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
ggsave(plot = g, filename = "Textures/Images/Stroma_WithNormal/Gexp/Stroma_gexp_top20_selectedpathways_barplot.png", width = max(5.5, 5+max(nchar(fgseaRes_res2$pathway))/23), height = max(2, (1+nrow(fgseaRes_res2)/5.5)), units = 'in', dpi = 300)




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
# png(paste0("/Users/oscarbruck/OneDrive - University of Helsinki/Tutkimus/Projekteja/MDS/Analysis/HE/Results/Volcanoplot_mIHC_singlemarker_", dico, ".png"),
# width = 9, height = 9, units = 'in', res = 300)
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
  geom_label_repel(data=pvalue_df1[pvalue_df1$pvalue<0.00001 & pvalue_df1$fold_change1 != 1,],
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
  geom_label_repel(data=pvalue_df1[pvalue_df1$pvalue<0.000001 & pvalue_df1$fold_change1 != 1,],
                   # label.padding = 0.5,
                   segment.size = 1,
                   aes(label=genes), size = 5, fontface = "bold")
ggsave(plot = g1, filename = "Textures/Images/Stroma_WithNormal/Gexp/Stroma_gexp_volcano_p000001.png", width = 7, height = 7, units = 'in', dpi = 300)
ggsave(plot = g2, filename = "Textures/Images/Stroma_WithNormal/Gexp/Stroma_gexp_volcano_p0000001.png", width = 7, height = 7, units = 'in', dpi = 300)







##################################
##### STROMA W AND WO NORMAL #####
##################################


# Read gexp results data
pvalue_df_w_normal <- read_xlsx("Textures/Images/Stroma_WithNormal/Gexp/gexp_genes.xlsx") %>%
  dplyr::mutate(FC = log(median1/median2))
colnames(pvalue_df_w_normal)[2:ncol(pvalue_df_w_normal)] <- paste0(colnames(pvalue_df_w_normal)[2:ncol(pvalue_df_w_normal)], "_w_normal")
pvalue_df_wo_normal <- read_xlsx("Textures/Images/Stroma/Gexp/gexp_genes.xlsx") %>%
  dplyr::mutate(FC = log(median1/median2))
colnames(pvalue_df_wo_normal)[2:ncol(pvalue_df_wo_normal)] <- paste0(colnames(pvalue_df_wo_normal)[2:ncol(pvalue_df_wo_normal)], "_wo_normal")

# Join
pvalue_df <- full_join(pvalue_df_w_normal, pvalue_df_wo_normal)

# Modify
pvalue_df <- pvalue_df %>%
  dplyr::mutate(
    # Significance = ifelse(p_adjusted_w_normal < 0.05 & p_adjusted_wo_normal < 0.05, "#4daf4a",
    #                       ifelse(p_adjusted_w_normal < 0.05, "#e41a1c",
    #                              ifelse(p_adjusted_wo_normal < 0.05, "#377eb8", "Black"))),
    Significance = ifelse(p_adjusted_w_normal < 0.05 & p_adjusted_wo_normal < 0.10, "N+&N-",
                          ifelse(p_adjusted_w_normal < 0.05, "N+",
                                 ifelse(p_adjusted_wo_normal < 0.10, "N-", "AdjP>0.10"))),
    Significance = factor(Significance, levels = c("N+", "N-", "N+&N-", "AdjP>0.10")),
    Sig = ifelse(Significance == "AdjP>0.10", 1, 21),
    Sig = factor(Sig)
  )

# Find 20 genes with lowest p value
a1 <- pvalue_df %>% 
  dplyr::filter(p_adjusted_w_normal < 0.05) %>%
  arrange(desc(FC_w_normal)) %>%
  slice(1:10)
a2 <- pvalue_df %>% 
  dplyr::filter(p_adjusted_wo_normal < 0.10) %>%
  arrange(desc(FC_wo_normal)) %>%
  slice(1:10)
a3 <- pvalue_df %>% 
  dplyr::filter(p_adjusted_w_normal < 0.05) %>%
  arrange(desc(-FC_w_normal)) %>%
  slice(1:10)
a4 <- pvalue_df %>% 
  dplyr::filter(p_adjusted_wo_normal < 0.10) %>%
  arrange(desc(-FC_wo_normal)) %>%
  slice(1:10)
pvalue_df <- pvalue_df %>%
  dplyr::mutate(
    Label = ifelse(genes %in% c(a1$genes, a2$genes, a3$genes, a4$genes), 1, 0)
  )

# Plot
g1 <- ggplot(pvalue_df, aes(FC_w_normal, FC_wo_normal, fill=Significance)) +  # , order=Significance
  geom_point(aes(size=Sig, shape = Sig), colour = "black") +
  # geom_label_repel(data=pvalue_df[c(pvalue_df$p_adjusted_w_normal<0.001 | pvalue_df$p_adjusted_wo_normal<0.001),],
  geom_label_repel(data=pvalue_df[pvalue_df$Label==1,],
                   # label.padding = 0.25,
                   segment.size = 1,
                   # nudge_x = 0.01,
                   # nudge_y = 0.01,
                   show.legend = FALSE,
                   min.segment.length = 0, # draw all line segments
                   # box.padding = 1,
                   aes(label=genes, color = Significance),
                   size = 5, fontface = "bold", segment.color = "black", fill = scales::alpha(c("white"), 0.85)) +
  geom_abline(intercept = 0, slope = -1, size = 1, linetype = 2) +
  scale_fill_brewer(palette = "Set1") +
  ylim(-0.15, 0.15) +
  xlim(-0.2, 0.2) +
  scale_color_brewer(palette = "Set1") +
  # scale_fill_manual(name="Significance"),
  # values=c("#4daf4a", "#377eb8", "#e41a1c", "Black"),
  # labels=c("With Normal p.adj<0.05", "Wo Normal p.adj<0.05", "W&Wo Normal p.adj<0.05", "p.adj>0.05")) +
  scale_size_manual(values=c(0.6, 4)) +
  scale_shape_manual(values=c(16, 21)) +
  labs(x="LOG10 fold change in gene expression", y="LOG10 fold change in gene expression") +
  guides(fill = guide_legend(title = "Sig", nrow=2, byrow=TRUE, override.aes = list(size = 6, shape = 21, fill = c("#e41a1c", "#377eb8", "#4daf4a", "Black"))),
         color = FALSE, size = FALSE, shape = FALSE) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.text.x = element_text(size = 15, face = "bold"),
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 15, face = "bold"),
        legend.direction = 'horizontal',
        # legend.box="vertical", 
        legend.position = 'bottom',
        legend.margin=margin(),
        legend.key = element_rect(size = 5),
        legend.key.size = unit(1, 'lines'))
ggsave(plot = g1, filename = "Textures/Images/Stroma/Stroma_w_vs_wo_normal_scatter.png", width = 6, height = 7, units = 'in', dpi = 300)

