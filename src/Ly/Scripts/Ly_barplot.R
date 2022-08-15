rm(list=ls())

print("Start Ly/Scripts/Ly_barplot.R")

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
tcga_kirc$`texture_cancer_%` <- 100*tcga_kirc$texture_cancer / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)

tcga_kirc <- tcga_kirc %>%
  dplyr::filter(`texture_cancer_%` > 5) %>%
  dplyr::mutate(tissue_source_site = gsub("-[[:print:]]{4}", "", gsub("TCGA-", "", tcga_id)))


############################# PLOT ##########################################################################################################


# Rename
colnames(tcga_kirc)[grep(pattern = "^inf_bin_lymphocytes_[[:print:]]*", colnames(tcga_kirc))] <- c("Ly_Blood", "Ly_Cancer", "Ly_Normal", "Ly_Stroma", "Ly_Other")
textures <- c("Ly_Blood", "Ly_Cancer", "Ly_Normal", "Ly_Stroma", "Ly_Other")

# Proportions
# Normalize in all textures (= values to equal 100%)
tcga_kirc$Ly_Blood1 <- 100*tcga_kirc$Ly_Blood / (tcga_kirc$Ly_Blood + tcga_kirc$Ly_Cancer + tcga_kirc$Ly_Normal + tcga_kirc$Ly_Stroma + tcga_kirc$Ly_Other)
tcga_kirc$Ly_Cancer1 <- 100*tcga_kirc$Ly_Cancer / (tcga_kirc$Ly_Blood + tcga_kirc$Ly_Cancer + tcga_kirc$Ly_Normal + tcga_kirc$Ly_Stroma + tcga_kirc$Ly_Other)
tcga_kirc$Ly_Normal1 <- 100*tcga_kirc$Ly_Normal / (tcga_kirc$Ly_Blood + tcga_kirc$Ly_Cancer + tcga_kirc$Ly_Normal + tcga_kirc$Ly_Stroma + tcga_kirc$Ly_Other)
tcga_kirc$Ly_Stroma1 <- 100*tcga_kirc$Ly_Stroma / (tcga_kirc$Ly_Blood + tcga_kirc$Ly_Cancer + tcga_kirc$Ly_Normal + tcga_kirc$Ly_Stroma + tcga_kirc$Ly_Other)
tcga_kirc$Ly_Other1 <- 100*tcga_kirc$Ly_Other / (tcga_kirc$Ly_Blood + tcga_kirc$Ly_Cancer + tcga_kirc$Ly_Normal + tcga_kirc$Ly_Stroma + tcga_kirc$Ly_Other)

# Replace
tcga_kirc <- tcga_kirc %>%
  dplyr::select(-c(Ly_Other, Ly_Cancer, Ly_Normal, Ly_Stroma, Ly_Blood)) %>%
  dplyr::rename(
    Ly_Other = Ly_Other1,
    Ly_Cancer = Ly_Cancer1,
    Ly_Normal = Ly_Normal1,
    Ly_Stroma = Ly_Stroma1,
    Ly_Blood = Ly_Blood1)

a <- 100/(tcga_kirc$Ly_Blood + tcga_kirc$Ly_Cancer + tcga_kirc$Ly_Normal + tcga_kirc$Ly_Stroma + tcga_kirc$Ly_Other)
tcga_kirc[textures] <- sapply(tcga_kirc[textures], function(x) x*a)


# Melt longer
tcga_kirc_long <- tcga_kirc %>%
  dplyr::select("tcga_id", "Ly_Blood", "Ly_Cancer", "Ly_Normal", "Ly_Stroma", "Ly_Other") %>%
  dplyr::rename("ID" = "tcga_id") %>%
  reshape2::melt() %>%
  dplyr::rename("Lymphocytes" = "variable") %>%
  arrange(desc(value))

# Order
tcga_kirc_long <- tcga_kirc_long %>%
  dplyr::mutate(orderid = ifelse(Lymphocytes=="Ly_Cancer", value, NA)) %>%
  dplyr::filter(!is.na(orderid)) %>%
  dplyr::select(ID, orderid) %>%
  dplyr::right_join(tcga_kirc_long) %>%
  dplyr::mutate(Lymphocytes = gsub("Ly_", "", Lymphocytes))

# Correlation plots
# pairs(tcga_kirc[textures])
# corr.test(tcga_kirc[textures], adjust = "BH", ci=F, method = "spearman")
# cor.test(x = tcga_kirc$Ly_Normal, y = tcga_kirc$Ly_Other, method = "spearman")

# Plot
g <- ggplot(tcga_kirc_long, aes(x = reorder(ID, -orderid), y = value, fill = Lymphocytes)) +
  geom_bar(stat = "identity") +
  labs(y="Lymphocyte proportion (%)") +
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
ggsave(plot = g, filename = "Ly/Images/TCGA_ly_barplot_relative.png", width = 18, height = 5, units = 'in', dpi = 300) #original pointsize = 12


g <- ggplot(tcga_kirc, aes(x = Ly_Normal, y = Ly_Cancer)) +
  geom_point(size=3, color = "black") +
  geom_smooth(method=lm) +
  labs(y="Lymphocytes in Cancer (%)", x="Lymphocytes in Normal (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, colour = "black"),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title.y = element_text(size=14, face="bold", colour = "black"),
        axis.title.x = element_text(size=14, face="bold", colour = "black"),
        axis.ticks.y = ,
        legend.position = "none") +
  scale_x_continuous(expand = c(0.02, 0.01), limits = c(0, NA)) +
  scale_y_continuous(expand = c(0.02, 0.01), limits = c(0, NA))
  # ylim(0, max(tcga_kirc$Ly_Cancer, na.rm=TRUE))
  # xlim(0, max(tcga_kirc$Ly_Normal, na.rm=TRUE))
ggsave(plot = g, filename = "Ly/Images/TCGA_ly_scatterplot_relative_cancer_normal.png", width = 5, height = 5, units = 'in', dpi = 300)



############################# TEXTURE-NORMALIZED LYMPHOCYTE PROPORTION ##########################################################################################################


# (Ly (relative) * absolute area of each textures) / total absolute sample area
tcga_kirc$Ly <- (tcga_kirc$Ly_Blood * tcga_kirc$texture_blood + tcga_kirc$Ly_Cancer * tcga_kirc$texture_cancer +
                   tcga_kirc$Ly_Stroma * tcga_kirc$texture_stroma + tcga_kirc$Ly_Normal * tcga_kirc$texture_normal +
                   tcga_kirc$Ly_Other * tcga_kirc$texture_other) /
  (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)

sapply(tcga_kirc["Ly"], quantile)

# TCGA clinical sources
clinical_sources <- read_xlsx("../data/TCGA_source_sites.xlsx") %>%
  dplyr::select(-StudyName, -BCR)

# Join Otso's data and TCGA centers
tcga_kirc <- tcga_kirc %>%
  dplyr::left_join(clinical_sources)

# tcga_kirc1 <- tcga_kirc %>%
#   dplyr::filter(`texture_normal_%` < 0.01)
# tcga_kirc2 <- tcga_kirc %>%
#   dplyr::filter(`texture_normal_%` >= 0.01)
# tcga_kirc <- tcga_kirc2


# Plot
g <- ggplot(tcga_kirc, aes(x=ClinicalCenter, y=Ly)) +
  # geom_rect(xmin=0, xmax=17, ymin=3.82, ymax=17.78, alpha=0.5,fill="grey90") +
  geom_jitter(size=3, width = 0.2, aes(fill=ClinicalCenter), shape = 21, color = "black") +
  annotate("rect", xmin=0, xmax=17, ymin=as.numeric(quantile(tcga_kirc$Ly, 0.25)), ymax=as.numeric(quantile(tcga_kirc$Ly, 0.75)), fill="grey30", alpha=0.15) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_hline(yintercept = median(tcga_kirc$Ly), linetype="dashed", color = "red", size = 1.5) +
  labs(x="Clinical Center", y=paste0("Lymphocyte proportion (%)")) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, colour = "black", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title=element_text(size=14, face="bold", colour = "black"),
        legend.position = "none",
        plot.margin = margin(0.1, 0.1, 0.1, 0.3, "in")) +
  # scale_fill_brewer(palette = c("Set1")) +
  # scale_fill_manual(values = pal1) +
  stat_compare_means(method = "kruskal.test",
                     # label.y = max(tcga_kirc[texture], na.rm = TRUE),
                     label.x = 6,
                     # label.y = a5,
                     size = 5)
ggsave(plot = g, filename = "Clinical_center/Images/TCGA_texture_normalized_ly_center.png",
       width = 8, height = 8, dpi = 300, units = 'in') #original pointsize = 12


# Plot2
tcga_kirc_tmp <- tcga_kirc %>%
  dplyr::filter(str_detect(ClinicalCenter, "Harvard|International|Anderson|MSKCC|Pittsburgh")) %>%
  group_by(ClinicalCenter) %>%
  mutate(Ly_z = factor(ntile(Ly, 2)),
         ClinicalCenter2=ifelse(str_detect(ClinicalCenter, "Harvard"), 1,
                                ifelse(str_detect(ClinicalCenter, "International"), 2,
                                       ifelse(str_detect(ClinicalCenter, "Anderson"), 3,
                                              ifelse(str_detect(ClinicalCenter, "MSKCC"), 4, 5)))),
         ClinicalCenter2=factor(ClinicalCenter2)) %>%
  ungroup()

g <- ggplot(tcga_kirc_tmp, aes(x=ClinicalCenter2, y=Ly)) +
  geom_jitter(size=3, width = 0.2, aes(fill=Ly_z), shape = 21, color = "black") +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  scale_fill_brewer(palette="Set1", direction = -1) +
  labs(x="Clinical Center", y=paste0("Lymphocyte proportion (%)")) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, colour = "black"),
        # axis.text.y = element_text(size=12, colour = "black"),
        # axis.title = element_text(size=14, face="bold", colour = "black"),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none",
        plot.margin = margin(0.1, 0.1, 0.1, 0.3, "in"))
ggsave(plot = g, filename = "Clinical_center/Images/TCGA_texture_normalized_ly_center2.png",
       width = 4, height = 3, dpi = 300, units = 'in') #original pointsize = 12



# Melt longer
tcga_kirc_long <- tcga_kirc %>%
  dplyr::select("tcga_id", "ClinicalCenter", "Ly", "Ly_Cancer")


# Wilcoxon
df <- NULL
for (variable1 in unique(tcga_kirc_long$ClinicalCenter)) {
  # Create tmp variable
  tcga_kirc_long$var1 <- ifelse(tcga_kirc_long$ClinicalCenter == variable1, 1, 0)
  # Median
  medianF = median(tcga_kirc_long[tcga_kirc_long$var1 == 0,]$Ly, na.rm = TRUE)
  medianT = median(tcga_kirc_long[tcga_kirc_long$var1 == 1,]$Ly, na.rm = TRUE)
  # Wilcoxon test
  wilcox_test <- wilcox.test(tcga_kirc_long$Ly ~ tcga_kirc_long$var1, correct=FALSE, exact=FALSE, estimate=TRUE, na.rm=TRUE)
  df_tmp <- data.frame(variable1, wilcox_test$p.value, medianT, medianF)
  # Collect to dataframe
  if (!exists("df")) {
    df <- df_tmp
  } else {
    df <- rbind(df, df_tmp)
  }
  # Remove var1
  tcga_kirc_long$var1 <- NULL
}

# Adjust p values
df$adj <- p.adjust(df$wilcox_test.p.value, method = "BH")

# Export
writexl::write_xlsx(df, "Clinical_center/Images/TCGA_texture_normalized_ly_center.xlsx")


# # Image analysis results
# tcga_kirc <- read_xlsx("../data/image_analysis_results_final.xlsx") %>%
#   dplyr::mutate(tissue_source_site = gsub("-[[:print:]]{4}", "", gsub("TCGA-", "", tcga_id))) %>%
#   dplyr::distinct(tissue_source_site)
# 
# # TCGA clinical sources
# clinical_sources <- read_xlsx("../data/TCGA_source_sites.xlsx") %>%
#   dplyr::select(-StudyName, -BCR) %>%
#   dplyr::filter(tissue_source_site %in% tcga_kirc$tissue_source_site)
# 
# # TCGA Texture-normalized lymphocytes by center
# df <- read_xlsx("Clinical_center/Images/TCGA_texture_normalized_ly_center.xlsx") %>%
#   dplyr::filter(wilcox_test.p.value < 0.05) %>%
#   dplyr::rename(ClinicalCenter = variable1) %>%
#   dplyr::left_join(clinical_sources)
# 
# 
# 
# # Normalize Ly
# tcga_kirc_long1 <- tcga_kirc_long %>%
#   dplyr::left_join(df %>% dplyr::rename(ClinicalCenter = variable1))
# tcga_kirc_long1$Ly_norm <- tcga_kirc_long1$Ly / tcga_kirc_long1$medianT
# 
# # Plot
# g <- ggplot(tcga_kirc_long1, aes_string(x="ClinicalCenter", y="Ly_norm")) +
#   geom_jitter(size=3, width = 0.2, aes(fill=ClinicalCenter), shape = 21, color = "black") +
#   geom_boxplot(outlier.shape = NA, alpha = 0.5) +
#   labs(x="Clinical Center", y=paste0("Texture-normalized proportion of lymphocytes (%)")) +
#   theme_bw() +
#   theme(axis.text.x = element_text(size=12, colour = "black", angle = 45, vjust = 1, hjust = 1),
#         axis.text.y = element_text(size=12, colour = "black"),
#         axis.title=element_text(size=14, face="bold", colour = "black"),
#         legend.position = "none",
#         plot.margin = margin(0.1, 0.1, 0.1, 0.3, "in")) +
#   # scale_fill_brewer(palette = c("Set1")) +
#   stat_compare_means(method = "kruskal.test",
#                      # label.y = max(tcga_kirc[texture], na.rm = TRUE),
#                      label.x = 6,
#                      # label.y = a5,
#                      size = 5)
# ggsave(plot = g, filename = "Clinical_center/Images/TCGA_texture_normalized_ly_center_after_normalization.png", width = 8, height = 8, units = 'in', dpi = 300)


############################# LOAD DATA ##########################################################################################################


# Otso's data
tcga_kirc <- read_xlsx("../data/image_analysis_results_final.xlsx")

# Normalize by removing empty
tcga_kirc$`texture_cancer_%` <- 100*tcga_kirc$texture_cancer / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)

tcga_kirc <- tcga_kirc %>%
  dplyr::filter(`texture_cancer_%` > 5) %>%
  dplyr::mutate(tissue_source_site = gsub("-[[:print:]]{4}", "", gsub("TCGA-", "", tcga_id)))



############################# PLOT ##########################################################################################################


# Rename
colnames(tcga_kirc)[grep(pattern = "^inf_bin_lymphocytes_[[:print:]]*", colnames(tcga_kirc))] <- c("Ly_Blood", "Ly_Cancer", "Ly_Normal", "Ly_Stroma", "Ly_Other")
# colnames(tcga_kirc)[grep(pattern = "^bin_lymphocytes_[[:print:]]*", colnames(tcga_kirc))] <- c("Ly_Blood", "Ly_Cancer", "Ly_Normal", "Ly_Stroma", "Ly_Other")
textures <- c("Ly_Blood", "Ly_Cancer", "Ly_Normal", "Ly_Stroma", "Ly_Other")


# Quantity

# # (Ly (relative) * absolute area of each textures) / total absolute sample area
# tcga_kirc$Ly <- (tcga_kirc$Ly_Blood + tcga_kirc$Ly_Cancer +
#                    tcga_kirc$Ly_Stroma+ tcga_kirc$Ly_Normal +
#                    tcga_kirc$Ly_Other) /
#   (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other) * 100
# 
# # Normalize in all textures (= values to equal 100%)
# # (Ly (relative) * absolute area of each textures) / total absolute sample area
# tcga_kirc$Ly_Blood <- tcga_kirc$Ly_Blood /
#   (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
# tcga_kirc$Ly_Cancer <- tcga_kirc$Ly_Cancer /
#   (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
# tcga_kirc$Ly_Stroma <- tcga_kirc$Ly_Stroma /
#   (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
# tcga_kirc$Ly_Normal <- tcga_kirc$Ly_Normal /
#   (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
# tcga_kirc$Ly_Other <- tcga_kirc$Ly_Other /
#   (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)

# Change to percents
tcga_kirc[textures] <- sapply(tcga_kirc[textures], function(x) x*100)  # 100 = 100%, 5 = 5 textures so max is 500% which looks weird in the plot
sapply(tcga_kirc[textures], quantile)

# Melt longer
tcga_kirc_long <- tcga_kirc %>%
  dplyr::select("tcga_id", "Ly_Blood", "Ly_Cancer", "Ly_Normal", "Ly_Stroma", "Ly_Other") %>%
  dplyr::rename("ID" = "tcga_id") %>%
  reshape2::melt() %>%
  dplyr::rename("Lymphocytes" = "variable") %>%
  arrange(desc(value))

# Order
tcga_kirc_long <- tcga_kirc_long %>%
  group_by(ID) %>%
  summarise(orderid = sum(value, na.rm = TRUE)) %>%
  dplyr::right_join(tcga_kirc_long) %>%
  dplyr::mutate(Lymphocytes = gsub("Ly_", "", Lymphocytes))

  
  
# Plot
g <- ggplot(tcga_kirc_long, aes(x = reorder(ID, -orderid), y = value, fill = Lymphocytes)) +
  geom_bar(stat = "identity") +
  labs(y="Lymphocyte proportion (%)") +
  scale_y_continuous(limits = c(0, max(ceiling(tcga_kirc_long$orderid)/10, 0)*10), expand = c(0, 0)) +
  # scale_fill_brewer(palette="Set1")
  scale_fill_manual(values=c("#4daf4a", "#e41a1c", "#a6cee3", "#984ea3", "#ffff33")) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title.y=element_text(size=14, face="bold", colour = "black"),
        axis.title.x=element_blank(),
        legend.title = element_text(size=14, face="bold", colour = "black"),
        legend.text = element_text(size=12, colour = "black"),
        legend.key = element_rect(color="black"),
        legend.position = "bottom")
ggsave(plot = g, filename = "Ly/Images/TCGA_ly_barplot_quantity.png", width = 18, height = 5, units = 'in', dpi = 300)


# Plot
## Plot parameters
a <- max(tcga_kirc_long$value, na.rm=TRUE)
## Adjust p values
pairwise.test = tcga_kirc_long %>%
  wilcox_test(value~Lymphocytes) %>%
  adjust_pvalue(method = 'BH')

g <- ggplot(tcga_kirc_long, aes(x = Lymphocytes, y = value)) +
  geom_jitter(size=5, width = 0.2, aes(fill=Lymphocytes), shape = 21, color = "black") +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
  labs(y="Lymphocyte proportion (%)", x="Texture") +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, colour = "black"),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title.y = element_text(size=14, face="bold", colour = "black"),
        axis.title.x = element_text(size=14, face="bold", colour = "black"),
        legend.position = "none") +
  scale_fill_brewer(palette = c("Set1")) +
  scale_y_continuous(breaks = c(0, 50, 100), labels = c(0, 50, 100)) +
  stat_compare_means(method = "kruskal.test",
                     label.y = 1.65*a,
                     size = 5) +
  stat_pvalue_manual(
    pairwise.test,
    label = "p.adj.signif",
    bracket.size = 1.0,
    size = 5,
    y.position = c(1.18*a, 1.36*a, 1.48*a,
                   1.54*a, 1.12*a, 1.30*a,
                   1.42*a, 1.06*a, 1.24*a,
                   1.00*a)
  )
ggsave(plot = g, filename = "Ly/Images/TCGA_ly_scatterplot_quantity.png", width = 8, height = 5.5, units = 'in', dpi = 300)
