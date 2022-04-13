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


# Otso's data
tcga_kirc <- read_xlsx("/Users/oscarbruck/OneDrive - University of Helsinki/RCC/Otso/Data/otso/raw_data.xlsx")

# Normalize by removing empty
tcga_kirc$`texture_cancer_%` <- 100*tcga_kirc$texture_cancer / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)

tcga_kirc <- tcga_kirc %>%
  dplyr::filter(`texture_cancer_%` > 5) %>%
  dplyr::filter(is.na(PoorQuality)) %>%
  dplyr::mutate(tissue_source_site = gsub("-[[:print:]]{4}", "", gsub("TCGA-", "", tcga_id)))



############################# PLOT ##########################################################################################################


# Rename
colnames(tcga_kirc)[grep(pattern = "^inf_bin_lymphocytes_[[:print:]]*", colnames(tcga_kirc))] <- c("Ly_Blood", "Ly_Cancer", "Ly_Normal", "Ly_Stroma", "Ly_Other")
textures <- c("Ly_Blood", "Ly_Cancer", "Ly_Normal", "Ly_Stroma", "Ly_Other")

# Proportions
# Normalize in all textures (= values to equal 100%)
tcga_kirc$Ly_Blood <- 100*tcga_kirc$Ly_Blood / (tcga_kirc$Ly_Blood + tcga_kirc$Ly_Cancer + tcga_kirc$Ly_Normal + tcga_kirc$Ly_Stroma + tcga_kirc$Ly_Other)
tcga_kirc$Ly_Cancer <- 100*tcga_kirc$Ly_Cancer / (tcga_kirc$Ly_Blood + tcga_kirc$Ly_Cancer + tcga_kirc$Ly_Normal + tcga_kirc$Ly_Stroma + tcga_kirc$Ly_Other)
tcga_kirc$Ly_Normal <- 100*tcga_kirc$Ly_Normal / (tcga_kirc$Ly_Blood + tcga_kirc$Ly_Cancer + tcga_kirc$Ly_Normal + tcga_kirc$Ly_Stroma + tcga_kirc$Ly_Other)
tcga_kirc$Ly_Stroma <- 100*tcga_kirc$Ly_Stroma / (tcga_kirc$Ly_Blood + tcga_kirc$Ly_Cancer + tcga_kirc$Ly_Normal + tcga_kirc$Ly_Stroma + tcga_kirc$Ly_Other)
tcga_kirc$Ly_Other <- 100*tcga_kirc$Ly_Other / (tcga_kirc$Ly_Blood + tcga_kirc$Ly_Cancer + tcga_kirc$Ly_Normal + tcga_kirc$Ly_Stroma + tcga_kirc$Ly_Other)

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
  dplyr::right_join(tcga_kirc_long)


# Plot
png("/Users/oscarbruck/OneDrive - University of Helsinki/RCC/Otso/Analysis/Ly/Images/TCGA_ly_barplot_relative.png", width = 18, height = 5, units = 'in', res = 300, pointsize = 12) #original pointsize = 12
ggplot(tcga_kirc_long, aes(x = reorder(ID, -orderid), y = value, fill = Lymphocytes)) +
  geom_bar(stat = "identity") +
  labs(y="Ly Density (%)") +
  scale_y_continuous(limits = c(0,100.1), expand = c(0, 0)) +
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
dev.off()



############################# TEXTURE-NORMALIZED LYMPHOCYTE PROPORTION ##########################################################################################################



# (Ly (relative) * absolute area of each textures) / total absolute sample area
tcga_kirc$Ly <- (tcga_kirc$Ly_Blood * tcga_kirc$texture_blood + tcga_kirc$Ly_Cancer * tcga_kirc$texture_cancer +
                   tcga_kirc$Ly_Stroma * tcga_kirc$texture_stroma + tcga_kirc$Ly_Normal * tcga_kirc$texture_normal +
                   tcga_kirc$Ly_Other * tcga_kirc$texture_other) /
  (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)

# TCGA clinical sources
clinical_sources <- read_xlsx("/Users/oscarbruck/OneDrive - University of Helsinki/RCC/Otso/Data/TCGA_source_sites.xlsx") %>%
  dplyr::select(-StudyName, -BCR)

# Join Otso's data and TCGA centers
tcga_kirc <- tcga_kirc %>%
  dplyr::left_join(clinical_sources)


# Plot
png("/Users/oscarbruck/OneDrive - University of Helsinki/RCC/Otso/Analysis/Clinical_center/Images/TCGA_texture_normalized_ly_center.png", width = 8, height = 8, units = 'in', res = 300, pointsize = 12) #original pointsize = 12
ggplot(data=tcga_kirc, aes(x=ClinicalCenter, y=Ly)) +
  # geom_rect(xmin=0, xmax=17, ymin=3.82, ymax=17.78, alpha=0.5,fill="grey90") +
  geom_jitter(size=3, width = 0.2, aes(fill=ClinicalCenter), shape = 21, color = "black") +
  annotate("rect", xmin=0, xmax=17, ymin=as.numeric(quantile(tcga_kirc$Ly, 0.25)), ymax=as.numeric(quantile(tcga_kirc$Ly, 0.75)), fill="grey30", alpha=0.15) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_hline(yintercept = median(tcga_kirc$Ly), linetype="dashed", color = "red", size = 1.5) +
  labs(x="Clinical Center", y=paste0("Proportion of lymphocytes (%)")) +
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
dev.off()


# Melt longer
tcga_kirc_long <- tcga_kirc %>%
  dplyr::select("tcga_id", "ClinicalCenter", "Ly")
#   reshape2::melt() %>%
#   dplyr::rename("Lymphocytes" = "variable") %>%
#   arrange(desc(value))
# 
# # Order
# tcga_kirc_long <- tcga_kirc_long %>%
#   group_by(tcga_id) %>%
#   summarise(orderid = sum(value, na.rm = TRUE)) %>%
#   dplyr::right_join(tcga_kirc_long)

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


# Export
writexl::write_xlsx(df, "/Users/oscarbruck/OneDrive - University of Helsinki/RCC/Otso/Analysis/Clinical_center/Images/TCGA_texture_normalized_ly_center.xlsx")

# Normalize Ly
tcga_kirc_long1 <- tcga_kirc_long %>%
  dplyr::left_join(df %>% dplyr::rename(ClinicalCenter = variable1))
tcga_kirc_long1$Ly_norm <- tcga_kirc_long1$Ly / tcga_kirc_long1$medianT

# Plot
png("/Users/oscarbruck/OneDrive - University of Helsinki/RCC/Otso/Analysis/Clinical_center/Images/TCGA_texture_normalized_ly_center_after_normalization.png", width = 8, height = 8, units = 'in', res = 300, pointsize = 12) #original pointsize = 12
ggplot(tcga_kirc_long1, aes_string(x="ClinicalCenter", y="Ly_norm")) +
  geom_jitter(size=3, width = 0.2, aes(fill=ClinicalCenter), shape = 21, color = "black") +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  labs(x="Clinical Center", y=paste0("Texture-normalized proportion of lymphocytes (%)")) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, colour = "black", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title=element_text(size=14, face="bold", colour = "black"),
        legend.position = "none",
        plot.margin = margin(0.1, 0.1, 0.1, 0.3, "in")) +
  # scale_fill_brewer(palette = c("Set1")) +
  stat_compare_means(method = "kruskal.test",
                     # label.y = max(tcga_kirc[texture], na.rm = TRUE),
                     label.x = 6,
                     # label.y = a5,
                     size = 5)
dev.off()

colnames(df2)[1] <- "tcga_id"


# tcga_kirc_long2 <- tcga_kirc_long1 %>% inner_join(df2)
# tcga_kirc_long2 <- tcga_kirc_long2 %>% dplyr::left_join(tcga_kirc %>% dplyr::select(tcga_id, starts_with("texture")))
# tcga_kirc_long2 <- tcga_kirc_long2 %>% dplyr::left_join(tcga_kirc %>% dplyr::select(tcga_id, starts_with("bin_lymphocytes_")))
# tcga_kirc_long2$Ly_c <- tcga_kirc_long2$Ly*(tcga_kirc_long2$texture_cancer)
# tcga_kirc_long2$Ly_norm_c <- tcga_kirc_long2$Ly_norm*(tcga_kirc_long2$texture_cancer)
# tcga_kirc_long2$Ly_norm_c <- tcga_kirc_long2$Ly*(tcga_kirc_long2$texture_cancer + tcga_kirc_long2$texture_normal + tcga_kirc_long2$texture_other + tcga_kirc_long2$texture_blood + tcga_kirc_long2$texture_stroma)
# tcga_kirc_long2$bin_ly <- (tcga_kirc_long2$bin_lymphocytes_blood + tcga_kirc_long2$bin_lymphocytes_cancer + tcga_kirc_long2$bin_lymphocytes_normal + tcga_kirc_long2$bin_lymphocytes_stroma + tcga_kirc_long2$bin_lymphocytes_other) / (tcga_kirc_long2$texture_cancer + tcga_kirc_long2$texture_normal + tcga_kirc_long2$texture_other + tcga_kirc_long2$texture_blood + tcga_kirc_long2$texture_stroma)
# tcga_kirc_long2$bin_ly <- tcga_kirc_long2$bin_lymphocytes_cancer/ tcga_kirc_long2$texture_cancer
# tcga_kirc_long2[,grep(x = colnames(tcga_kirc_long2), pattern = "GEXP")] <- sapply(tcga_kirc_long2[,grep(x = colnames(tcga_kirc_long2), pattern = "GEXP")], as.numeric)
# cyt <- tcga_kirc_long2 %>% dplyr::select(`N:GEXP:GZMB`,`N:GEXP:PRF1`)
# cyt$cyt_score <- exp(log(rowMeans(cyt)))
# tcga_kirc_long2 <- tcga_kirc_long2 %>% left_join(cyt)
# tcga_kirc_long3 <- tcga_kirc_long2 %>% dplyr::filter(!ClinicalCenter %in% c("Fox Chase", "MD Anderson Cancer Center"))
# plot(log(tcga_kirc_long3$bin_ly), log(tcga_kirc_long3$`N:GEXP:GZMB`))
# cor.test(log(tcga_kirc_long3$bin_ly), log(tcga_kirc_long3$`N:GEXP:GZMB`), method = "spearman")
# plot(log(tcga_kirc_long3$bin_ly), log(tcga_kirc_long3$cyt_score))
# cor.test(log(tcga_kirc_long3$bin_ly), log(tcga_kirc_long3$cyt_score), method = "spearman")



############################# LOAD DATA ##########################################################################################################


# Otso's data
tcga_kirc <- read_xlsx("/Users/oscarbruck/OneDrive - University of Helsinki/RCC/Otso/Data/otso/raw_data.xlsx")

# Normalize by removing empty
tcga_kirc$`texture_cancer_%` <- 100*tcga_kirc$texture_cancer / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)

tcga_kirc <- tcga_kirc %>%
  dplyr::filter(`texture_cancer_%` > 5) %>%
  dplyr::filter(is.na(PoorQuality)) %>%
  dplyr::mutate(tissue_source_site = gsub("-[[:print:]]{4}", "", gsub("TCGA-", "", tcga_id)))



############################# PLOT ##########################################################################################################


# Rename
colnames(tcga_kirc)[grep(pattern = "^inf_bin_lymphocytes_[[:print:]]*", colnames(tcga_kirc))] <- c("Ly_Blood", "Ly_Cancer", "Ly_Normal", "Ly_Stroma", "Ly_Other")
textures <- c("Ly_Blood", "Ly_Cancer", "Ly_Normal", "Ly_Stroma", "Ly_Other")


# Quantity
# Normalize in all textures (= values to equal 100%)
tcga_kirc[textures] <- sapply(tcga_kirc[textures], function(x) x*100/5)  # 100 = 100%, 5 = 5 textures so max is 500% which looks weird in the plot



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
  dplyr::right_join(tcga_kirc_long)


# Plot
png("/Users/oscarbruck/OneDrive - University of Helsinki/RCC/Otso/Analysis/Ly/Images/TCGA_ly_barplot_quantity.png", width = 18, height = 5, units = 'in', res = 300, pointsize = 12) #original pointsize = 12
ggplot(tcga_kirc_long, aes(x = reorder(ID, -orderid), y = value, fill = Lymphocytes)) +
  geom_bar(stat = "identity") +
  labs(y="Ly Density (%)") +
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
dev.off()