rm(list=ls())

print("Start Textures/Scripts/Margin/Margin_ly_barplot.R")

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
tcga_kirc$`texture_normal_%` <- 100*tcga_kirc$texture_normal / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)


tcga_kirc <- tcga_kirc %>%
  dplyr::filter(`texture_cancer_%` > 5) %>%
  dplyr::filter(`texture_normal_%` >= 1) %>%
  dplyr::mutate(tissue_source_site = gsub("-[[:print:]]{4}", "", gsub("TCGA-", "", tcga_id)))



############################# PLOT ##########################################################################################################


# Rename
colnames(tcga_kirc)[grep(pattern = "inf_margin_bin_lymphocytes_[[:print:]]", colnames(tcga_kirc))] <- c("Ly_Blood", "Ly_Normal", "Ly_Stroma", "Ly_Other")
textures <- c("Ly_Blood", "Ly_Normal", "Ly_Stroma", "Ly_Other")


# Proportions
# Normalize in all textures (= values to equal 100%)
tcga_kirc$Ly_Blood1 <- 100*tcga_kirc$Ly_Blood / (tcga_kirc$Ly_Blood + tcga_kirc$Ly_Normal + tcga_kirc$Ly_Stroma + tcga_kirc$Ly_Other)
tcga_kirc$Ly_Normal1 <- 100*tcga_kirc$Ly_Normal / (tcga_kirc$Ly_Blood + tcga_kirc$Ly_Normal + tcga_kirc$Ly_Stroma + tcga_kirc$Ly_Other)
tcga_kirc$Ly_Stroma1 <- 100*tcga_kirc$Ly_Stroma / (tcga_kirc$Ly_Blood + tcga_kirc$Ly_Normal + tcga_kirc$Ly_Stroma + tcga_kirc$Ly_Other)
tcga_kirc$Ly_Other1 <- 100*tcga_kirc$Ly_Other / (tcga_kirc$Ly_Blood + tcga_kirc$Ly_Normal + tcga_kirc$Ly_Stroma + tcga_kirc$Ly_Other)

# Replace
tcga_kirc <- tcga_kirc %>%
  dplyr::select(-c(Ly_Other, Ly_Normal, Ly_Stroma, Ly_Blood)) %>%
  dplyr::rename(
    Ly_Other = Ly_Other1,
    Ly_Normal = Ly_Normal1,
    Ly_Stroma = Ly_Stroma1,
    Ly_Blood = Ly_Blood1)

a <- 100/(tcga_kirc$Ly_Blood + tcga_kirc$Ly_Normal + tcga_kirc$Ly_Stroma + tcga_kirc$Ly_Other)
tcga_kirc[textures] <- sapply(tcga_kirc[textures], function(x) x*a)


# Melt longer
tcga_kirc_long <- tcga_kirc %>%
  dplyr::select("tcga_id", "Ly_Blood", "Ly_Normal", "Ly_Stroma", "Ly_Other") %>%
  dplyr::rename("ID" = "tcga_id") %>%
  reshape2::melt() %>%
  dplyr::rename("Texture" = "variable") %>%
  arrange(desc(value))

# Order
tcga_kirc_long <- tcga_kirc_long %>%
  dplyr::mutate(orderid = ifelse(Texture=="Ly_Stroma", value, NA)) %>%
  dplyr::filter(!is.na(orderid)) %>%
  dplyr::select(ID, orderid) %>%
  dplyr::right_join(tcga_kirc_long) %>%
  dplyr::mutate(Texture = gsub("Ly_", "", Texture))

# Plot
g <- ggplot(tcga_kirc_long, aes(x = reorder(ID, -orderid), y = value, fill = Texture)) +
  geom_bar(stat = "identity") +
  labs(y="Proportion (%)") +
  scale_y_continuous(limits = c(0,100.1), expand = c(0, 0)) +
  # scale_fill_brewer(palette="Set1")
  scale_fill_manual(values=c("#e41a1c","#4daf4a","#984ea3","#ff7f00")) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title.y=element_text(size=14, face="bold", colour = "black"),
        axis.title.x=element_blank(),
        legend.title = element_text(size=14, face="bold", colour = "black"),
        legend.text = element_text(size=12, colour = "black"),
        legend.key = element_rect(color="black"),
        legend.position = "bottom")
ggsave(plot = g, filename = "Textures/Images/Margin/TCGA_margin_ly_barplot_relative.png", width = 18, height = 5, units = 'in', dpi = 300)


g <- ggplot(tcga_kirc, aes(x = Ly_Normal, y = Ly_Stroma)) +
  geom_point(size=3, color = "black") +
  geom_smooth(method=lm) +
  labs(y="Lymphocytes in Stroma (%)", x="Lymphocytes in Normal (%)") +
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
ggsave(plot = g, filename = "Textures/Images/Margin/TCGA_ly_scatterplot_relative_cancer_normal.png", width = 5, height = 5, units = 'in', dpi = 300)



# Adjust p values
pairwise.test = tcga_kirc_long %>% 
  wilcox_test(value~Texture, paired = TRUE) %>% 
  adjust_pvalue(method = 'BH') %>%
  # mutate(p.adj = round(p.adj, 2))
  mutate(p.adj = ifelse(p.adj < 0.001, "***",
                        ifelse(p.adj < 0.01, "**",
                               ifelse(p.adj < 0.05, "*", round(p.adj, 2)))))

a <- max(tcga_kirc_long$value)

g <- ggplot(tcga_kirc_long, aes(x = Texture, y = value)) +
  geom_jitter(size=5, width = 0.2, aes(fill=Texture), shape = 21, color = "black") +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
  labs(y="Proportion (%)") +
  # scale_y_continuous(limits = c(0,100), expand = c(0, 0)) +
  # scale_fill_brewer(palette="Set1")
  scale_fill_manual(values=c("#e41a1c","#4daf4a","#984ea3","#ff7f00")) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, colour = "black"),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title=element_text(size=14, face="bold", colour = "black"),
        legend.position = "none") +
  # scale_fill_brewer(palette = c("Set1")) +
  # scale_x_discrete(labels=c("0" = "Low", "1" = "High")) +
  stat_compare_means(method = "kruskal.test",
                     label.y = 1.45*a,
                     size = 6) +
  stat_pvalue_manual(
    pairwise.test,
    label = "p.adj",
    bracket.size = 1.5,
    size = 5,
    # y.position = c(1.12*a, 1.24*a, 1.3*a, 1.0*a, 1.18*a, 1.06*a)
    y.position = c(1.14*a, 1.28*a, 1.35*a, 1.0*a, 1.21*a, 1.07*a)
  )
ggsave(plot = g, filename = "Textures/Images/Margin/TCGA_margin_ly_scatterplot.png", width = 6, height = 6, dpi = 300)



############################# MARGIN VS. SAMPLE TEXTURE ########################################################################


# Image data
tcga_kirc <- read_xlsx("../data/image_analysis_results_final.xlsx")

# Normalize by removing empty
tcga_kirc$`texture_cancer_%` <- 100*tcga_kirc$texture_cancer / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`texture_normal_%` <- 100*tcga_kirc$texture_normal / (tcga_kirc$texture_blood + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
# tcga_kirc$inf_bin_lymphocytes_blood <- 100*tcga_kirc$inf_bin_lymphocytes_blood / (tcga_kirc$inf_bin_lymphocytes_blood + tcga_kirc$inf_bin_lymphocytes_normal + tcga_kirc$inf_bin_lymphocytes_stroma + tcga_kirc$inf_bin_lymphocytes_other)
# tcga_kirc$inf_bin_lymphocytes_normal <- 100*tcga_kirc$inf_bin_lymphocytes_normal / (tcga_kirc$inf_bin_lymphocytes_blood + tcga_kirc$inf_bin_lymphocytes_normal + tcga_kirc$inf_bin_lymphocytes_stroma + tcga_kirc$inf_bin_lymphocytes_other)
# tcga_kirc$inf_bin_lymphocytes_stroma <- 100*tcga_kirc$inf_bin_lymphocytes_stroma / (tcga_kirc$inf_bin_lymphocytes_blood + tcga_kirc$inf_bin_lymphocytes_normal + tcga_kirc$inf_bin_lymphocytes_stroma + tcga_kirc$inf_bin_lymphocytes_other)
# tcga_kirc$inf_bin_lymphocytes_other <- 100*tcga_kirc$inf_bin_lymphocytes_other / (tcga_kirc$inf_bin_lymphocytes_blood + tcga_kirc$inf_bin_lymphocytes_normal + tcga_kirc$inf_bin_lymphocytes_stroma + tcga_kirc$inf_bin_lymphocytes_other)
# tcga_kirc$inf_margin_bin_lymphocytes_blood <- 100*tcga_kirc$inf_margin_bin_lymphocytes_blood / (tcga_kirc$inf_margin_bin_lymphocytes_blood + tcga_kirc$inf_margin_bin_lymphocytes_normal + tcga_kirc$inf_margin_bin_lymphocytes_stroma + tcga_kirc$inf_margin_bin_lymphocytes_other)
# tcga_kirc$inf_margin_bin_lymphocytes_normal <- 100*tcga_kirc$inf_margin_bin_lymphocytes_normal / (tcga_kirc$inf_margin_bin_lymphocytes_blood + tcga_kirc$inf_margin_bin_lymphocytes_normal + tcga_kirc$inf_margin_bin_lymphocytes_stroma + tcga_kirc$inf_margin_bin_lymphocytes_other)
# tcga_kirc$inf_margin_bin_lymphocytes_stroma <- 100*tcga_kirc$inf_margin_bin_lymphocytes_stroma / (tcga_kirc$inf_margin_bin_lymphocytes_blood + tcga_kirc$inf_margin_bin_lymphocytes_normal + tcga_kirc$inf_margin_bin_lymphocytes_stroma + tcga_kirc$inf_margin_bin_lymphocytes_other)
# tcga_kirc$inf_margin_bin_lymphocytes_other <- 100*tcga_kirc$inf_margin_bin_lymphocytes_other / (tcga_kirc$inf_margin_bin_lymphocytes_blood + tcga_kirc$inf_margin_bin_lymphocytes_normal + tcga_kirc$inf_margin_bin_lymphocytes_stroma + tcga_kirc$inf_margin_bin_lymphocytes_other)
# tcga_kirc$inf_non_margin_bin_lymphocytes_blood <- 100*tcga_kirc$inf_non_margin_bin_lymphocytes_blood / (tcga_kirc$inf_non_margin_bin_lymphocytes_blood + tcga_kirc$inf_non_margin_bin_lymphocytes_normal + tcga_kirc$inf_non_margin_bin_lymphocytes_stroma + tcga_kirc$inf_non_margin_bin_lymphocytes_other)
# tcga_kirc$inf_non_margin_bin_lymphocytes_normal <- 100*tcga_kirc$inf_non_margin_bin_lymphocytes_normal / (tcga_kirc$inf_non_margin_bin_lymphocytes_blood + tcga_kirc$inf_non_margin_bin_lymphocytes_normal + tcga_kirc$inf_non_margin_bin_lymphocytes_stroma + tcga_kirc$inf_non_margin_bin_lymphocytes_other)
# tcga_kirc$inf_non_margin_bin_lymphocytes_stroma <- 100*tcga_kirc$inf_non_margin_bin_lymphocytes_stroma / (tcga_kirc$inf_non_margin_bin_lymphocytes_blood + tcga_kirc$inf_non_margin_bin_lymphocytes_normal + tcga_kirc$inf_non_margin_bin_lymphocytes_stroma + tcga_kirc$inf_non_margin_bin_lymphocytes_other)
# tcga_kirc$inf_non_margin_bin_lymphocytes_other <- 100*tcga_kirc$inf_non_margin_bin_lymphocytes_other / (tcga_kirc$inf_non_margin_bin_lymphocytes_blood + tcga_kirc$inf_non_margin_bin_lymphocytes_normal + tcga_kirc$inf_non_margin_bin_lymphocytes_stroma + tcga_kirc$inf_non_margin_bin_lymphocytes_other)


tcga_kirc <- tcga_kirc %>%
  dplyr::filter(`texture_cancer_%` > 5) %>%
  dplyr::filter(`texture_normal_%` > 1) %>%
  dplyr::mutate(tissue_source_site = gsub("-[[:print:]]{4}", "", gsub("TCGA-", "", tcga_id)))


############################# PLOT ##########################################################################################################


# Rename
colnames(tcga_kirc)[grep(pattern = "^inf_bin_lymphocytes_", colnames(tcga_kirc))] <- c("Blood", "Cancer", "Normal", "Stroma", "Other")
colnames(tcga_kirc)[grep(pattern = "^inf_margin_bin_lymphocytes_", colnames(tcga_kirc))] <- c("Blood_Margin", "Normal_Margin", "Stroma_Margin", "Other_Margin")
colnames(tcga_kirc)[grep(pattern = "^inf_non_margin_bin_lymphocytes_", colnames(tcga_kirc))] <- c("Blood_NonMargin", "Cancer_NonMargin", "Normal_NonMargin", "Stroma_NonMargin", "Other_NonMargin")
textures <- c("Blood", "Normal", "Stroma", "Other", "Blood_Margin", "Normal_Margin", "Stroma_Margin", "Other_Margin")

# Multiply to 100%
tcga_kirc[c("Blood", "Cancer", "Normal", "Stroma", "Other", "Blood_Margin", "Normal_Margin", "Stroma_Margin", "Other_Margin", "Blood_NonMargin", "Cancer_NonMargin", "Normal_NonMargin", "Stroma_NonMargin", "Other_NonMargin")] <- sapply(tcga_kirc[c("Blood", "Cancer", "Normal", "Stroma", "Other", "Blood_Margin", "Normal_Margin", "Stroma_Margin", "Other_Margin", "Blood_NonMargin", "Cancer_NonMargin", "Normal_NonMargin", "Stroma_NonMargin", "Other_NonMargin")], function(x) x*100)


##### MARGIN VS WHOLE SLIDE #####


# Melt longer
tcga_kirc_long <- tcga_kirc %>%
  dplyr::select("tcga_id", "Blood", "Normal", "Stroma", "Other", "Blood_Margin", "Normal_Margin", "Stroma_Margin", "Other_Margin") %>%
  dplyr::rename("ID" = "tcga_id") %>%
  reshape2::melt() %>%
  dplyr::rename("Texture" = "variable") %>%
  arrange(desc(value))

# Prepare
tcga_kirc_long <- tcga_kirc_long %>%
  dplyr::mutate(orderid = ifelse(Texture=="Stroma", value, NA)) %>%
  dplyr::filter(!is.na(orderid)) %>%
  dplyr::select(ID, orderid) %>%
  dplyr::right_join(tcga_kirc_long) %>%
  dplyr::mutate(
    Class = gsub("_Margin", "", Texture),
    Region = ifelse(str_detect(Texture, "Margin"), "Margin", "Whole"),
    Region1 = ifelse(str_detect(Texture, "Margin"), 1, 0)
  )

# Summarise
tcga_kirc_long1 <- tcga_kirc_long %>%
  group_by(Texture, Class, Region, Region1) %>%
  summarise(value = median(value, na.rm=TRUE),
            iqr = IQR(value, na.rm=TRUE)) %>%
  ungroup()

# Statistics
ptest = tcga_kirc_long %>%
  group_by(Class) %>%
  wilcox_test(value~Region, paired = FALSE) %>%
  add_xy_position(x = "supp") %>%
  ungroup()
ptest

# Plot
g <- ggplot() +
  geom_bar(data = tcga_kirc_long1, aes(x = reorder(Texture, -Region1), y = value, fill = Region), stat = "identity", color="black", position=position_dodge()) +
  geom_errorbar(data = tcga_kirc_long1, aes(x = reorder(Texture, -Region1), ymin=value-iqr, ymax=value+iqr), width=0.2, size = 2,
                position=position_dodge(.9)) +
  labs(y="Proportion (%)") +
  scale_y_continuous(limits = c(0, ceiling(max(tcga_kirc_long1$value, na.rm=TRUE)/10)*10+10 ), expand = c(0, 0)) +
  scale_fill_brewer(palette="Set1") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title.y=element_text(size=14, face="bold", colour = "black"),
        axis.title.x=element_blank(),
        legend.title = element_text(size=14, face="bold", colour = "black"),
        legend.text = element_text(size=12, colour = "black"),
        legend.key = element_rect(color="black"),
        legend.position = "bottom") +
  facet_wrap("Class", scales = "free_x") +
  scale_x_discrete(breaks = NULL)
ggsave(plot = g, filename = "Textures/Images/Margin/TCGA_margin_ly_sample_vs_margin_barplot.png", width = 5, height = 5.5, units = 'in', dpi = 300)



##### MARGIN VS NONMARGIN #####
textures2 <- c("Blood_Margin", "Normal_Margin", "Stroma_Margin", "Other_Margin", "Blood_NonMargin", "Normal_NonMargin", "Stroma_NonMargin", "Other_NonMargin")


# Melt longer
tcga_kirc_long <- tcga_kirc %>%
  dplyr::select("tcga_id", "Blood_Margin", "Normal_Margin", "Stroma_Margin", "Other_Margin", "Blood_NonMargin", "Normal_NonMargin", "Stroma_NonMargin", "Other_NonMargin") %>%
  dplyr::rename("ID" = "tcga_id") %>%
  reshape2::melt() %>%
  dplyr::rename("Texture" = "variable") %>%
  arrange(desc(value))

# Prepare
tcga_kirc_long <- tcga_kirc_long %>%
  # dplyr::mutate(orderid = ifelse(Texture=="Stroma", value, NA)) %>%
  # dplyr::filter(!is.na(orderid)) %>%
  # dplyr::select(ID, orderid) %>%
  dplyr::right_join(tcga_kirc_long) %>%
  dplyr::mutate(
    Class = gsub("_Margin", "", gsub("_NonMargin", "", Texture)),
    Region = ifelse(str_detect(Texture, "NonMargin"), "NonMargin", "Margin"),
    Region1 = ifelse(str_detect(Texture, "Margin"), 1, 0)
  )

# Summarise
tcga_kirc_long1 <- tcga_kirc_long %>%
  group_by(Texture, Class, Region, Region1) %>%
  summarise(value = median(value, na.rm=TRUE),
            iqr = IQR(value, na.rm=TRUE)) %>%
  ungroup()

# Statistics
ptest = tcga_kirc_long %>%
  group_by(Class) %>%
  wilcox_test(value~Region, paired = FALSE) %>%
  add_xy_position(x = "supp") %>%
  ungroup()
ptest

# Plot
g <- ggplot() +
  geom_bar(data = tcga_kirc_long1, aes(x = reorder(Texture, -Region1), y = value, fill = Region), stat = "identity", color="black", position=position_dodge()) +
  geom_errorbar(data = tcga_kirc_long1, aes(x = reorder(Texture, -Region1), ymin=value-iqr, ymax=value+iqr), width=0.2, size = 2,
                position=position_dodge(.9)) +
  labs(y="Proportion (%)") +
  scale_y_continuous(limits = c(0, ceiling(max(tcga_kirc_long1$value, na.rm=TRUE)/10)*10+10 ), expand = c(0, 0)) +
  scale_fill_brewer(palette="Set1") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title.y=element_text(size=14, face="bold", colour = "black"),
        axis.title.x=element_blank(),
        legend.title = element_text(size=14, face="bold", colour = "black"),
        legend.text = element_text(size=12, colour = "black"),
        legend.key = element_rect(color="black"),
        legend.position = "bottom") +
  facet_wrap("Class", scales = "free_x") +
  scale_x_discrete(breaks = NULL)
ggsave(plot = g, filename = "Textures/Images/Margin/TCGA_margin_ly_margin_vs_nonmargin_barplot.png", width = 5, height = 5.5, units = 'in', dpi = 300)
