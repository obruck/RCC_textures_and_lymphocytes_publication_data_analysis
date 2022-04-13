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
tcga_kirc$`margin_texture_blood_%` <- 100*tcga_kirc$margin_texture_blood / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)
tcga_kirc$`margin_texture_cancer_%` <- 100*tcga_kirc$margin_texture_cancer / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)
tcga_kirc$`margin_texture_normal_%` <- 100*tcga_kirc$margin_texture_normal / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)
tcga_kirc$`margin_texture_stroma_%` <- 100*tcga_kirc$margin_texture_stroma / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)
tcga_kirc$`margin_texture_other_%` <- 100*tcga_kirc$margin_texture_other / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)


tcga_kirc <- tcga_kirc %>%
  dplyr::filter(`texture_cancer_%` > 5) %>%
  dplyr::filter(`margin_texture_normal_%` > 1) %>%
  dplyr::filter(is.na(PoorQuality)) %>%
  dplyr::mutate(tissue_source_site = gsub("-[[:print:]]{4}", "", gsub("TCGA-", "", tcga_id)))



############################# PLOT ##########################################################################################################


# Rename
colnames(tcga_kirc)[grep(pattern = "^margin_texture_[[:print:]]*%", colnames(tcga_kirc))] <- c("Blood", "Cancer", "Normal", "Stroma", "Other")
textures <- c("Blood", "Normal", "Stroma", "Other")

# Factor IDs
# tcga_kirc <- tcga_kirc %>% arrange(Cancer) %>% dplyr::mutate(tcga_id = factor(tcga_id))
# tcga_kirc <- tcga_kirc %>% mutate(tcga_id = factor(Cancer, levels = Cancer))

# Melt longer
tcga_kirc_long <- tcga_kirc %>%
  dplyr::select("tcga_id", "Blood", "Normal", "Stroma", "Other") %>%
  dplyr::rename("ID" = "tcga_id") %>%
  reshape2::melt() %>%
  dplyr::rename("Texture" = "variable") %>%
  arrange(desc(value))

# Order
tcga_kirc_long <- tcga_kirc_long %>%
  dplyr::mutate(orderid = ifelse(Texture=="Stroma", value, NA)) %>%
  dplyr::filter(!is.na(orderid)) %>%
  dplyr::select(ID, orderid) %>%
  dplyr::right_join(tcga_kirc_long)

# Plot
png("/Users/oscarbruck/OneDrive - University of Helsinki/RCC/Otso/Analysis/Textures/Images/Margin/TCGA_margin_texture_barplot.png", width = 18, height = 5, units = 'in', res = 300, pointsize = 12) #original pointsize = 12
ggplot(tcga_kirc_long, aes(x = reorder(ID, -orderid), y = value, fill = Texture)) +
  geom_bar(stat = "identity") +
  labs(y="Proportion (%)") +
  scale_y_continuous(limits = c(0,100), expand = c(0, 0)) +
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

# Adjust p values
pairwise.test = tcga_kirc_long %>% 
  wilcox_test(value~Texture, paired = TRUE) %>% 
  adjust_pvalue(method = 'BH') %>%
  # mutate(p.adj = round(p.adj, 2))
  mutate(p.adj = ifelse(p.adj < 0.001, "***",
                        ifelse(p.adj < 0.01, "**",
                               ifelse(p.adj < 0.05, "*", round(p.adj, 2)))))

a <- max(tcga_kirc_long$value)

png("/Users/oscarbruck/OneDrive - University of Helsinki/RCC/Otso/Analysis/Textures/Images/Margin/TCGA_margin_texture_scatterplot.png", width = 6, height = 6, units = 'in', res = 300, pointsize = 12) #original pointsize = 12
ggplot(tcga_kirc_long, aes(x = Texture, y = value)) +
  geom_jitter(size=5, width = 0.2, aes(fill=Texture), shape = 21, color = "black") +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
  labs(y="Proportion (%)") +
  # scale_y_continuous(limits = c(0,100), expand = c(0, 0)) +
  # scale_fill_brewer(palette="Set1")
  scale_fill_manual(values=c("#4daf4a", "#e41a1c", "#a6cee3", "#984ea3", "#ffff33")) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, colour = "black"),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title=element_text(size=14, face="bold", colour = "black"),
        legend.position = "none") +
  scale_fill_brewer(palette = c("Set1")) +
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
dev.off()