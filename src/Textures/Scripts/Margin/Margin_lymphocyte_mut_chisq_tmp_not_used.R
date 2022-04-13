rm(list=ls())

# Load libraries
library(readxl)
library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(RColorBrewer)



############################# LOAD DATA ##########################################################################################################


# Otso's data
tcga_kirc <- read_xlsx("/Users/oscarbruck/OneDrive - University of Helsinki/RCC/Otso/Data/otso/raw_data.xlsx")

# Normalize by removing empty
tcga_kirc$`texture_cancer_%` <- 100*tcga_kirc$texture_cancer / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`margin_texture_normal_%` <- 100*tcga_kirc$margin_texture_normal / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_cancer + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)
tcga_kirc$`margin_texture_blood_%` <- 100*tcga_kirc$margin_texture_blood / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_cancer + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)
tcga_kirc$`margin_texture_cancer_%` <- 100*tcga_kirc$margin_texture_cancer / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_cancer + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)
tcga_kirc$`margin_texture_normal_%` <- 100*tcga_kirc$margin_texture_normal / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_cancer + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)
tcga_kirc$`margin_texture_stroma_%` <- 100*tcga_kirc$margin_texture_stroma / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_cancer + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)
tcga_kirc$`margin_texture_other_%` <- 100*tcga_kirc$margin_texture_other / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_cancer + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)


tcga_kirc <- tcga_kirc %>%
  dplyr::filter(`texture_cancer_%` > 5) %>%
  dplyr::filter(`margin_texture_normal_%` > 1) %>%
  dplyr::filter(is.na(PoorQuality)) %>%
  dplyr::mutate(tissue_source_site = gsub("-[[:print:]]{4}", "", gsub("TCGA-", "", tcga_id)))


############################# PLOT ##########################################################################################################


# Rename
colnames(tcga_kirc)[grep(pattern = "^inf_margin_bin_lymphocytes", colnames(tcga_kirc))] <- c("Margin_Blood", "Margin_Normal", "Margin_Stroma", "Margin_Other")
colnames(tcga_kirc)[grep(pattern = "^inf_non_margin_bin_lymphocytes", colnames(tcga_kirc))] <- c("NonMargin_Blood", "NonMargin_Cancer", "NonMargin_Normal", "NonMargin_Stroma", "NonMargin_Other")
margin <- c("Margin_Blood", "Margin_Normal", "Margin_Stroma", "Margin_Other")
nonmargin <- c("NonMargin_Blood", "NonMargin_Cancer", "NonMargin_Normal", "NonMargin_Stroma", "NonMargin_Other")


# To percent
tcga_kirc[margin] <- sapply(tcga_kirc[margin], function(x) 100*x)
tcga_kirc[nonmargin] <- sapply(tcga_kirc[nonmargin], function(x) 100*x)

# Factor IDs
# tcga_kirc <- tcga_kirc %>% arrange(Cancer) %>% dplyr::mutate(tcga_id = factor(tcga_id))
# tcga_kirc <- tcga_kirc %>% mutate(tcga_id = factor(Cancer, levels = Cancer))

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


# Adjust p values
pairwise.test = tcga_kirc_long %>% 
  group_by(Texture) %>%
  wilcox_test(value~Margin, paired = TRUE) %>% 
  adjust_pvalue(method = 'BH') %>%
  # mutate(p.adj = round(p.adj, 2))
  mutate(p.adj = ifelse(p.adj < 0.001, "***",
                        ifelse(p.adj < 0.01, "**",
                               ifelse(p.adj < 0.05, "*", round(p.adj, 2)))))


# Plot
png("/Users/oscarbruck/OneDrive - University of Helsinki/RCC/Otso/Analysis/Textures/Images/Margin/TCGA_margin_Ly_scatterplot.png", width = 6, height = 6, units = 'in', res = 300, pointsize = 12) #original pointsize = 12
ggplot(tcga_kirc_long, aes(x = Margin, y = value, group = variable)) +
  geom_point(size=5, aes(fill=Margin), shape = 21, color = "black", position=position_jitterdodge()) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
  labs(y="Proportion (%)", x="Textures") +
  scale_y_continuous(limits = c(0,100), expand = c(0, 0)) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, colour = "black"),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title=element_text(size=14, face="bold", colour = "black"),
        legend.position = "none") +
  scale_fill_brewer(palette = c("Set1")) +
  facet_wrap("Texture") +
  stat_pvalue_manual(
    pairwise.test,
    label = "p.adj",
    bracket.size = 1.5,
    size = 5,
    y.position = 90
  )
dev.off()





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
tcga_kirc_wide
# Join margin and non-margin dataframes
tcga_kirc_wide <- tcga_kirc_margin %>%
  full_join(tcga_kirc_nonmargin) %>%
  mutate(FC = ifelse(Margin>NonMargin, "Margin", ifelse(NonMargin>Margin, "NonMargin", NA)))
  # mutate(FC = ifelse(Margin>NonMargin, "Margin", ifelse(NonMargin>Margin, "NonMargin", NA)))

# Read mutation data
tcga_kirc <- read_xlsx("/Users/oscarbruck/OneDrive - University of Helsinki/RCC/Otso/Data/otso/data_with_ak.xlsx") %>%
  janitor::clean_names() %>%
  dplyr::rename(ID = tcga_id)
## Join
tcga_kirc <- tcga_kirc_wide %>% inner_join(tcga_kirc)


# Function to find the colnumber
findcolnumber <- function(df, thecolumnname){
  which(colnames(df) == deparse(substitute(thecolumnname)))
}

# Select variables
tcga_kirc <- tcga_kirc %>% dplyr::select(ID, Margin, Texture, FC, purity, ploidy, contains("cluster"), contains("mutation"), texture_blood:ncol(tcga_kirc))

# Modify data
tcga_kirc <- tcga_kirc %>% 
  dplyr::mutate_at(.vars = c("purity", "ploidy"), as.numeric)

tcga_kirc <- tcga_kirc %>% 
  dplyr::mutate_at(.vars = c("m_rna_cluster", "mi_rna_cluster", "methylation_cluster"), function(x) gsub(pattern = "Unavailable", "", x)) %>%
  dplyr::mutate(m_rna_cluster = ifelse(m_rna_cluster %in% c(3,4,5), "Other", m_rna_cluster),
                mi_rna_cluster = ifelse(mi_rna_cluster %in% c(2,5,6), "Other", mi_rna_cluster),
                m_rna_cluster = ifelse(m_rna_cluster == "", NA, m_rna_cluster),
                mi_rna_cluster = ifelse(mi_rna_cluster == "", NA, mi_rna_cluster),
                methylation_cluster = ifelse(methylation_cluster == "", NA, methylation_cluster)) %>%
  dplyr::mutate_at(.vars = c("m_rna_cluster", "mi_rna_cluster", "methylation_cluster"), as.factor) %>%
  dplyr::mutate_at(.vars = vars(ends_with("_mutation")), function(x) ifelse(is.na(x), "WT",
                                                                            ifelse(x == "Unavailable", NA, "MUT"))) %>%
  dplyr::mutate(mutations_total = ifelse(mutations_total == "Unavailable", NA, as.numeric(as.character(mutations_total))),
                mutations_total = ifelse(mutations_total > median(mutations_total, na.rm=TRUE), "High", "Low")) %>%
  dplyr::mutate(ploidy = ifelse(ploidy == "Unavailable", NA, as.numeric(as.character(ploidy))),
                ploidy = ifelse(ploidy > median(ploidy, na.rm=TRUE), "High", "Low"))


## Multiple comparison (two-sided, unpaired, Mann-Whitney U test)
multiple_t_tests_p_value <- lapply(tcga_kirc[c(str_detect(colnames(tcga_kirc), "_mutation") | colnames(tcga_kirc)=="mutations_total" | colnames(tcga_kirc)=="ploidy")],
                                   function(x) chisq.test(tcga_kirc$FC, x))
### P-values can be extracted from the result object
pvalue <- data.frame(p.value = sapply(multiple_t_tests_p_value, getElement, name = "p.value"))
### Observed counts
observed_counts = lapply(multiple_t_tests_p_value, getElement, name = "observed")
### Create a matrix and dataframe of the p-values
pvalue_df2 <- pvalue %>% data.frame() %>% mutate(
  #### Add the p values to a new dataframe
  p_adjusted = p.adjust(p = as.matrix(pvalue), method = "BH")
)
pvalue_df2$variable <- names(tcga_kirc[c(str_detect(colnames(tcga_kirc), "_mutation") | colnames(tcga_kirc)=="mutations_total" | colnames(tcga_kirc)=="ploidy")])
pvalue_df2 <- pvalue_df2 %>%
  arrange(p.value)

tmp1 <- NULL
for (var1 in 1:length(pvalue_df2$variable)) {
  i = pvalue_df2$variable[var1]
  print(i)
  tmp <- data.frame(cbind(i, observed_counts[[i]][1], observed_counts[[i]][2], observed_counts[[i]][3], observed_counts[[i]][4]))
  colnames(tmp) <- c("variable", "MUT_Margin", "WT_Margin", "MUT_NonMargin", "WT_NonMargin")
  tmp <- as.data.frame(tmp)
  tmp$variable <- as.character(tmp$variable)
  
  # Collect
  if (exists("tmp1")) {
    tmp1 <- tmp1 %>% rbind(tmp)
  } else {
    tmp1 <- tmp
  }
}

# Join
pvalue_df2 <- pvalue_df2 %>% dplyr::left_join(tmp1)

