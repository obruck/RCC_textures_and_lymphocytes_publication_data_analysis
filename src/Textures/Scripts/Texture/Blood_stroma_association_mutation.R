rm(list = ls())

print("Start Textures/Scripts/Texture/Blood_stroma_association_mutation.R")

# Load libraries
library(readxl)
library(tidyverse)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(rstatix)
library(GSA)
library(ggpubr)


############################# LOAD DATA ##########################################################################################################


# Image data
tcga_kirc <- read_xlsx("../data/mutations_final.xlsx")


# Normalize by removing empty
tcga_kirc$`texture_blood_%` <- 100*tcga_kirc$texture_blood / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`texture_cancer_%` <- 100*tcga_kirc$texture_cancer / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`texture_normal_%` <- 100*tcga_kirc$texture_normal / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`texture_stroma_%` <- 100*tcga_kirc$texture_stroma / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`texture_other_%` <- 100*tcga_kirc$texture_other / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)

tcga_kirc <- tcga_kirc %>%
  dplyr::filter(`texture_cancer_%` > 5) %>%
  janitor::clean_names() %>%
  dplyr::mutate(tissue_source_site = gsub("-[[:print:]]{4}", "", gsub("TCGA-", "", tcga_id)))


############################# PLOT ##########################################################################################################


# Rename
colnames(tcga_kirc)[grep(pattern = "^texture_[[:print:]]*percent", colnames(tcga_kirc))] <- c("Blood", "Cancer", "Normal", "Stroma", "Other")
textures <- c("Blood", "Cancer", "Normal", "Stroma", "Other")


################################
##### BLOOD WITHOUT NORMAL #####
################################


# Function to find the colnumber
findcolnumber <- function(df, thecolumnname){
  which(colnames(df) == deparse(substitute(thecolumnname)))
}

# Filter only samples with <1% normal tissue
tcga_kirc <- tcga_kirc %>% dplyr::filter(Normal < 1)


# Select variables
tcga_kirc <- tcga_kirc %>% dplyr::select(tcga_id, purity, ploidy, contains("cluster"), contains("mutation"), texture_blood:ncol(tcga_kirc))

# Modify data
tcga_kirc <- tcga_kirc %>% 
  dplyr::mutate_at(.vars = c("purity", "ploidy"), as.numeric)

table(tcga_kirc$m_rna_cluster)
tcga_kirc <- tcga_kirc %>% 
  dplyr::mutate_at(.vars = c("m_rna_cluster", "mi_rna_cluster", "methylation_cluster"), function(x) gsub(pattern = "Unavailable", "", x)) %>%
  dplyr::mutate(mi_rna_cluster = ifelse(mi_rna_cluster %in% c(2,5,6), "Other", mi_rna_cluster),
                m_rna_cluster = ifelse(m_rna_cluster %in% c(3,4,5), "Other", m_rna_cluster),
                m_rna_cluster = ifelse(m_rna_cluster %in% c("", "n/a"), NA, m_rna_cluster),
                mi_rna_cluster = ifelse(mi_rna_cluster %in% c("", "n/a"), NA, mi_rna_cluster),
                methylation_cluster = ifelse(methylation_cluster %in% c("", "n/a"), NA, methylation_cluster)) %>%
  dplyr::mutate_at(.vars = c("m_rna_cluster", "mi_rna_cluster", "methylation_cluster"), as.factor) %>%
  dplyr::mutate_at(.vars = vars(ends_with("_mutation")), function(x) ifelse(is.na(x), "WT",
                                                                      ifelse(x == "Unavailable", NA, "MUT"))) %>%
  dplyr::mutate(mutations_total = ifelse(mutations_total == "Unavailable", NA, as.numeric(as.character(mutations_total))),
                mutations_total = ifelse(mutations_total > median(mutations_total, na.rm=TRUE), "High", "Low")) %>%
  dplyr::mutate(ploidy = ifelse(ploidy == "Unavailable", NA, as.numeric(as.character(ploidy))),
                ploidy = ifelse(ploidy > median(ploidy, na.rm=TRUE), "High", "Low"))


# Find variables with 2 values excluding NA
unique_values = sapply(tcga_kirc[c(str_detect(colnames(tcga_kirc), "_mutation") | colnames(tcga_kirc)=="mutations_total" | colnames(tcga_kirc)=="ploidy")], function(x) length(unique(x[!is.na(x)])))
unique_values = names(unique_values[unique_values>1])

# Multiple comparison (two-sided, unpaired, Mann-Whitney U test)
# gexp_data <- tcga_kirc[(findcolnumber(tcga_kirc, tissue_source_site)+1):ncol(tcga_kirc)]
multiple_t_tests_p_value <- lapply(tcga_kirc[unique_values],
                                   function(x) wilcox.test(tcga_kirc$Blood ~ x, correct=FALSE, exact=FALSE, estimate=TRUE, na.rm=TRUE))
## P-values can be extracted from the result object
pvalue <- data.frame(p.value = sapply(multiple_t_tests_p_value, getElement, name = "p.value"))
## Create a matrix and dataframe of the p-values
pvalue_mat <- as.matrix(pvalue)
pvalue_df <- data.frame(pvalue)
## Add the p values to a new dataframe
pvalue_df$p_adjusted <- p.adjust(p = pvalue_mat, method = "BH")
## Add also the t values, 95%CI to the same dataframe
pvalue_df$median1 <- unlist(lapply(tcga_kirc[unique_values], function(x) median(tcga_kirc[x %in% c("High", "MUT"),]$Blood, na.rm=TRUE)))
pvalue_df$median2 <- unlist(lapply(tcga_kirc[unique_values], function(x) median(tcga_kirc[x %in% c("Low", "WT"),]$Blood, na.rm=TRUE)))
## Rownames to column
pvalue_df <- pvalue_df %>%
  rownames_to_column() %>%
  rename(genes = rowname,
         pvalue = p.value
  ) %>%
  arrange(pvalue)
rownames(pvalue_df) = NULL


# Rename variables
# pvalue_df <- pvalue_df %>%
#   mutate(genes = gsub("_mutation", "", genes),
#          genes = toupper(genes)
#   )


# Export data
writexl::write_xlsx(pvalue_df, "Textures/Images/Blood/Mut/Blood_mut.xlsx")
tcga_kirc0 <- tcga_kirc


# Plot parameters
a <- max(tcga_kirc$Blood, na.rm=TRUE)


for (two1 in c(pvalue_df$genes, "methylation_cluster")) {
  
  # Assign value
  tcga_kirc0$two1 <- tcga_kirc0[[two1]]
  
  # Filter NAs
  tcga_kirc1 <- tcga_kirc0 %>%
    dplyr::filter(!is.na(two1)) %>%
    dplyr::mutate(two1 = as.factor(two1))
  
  # Adjust p values
  pairwise.test = tcga_kirc1 %>%
    wilcox_test(Blood~two1) %>%
    adjust_pvalue(method = 'BH') %>%
    mutate(p.adj = ifelse(p.adj < 0.001, "***",
                      ifelse(p.adj < 0.01, "**",
                             ifelse(p.adj < 0.05, "*", round(p.adj, 2))))) %>%
    mutate(p = ifelse(p < 0.001, "***",
                      ifelse(p < 0.01, "**",
                             ifelse(p < 0.05, "*", round(p, 2)))))
  
  
  # Plot
  g <- ggplot(tcga_kirc1, aes(x = two1, y = Blood)) +
    geom_jitter(size=5, width = 0.2, aes(fill=two1), shape = 21, color = "black") +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
    labs(y="Blood proportion (%)", x=toupper(gsub("_mutation", "", two1))) +
    theme_bw() +
    theme(axis.text.x = element_text(size=12, colour = "black"),
          axis.text.y = element_text(size=12, colour = "black"),
          axis.title=element_text(size=14, face="bold", colour = "black"),
          legend.position = "none") +
    scale_fill_brewer(palette = c("Set1")) +
    # scale_x_discrete(labels=c("0" = "Low", "1" = "High")) +
    # stat_compare_means(method = "wilcoxon.test",
    #                    label.y = 1.0*a,
    #                    label = "p.signif",
    #                    bracket.size = 1,
    #                    size = 5) +
    stat_pvalue_manual(
      pairwise.test,
      label = "p",
      bracket.size = 1,
      size = 5,
      # y.position = c(1.12*a, 1.24*a, 1.3*a, 1.0*a, 1.18*a, 1.06*a)
      y.position = 1.0*a
    )
  ggsave(plot = g, filename = paste0("Textures/Images/Blood/Mut/Blood_", two1, "_mut.png"), width = 5, height = 5, units = 'in', dpi = 300)
  
}


for (two1 in c("m_rna_cluster")) {
  
  # Assign value
  tcga_kirc0$two1 <- tcga_kirc0[[two1]]
  
  # Filter NAs
  tcga_kirc1 <- tcga_kirc0 %>%
    dplyr::filter(!is.na(two1)) %>%
    dplyr::mutate(two1 = as.factor(two1))
  
  # Adjust p values
  pairwise.test = tcga_kirc1 %>%
    wilcox_test(Blood~two1) %>%
    adjust_pvalue(method = 'BH') %>%
    mutate(p.adj = ifelse(p.adj < 0.001, "***",
                          ifelse(p.adj < 0.01, "**",
                                 ifelse(p.adj < 0.05, "*", round(p.adj, 2))))) %>%
    mutate(p = ifelse(p < 0.001, "***",
                      ifelse(p < 0.01, "**",
                             ifelse(p < 0.05, "*", round(p, 2)))))
  
  
  # Plot
  g <- ggplot(tcga_kirc1, aes(x = two1, y = Blood)) +
    geom_jitter(size=5, width = 0.2, aes(fill=two1), shape = 21, color = "black") +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
    labs(y="Blood proportion (%)", x=gsub("_", " ", gsub("methylation", "Methylation", gsub("m_rna", "mRNA", two1)))) +
    theme_bw() +
    theme(axis.text.x = element_text(size=12, colour = "black"),
          axis.text.y = element_text(size=12, colour = "black"),
          axis.title=element_text(size=14, face="bold", colour = "black"),
          legend.position = "none") +
    scale_fill_brewer(palette = c("Set1")) +
    stat_compare_means(method = "kruskal.test",
                       label.y = 1.25*a,
                       size = 6) +
    stat_pvalue_manual(
      pairwise.test,
      label = "p.adj",
      bracket.size = 1,
      size = 5,
      y.position = c(1.0*a, 1.12*a, 1.06*a)
    )
  ggsave(plot = g, filename = paste0("Textures/Images/Blood/Mut/Blood_", two1, "_cluster.png"), width = 5, height = 5, units = 'in', dpi = 300)
  
}


for (two1 in c("mi_rna_cluster")) {
  
  # Assign value
  tcga_kirc0$two1 <- tcga_kirc0[[two1]]
  
  # Filter NAs
  tcga_kirc1 <- tcga_kirc0 %>%
    dplyr::filter(!is.na(two1)) %>%
    dplyr::filter(!(two1 == "n/a")) %>%
    dplyr::mutate(two1 = factor(two1))
  
  # Adjust p values
  pairwise.test = tcga_kirc1 %>%
    wilcox_test(Blood~two1) %>%
    adjust_pvalue(method = 'BH') %>%
    mutate(p.adj = ifelse(p.adj < 0.001, "***",
                          ifelse(p.adj < 0.01, "**",
                                 ifelse(p.adj < 0.05, "*", round(p.adj, 2))))) %>%
    mutate(p = ifelse(p < 0.001, "***",
                      ifelse(p < 0.01, "**",
                             ifelse(p < 0.05, "*", round(p, 2)))))
  
  
  # Plot
  g <- ggplot(tcga_kirc1, aes(x = two1, y = Blood)) +
    geom_jitter(size=5, width = 0.2, aes(fill=two1), shape = 21, color = "black") +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
    labs(y="Blood proportion (%)", x=gsub("_", " ", gsub("mi_rna", "miRNA", two1))) +
    theme_bw() +
    theme(axis.text.x = element_text(size=12, colour = "black"),
          axis.text.y = element_text(size=12, colour = "black"),
          axis.title=element_text(size=14, face="bold", colour = "black"),
          legend.position = "none") +
    scale_fill_brewer(palette = c("Set1")) +
    stat_compare_means(method = "kruskal.test",
                       label.y = 1.45*a,
                       size = 6) +
    stat_pvalue_manual(
      pairwise.test,
      label = "p.adj",
      bracket.size = 1,
      size = 5,
      # y.position = c(1.12*a, 1.24*a, 1.3*a, 1.0*a, 1.18*a, 1.06*a)
      y.position = c(1.14*a, 1.28*a, 1.35*a, 1.0*a, 1.21*a, 1.07*a)
    )
  ggsave(plot = g, filename = paste0("Textures/Images/Blood/Mut/Blood_", two1, "_cluster.png"), width = 5, height = 5, units = 'in', dpi = 300)
  
}




###########################
##### BLOOD WO NORMAL #####
###########################


############################# LOAD DATA ##########################################################################################################


# Image data
tcga_kirc <- read_xlsx("../data/mutations_final.xlsx")


# Normalize by removing empty
tcga_kirc$`texture_blood_%` <- 100*tcga_kirc$texture_blood / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`texture_cancer_%` <- 100*tcga_kirc$texture_cancer / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`texture_normal_%` <- 100*tcga_kirc$texture_normal / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`texture_stroma_%` <- 100*tcga_kirc$texture_stroma / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`texture_other_%` <- 100*tcga_kirc$texture_other / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)

tcga_kirc <- tcga_kirc %>%
  dplyr::filter(`texture_cancer_%` > 5) %>%
  janitor::clean_names() %>%
  dplyr::mutate(tissue_source_site = gsub("-[[:print:]]{4}", "", gsub("TCGA-", "", tcga_id)))


############################# PLOT ##########################################################################################################


# Rename
colnames(tcga_kirc)[grep(pattern = "^texture_[[:print:]]*percent", colnames(tcga_kirc))] <- c("Blood", "Cancer", "Normal", "Stroma", "Other")
textures <- c("Blood", "Cancer", "Normal", "Stroma", "Other")


# Function to find the colnumber
findcolnumber <- function(df, thecolumnname){
  which(colnames(df) == deparse(substitute(thecolumnname)))
}

# Filter only samples with <1% normal tissue
tcga_kirc <- tcga_kirc %>% dplyr::filter(Normal >= 1)


# Select variables
tcga_kirc <- tcga_kirc %>% dplyr::select(tcga_id, purity, ploidy, contains("cluster"), contains("mutation"), texture_blood:ncol(tcga_kirc))

# Modify data
tcga_kirc <- tcga_kirc %>% 
  dplyr::mutate_at(.vars = c("purity", "ploidy"), as.numeric)


tcga_kirc <- tcga_kirc %>% 
  dplyr::mutate_at(.vars = c("m_rna_cluster", "mi_rna_cluster", "methylation_cluster"), function(x) gsub(pattern = "Unavailable", "", x)) %>%
  dplyr::mutate(mi_rna_cluster = ifelse(mi_rna_cluster %in% c(2,5,6), "Other", mi_rna_cluster),
                m_rna_cluster = ifelse(m_rna_cluster %in% c(3,4,5), "Other", m_rna_cluster),
                m_rna_cluster = ifelse(m_rna_cluster %in% c("", "n/a"), NA, m_rna_cluster),
                mi_rna_cluster = ifelse(mi_rna_cluster %in% c("", "n/a"), NA, mi_rna_cluster),
                methylation_cluster = ifelse(methylation_cluster %in% c("", "n/a"), NA, methylation_cluster)) %>%
  dplyr::mutate_at(.vars = c("m_rna_cluster", "mi_rna_cluster", "methylation_cluster"), as.factor) %>%
  dplyr::mutate_at(.vars = vars(ends_with("_mutation")), function(x) ifelse(is.na(x), "WT",
                                                                            ifelse(x == "Unavailable", NA, "MUT"))) %>%
  dplyr::mutate(mutations_total = ifelse(mutations_total == "Unavailable", NA, as.numeric(as.character(mutations_total))),
                mutations_total = ifelse(mutations_total > median(mutations_total, na.rm=TRUE), "High", "Low")) %>%
  dplyr::mutate(ploidy = ifelse(ploidy == "Unavailable", NA, as.numeric(as.character(ploidy))),
                ploidy = ifelse(ploidy > median(ploidy, na.rm=TRUE), "High", "Low"))

# a <- data.frame(tcga_kirc[c(str_detect(colnames(tcga_kirc), "_mutation") | colnames(tcga_kirc)=="mutations_total" | colnames(tcga_kirc)=="ploidy")])
# b <- sapply(a, table)

# Find variables with 2 values excluding NA
unique_values_2 = sapply(tcga_kirc[c((str_detect(colnames(tcga_kirc), "_mutation") & !colnames(tcga_kirc)=="met_mutation") | colnames(tcga_kirc)=="mutations_total" | colnames(tcga_kirc)=="ploidy")], function(x) length(unique(x[!is.na(x)])))
unique_values_2 = names(unique_values_2[unique_values_2==2])

# Multiple comparison (two-sided, unpaired, Mann-Whitney U test)
# gexp_data <- tcga_kirc[(findcolnumber(tcga_kirc, tissue_source_site)+1):ncol(tcga_kirc)]
multiple_t_tests_p_value <- lapply(tcga_kirc[unique_values_2],
                                   function(x) wilcox.test(tcga_kirc$Blood ~ x, correct=FALSE, exact=FALSE, estimate=TRUE, na.rm=TRUE))
## P-values can be extracted from the result object
pvalue <- data.frame(p.value = sapply(multiple_t_tests_p_value, getElement, name = "p.value"))
## Create a matrix and dataframe of the p-values
pvalue_mat <- as.matrix(pvalue)
pvalue_df <- data.frame(pvalue)
## Add the p values to a new dataframe
pvalue_df$p_adjusted <- p.adjust(p = pvalue_mat, method = "BH")
## Add also the t values, 95%CI to the same dataframe
pvalue_df$median1 <- unlist(lapply(tcga_kirc[unique_values_2], function(x) median(tcga_kirc[x %in% c("High", "MUT"),]$Blood, na.rm=TRUE)))
pvalue_df$median2 <- unlist(lapply(tcga_kirc[unique_values_2], function(x) median(tcga_kirc[x %in% c("Low", "WT"),]$Blood, na.rm=TRUE)))
## Rownames to column
pvalue_df <- pvalue_df %>%
  rownames_to_column() %>%
  rename(genes = rowname,
         pvalue = p.value
  ) %>%
  arrange(pvalue)
rownames(pvalue_df) = NULL


# Rename variables
# pvalue_df <- pvalue_df %>%
#   mutate(genes = gsub("_mutation", "", genes),
#          genes = toupper(genes)
#   )


# Export data
writexl::write_xlsx(pvalue_df, "Textures/Images/Blood_WithNormal/Mut/Blood_mut.xlsx")
tcga_kirc0 <- tcga_kirc


# Plot parameters
a <- max(tcga_kirc$Blood, na.rm=TRUE)


for (two1 in c(pvalue_df$genes)) {
  
  # Assign value
  tcga_kirc0$two1 <- tcga_kirc0[[two1]]
  
  # Filter NAs
  tcga_kirc1 <- tcga_kirc0 %>%
    dplyr::filter(!is.na(two1)) %>%
    dplyr::mutate(two1 = as.factor(two1))
  
  # Adjust p values
  pairwise.test = tcga_kirc1 %>%
    wilcox_test(Blood~two1) %>%
    adjust_pvalue(method = 'BH') %>%
    mutate(p.adj = ifelse(p.adj < 0.001, "***",
                          ifelse(p.adj < 0.01, "**",
                                 ifelse(p.adj < 0.05, "*", round(p.adj, 2))))) %>%
    mutate(p = ifelse(p < 0.001, "***",
                      ifelse(p < 0.01, "**",
                             ifelse(p < 0.05, "*", round(p, 2)))))
  
  
  # Plot
  g <- ggplot(tcga_kirc1, aes(x = two1, y = Blood)) +
    geom_jitter(size=5, width = 0.2, aes(fill=two1), shape = 21, color = "black") +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
    labs(y="Blood proportion (%)", x=toupper(gsub("_mutation", "", two1))) +
    theme_bw() +
    theme(axis.text.x = element_text(size=12, colour = "black"),
          axis.text.y = element_text(size=12, colour = "black"),
          axis.title=element_text(size=14, face="bold", colour = "black"),
          legend.position = "none") +
    scale_fill_brewer(palette = c("Set1")) +
    # scale_x_discrete(labels=c("0" = "Low", "1" = "High")) +
    # stat_compare_means(method = "wilcoxon.test",
    #                    label.y = 1.0*a,
    #                    label = "p.signif",
    #                    bracket.size = 1,
    #                    size = 5) +
    stat_pvalue_manual(
      pairwise.test,
      label = "p",
      bracket.size = 1,
      size = 5,
      # y.position = c(1.12*a, 1.24*a, 1.3*a, 1.0*a, 1.18*a, 1.06*a)
      y.position = 1.0*a
    )
  ggsave(plot = g, filename = paste0("Textures/Images/Blood_WithNormal/Mut/Blood_", two1, "_mut.png"), width = 5, height = 5, units = 'in', dpi = 300)
  
}


for (two1 in c("m_rna_cluster", "methylation_cluster")) {
  
  # Assign value
  tcga_kirc0$two1 <- tcga_kirc0[[two1]]
  
  # Filter NAs
  tcga_kirc1 <- tcga_kirc0 %>%
    dplyr::filter(!is.na(two1)) %>%
    dplyr::mutate(two1 = as.factor(two1))
  
  # Adjust p values
  pairwise.test = tcga_kirc1 %>%
    wilcox_test(Blood~two1) %>%
    adjust_pvalue(method = 'BH') %>%
    mutate(p.adj = ifelse(p.adj < 0.001, "***",
                          ifelse(p.adj < 0.01, "**",
                                 ifelse(p.adj < 0.05, "*", round(p.adj, 2))))) %>%
    mutate(p = ifelse(p < 0.001, "***",
                      ifelse(p < 0.01, "**",
                             ifelse(p < 0.05, "*", round(p, 2)))))
  
  
  # Plot
  g <- ggplot(tcga_kirc1, aes(x = two1, y = Blood)) +
    geom_jitter(size=5, width = 0.2, aes(fill=two1), shape = 21, color = "black") +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
    labs(y="Blood proportion (%)", x=gsub("_", " ", gsub("methylation", "Methylation", gsub("m_rna", "mRNA", two1)))) +
    theme_bw() +
    theme(axis.text.x = element_text(size=12, colour = "black"),
          axis.text.y = element_text(size=12, colour = "black"),
          axis.title=element_text(size=14, face="bold", colour = "black"),
          legend.position = "none") +
    scale_fill_brewer(palette = c("Set1")) +
    stat_compare_means(method = "kruskal.test",
                       label.y = 1.25*a,
                       size = 6) +
    stat_pvalue_manual(
      pairwise.test,
      label = "p.adj",
      bracket.size = 1,
      size = 5,
      y.position = c(1.0*a, 1.12*a, 1.06*a)
    )
  ggsave(plot = g, filename = paste0("Textures/Images/Blood_WithNormal/Mut/Blood_", two1, "_cluster.png"), width = 5, height = 5, units = 'in', dpi = 300)
  
}


for (two1 in c("mi_rna_cluster")) {
  
  # Assign value
  tcga_kirc0$two1 <- tcga_kirc0[[two1]]
  
  # Filter NAs
  tcga_kirc1 <- tcga_kirc0 %>%
    dplyr::filter(!is.na(two1)) %>%
    dplyr::filter(!(two1 == "n/a")) %>%
    dplyr::mutate(two1 = factor(two1))
  
  # Adjust p values
  pairwise.test = tcga_kirc1 %>%
    wilcox_test(Blood~two1) %>%
    adjust_pvalue(method = 'BH') %>%
    mutate(p.adj = ifelse(p.adj < 0.001, "***",
                          ifelse(p.adj < 0.01, "**",
                                 ifelse(p.adj < 0.05, "*", round(p.adj, 2))))) %>%
    mutate(p = ifelse(p < 0.001, "***",
                      ifelse(p < 0.01, "**",
                             ifelse(p < 0.05, "*", round(p, 2)))))
  
  
  # Plot
  g <- ggplot(tcga_kirc1, aes(x = two1, y = Blood)) +
    geom_jitter(size=5, width = 0.2, aes(fill=two1), shape = 21, color = "black") +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
    labs(y="Blood proportion (%)", x=gsub("_", " ", gsub("mi_rna", "miRNA", two1))) +
    theme_bw() +
    theme(axis.text.x = element_text(size=12, colour = "black"),
          axis.text.y = element_text(size=12, colour = "black"),
          axis.title=element_text(size=14, face="bold", colour = "black"),
          legend.position = "none") +
    scale_fill_brewer(palette = c("Set1")) +
    stat_compare_means(method = "kruskal.test",
                       label.y = 1.45*a,
                       size = 6) +
    stat_pvalue_manual(
      pairwise.test,
      label = "p.adj",
      bracket.size = 1,
      size = 5,
      # y.position = c(1.12*a, 1.24*a, 1.3*a, 1.0*a, 1.18*a, 1.06*a)
      y.position = c(1.14*a, 1.28*a, 1.35*a, 1.0*a, 1.21*a, 1.07*a)
    )
  ggsave(plot = g, filename = paste0("Textures/Images/Blood_WithNormal/Mut/Blood_", two1, "_cluster.png"), width = 5, height = 5, units = 'in', dpi = 300)
  
}




###########################€
##### STROMA WO NORMAL #####
############################


############################# LOAD DATA ##########################################################################################################


# Image data
tcga_kirc <- read_xlsx("../data/mutations_final.xlsx")


# Normalize by removing empty
tcga_kirc$`texture_blood_%` <- 100*tcga_kirc$texture_blood / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`texture_cancer_%` <- 100*tcga_kirc$texture_cancer / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`texture_normal_%` <- 100*tcga_kirc$texture_normal / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`texture_stroma_%` <- 100*tcga_kirc$texture_stroma / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`texture_other_%` <- 100*tcga_kirc$texture_other / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)

tcga_kirc <- tcga_kirc %>%
  dplyr::filter(`texture_cancer_%` > 5) %>%
  janitor::clean_names() %>%
  dplyr::mutate(tissue_source_site = gsub("-[[:print:]]{4}", "", gsub("TCGA-", "", tcga_id)))


############################# PLOT ##########################################################################################################


# Rename
colnames(tcga_kirc)[grep(pattern = "^texture_[[:print:]]*percent", colnames(tcga_kirc))] <- c("Blood", "Cancer", "Normal", "Stroma", "Other")
textures <- c("Blood", "Cancer", "Normal", "Stroma", "Other")


# Function to find the colnumber
findcolnumber <- function(df, thecolumnname){
  which(colnames(df) == deparse(substitute(thecolumnname)))
}

# Filter only samples with <1% normal tissue
tcga_kirc <- tcga_kirc %>% dplyr::filter(Normal < 1)


# Select variables
tcga_kirc <- tcga_kirc %>% dplyr::select(tcga_id, purity, ploidy, contains("cluster"), contains("mutation"), texture_blood:ncol(tcga_kirc))

# Modify data
tcga_kirc <- tcga_kirc %>% 
  dplyr::mutate_at(.vars = c("purity", "ploidy"), as.numeric)


tcga_kirc <- tcga_kirc %>% 
  dplyr::mutate_at(.vars = c("m_rna_cluster", "mi_rna_cluster", "methylation_cluster"), function(x) gsub(pattern = "Unavailable", "", x)) %>%
  dplyr::mutate(mi_rna_cluster = ifelse(mi_rna_cluster %in% c(2,5,6), "Other", mi_rna_cluster),
                m_rna_cluster = ifelse(m_rna_cluster %in% c(3,4,5), "Other", m_rna_cluster),
                m_rna_cluster = ifelse(m_rna_cluster %in% c("", "n/a"), NA, m_rna_cluster),
                mi_rna_cluster = ifelse(mi_rna_cluster %in% c("", "n/a"), NA, mi_rna_cluster),
                methylation_cluster = ifelse(methylation_cluster %in% c("", "n/a"), NA, methylation_cluster)) %>%
  dplyr::mutate_at(.vars = c("m_rna_cluster", "mi_rna_cluster", "methylation_cluster"), as.factor) %>%
  dplyr::mutate_at(.vars = vars(ends_with("_mutation")), function(x) ifelse(is.na(x), "WT",
                                                                            ifelse(x == "Unavailable", NA, "MUT"))) %>%
  dplyr::mutate(mutations_total = ifelse(mutations_total == "Unavailable", NA, as.numeric(as.character(mutations_total))),
                mutations_total = ifelse(mutations_total > median(mutations_total, na.rm=TRUE), "High", "Low")) %>%
  dplyr::mutate(ploidy = ifelse(ploidy == "Unavailable", NA, as.numeric(as.character(ploidy))),
                ploidy = ifelse(ploidy > median(ploidy, na.rm=TRUE), "High", "Low"))

# a <- data.frame(tcga_kirc[c(str_detect(colnames(tcga_kirc), "_mutation") | colnames(tcga_kirc)=="mutations_total" | colnames(tcga_kirc)=="ploidy")])
# b <- sapply(a, table)

# Find variables with 2 values excluding NA
unique_values_2 = sapply(tcga_kirc[c((str_detect(colnames(tcga_kirc), "_mutation") & !colnames(tcga_kirc)=="met_mutation") | colnames(tcga_kirc)=="mutations_total" | colnames(tcga_kirc)=="ploidy")], function(x) length(unique(x[!is.na(x)])))
unique_values_2 = names(unique_values_2[unique_values_2==2])

# Multiple comparison (two-sided, unpaired, Mann-Whitney U test)
# gexp_data <- tcga_kirc[(findcolnumber(tcga_kirc, tissue_source_site)+1):ncol(tcga_kirc)]
multiple_t_tests_p_value <- lapply(tcga_kirc[unique_values_2],
                                   function(x) wilcox.test(tcga_kirc$Stroma ~ x, correct=FALSE, exact=FALSE, estimate=TRUE, na.rm=TRUE))
## P-values can be extracted from the result object
pvalue <- data.frame(p.value = sapply(multiple_t_tests_p_value, getElement, name = "p.value"))
## Create a matrix and dataframe of the p-values
pvalue_mat <- as.matrix(pvalue)
pvalue_df <- data.frame(pvalue)
## Add the p values to a new dataframe
pvalue_df$p_adjusted <- p.adjust(p = pvalue_mat, method = "BH")
## Add also the t values, 95%CI to the same dataframe
pvalue_df$median1 <- unlist(lapply(tcga_kirc[unique_values_2], function(x) median(tcga_kirc[x %in% c("High", "MUT"),]$Stroma, na.rm=TRUE)))
pvalue_df$median2 <- unlist(lapply(tcga_kirc[unique_values_2], function(x) median(tcga_kirc[x %in% c("Low", "WT"),]$Stroma, na.rm=TRUE)))
## Rownames to column
pvalue_df <- pvalue_df %>%
  rownames_to_column() %>%
  rename(genes = rowname,
         pvalue = p.value
  ) %>%
  arrange(pvalue)
rownames(pvalue_df) = NULL


# Rename variables
# pvalue_df <- pvalue_df %>%
#   mutate(genes = gsub("_mutation", "", genes),
#          genes = toupper(genes)
#   )


# Export data
writexl::write_xlsx(pvalue_df, "Textures/Images/Stroma/Mut/Stroma_mut.xlsx")
tcga_kirc0 <- tcga_kirc


# Plot parameters
a <- max(tcga_kirc$Stroma, na.rm=TRUE)


for (two1 in c(pvalue_df$genes, "methylation_cluster")) {
  
  # Assign value
  tcga_kirc0$two1 <- tcga_kirc0[[two1]]
  
  # Filter NAs
  tcga_kirc1 <- tcga_kirc0 %>%
    dplyr::filter(!is.na(two1)) %>%
    dplyr::mutate(two1 = as.factor(two1))
  
  # Adjust p values
  pairwise.test = tcga_kirc1 %>%
    wilcox_test(Stroma~two1) %>%
    adjust_pvalue(method = 'BH') %>%
    mutate(p.adj = ifelse(p.adj < 0.001, "***",
                          ifelse(p.adj < 0.01, "**",
                                 ifelse(p.adj < 0.05, "*", round(p.adj, 2))))) %>%
    mutate(p = ifelse(p < 0.001, "***",
                      ifelse(p < 0.01, "**",
                             ifelse(p < 0.05, "*", round(p, 2)))))
  
  
  # Plot
  g <- ggplot(tcga_kirc1, aes(x = two1, y = Stroma)) +
    geom_jitter(size=5, width = 0.2, aes(fill=two1), shape = 21, color = "black") +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
    labs(y="Stroma proportion (%)", x=toupper(gsub("_mutation", "", two1))) +
    theme_bw() +
    theme(axis.text.x = element_text(size=12, colour = "black"),
          axis.text.y = element_text(size=12, colour = "black"),
          axis.title=element_text(size=14, face="bold", colour = "black"),
          legend.position = "none") +
    scale_fill_brewer(palette = c("Set1")) +
    # scale_x_discrete(labels=c("0" = "Low", "1" = "High")) +
    # stat_compare_means(method = "wilcoxon.test",
    #                    label.y = 1.0*a,
    #                    label = "p.signif",
    #                    bracket.size = 1,
    #                    size = 5) +
    stat_pvalue_manual(
      pairwise.test,
      label = "p",
      bracket.size = 1,
      size = 5,
      # y.position = c(1.12*a, 1.24*a, 1.3*a, 1.0*a, 1.18*a, 1.06*a)
      y.position = 1.0*a
    )
  ggsave(plot = g, filename = paste0("Textures/Images/Stroma/Mut/Stroma_", two1, "_mut.png"), width = 5, height = 5, units = 'in', dpi = 300)
  
}


for (two1 in c("m_rna_cluster")) {
  
  # Assign value
  tcga_kirc0$two1 <- tcga_kirc0[[two1]]
  
  # Filter NAs
  tcga_kirc1 <- tcga_kirc0 %>%
    dplyr::filter(!is.na(two1)) %>%
    dplyr::mutate(two1 = as.factor(two1))
  
  # Adjust p values
  pairwise.test = tcga_kirc1 %>%
    wilcox_test(Stroma~two1) %>%
    adjust_pvalue(method = 'BH') %>%
    mutate(p.adj = ifelse(p.adj < 0.001, "***",
                          ifelse(p.adj < 0.01, "**",
                                 ifelse(p.adj < 0.05, "*", round(p.adj, 2))))) %>%
    mutate(p = ifelse(p < 0.001, "***",
                      ifelse(p < 0.01, "**",
                             ifelse(p < 0.05, "*", round(p, 2)))))
  
  
  # Plot
  g <- ggplot(tcga_kirc1, aes(x = two1, y = Stroma)) +
    geom_jitter(size=5, width = 0.2, aes(fill=two1), shape = 21, color = "black") +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
    labs(y="Stroma proportion (%)", x=gsub("_", " ", gsub("methylation", "Methylation", gsub("m_rna", "mRNA", two1)))) +
    theme_bw() +
    theme(axis.text.x = element_text(size=12, colour = "black"),
          axis.text.y = element_text(size=12, colour = "black"),
          axis.title=element_text(size=14, face="bold", colour = "black"),
          legend.position = "none") +
    scale_fill_brewer(palette = c("Set1")) +
    stat_compare_means(method = "kruskal.test",
                       label.y = 1.25*a,
                       size = 6) +
    stat_pvalue_manual(
      pairwise.test,
      label = "p.adj",
      bracket.size = 1,
      size = 5,
      y.position = c(1.0*a, 1.12*a, 1.06*a)
    )
  ggsave(plot = g, filename = paste0("Textures/Images/Stroma/Mut/Stroma_", two1, "_cluster.png"), width = 5, height = 5, units = 'in', dpi = 300)
  
}


for (two1 in c("mi_rna_cluster")) {
  
  # Assign value
  tcga_kirc0$two1 <- tcga_kirc0[[two1]]
  
  # Filter NAs
  tcga_kirc1 <- tcga_kirc0 %>%
    dplyr::filter(!is.na(two1)) %>%
    dplyr::filter(!(two1 == "n/a")) %>%
    dplyr::mutate(two1 = factor(two1))
  
  # Adjust p values
  pairwise.test = tcga_kirc1 %>%
    wilcox_test(Stroma~two1) %>%
    adjust_pvalue(method = 'BH') %>%
    mutate(p.adj = ifelse(p.adj < 0.001, "***",
                          ifelse(p.adj < 0.01, "**",
                                 ifelse(p.adj < 0.05, "*", round(p.adj, 2))))) %>%
    mutate(p = ifelse(p < 0.001, "***",
                      ifelse(p < 0.01, "**",
                             ifelse(p < 0.05, "*", round(p, 2)))))
  
  
  # Plot
  g <- ggplot(tcga_kirc1, aes(x = two1, y = Stroma)) +
    geom_jitter(size=5, width = 0.2, aes(fill=two1), shape = 21, color = "black") +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
    labs(y="Stroma proportion (%)", x=gsub("_", " ", gsub("mi_rna", "miRNA", two1))) +
    theme_bw() +
    theme(axis.text.x = element_text(size=12, colour = "black"),
          axis.text.y = element_text(size=12, colour = "black"),
          axis.title=element_text(size=14, face="bold", colour = "black"),
          legend.position = "none") +
    scale_fill_brewer(palette = c("Set1")) +
    stat_compare_means(method = "kruskal.test",
                       label.y = 1.45*a,
                       size = 6) +
    stat_pvalue_manual(
      pairwise.test,
      label = "p.adj",
      bracket.size = 1,
      size = 5,
      # y.position = c(1.12*a, 1.24*a, 1.3*a, 1.0*a, 1.18*a, 1.06*a)
      y.position = c(1.14*a, 1.28*a, 1.35*a, 1.0*a, 1.21*a, 1.07*a)
    )
  ggsave(plot = g, filename = paste0("Textures/Images/Stroma/Mut/Stroma_", two1, "_cluster.png"), width = 5, height = 5, units = 'in', dpi = 300)
  
}







##############################
##### STROMA WITH NORMAL #####
##############################


############################# LOAD DATA ##########################################################################################################


# Image data
tcga_kirc <- read_xlsx("../data/mutations_final.xlsx")


# Normalize by removing empty
tcga_kirc$`texture_blood_%` <- 100*tcga_kirc$texture_blood / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`texture_cancer_%` <- 100*tcga_kirc$texture_cancer / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`texture_normal_%` <- 100*tcga_kirc$texture_normal / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`texture_stroma_%` <- 100*tcga_kirc$texture_stroma / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`texture_other_%` <- 100*tcga_kirc$texture_other / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)

tcga_kirc <- tcga_kirc %>%
  dplyr::filter(`texture_cancer_%` > 5) %>%
  janitor::clean_names() %>%
  dplyr::mutate(tissue_source_site = gsub("-[[:print:]]{4}", "", gsub("TCGA-", "", tcga_id)))


############################# PLOT ##########################################################################################################


# Rename
colnames(tcga_kirc)[grep(pattern = "^texture_[[:print:]]*percent", colnames(tcga_kirc))] <- c("Blood", "Cancer", "Normal", "Stroma", "Other")
textures <- c("Blood", "Cancer", "Normal", "Stroma", "Other")


# Function to find the colnumber
findcolnumber <- function(df, thecolumnname){
  which(colnames(df) == deparse(substitute(thecolumnname)))
}

# Filter only samples with <1% normal tissue
tcga_kirc <- tcga_kirc %>% dplyr::filter(Normal >= 1)


# Select variables
tcga_kirc <- tcga_kirc %>% dplyr::select(tcga_id, purity, ploidy, contains("cluster"), contains("mutation"), texture_blood:ncol(tcga_kirc))

# Modify data
tcga_kirc <- tcga_kirc %>% 
  dplyr::mutate_at(.vars = c("purity", "ploidy"), as.numeric)


tcga_kirc <- tcga_kirc %>% 
  dplyr::mutate_at(.vars = c("m_rna_cluster", "mi_rna_cluster", "methylation_cluster"), function(x) gsub(pattern = "Unavailable", "", x)) %>%
  dplyr::mutate(mi_rna_cluster = ifelse(mi_rna_cluster %in% c(2,5,6), "Other", mi_rna_cluster),
                m_rna_cluster = ifelse(m_rna_cluster %in% c(3,4,5), "Other", m_rna_cluster),
                m_rna_cluster = ifelse(m_rna_cluster %in% c("", "n/a"), NA, m_rna_cluster),
                mi_rna_cluster = ifelse(mi_rna_cluster %in% c("", "n/a"), NA, mi_rna_cluster),
                methylation_cluster = ifelse(methylation_cluster %in% c("", "n/a"), NA, methylation_cluster)) %>%
  dplyr::mutate_at(.vars = c("m_rna_cluster", "mi_rna_cluster", "methylation_cluster"), as.factor) %>%
  dplyr::mutate_at(.vars = vars(ends_with("_mutation")), function(x) ifelse(is.na(x), "WT",
                                                                            ifelse(x == "Unavailable", NA, "MUT"))) %>%
  dplyr::mutate(mutations_total = ifelse(mutations_total == "Unavailable", NA, as.numeric(as.character(mutations_total))),
                mutations_total = ifelse(mutations_total > median(mutations_total, na.rm=TRUE), "High", "Low")) %>%
  dplyr::mutate(ploidy = ifelse(ploidy == "Unavailable", NA, as.numeric(as.character(ploidy))),
                ploidy = ifelse(ploidy > median(ploidy, na.rm=TRUE), "High", "Low"))


# Multiple comparison (two-sided, unpaired, Mann-Whitney U test)
# gexp_data <- tcga_kirc[(findcolnumber(tcga_kirc, tissue_source_site)+1):ncol(tcga_kirc)]
multiple_t_tests_p_value <- lapply(tcga_kirc[c((str_detect(colnames(tcga_kirc), "_mutation") & !colnames(tcga_kirc)=="met_mutation") | colnames(tcga_kirc)=="mutations_total" | colnames(tcga_kirc)=="ploidy")],
                                   function(x) wilcox.test(tcga_kirc$Stroma ~ x, correct=FALSE, exact=FALSE, estimate=TRUE, na.rm=TRUE))
## P-values can be extracted from the result object
pvalue <- data.frame(p.value = sapply(multiple_t_tests_p_value, getElement, name = "p.value"))
## Create a matrix and dataframe of the p-values
pvalue_mat <- as.matrix(pvalue)
pvalue_df <- data.frame(pvalue)
## Add the p values to a new dataframe
pvalue_df$p_adjusted <- p.adjust(p = pvalue_mat, method = "BH")
## Add also the t values, 95%CI to the same dataframe
pvalue_df$median1 <- unlist(lapply(tcga_kirc[c((str_detect(colnames(tcga_kirc), "_mutation") & !colnames(tcga_kirc)=="met_mutation") | colnames(tcga_kirc)=="mutations_total" | colnames(tcga_kirc)=="ploidy")], function(x) median(tcga_kirc[x %in% c("High", "MUT"),]$Stroma, na.rm=TRUE)))
pvalue_df$median2 <- unlist(lapply(tcga_kirc[c((str_detect(colnames(tcga_kirc), "_mutation") & !colnames(tcga_kirc)=="met_mutation") | colnames(tcga_kirc)=="mutations_total" | colnames(tcga_kirc)=="ploidy")], function(x) median(tcga_kirc[x %in% c("Low", "WT"),]$Stroma, na.rm=TRUE)))
## Rownames to column
pvalue_df <- pvalue_df %>%
  rownames_to_column() %>%
  rename(genes = rowname,
         pvalue = p.value
  ) %>%
  arrange(pvalue)
rownames(pvalue_df) = NULL


# Rename variables
# pvalue_df <- pvalue_df %>%
#   mutate(genes = gsub("_mutation", "", genes),
#          genes = toupper(genes)
#   )


# Export data
writexl::write_xlsx(pvalue_df, "Textures/Images/Stroma_WithNormal/Mut/Stroma_mut.xlsx")
tcga_kirc0 <- tcga_kirc


# Plot parameters
a <- max(tcga_kirc$Stroma, na.rm=TRUE)


for (two1 in pvalue_df$genes) {
  
  # Assign value
  tcga_kirc0$two1 <- tcga_kirc0[[two1]]
  
  # Filter NAs
  tcga_kirc1 <- tcga_kirc0 %>%
    dplyr::filter(!is.na(two1)) %>%
    dplyr::mutate(two1 = as.factor(two1))
  
  # Adjust p values
  pairwise.test = tcga_kirc1 %>%
    wilcox_test(Stroma~two1) %>%
    adjust_pvalue(method = 'BH') %>%
    mutate(p.adj = ifelse(p.adj < 0.001, "***",
                          ifelse(p.adj < 0.01, "**",
                                 ifelse(p.adj < 0.05, "*", round(p.adj, 2))))) %>%
    mutate(p = ifelse(p < 0.001, "***",
                      ifelse(p < 0.01, "**",
                             ifelse(p < 0.05, "*", round(p, 2)))))
  
  
  # Plot
  g <- ggplot(tcga_kirc1, aes(x = two1, y = Stroma)) +
    geom_jitter(size=5, width = 0.2, aes(fill=two1), shape = 21, color = "black") +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
    labs(y="Stroma proportion (%)", x=toupper(gsub("_mutation", "", two1))) +
    theme_bw() +
    theme(axis.text.x = element_text(size=12, colour = "black"),
          axis.text.y = element_text(size=12, colour = "black"),
          axis.title=element_text(size=14, face="bold", colour = "black"),
          legend.position = "none") +
    scale_fill_brewer(palette = c("Set1")) +
    # scale_x_discrete(labels=c("0" = "Low", "1" = "High")) +
    # stat_compare_means(method = "wilcoxon.test",
    #                    label.y = 1.0*a,
    #                    label = "p.signif",
    #                    bracket.size = 1,
    #                    size = 5) +
    stat_pvalue_manual(
      pairwise.test,
      label = "p",
      bracket.size = 1,
      size = 5,
      # y.position = c(1.12*a, 1.24*a, 1.3*a, 1.0*a, 1.18*a, 1.06*a)
      y.position = 1.0*a
    )
  ggsave(plot = g, filename = paste0("Textures/Images/Stroma_WithNormal/Mut/Stroma_", two1, "_mut.png"), width = 5, height = 5, units = 'in', dpi = 300)
  
}


for (two1 in c("m_rna_cluster", "methylation_cluster")) {
  
  # Assign value
  tcga_kirc0$two1 <- tcga_kirc0[[two1]]
  
  # Filter NAs
  tcga_kirc1 <- tcga_kirc0 %>%
    dplyr::filter(!is.na(two1)) %>%
    dplyr::mutate(two1 = as.factor(two1))
  
  # Adjust p values
  pairwise.test = tcga_kirc1 %>%
    wilcox_test(Stroma~two1) %>%
    adjust_pvalue(method = 'BH') %>%
    mutate(p.adj = ifelse(p.adj < 0.001, "***",
                          ifelse(p.adj < 0.01, "**",
                                 ifelse(p.adj < 0.05, "*", round(p.adj, 2))))) %>%
    mutate(p = ifelse(p < 0.001, "***",
                      ifelse(p < 0.01, "**",
                             ifelse(p < 0.05, "*", round(p, 2)))))
  
  
  # Plot
  g <- ggplot(tcga_kirc1, aes(x = two1, y = Stroma)) +
    geom_jitter(size=5, width = 0.2, aes(fill=two1), shape = 21, color = "black") +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
    labs(y="Stroma proportion (%)", x=gsub("_", " ", gsub("methylation", "Methylation", gsub("m_rna", "mRNA", two1)))) +
    theme_bw() +
    theme(axis.text.x = element_text(size=12, colour = "black"),
          axis.text.y = element_text(size=12, colour = "black"),
          axis.title=element_text(size=14, face="bold", colour = "black"),
          legend.position = "none") +
    scale_fill_brewer(palette = c("Set1")) +
    stat_compare_means(method = "kruskal.test",
                       label.y = 1.25*a,
                       size = 6) +
    stat_pvalue_manual(
      pairwise.test,
      label = "p.adj",
      bracket.size = 1,
      size = 5,
      y.position = c(1.0*a, 1.12*a, 1.06*a)
    )
  ggsave(plot = g, filename = paste0("Textures/Images/Stroma_WithNormal/Mut/Stroma_", two1, "_cluster.png"), width = 5, height = 5, units = 'in', dpi = 300)
  
}


for (two1 in c("mi_rna_cluster")) {
  
  # Assign value
  tcga_kirc0$two1 <- tcga_kirc0[[two1]]
  
  # Filter NAs
  tcga_kirc1 <- tcga_kirc0 %>%
    dplyr::filter(!is.na(two1)) %>%
    dplyr::filter(!(two1 == "n/a")) %>%
    dplyr::mutate(two1 = factor(two1))
  
  # Adjust p values
  pairwise.test = tcga_kirc1 %>%
    wilcox_test(Stroma~two1) %>%
    adjust_pvalue(method = 'BH') %>%
    mutate(p.adj = ifelse(p.adj < 0.001, "***",
                          ifelse(p.adj < 0.01, "**",
                                 ifelse(p.adj < 0.05, "*", round(p.adj, 2))))) %>%
    mutate(p = ifelse(p < 0.001, "***",
                      ifelse(p < 0.01, "**",
                             ifelse(p < 0.05, "*", round(p, 2)))))
  
  
  # Plot
  g <- ggplot(tcga_kirc1, aes(x = two1, y = Stroma)) +
    geom_jitter(size=5, width = 0.2, aes(fill=two1), shape = 21, color = "black") +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
    labs(y="Stroma proportion (%)", x=gsub("_", " ", gsub("mi_rna", "miRNA", two1))) +
    theme_bw() +
    theme(axis.text.x = element_text(size=12, colour = "black"),
          axis.text.y = element_text(size=12, colour = "black"),
          axis.title=element_text(size=14, face="bold", colour = "black"),
          legend.position = "none") +
    scale_fill_brewer(palette = c("Set1")) +
    stat_compare_means(method = "kruskal.test",
                       label.y = 1.45*a,
                       size = 6) +
    stat_pvalue_manual(
      pairwise.test,
      label = "p.adj",
      bracket.size = 1,
      size = 5,
      # y.position = c(1.12*a, 1.24*a, 1.3*a, 1.0*a, 1.18*a, 1.06*a)
      y.position = c(1.14*a, 1.28*a, 1.35*a, 1.0*a, 1.21*a, 1.07*a)
    )
  ggsave(plot = g, filename = paste0("Textures/Images/Stroma_WithNormal/Mut/Stroma_", two1, "_cluster.png"), width = 5, height = 5, units = 'in', dpi = 300)
  
}
