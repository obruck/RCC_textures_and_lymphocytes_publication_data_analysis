rm(list=ls())

print("Start Ly/Scripts/Ly_mut.R")

# Load libraries
library(readxl)
library(tidyverse)
library(data.table, warn.conflicts = TRUE)
library(ggplot2)
library(ggpubr)
library(janitor)
library(rstatix)
library(RColorBrewer)
library(fastDummies)
library(survival)
library(survminer)
library(sjPlot)


# Set theme
set_theme(theme_bw(), axis.title.color = "black", axis.title.size = 1.1, axis.textcolor = "black", axis.textsize = 1.05, geom.label.size = 3.5, legend.pos="bottom")


############################# LOAD DATA ##########################################################################################################


# Image data
tcga_kirc <- read_xlsx("../data/mutations_final.xlsx") %>%
  janitor::clean_names()

# Normalize by removing empty
tcga_kirc$texture_cancer_percent <- 100*tcga_kirc$texture_cancer / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$texture_normal_percent <- 100*tcga_kirc$texture_normal / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$inf_bin_lymphocytes_total <- 100*(tcga_kirc$bin_lymphocytes_blood + tcga_kirc$bin_lymphocytes_normal + tcga_kirc$bin_lymphocytes_cancer + tcga_kirc$bin_lymphocytes_stroma + tcga_kirc$bin_lymphocytes_other) / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)

tcga_kirc <- tcga_kirc %>%
  dplyr::filter(texture_cancer_percent > 5) %>%
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
# tcga_kirc <- tcga_kirc %>% dplyr::filter(!tissue_source_site %in% c("BP", "AK") )


# Save
tcga_kirc0 <- tcga_kirc


############################# MODIFY CLINICAL DATA ##########################################################################################################


# Rename
colnames(tcga_kirc)[grep(pattern = "^inf_bin_lymphocytes_[[:print:]]*", colnames(tcga_kirc))] <- c("Blood", "Cancer", "Normal", "Stroma", "Other", "Total")
textures <- c("Blood", "Cancer", "Normal", "Stroma", "Other", "Total")


############################# PLOT ##########################################################################################################


# To percent
tcga_kirc[textures] <- sapply(tcga_kirc[textures], function(x) 100*x)


# Function to find the colnumber
findcolnumber <- function(df, thecolumnname){
  which(colnames(df) == deparse(substitute(thecolumnname)))
}


# Select variables
tcga_kirc <- tcga_kirc %>%
  dplyr::select(tcga_id, !!textures, purity, ploidy, contains("cluster"), contains("mutation"), tissue_source_site)

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
  dplyr::mutate(mi_rna_cluster_1 = ifelse(mi_rna_cluster == 1, 1, 0),
                mi_rna_cluster_3 = ifelse(mi_rna_cluster == 3, 1, 0),
                mi_rna_cluster_4 = ifelse(mi_rna_cluster == 4, 1, 0),
                mi_rna_cluster_other = ifelse(mi_rna_cluster == "Other", 1, 0),
                m_rna_cluster_1 = ifelse(m_rna_cluster == 1, 1, 0),
                m_rna_cluster_2 = ifelse(m_rna_cluster == 2, 1, 0),
                m_rna_cluster_other = ifelse(m_rna_cluster == "Other", 1, 0),
                methylation_cluster_1 = ifelse(methylation_cluster == 1, 1, 0),
                methylation_cluster_2 = ifelse(methylation_cluster == 2, 1, 0),
                methylation_cluster_3 = ifelse(methylation_cluster == 3, 1, 0)) %>%
  dplyr::mutate_at(.vars = c("m_rna_cluster", "mi_rna_cluster", "methylation_cluster"), as.factor) %>%
  dplyr::mutate_at(.vars = vars(ends_with("_mutation")), function(x) ifelse(is.na(x), "WT",
                                                                            ifelse(x == "Unavailable", NA, "MUT"))) %>%
  dplyr::mutate(mutations_total = ifelse(mutations_total == "Unavailable", NA, as.numeric(as.character(mutations_total))),
                mutations_total = ifelse(mutations_total > median(mutations_total, na.rm=TRUE), "High", "Low")) %>%
  dplyr::mutate(ploidy = ifelse(ploidy == "Unavailable", NA, as.numeric(as.character(ploidy))),
                ploidy = ifelse(ploidy > median(ploidy, na.rm=TRUE), "High", "Low"))


# Loop over all textures
## Save data
tcga_kirc0 <- tcga_kirc

for (texture1 in textures) {
  # texture1 = textures[3]
  print(texture1)
  
  # Reset data
  if (texture1 == "Normal") {
    tcga_kirc <- tcga_kirc0 %>%
      dplyr::filter(Normal > 1)
  } else {
    tcga_kirc <- tcga_kirc0
  }
  
  # Divide patients into low and high by clinical centers
  tcga_kirc$FC <- tcga_kirc[[texture1]]
  tcga_kirc <- tcga_kirc %>%
    group_by(tissue_source_site) %>%
    mutate(FC = (FC > median(FC, na.rm = TRUE)) + 1) %>%
    mutate(FC = ifelse(FC == 2, "High", "Low")) %>%
    ungroup()
  
  # Variables with 2 values
  unique_values_2 = sapply(tcga_kirc[c(str_detect(colnames(tcga_kirc), "_mutation") | str_detect(colnames(tcga_kirc), "total") | str_detect(colnames(tcga_kirc), "ploidy") | str_detect(colnames(tcga_kirc), "cluster_"))], function(x) length(unique(x[!is.na(x)])))
  unique_values_2 = names(unique_values_2[unique_values_2==2])
  
  # Multiple comparison (chisq test)
  # multiple_t_tests_p_value <- lapply(tcga_kirc[unique_values_2],
  #                                    function(x) wilcox.test(tcga_kirc$FC ~ x, correct=FALSE, exact=FALSE, estimate=TRUE, na.rm=TRUE))
  multiple_t_tests_p_value <- lapply(tcga_kirc[unique_values_2],
                                     function(x) chisq.test(x, tcga_kirc$FC))
  ## P-values can be extracted from the result object
  pvalue <- data.frame(p.value = sapply(multiple_t_tests_p_value, getElement, name = "p.value"))
  ### Observed counts
  observed_counts = lapply(multiple_t_tests_p_value, getElement, name = "observed")
  ## Create a matrix and dataframe of the p-values
  pvalue_mat <- as.matrix(pvalue)
  pvalue_df <- data.frame(pvalue)
  ## Add the p values to a new dataframe
  pvalue_df$p_adjusted <- p.adjust(p = pvalue_mat, method = "BH")
  ## Add also the t values, 95%CI to the same dataframe
  pvalue_df$median1 <- NA; pvalue_df$median2 <- NA;
  for (i in 1:length(observed_counts)) {
    pvalue_df[i,]$median1 <- observed_counts[[i]][2] / (observed_counts[[i]][1] + observed_counts[[i]][2])
    pvalue_df[i,]$median2 <- observed_counts[[i]][4] / (observed_counts[[i]][3] + observed_counts[[i]][4])
  }
  ## Rownames to column
  pvalue_df <- pvalue_df %>%
    rownames_to_column() %>%
    rename(genes = rowname,
           pvalue = p.value
    ) %>%
    arrange(pvalue)
  rownames(pvalue_df) = NULL
  
  
  
  # Export data
  writexl::write_xlsx(pvalue_df, paste0("Ly/Images/Mut/Ly_", texture1, "_mut.xlsx"))
  
  
  # Plot parameters
  # a <- max(tcga_kirc$FC, na.rm=TRUE)
  
  
  for (two1 in pvalue_df$genes) {
    
    # Reset data
    tcga_kirc1 <- tcga_kirc
    
    # Assign value
    tcga_kirc1$two1 <- ifelse(tcga_kirc1[[two1]] %in% c("High", "MUT", 1), TRUE,
                              ifelse(tcga_kirc1[[two1]] %in% c("Low", "WT", 0), FALSE, NA))
    
    # Filter NAs
    tcga_kirc1 <- tcga_kirc1 %>%
      dplyr::filter(!is.na(two1)) %>%
      dplyr::mutate(two1 = as.factor(two1),
                    FC = factor(FC, levels = c("Low", "High")))
    
    # Barplot
    g <- sjp.xtab(
      tcga_kirc1$FC,
      tcga_kirc1$two1,
      show.total = FALSE,
      margin = "row",
      # bar.pos = "stack",
      bar.pos = "dodge",
      show.summary = TRUE,
      geom.colors=c("#1B9E77", "#D95F02"),
      axis.titles=c(paste0("Lymphocyte proportion in ", texture1), toupper(gsub("_mutation", "", two1))),
      legend.title="")
    
    ggsave(plot = g, filename = paste0("Ly/Images/Mut/Ly_", two1, "_", texture1, "_mut.png"), width = 4, height = 4.5, units = 'in', dpi = 300)
    
  }
  
}

