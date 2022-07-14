rm(list=ls())

print("Start Textures/Scripts/Margin/Margin_texture_mut.R")

# Load libraries
library(readxl)
library(tidyverse)
library(data.table)
library(ggrepel)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(RColorBrewer)


############################# LOAD DATA ##########################################################################################################


# Image data
tcga_kirc <- read_xlsx("../data/image_analysis_results_final.xlsx")

# Normalize by removing empty
tcga_kirc$`texture_cancer_%` <- 100*tcga_kirc$texture_cancer / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`texture_normal_%` <- 100*tcga_kirc$texture_normal / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`margin_texture_blood_%` <- 100*tcga_kirc$margin_texture_blood / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)
tcga_kirc$`margin_texture_cancer_%` <- 100*tcga_kirc$margin_texture_cancer / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)
tcga_kirc$`margin_texture_normal_%` <- 100*tcga_kirc$margin_texture_normal / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)
tcga_kirc$`margin_texture_stroma_%` <- 100*tcga_kirc$margin_texture_stroma / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)
tcga_kirc$`margin_texture_other_%` <- 100*tcga_kirc$margin_texture_other / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)
tcga_kirc$`non_margin_texture_blood_%` <- 100*tcga_kirc$non_margin_texture_blood / (tcga_kirc$non_margin_texture_blood + tcga_kirc$non_margin_texture_normal + tcga_kirc$non_margin_texture_stroma + tcga_kirc$non_margin_texture_other)
tcga_kirc$`non_margin_texture_cancer_%` <- 100*tcga_kirc$non_margin_texture_cancer / (tcga_kirc$non_margin_texture_blood + tcga_kirc$non_margin_texture_cancer + tcga_kirc$non_margin_texture_normal + tcga_kirc$non_margin_texture_stroma + tcga_kirc$non_margin_texture_other)
tcga_kirc$`non_margin_texture_normal_%` <- 100*tcga_kirc$non_margin_texture_normal / (tcga_kirc$non_margin_texture_blood + tcga_kirc$non_margin_texture_normal + tcga_kirc$non_margin_texture_stroma + tcga_kirc$non_margin_texture_other)
tcga_kirc$`non_margin_texture_stroma_%` <- 100*tcga_kirc$non_margin_texture_stroma / (tcga_kirc$non_margin_texture_blood + tcga_kirc$non_margin_texture_normal + tcga_kirc$non_margin_texture_stroma + tcga_kirc$non_margin_texture_other)
tcga_kirc$`non_margin_texture_other_%` <- 100*tcga_kirc$non_margin_texture_other / (tcga_kirc$non_margin_texture_blood + tcga_kirc$non_margin_texture_normal + tcga_kirc$non_margin_texture_stroma + tcga_kirc$non_margin_texture_other)


# Filter
tcga_kirc <- tcga_kirc %>%
  dplyr::filter(`texture_cancer_%` > 5) %>%
  dplyr::filter(`texture_normal_%` > 1) %>%
  dplyr::mutate(tissue_source_site = gsub("-[[:print:]]{4}", "", gsub("TCGA-", "", tcga_id)))

# Save data
tcga_kirc0 <- tcga_kirc


############################# PLOT ##########################################################################################################


# Rename
colnames(tcga_kirc)[grep(pattern = "[[:print:]]{0,4}margin_texture_[[:print:]]*%", colnames(tcga_kirc))] <- c("Margin_Blood", "NonMargin_Blood", "Margin_Cancer", "NonMargin_Cancer", "Margin_Normal", "NonMargin_Normal", "Margin_Stroma", "NonMargin_Stroma", "Margin_Other", "NonMargin_Other")
margin <- c("Margin_Blood", "Margin_Normal", "Margin_Stroma", "Margin_Other")
nonmargin <- c("NonMargin_Blood", "NonMargin_Normal", "NonMargin_Stroma", "NonMargin_Other")


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
g <- ggplot(tcga_kirc_long, aes(x = Margin, y = value, group = variable)) +
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
ggsave(plot = g, filename = "Textures/Images/Margin/TCGA_margin_Texture_scatterplot.png", width = 6, height = 6, units = 'in', dpi = 300)






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

# Read mutation data
tcga_kirc <- read_xlsx("../data/mutations_final.xlsx") %>%
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


# Join normal texture
tcga_kirc <- tcga_kirc %>% dplyr::left_join(tcga_kirc0 %>% dplyr::select(ID = tcga_id, `texture_normal_%`))


# Loop over all textures
## Save data
tcga_kirc0 <- tcga_kirc

for (texture1 in unique(tcga_kirc0$Texture) ) {
  # texture1 = unique(tcga_kirc0$Texture)[3]
  print(texture1)
  
  # Reset data
  if (texture1 == "Normal") {
    
    tcga_kirc <- tcga_kirc0 %>%
      dplyr::filter(`texture_normal_%` > 1) %>%
      dplyr::filter(Texture %in% texture1)
    
  } else {
  tcga_kirc <- tcga_kirc0 %>%
    dplyr::filter(Texture %in% texture1)
  }
  
  # Multiple comparison (two-sided, unpaired, Mann-Whitney U test)
  # gexp_data <- tcga_kirc[(findcolnumber(tcga_kirc, tissue_source_site)+1):ncol(tcga_kirc)]
  multiple_t_tests_p_value <- lapply(tcga_kirc[c(str_detect(colnames(tcga_kirc), "_mutation") | colnames(tcga_kirc)=="mutations_total" | colnames(tcga_kirc)=="ploidy")],
                                     function(x) wilcox.test(tcga_kirc$FC ~ x, correct=FALSE, exact=FALSE, estimate=TRUE, na.rm=TRUE))
  ## P-values can be extracted from the result object
  pvalue <- data.frame(p.value = sapply(multiple_t_tests_p_value, getElement, name = "p.value"))
  ## Create a matrix and dataframe of the p-values
  pvalue_mat <- as.matrix(pvalue)
  pvalue_df <- data.frame(pvalue)
  ## Add the p values to a new dataframe
  pvalue_df$p_adjusted <- p.adjust(p = pvalue_mat, method = "BH")
  ## Add also the t values, 95%CI to the same dataframe
  pvalue_df$median1 <- unlist(lapply(tcga_kirc[c(str_detect(colnames(tcga_kirc), "_mutation") | colnames(tcga_kirc)=="mutations_total" | colnames(tcga_kirc)=="ploidy")], function(x) median(tcga_kirc[x %in% c("High", "MUT"),]$FC, na.rm=TRUE)))
  pvalue_df$median2 <- unlist(lapply(tcga_kirc[c(str_detect(colnames(tcga_kirc), "_mutation") | colnames(tcga_kirc)=="mutations_total" | colnames(tcga_kirc)=="ploidy")], function(x) median(tcga_kirc[x %in% c("Low", "WT"),]$FC, na.rm=TRUE)))
  ## Extract frequencies
  pvalue_df1 <- unlist(lapply(tcga_kirc[c(str_detect(colnames(tcga_kirc), "_mutation") | colnames(tcga_kirc)=="mutations_total" | colnames(tcga_kirc)=="ploidy")], function(x) table(x %in% c("High", "MUT"))))
  pvalue_df$freq_TRUE <- pvalue_df1[grep(pattern = "TRUE", x = names(pvalue_df1))]
  pvalue_df$freq_FALSE <- pvalue_df1[grep(pattern = "FALSE", x = names(pvalue_df1))]
  pvalue_df$prop = pvalue_df$freq_TRUE / (pvalue_df$freq_TRUE + pvalue_df$freq_FALSE)
  ## Rownames to column
  pvalue_df <- pvalue_df %>%
    rownames_to_column() %>%
    rename(genes = rowname,
           pvalue = p.value
    ) %>%
    arrange(pvalue)
  rownames(pvalue_df) = NULL
  
  
  
  # Export data
  writexl::write_xlsx(pvalue_df, paste0("Textures/Images/Margin/Margin_vs_non_margin/Texture/Mut/Ly_in_margin_vs_nonmargin_", texture1, "_mut.xlsx"))
  
  
  # Plot parameters
  a <- max(tcga_kirc$FC, na.rm=TRUE)
  
  
  for (two1 in pvalue_df$genes) {
    
    
    # Reset data
    tcga_kirc1 <- tcga_kirc
    
    # Assign value
    # two1 = pvalue_df$genes[1]
    tcga_kirc1$two1 <- tcga_kirc1[[two1]]
    
    # Filter NAs
    tcga_kirc1 <- tcga_kirc1 %>%
      dplyr::filter(!is.na(two1)) %>%
      dplyr::mutate(two1 = as.factor(two1))
    
    
    # Adjust p values
    pairwise.test = tcga_kirc1 %>%
      wilcox_test(FC~two1) %>%
      adjust_pvalue(method = 'BH') %>%
      # mutate(p.adj = round(p.adj, 2))
      mutate(p = ifelse(p < 0.001, "***",
                        ifelse(p < 0.01, "**",
                               ifelse(p < 0.05, "*", round(p, 2)))))
    
    
    # Plot
    g <- ggplot(tcga_kirc1, aes(x = two1, y = FC)) +
      geom_jitter(size=5, width = 0.2, aes(fill=two1), shape = 21, color = "black") +
      geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
      labs(y=paste0(texture1, " texture proportion fold change (margin vs. non-margin)"), x=toupper(gsub("_mutation", "", two1))) +
      theme_bw() +
      theme(axis.text.x = element_text(size=12, colour = "black"),
            axis.text.y = element_text(size=12, colour = "black"),
            axis.title.y = element_text(size=10, face="bold", colour = "black"),
            axis.title.x = element_text(size=14, face="bold", colour = "black"),
            legend.position = "none") +
      scale_fill_brewer(palette = c("Set1")) +
      stat_pvalue_manual(
        pairwise.test,
        label = "p",
        bracket.size = 1.5,
        size = 5,
        # y.position = c(1.12*a, 1.24*a, 1.3*a, 1.0*a, 1.18*a, 1.06*a)
        y.position = 1.05*a
      )
    g
    ggsave(plot = g, filename = paste0("Textures/Images/Margin/Margin_vs_non_margin/Texture/Mut/Ly_in_margin_vs_nonmargin_", two1, "_", texture1, "_mut.png"), width = 5, height = 5, units = 'in', dpi = 300)
    
  }
  
}



################################ Prepare linear scatter plots for mutations ################################################################


# Load data
rm(pvalue_df, pvalue_df1)
for (texture1 in gsub("Margin_", "", margin)) {
  pvalue_df1 <- readxl::read_xlsx(paste0("Textures/Images/Margin/Margin_vs_non_margin/Texture/Mut/Ly_in_margin_vs_nonmargin_", texture1, "_mut.xlsx")) %>%
    mutate(Pvalue = ifelse(pvalue < 0.05, "p<0.05",
                        ifelse(pvalue < 0.10, "p<0.10", NA)),
           Texture = texture1) %>%
    dplyr::filter(!is.na(Pvalue))
  if (exists("pvalue_df")) {
    pvalue_df <- rbind(pvalue_df, pvalue_df1)
  } else {
    pvalue_df <- pvalue_df1
  }
}

# Modify data
pvalue_df <- pvalue_df %>%
  dplyr::mutate(
    genes = gsub("_mutation", "", genes),
    genes = toupper(genes),
    prop = prop * 100
  )


g <- ggplot(pvalue_df, aes(x = median1, y = median2)) +
  geom_point(aes(shape = Pvalue, fill = Texture, size = prop), color = "black") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x="Median margin:nonmargin texture ratio in mutated samples", y="Median margin:nonmargin texture ratio in wild-type samples") +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 3)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3)) +
  theme(axis.text.x = element_text(size=12, colour = "black"),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title.y = element_text(size=11, face="bold", colour = "black"),
        axis.title.x = element_text(size=11, face="bold", colour = "black"),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.key.size = unit(0.01, "cm"),
        legend.margin = margin(),
        legend.spacing.y = unit(0.1, "cm")) +
  scale_fill_manual(values = c("#e41a1c", "#4daf4a", "#984ea3", "#ff7f00")) +
  scale_shape_manual(values = c(21, 24)) +
  scale_size(breaks = c(10, 30, 50), range = c(1, 10)) +
  guides(fill = guide_legend(override.aes=list(shape=21)),
         size = guide_legend("Proportion (%)")) +
  geom_label_repel(segment.size = 1,
                   aes(label=genes), size = 4, fontface = "bold")
g
ggsave(plot = g, filename = "Textures/Images/Margin/Margin_vs_non_margin/Texture/Mut/Ly_in_margin_vs_nonmargin_mut_only_significant.png", width = 5.0, height = 5.5, units = 'in', dpi = 300)
