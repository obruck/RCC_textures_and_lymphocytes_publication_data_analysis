rm(list=ls())

# Load libraries
library(readxl)
library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(RColorBrewer)
library(fastDummies)


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
tcga_kirc$`non_margin_texture_normal_%` <- 100*tcga_kirc$non_margin_texture_normal / (tcga_kirc$non_margin_texture_blood + tcga_kirc$non_margin_texture_cancer + tcga_kirc$non_margin_texture_normal + tcga_kirc$non_margin_texture_stroma + tcga_kirc$non_margin_texture_other)
tcga_kirc$`non_margin_texture_blood_%` <- 100*tcga_kirc$non_margin_texture_blood / (tcga_kirc$non_margin_texture_blood + tcga_kirc$non_margin_texture_cancer + tcga_kirc$non_margin_texture_normal + tcga_kirc$non_margin_texture_stroma + tcga_kirc$non_margin_texture_other)
tcga_kirc$`non_margin_texture_cancer_%` <- 100*tcga_kirc$non_margin_texture_cancer / (tcga_kirc$non_margin_texture_blood + tcga_kirc$non_margin_texture_cancer + tcga_kirc$non_margin_texture_normal + tcga_kirc$non_margin_texture_stroma + tcga_kirc$non_margin_texture_other)
tcga_kirc$`non_margin_texture_normal_%` <- 100*tcga_kirc$non_margin_texture_normal / (tcga_kirc$non_margin_texture_blood + tcga_kirc$non_margin_texture_cancer + tcga_kirc$non_margin_texture_normal + tcga_kirc$non_margin_texture_stroma + tcga_kirc$non_margin_texture_other)
tcga_kirc$`non_margin_texture_stroma_%` <- 100*tcga_kirc$non_margin_texture_stroma / (tcga_kirc$non_margin_texture_blood + tcga_kirc$non_margin_texture_cancer + tcga_kirc$non_margin_texture_normal + tcga_kirc$non_margin_texture_stroma + tcga_kirc$non_margin_texture_other)
tcga_kirc$`non_margin_texture_other_%` <- 100*tcga_kirc$non_margin_texture_other / (tcga_kirc$non_margin_texture_blood + tcga_kirc$non_margin_texture_cancer + tcga_kirc$non_margin_texture_normal + tcga_kirc$non_margin_texture_stroma + tcga_kirc$non_margin_texture_other)

tcga_kirc <- tcga_kirc %>%
  dplyr::filter(`texture_cancer_%` > 5) %>%
  dplyr::filter(`margin_texture_normal_%` > 1) %>%
  dplyr::filter(is.na(PoorQuality)) %>%
  dplyr::mutate(tissue_source_site = gsub("-[[:print:]]{4}", "", gsub("TCGA-", "", tcga_id)))


############################# PLOT ##########################################################################################################


# Rename
colnames(tcga_kirc)[grep(pattern = "[[:print:]]{0,4}margin_texture_[[:print:]]*%", colnames(tcga_kirc))] <- c("Margin_Blood", "NonMargin_Blood", "Margin_Cancer", "NonMargin_Cancer", "Margin_Normal", "NonMargin_Normal", "Margin_Stroma", "NonMargin_Stroma", "Margin_Other", "NonMargin_Other")
textures <- c("Margin_Blood", "NonMargin_Blood", "Margin_Normal", "NonMargin_Normal", "Margin_Stroma", "NonMargin_Stroma", "Margin_Other", "NonMargin_Other")
# colnames(tcga_kirc)[grep(pattern = "^inf_margin_bin_lymphocytes", colnames(tcga_kirc))] <- c("Margin_Blood", "Margin_Normal", "Margin_Stroma", "Margin_Other")
# colnames(tcga_kirc)[grep(pattern = "^inf_non_margin_bin_lymphocytes", colnames(tcga_kirc))] <- c("NonMargin_Blood", "NonMargin_Cancer", "NonMargin_Normal", "NonMargin_Stroma", "NonMargin_Other")
# margin <- c("Margin_Blood", "Margin_Normal", "Margin_Stroma", "Margin_Other")
# nonmargin <- c("NonMargin_Blood", "NonMargin_Normal", "NonMargin_Stroma", "NonMargin_Other")
# textures <- c(margin, nonrmargin)


# Petri's data
df <- readRDS("/Users/oscarbruck/OneDrive - University of Helsinki/RCC/Otso/Data/petri/KIRC.fm.rds")
# Filter patients
df <- df %>%
  dplyr::select(one_of(tcga_kirc$tcga_id))

# Clean clinical data
df1 <- df[grep(pattern = "CLIN:", x = rownames(df)),]
df1 <- df1 %>% rownames_to_column()
df1 <- df1[apply(df1[,-1], 1, function(x) !all(x==0)),]
df1 <- df1[apply(df1[,-1], 1, function(x) !all(x==1)),]
df1 <- df1 %>% filter_all(any_vars(!is.na(.)))
df2 <- df1 %>% 
  pivot_longer(!rowname, names_to = "col1", values_to = "col2") %>% 
  pivot_wider(names_from = "rowname", values_from = "col2") %>%
  dplyr::rename(tcga_id = col1)

# Join
tcga_kirc <- tcga_kirc %>%
  dplyr::left_join(df2)

# Multiple binary to one categorical
# colnames(tcga_kirc[colnames(tcga_kirc)[grep(pattern = "CLIN:pathologic_T", colnames(tcga_kirc))]])
## CLIN:pathologic_T
cols1 = grep(pattern = "CLIN:pathologic_T", colnames(tcga_kirc))
mc = max.col(tcga_kirc[cols1], ties.method = "first")
tcga_kirc$Tumor_size <- gsub("B:CLIN:pathologic_T_", "", colnames(tcga_kirc[cols1])[mc])
tcga_kirc$Tumor_size <- substr(tcga_kirc$Tumor_size, 1, 2)

## CLIN:pathologic_N
cols1 = grep(pattern = "CLIN:pathologic_N", colnames(tcga_kirc))
mc = max.col(tcga_kirc[cols1], ties.method = "first")
tcga_kirc$Lymph_node_metastasis <- gsub("B:CLIN:pathologic_N_", "", colnames(tcga_kirc[cols1])[mc])

## CLIN:pathologic_M
cols1 = grep(pattern = "CLIN:pathologic_M", colnames(tcga_kirc))
mc = max.col(tcga_kirc[cols1], ties.method = "first")
tcga_kirc$Organ_metastasis <- gsub("B:CLIN:pathologic_M_", "", colnames(tcga_kirc[cols1])[mc])

## CLIN:clinical_M
cols1 = grep(pattern = "CLIN:clinical_M", colnames(tcga_kirc))
mc = max.col(tcga_kirc[cols1], ties.method = "first")
tcga_kirc$Organ_metastasis_ <- gsub("B:CLIN:clinical_M_", "", colnames(tcga_kirc[cols1])[mc])

## CLIN:white_cell_count_result
cols1 = grep(pattern = "CLIN:white_cell_count_result", colnames(tcga_kirc))
mc = max.col(tcga_kirc[cols1], ties.method = "first")
tcga_kirc$WBC <- gsub("B:CLIN:white_cell_count_result_", "", colnames(tcga_kirc[cols1])[mc])
tcga_kirc$WBC <- ifelse(tcga_kirc$WBC == "[Unknown]", NA, tcga_kirc$WBC)
tcga_kirc <- tcga_kirc %>% mutate(WBC = factor(WBC, levels = c("Low", "Normal", "Elevated", NA)))

## CLIN:gender
cols1 = grep(pattern = "CLIN:gender", colnames(tcga_kirc))
mc = max.col(tcga_kirc[cols1], ties.method = "first")
tcga_kirc$gender <- gsub("B:CLIN:gender_", "", colnames(tcga_kirc[cols1])[mc])

## CLIN:vital_status
cols1 = grep(pattern = "CLIN:vital_status", colnames(tcga_kirc))
mc = max.col(tcga_kirc[cols1], ties.method = "first")
tcga_kirc$vital_status <- gsub("B:CLIN:vital_status_", "", colnames(tcga_kirc[cols1])[mc])

## CLIN:race
cols1 = grep(pattern = "CLIN:race", colnames(tcga_kirc))
mc = max.col(tcga_kirc[cols1], ties.method = "first")
tcga_kirc$race <- gsub("B:CLIN:race_", "", colnames(tcga_kirc[cols1])[mc])
tcga_kirc$race <- ifelse(tcga_kirc$`B:CLIN:ethnicity_HISPANIC OR LATINO` == 1, "HISPANIC OR LATINO", tcga_kirc$race)
tcga_kirc$race <- ifelse(tcga_kirc$race == "BLACK OR AFRICAN AMERICAN", "AFRICAN AMERICAN",
                         ifelse(tcga_kirc$race == "HISPANIC OR LATINO", "HISPANIC", tcga_kirc$race))
tcga_kirc$race <- str_to_sentence(tcga_kirc$race)

## CLIN:pathologic_stage
cols1 = grep(pattern = "CLIN:pathologic_stage", colnames(tcga_kirc))
mc = max.col(tcga_kirc[cols1], ties.method = "first")
tcga_kirc$stage <- gsub("B:CLIN:pathologic_stage_", "", colnames(tcga_kirc[cols1])[mc])

## CLIN:primary_therapy_outcome_success
cols1 = grep(pattern = "CLIN:primary_therapy_outcome_success", colnames(tcga_kirc))
mc = max.col(tcga_kirc[cols1], ties.method = "first")
tcga_kirc$therapy_outcome <- gsub("B:CLIN:primary_therapy_outcome_success_", "", colnames(tcga_kirc[cols1])[mc])
tcga_kirc$therapy_outcome <- ifelse(tcga_kirc$therapy_outcome == "[Unknown]", NA, tcga_kirc$therapy_outcome)
tcga_kirc <- tcga_kirc %>% mutate(therapy_outcome = factor(therapy_outcome, levels = c("Complete Remission/Response", "Stable Disease", "Progressive Disease", NA)))

## CLIN:eastern_cancer_oncology_group
cols1 = grep(pattern = "CLIN:eastern_cancer_oncology_group", colnames(tcga_kirc))
mc = max.col(tcga_kirc[cols1], ties.method = "first")
tcga_kirc$ECOG <- gsub("B:CLIN:eastern_cancer_oncology_group_", "", colnames(tcga_kirc[cols1])[mc])
tcga_kirc$ECOG <- ifelse(tcga_kirc$ECOG == "[Unknown]", NA, tcga_kirc$ECOG)
tcga_kirc$ECOG <- ifelse(tcga_kirc$ECOG == "2", NA, tcga_kirc$ECOG)
# tcga_kirc$ECOG <- as.factor(tcga_kirc$ECOG))

## CLIN:karnofsky_performance_score
cols1 = grep(pattern = "CLIN:karnofsky_performance_score", colnames(tcga_kirc))
mc = max.col(tcga_kirc[cols1], ties.method = "first")
tcga_kirc$Karnofsky <- gsub("B:CLIN:karnofsky_performance_score_", "", colnames(tcga_kirc[cols1])[mc])
tcga_kirc$Karnofsky <- ifelse(tcga_kirc$Karnofsky == "[Unknown]", NA, tcga_kirc$Karnofsky)
tcga_kirc$Karnofsky <- as.numeric(as.character(tcga_kirc$Karnofsky))
tcga_kirc$Karnofsky <- ifelse(tcga_kirc$Karnofsky > 85, ">80", "≤80")
tcga_kirc <- tcga_kirc %>% mutate(Karnofsky = factor(Karnofsky, levels = c("≤80", ">80")))

## CLIN:neoplasm_histologic_grade
cols1 = grep(pattern = "CLIN:neoplasm_histologic_grade", colnames(tcga_kirc))
mc = max.col(tcga_kirc[cols1], ties.method = "first")
tcga_kirc$grade <- gsub("B:CLIN:neoplasm_histologic_grade_", "", colnames(tcga_kirc[cols1])[mc])
tcga_kirc$grade <- ifelse(tcga_kirc$grade == "GX", NA, tcga_kirc$grade)
tcga_kirc <- tcga_kirc %>% mutate(grade = factor(grade, levels = c("G1", "G2", "G3", "G4")))

## CLIN:tobacco_smoking_history
cols1 = grep(pattern = "CLIN:tobacco_smoking_history", colnames(tcga_kirc))
mc = max.col(tcga_kirc[cols1], ties.method = "first")
tcga_kirc$smoking_status <- gsub("B:CLIN:tobacco_smoking_history_", "", colnames(tcga_kirc[cols1])[mc])
tcga_kirc$smoking_status <- ifelse(tcga_kirc$smoking_status == "[Unknown]", NA, tcga_kirc$smoking_status)
tcga_kirc$smoking_status <- ifelse(str_detect(tcga_kirc$smoking_status, "eformed"), "Ex-smoker", tcga_kirc$smoking_status)
tcga_kirc$smoking_status <- gsub("Lifelong ", "", tcga_kirc$smoking_status)
tcga_kirc <- tcga_kirc %>% mutate(smoking_status = factor(smoking_status, levels = c("Non-smoker", "Current smoker", "Ex-smoker", NA)))

## CLIN:lactate_dehydrogenase
cols1 = grep(pattern = "CLIN:lactate_dehydrogenase", colnames(tcga_kirc))
mc = max.col(tcga_kirc[cols1], ties.method = "first")
tcga_kirc$LDH <- gsub("B:CLIN:lactate_dehydrogenase_result_", "", colnames(tcga_kirc[cols1])[mc])
tcga_kirc$LDH <- ifelse(tcga_kirc$LDH == "[Unknown]", NA, tcga_kirc$LDH)
tcga_kirc <- tcga_kirc %>% mutate(LDH = factor(LDH, levels = c("Normal", "Elevated", NA)))

## CLIN:hemoglobin_result
cols1 = grep(pattern = "CLIN:hemoglobin_result", colnames(tcga_kirc))
mc = max.col(tcga_kirc[cols1], ties.method = "first")
tcga_kirc$HB <- gsub("B:CLIN:hemoglobin_result_", "", colnames(tcga_kirc[cols1])[mc])
tcga_kirc$HB <- ifelse(tcga_kirc$HB == "[Unknown]", NA, tcga_kirc$HB)
tcga_kirc <- tcga_kirc %>% mutate(HB = factor(HB, levels = c("Low", "Normal", "Elevated", NA)))

## CLIN:serum_calcium_result
cols1 = grep(pattern = "CLIN:serum_calcium_result", colnames(tcga_kirc))
mc = max.col(tcga_kirc[cols1], ties.method = "first")
tcga_kirc$calcium <- gsub("B:CLIN:serum_calcium_result_", "", colnames(tcga_kirc[cols1])[mc])
tcga_kirc$calcium <- ifelse(tcga_kirc$calcium == "[Unknown]", NA, tcga_kirc$calcium)
tcga_kirc <- tcga_kirc %>% mutate(calcium = factor(calcium, levels = c("Low", "Normal", "Elevated", NA)))

## CLIN:platelet_qualitative
cols1 = grep(pattern = "CLIN:platelet_qualitative", colnames(tcga_kirc))
mc = max.col(tcga_kirc[cols1], ties.method = "first")
tcga_kirc$PLT <- gsub("B:CLIN:platelet_qualitative_result_", "", colnames(tcga_kirc[cols1])[mc])
tcga_kirc$PLT <- ifelse(tcga_kirc$PLT == "[Unknown]", NA, tcga_kirc$PLT)
tcga_kirc <- tcga_kirc %>% mutate(PLT = factor(PLT, levels = c("Low", "Normal", "Elevated", NA)))

# Summary
sapply(tcga_kirc[, c(which( colnames(tcga_kirc)=="Tumor_size" ) : which( colnames(tcga_kirc)=="PLT" ) )], function(x) length(unique(x)))
sapply(tcga_kirc[, c(which( colnames(tcga_kirc)=="Tumor_size" ) : which( colnames(tcga_kirc)=="PLT" ) )], function(x) unique(x))
twos = c("gender", "vital_status", "LDH", "ECOG", "Karnofsky")
threes = c("Lymph_node_metastasis", "Organ_metastasis", "Organ_metastasis_", "WBC", "therapy_outcome", "smoking_status", "HB", "calcium", "PLT")
fours = c("Tumor_size", "race", "stage", "grade")


############################# PLOT ##########################################################################################################



# To percent
# tcga_kirc[margin] <- sapply(tcga_kirc[margin], function(x) 100*x)
# tcga_kirc[nonmargin] <- sapply(tcga_kirc[nonmargin], function(x) 100*x)



## Variables to hotencode
to_hotencode <- tcga_kirc[c(twos, threes, fours)]
### k-cluster make one hot encoder
to_hotencode_var <- dummy_cols(tcga_kirc[colnames(to_hotencode)]) %>%
  dplyr::select(-ends_with("_NA"))
to_hotencode_var <- to_hotencode_var %>% dplyr::select(-(1:ncol(to_hotencode)))
### Replace NA with 0
# for(i in 1:ncol(to_hotencode_var)){
#   to_hotencode_var[is.na(to_hotencode_var[,i]), i] = 0
# }
### Remove duplicate binary variables
to_hotencode_var <- to_hotencode_var %>% dplyr::select(-c(gender_FEMALE, vital_status_Alive, LDH_Normal, ECOG_0, `Karnofsky_≤80`))
### Join
tcga_kirc <- cbind(tcga_kirc, to_hotencode_var)


# Loop over margin textures
for (margin_var in c(textures)) {
  # margin_var = c(margin, nonmargin)[1]
  tcga_kirc$margin_var1 = tcga_kirc[[margin_var]]
  
  # Multiple comparison (two-sided, unpaired, Mann-Whitney U test)
  # gexp_data <- tcga_kirc[(findcolnumber(tcga_kirc, tissue_source_site)+1):ncol(tcga_kirc)]
  multiple_t_tests_p_value <- lapply(tcga_kirc[colnames(to_hotencode_var)],
                                     function(x) wilcox.test(tcga_kirc$margin_var1 ~ x, correct=FALSE, exact=FALSE, estimate=TRUE, na.rm=TRUE))
  ## P-values can be extracted from the result object
  pvalue <- data.frame(p.value = sapply(multiple_t_tests_p_value, getElement, name = "p.value"))
  ## Create a matrix and dataframe of the p-values
  pvalue_mat <- as.matrix(pvalue)
  pvalue_df <- data.frame(pvalue)
  ## Add the p values to a new dataframe
  pvalue_df$p_adjusted <- p.adjust(p = pvalue_mat, method = "BH")
  ## Add also the t values, 95%CI to the same dataframe
  pvalue_df$median1 <- unlist(lapply(tcga_kirc[colnames(to_hotencode_var)], function(x) median(tcga_kirc[x == 1,]$margin_var1, na.rm=TRUE)))
  pvalue_df$median0 <- unlist(lapply(tcga_kirc[colnames(to_hotencode_var)], function(x) median(tcga_kirc[x == 0,]$margin_var1, na.rm=TRUE)))
  ## Rownames to column
  pvalue_df <- pvalue_df %>%
    rownames_to_column() %>%
    rename(genes = rowname,
           pvalue = p.value
    ) %>%
    arrange(pvalue) %>%
    dplyr::mutate(margin_var = margin_var)
  rownames(pvalue_df) = NULL
  ## Join
  if (!exists("pvalue_df0")) {
    pvalue_df0 <- pvalue_df
  } else {
    pvalue_df0 <- rbind(pvalue_df0, pvalue_df)
  }
  
}
pvalue_df <- pvalue_df0; rm(pvalue_df0)

# Adjust p value
pvalue_df$p_adjusted1 <- p.adjust(pvalue_df$pvalue, method = "BH")

# Filter variables with p_adjusted1 < 0.10
vars <- pvalue_df %>%
  dplyr::filter(p_adjusted1 < 0.10) %>%
  dplyr::select(genes) %>%
  dplyr::mutate(genes = gsub("_[[:print:]]*", "", genes)) %>%
  dplyr::distinct()
tmp <- data.frame(twos); tmp <- rbind(tmp, data.frame(threes) %>% rename(twos = threes)); tmp <- rbind(tmp, data.frame(fours) %>% rename(twos = fours))
for (i in vars$genes) {
  a <- tmp %>% dplyr::filter(str_detect(twos, i))
  ## Join
  if (!exists("b")) {
    b <- data.frame(a)
  } else {
    b <- rbind(b, a)
  }
}
a <- b; rm(b)
b <- as.character(a$twos)

# Melt longer
tcga_kirc_long <- tcga_kirc %>%
  dplyr::select("tcga_id", "Margin_Blood", "Margin_Normal", "Margin_Stroma", "Margin_Other", "NonMargin_Blood", "NonMargin_Normal", "NonMargin_Stroma", "NonMargin_Other", b) %>%
  dplyr::rename("ID" = "tcga_id") %>%
  reshape2::melt(c("ID", "Margin_Blood", "Margin_Normal", "Margin_Stroma", "Margin_Other", "NonMargin_Blood", "NonMargin_Normal", "NonMargin_Stroma", "NonMargin_Other")) %>%
  arrange(desc(value))
# dplyr::mutate(Margin = gsub("_[[:print:]]*", "", variable),
#               Texture = gsub("[[:print:]]*_", "", variable),
#               Margin = factor(Margin),
#               Texture = factor(Texture),
#               variable = factor(variable))


# Adjust p values
pairwise.test = tcga_kirc_long %>% 
  group_by(variable) %>%
  wilcox_test(Margin_Blood~value) %>% 
  adjust_pvalue(method = 'BH') %>%
  # mutate(p.adj = round(p.adj, 2))
  mutate(p.adj = ifelse(p.adj < 0.001, "***",
                        ifelse(p.adj < 0.01, "**",
                               ifelse(p.adj < 0.05, "*", round(p.adj, 2)))))


# Plot
png("/Users/oscarbruck/OneDrive - University of Helsinki/RCC/Otso/Analysis/Textures/Images/Margin/TCGA_margin_Ly_scatterplot.png", width = 6, height = 6, units = 'in', res = 300, pointsize = 12) #original pointsize = 12
tcga_kirc_long %>% 
  dplyr::filter(variable == "stage") %>%
  ggplot(aes(x = value, y = Margin_Blood)) +
  geom_point(size=5, aes(fill=value), shape = 21, color = "black", position=position_jitterdodge()) +
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

# Join margin and non-margin dataframes
tcga_kirc_wide <- tcga_kirc_margin %>%
  full_join(tcga_kirc_nonmargin) %>%
  mutate(FC = Margin/NonMargin,
         FC = ifelse(FC == "Inf", 10,
                     ifelse(FC > 10, 10,
                            ifelse(FC < 0.1, 0.1, FC))))

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



# Loop over all textures
## Save data
tcga_kirc0 <- tcga_kirc

for (texture1 in c("All", as.list(unique(tcga_kirc0$Texture))) ) {
  print(texture1)
  
  # Reset data
  if (texture1 == "All") {
    tcga_kirc <- tcga_kirc0 %>%
      dplyr::filter(Texture %in% unique(tcga_kirc0$Texture))
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
  ## Rownames to column
  pvalue_df <- pvalue_df %>%
    rownames_to_column() %>%
    rename(genes = rowname,
           pvalue = p.value
    ) %>%
    arrange(pvalue)
  rownames(pvalue_df) = NULL
  
  
  
  # Export data
  writexl::write_xlsx(pvalue_df, paste0("/Users/oscarbruck/OneDrive - University of Helsinki/RCC/Otso/Analysis/Textures/Images/Margin/Margin_vs_non_margin/Mut/Ly_in_margin_vs_nonmargin_", texture1, "_mut.xlsx"))
  
  
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
    
    
    
    
    if (length(unique(tcga_kirc1$two1)) == 2) {
      
      # Plot
      g <- ggplot(tcga_kirc1, aes(x = two1, y = FC)) +
        geom_jitter(size=5, width = 0.2, aes(fill=two1), shape = 21, color = "black") +
        geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
        labs(y=paste0("Lymphocyte proportion fold change (", texture1, " margin vs. non-margin)"), x=toupper(gsub("_mutation", "", two1))) +
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
      ggsave(plot = g, filename = paste0("/Users/oscarbruck/OneDrive - University of Helsinki/RCC/Otso/Analysis/Textures/Images/Margin/Margin_vs_non_margin/Clinical/Ly_in_margin_vs_nonmargin_", two1, "_", texture1, "_clin.png"), width = 5, height = 5, units = 'in', dpi = 300)
      
    } else if (length(unique(tcga_kirc1$two1)) == 3) {
      
      # Plot
      g <- ggplot(tcga_kirc1, aes(x = two1, y = FC)) +
        geom_jitter(size=5, width = 0.2, aes(fill=two1), shape = 21, color = "black") +
        geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
        labs(y=paste0("Lymphocyte proportion fold change (", texture1, " margin vs. non-margin)"), x=toupper(gsub("_mutation", "", two1))) +
        theme_bw() +
        theme(axis.text.x = element_text(size=12, colour = "black"),
              axis.text.y = element_text(size=12, colour = "black"),
              axis.title.y = element_text(size=10, face="bold", colour = "black"),
              axis.title.x = element_text(size=14, face="bold", colour = "black"),
              legend.position = "none") +
        scale_fill_brewer(palette = c("Set1")) +
        stat_compare_means(method = "kruskal.test",
                           label.y = 1.25*a,
                           size = 6) +
        stat_pvalue_manual(
          pairwise.test,
          label = "p.adj",
          bracket.size = 1.5,
          size = 5,
          # y.position = c(1.12*a, 1.24*a, 1.3*a, 1.0*a, 1.18*a, 1.06*a)
          y.position = c(1.0*a, 1.12*a, 1.06*a)
        )
      ggsave(plot = g, filename = paste0("/Users/oscarbruck/OneDrive - University of Helsinki/RCC/Otso/Analysis/Textures/Images/Margin/Margin_vs_non_margin/Clinical/Ly_in_margin_vs_nonmargin_", two1, "_", texture1, "_clin.png"), width = 5, height = 5, units = 'in', dpi = 300)
      
    } else if (length(unique(tcga_kirc1$two1)) == 3) {
      
      # Plot
      g <- ggplot(tcga_kirc1, aes(x = two1, y = FC)) +
        geom_jitter(size=5, width = 0.2, aes(fill=two1), shape = 21, color = "black") +
        geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
        labs(y=paste0("Lymphocyte proportion fold change (", texture1, " margin vs. non-margin)"), x=toupper(gsub("_mutation", "", two1))) +
        theme_bw() +
        theme(axis.text.x = element_text(size=12, colour = "black"),
              axis.text.y = element_text(size=12, colour = "black"),
              axis.title.y = element_text(size=10, face="bold", colour = "black"),
              axis.title.x = element_text(size=14, face="bold", colour = "black"),
              legend.position = "none") +
        scale_fill_brewer(palette = c("Set1")) +
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
      ggsave(plot = g, filename = paste0("/Users/oscarbruck/OneDrive - University of Helsinki/RCC/Otso/Analysis/Textures/Images/Margin/Margin_vs_non_margin/Clinical/Ly_in_margin_vs_nonmargin_", two1, "_", texture1, "_clin.png"), width = 5, height = 5, units = 'in', dpi = 300)
      
      
    }
    
  }
  
}
