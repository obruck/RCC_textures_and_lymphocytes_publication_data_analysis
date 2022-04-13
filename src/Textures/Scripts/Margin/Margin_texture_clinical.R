rm(list=ls())

print("Start Textures/Scripts/Margin/Margin_texture_clinical.R")

# Load libraries
library(readxl)
library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(RColorBrewer)
library(fastDummies)
library(survival)
library(survminer)


############################# LOAD DATA ##########################################################################################################


# Image data
tcga_kirc <- read_xlsx("../data/image_analysis_results_final.xlsx")

# Normalize by removing empty
tcga_kirc$`texture_cancer_%` <- 100*tcga_kirc$texture_cancer / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`texture_normal_%` <- 100*tcga_kirc$texture_normal / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`margin_texture_normal_%` <- 100*tcga_kirc$margin_texture_normal / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)
tcga_kirc$`margin_texture_blood_%` <- 100*tcga_kirc$margin_texture_blood / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)
tcga_kirc$`margin_texture_cancer_%` <- 100*tcga_kirc$margin_texture_cancer / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)
tcga_kirc$`margin_texture_stroma_%` <- 100*tcga_kirc$margin_texture_stroma / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)
tcga_kirc$`margin_texture_other_%` <- 100*tcga_kirc$margin_texture_other / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)
tcga_kirc$`non_margin_texture_normal_%` <- 100*tcga_kirc$non_margin_texture_normal / (tcga_kirc$non_margin_texture_blood + tcga_kirc$non_margin_texture_normal + tcga_kirc$non_margin_texture_stroma + tcga_kirc$non_margin_texture_other)
tcga_kirc$`non_margin_texture_blood_%` <- 100*tcga_kirc$non_margin_texture_blood / (tcga_kirc$non_margin_texture_blood + tcga_kirc$non_margin_texture_normal + tcga_kirc$non_margin_texture_stroma + tcga_kirc$non_margin_texture_other)
tcga_kirc$`non_margin_texture_cancer_%` <- 100*tcga_kirc$non_margin_texture_cancer / (tcga_kirc$non_margin_texture_blood + tcga_kirc$non_margin_texture_normal + tcga_kirc$non_margin_texture_stroma + tcga_kirc$non_margin_texture_other)
tcga_kirc$`non_margin_texture_normal_%` <- 100*tcga_kirc$non_margin_texture_normal / (tcga_kirc$non_margin_texture_blood + tcga_kirc$non_margin_texture_normal + tcga_kirc$non_margin_texture_stroma + tcga_kirc$non_margin_texture_other)
tcga_kirc$`non_margin_texture_stroma_%` <- 100*tcga_kirc$non_margin_texture_stroma / (tcga_kirc$non_margin_texture_blood + tcga_kirc$non_margin_texture_normal + tcga_kirc$non_margin_texture_stroma + tcga_kirc$non_margin_texture_other)
tcga_kirc$`non_margin_texture_other_%` <- 100*tcga_kirc$non_margin_texture_other / (tcga_kirc$non_margin_texture_blood + tcga_kirc$non_margin_texture_normal + tcga_kirc$non_margin_texture_stroma + tcga_kirc$non_margin_texture_other)

tcga_kirc <- tcga_kirc %>%
  dplyr::filter(`texture_cancer_%` > 5) %>%
  dplyr::mutate(tissue_source_site = gsub("-[[:print:]]{4}", "", gsub("TCGA-", "", tcga_id)))

# Save
tcga_kirc0 <- tcga_kirc


############################# MODIFY CLINICAL DATA ##########################################################################################################


# Rename
colnames(tcga_kirc)[grep(pattern = "[[:print:]]{0,4}margin_texture_[[:print:]]*%", colnames(tcga_kirc))] <- c("Margin_Blood", "NonMargin_Blood", "Margin_Cancer", "NonMargin_Cancer", "Margin_Normal", "NonMargin_Normal", "Margin_Stroma", "NonMargin_Stroma", "Margin_Other", "NonMargin_Other")
margin <- c("Margin_Blood", "Margin_Normal", "Margin_Stroma", "Margin_Other")
nonmargin <- c("NonMargin_Blood", "NonMargin_Normal", "NonMargin_Stroma", "NonMargin_Other")
textures <- c(margin, nonmargin)


# Petri's data
df <- readRDS("../data/clinical_transcriptome.rds")
# Filter patients
df <- df %>%
  dplyr::select(one_of(tcga_kirc$tcga_id))

# Clean clinical data
df1 <- df[grep(pattern = "CLIN:", x = rownames(df)),]
df1 <- df1 %>% rownames_to_column()
df1 <- df1[apply(df1[,-1], 1, function(x) !all(x==0)),]
df1 <- df1[apply(df1[,-1], 1, function(x) !all(x==1)),]
# df1 <- df1 %>% filter_all(any_vars(!is.na(.)))
df1 <- df1 %>% dplyr::filter(!is.na(rowname))
df2 <- df1 %>% 
  pivot_longer(!rowname, names_to = "col1", values_to = "col2") %>% 
  pivot_wider(names_from = "rowname", values_from = "col2") %>%
  dplyr::rename(tcga_id = col1)
## Join
tcga_kirc <- tcga_kirc %>%
  dplyr::left_join(df2)

# Clean sample data
df1 <- df[grep(pattern = "N:SAMP:", x = rownames(df)),]
df1 <- df1 %>% rownames_to_column()
df1 <- df1[apply(df1[,-1], 1, function(x) !all(x==0)),]
df1 <- df1[apply(df1[,-1], 1, function(x) !all(x==1)),]
# df1 <- df1 %>% filter_all(any_vars(!is.na(.)))
df1 <- df1 %>% dplyr::filter(!is.na(rowname))
df2 <- df1 %>% 
  pivot_longer(!rowname, names_to = "col1", values_to = "col2") %>% 
  pivot_wider(names_from = "rowname", values_from = "col2") %>%
  dplyr::rename(tcga_id = col1)
## Convert numeric data to binary data
df2[2:ncol(df2)] <- sapply(df2[2:ncol(df2)], function(x) ifelse(x>median(x, na.rm=TRUE), 1, 0))
colnames(df2) <- gsub("N:SAMP:", "", colnames(df2))
## Join
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

## CLIN:age_at_initial_pathologic_diagnosis
cols1 = grep(pattern = "CLIN:age_at_initial_pathologic_diagnosis", colnames(tcga_kirc))
mc = max.col(tcga_kirc[cols1], ties.method = "first")
tcga_kirc$Age <- tcga_kirc[[cols1]]
tcga_kirc$Age <- ifelse(tcga_kirc$Age > median(tcga_kirc$Age, na.rm=TRUE), "High", "Low")

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
twos = c("gender", "vital_status", "LDH", "ECOG", "Karnofsky", "Age", colnames(df2)[2:ncol(df2)])
threes = c("Lymph_node_metastasis", "Organ_metastasis", "Organ_metastasis_", "WBC", "therapy_outcome", "smoking_status", "HB", "calcium", "PLT")
fours = c("Tumor_size", "race", "stage", "grade")


# Convert variables to binary format
## Variables to hotencode
to_hotencode <- tcga_kirc[c(twos, threes, fours)]
### k-cluster make one hot encoder
to_hotencode_var <- dummy_cols(tcga_kirc[colnames(to_hotencode)]) %>%
  dplyr::select(-ends_with("_NA"))
### Remove duplicate binary variables
to_hotencode_var <- to_hotencode_var %>% dplyr::select(-c(Age_Low, gender_FEMALE, vital_status_Alive, LDH_Normal, ECOG_0, `Karnofsky_≤80`))
### Join
to_hotencode_var <- cbind(tcga_kirc %>% dplyr::select(tcga_id) %>% dplyr::rename(ID = tcga_id), to_hotencode_var)


############################# PLOT ##########################################################################################################


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

## Join
tcga_kirc <- tcga_kirc_wide %>% inner_join(to_hotencode_var)
tcga_kirc <- tcga_kirc %>% dplyr::left_join(tcga_kirc0 %>% dplyr::select(ID = tcga_id, `texture_normal_%`))



# Arrange columns
tcga_kirc <- tcga_kirc %>%
  dplyr::relocate("neoantigen_num", .after = "TGF-beta_Response") %>%
  dplyr::relocate("indel_num", "immunogenic_indel_num", "numberOfImmunogenicMutation", "Genome_doublings", .after = "Indel_Neoantigens") %>%
  dplyr::rename(Neoantigens = neoantigen_num,
                ImmunogenicMutations = numberOfImmunogenicMutation,
                Indels = indel_num,
                Immunogenic_indels = immunogenic_indel_num
  )


# Loop over all textures
## Save data
tcga_kirc0 <- tcga_kirc

for (texture1 in unique(tcga_kirc0$Texture) ) {
  # texture1 <- "Stroma"
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
  
  # Find column number
  findcolnumber <- function(df, thecolumnname){
    which(colnames(df) == deparse(substitute(thecolumnname)))
  }
  
  # Find variables with 2 values excluding NA
  unique_values_2 = sapply(tcga_kirc[(findcolnumber(tcga_kirc, Leukocyte_Fraction)):(ncol(tcga_kirc))], function(x) length(unique(x[!is.na(x)])))
  unique_values_2 = names(unique_values_2[unique_values_2==2])
  ## Remove character variables
  nums <- unlist(lapply(tcga_kirc[unique_values_2], is.numeric))
  unique_values_2 <- unique_values_2[nums]
  ## Remove uninteresting variables
  unique_values_2 <- unique_values_2[-which(unique_values_2 %in% c("Wound_Healing", "Lymphocyte_Infiltration_Signature_Score", "Macrophage_Regulation", "Lymphocytes", "Number_of_Segments", "Fraction_Altered",
                                                                   "BCR_Evenness", "BCR_Shannon", "TCR_Evenness", "TCR_Shannon", "CTA_Score", "Aneuploidy_Score",
                                                                   "OS_Time", "PFI_Time", "Macrophages_M0", "Mast_Cells", "Mast_Cells_Activated", "Mast_Cells_Resting",
                                                                   "T_Cells_Follicular_Helper", "T_Cells_CD4_Memory_Resting", "B_Cells_Naive",
                                                                   "Dendritic_Cells_Resting", "Dendritic_Cells_Activated", "Neutrophils_1", "Eosinophils_1",
                                                                   "HBV", "HCV", "HHV", "HIV", "HPV", "HTLV", "MCV",
                                                                   "therapy_outcome_Complete Remission/Response", "therapy_outcome_Stable Disease", "therapy_outcome_Progressive Disease",
                                                                   "purity", "ploidy","Coverage_for_80_power",
                                                                   "Cancer_DNA_fraction", "Subclonal_genome_fraction", "numberOfNonSynonymousSNP",
                                                                   "numberOfPeptideTested", "numberOfBindingPMHC", "numberOfBindingExpressedPMHC",
                                                                   "saimiriine_herpesvirus", "AS", "ASprime", "LOH_n_seg", "LOH_frac_altered", "Nrf2_Score"))]
  
  
  # Multiple comparison (two-sided, unpaired, Mann-Whitney U test)
  # multiple_t_tests_p_value <- lapply(tcga_kirc[(findcolnumber(tcga_kirc, Age)+1):(ncol(tcga_kirc)-1)],
  multiple_t_tests_p_value <- lapply(tcga_kirc[unique_values_2],
                                     function(x) wilcox.test(tcga_kirc$FC ~ x, correct=FALSE, exact=FALSE, estimate=TRUE, na.rm=TRUE))
  ## P-values can be extracted from the result object
  pvalue <- data.frame(p.value = sapply(multiple_t_tests_p_value, getElement, name = "p.value"))
  ## Create a matrix and dataframe of the p-values
  pvalue_mat <- as.matrix(pvalue)
  pvalue_df <- data.frame(pvalue)
  ## Add the p values to a new dataframe
  pvalue_df$p_adjusted <- p.adjust(p = pvalue_mat, method = "BH")
  ## Add also the t values, 95%CI to the same dataframe
  pvalue_df$median1 <- unlist(lapply(tcga_kirc[unique_values_2], function(x) median(tcga_kirc[x==1,]$FC, na.rm=TRUE)))
  pvalue_df$median2 <- unlist(lapply(tcga_kirc[unique_values_2], function(x) median(tcga_kirc[x==0,]$FC, na.rm=TRUE)))
  ## Rownames to column
  pvalue_df <- pvalue_df %>%
    rownames_to_column() %>%
    rename(genes = rowname,
           pvalue = p.value
    ) %>%
    arrange(pvalue)
  rownames(pvalue_df) = NULL
  
  
  
  # Export data
  writexl::write_xlsx(pvalue_df, paste0("Textures/Images/Margin/Margin_vs_non_margin/Texture/Clinical/Ly_in_margin_vs_nonmargin_", texture1, "_clinical.xlsx"))
  
  # Filter variables with p_adjusted1 < 0.10
  vars <- pvalue_df %>%
    dplyr::filter(pvalue < 0.10) %>%
    dplyr::select(genes) %>%
    dplyr::mutate(genes = gsub("_[[:print:]]*", "", genes)) %>%
    dplyr::distinct()
  # tmp <- data.frame(twos); tmp <- rbind(tmp, data.frame(threes) %>% rename(twos = threes)); tmp <- rbind(tmp, data.frame(fours) %>% rename(twos = fours))
  tmp = data.frame(twos = unique_values_2)
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
  
  
  # Plot parameters
  a <- max(tcga_kirc$FC, na.rm=TRUE)
  
  
  for (two1 in b) {
    # two1 = b[1]
    
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
      ggsave(plot = g, filename = paste0("Textures/Images/Margin/Margin_vs_non_margin/Texture/Clinical/Ly_in_margin_vs_nonmargin_", two1, "_", texture1, "_clin.png"), width = 5, height = 5, units = 'in', dpi = 300)
      
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
                           label.y = 1.30*a,
                           size = 6) +
        stat_pvalue_manual(
          pairwise.test,
          label = "p",
          bracket.size = 1.5,
          size = 5,
          # y.position = c(1.12*a, 1.24*a, 1.3*a, 1.0*a, 1.18*a, 1.06*a)
          y.position = c(1.0*a, 1.15*a, 1.075*a)
        )
      ggsave(plot = g, filename = paste0("Textures/Images/Margin/Margin_vs_non_margin/Texture/Clinical/Ly_in_margin_vs_nonmargin_", two1, "_", texture1, "_clin.png"), width = 5, height = 5, units = 'in', dpi = 300)
      
    } else if (length(unique(tcga_kirc1$two1)) == 4) {
      
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
                           label.y = 1.50*a,
                           size = 6) +
        stat_pvalue_manual(
          pairwise.test,
          label = "p",
          bracket.size = 1.5,
          size = 5,
          # y.position = c(1.12*a, 1.24*a, 1.3*a, 1.0*a, 1.18*a, 1.06*a)
          y.position = c(1.16*a, 1.32*a, 1.40*a, 1.0*a, 1.24*a, 1.08*a)
        )
      ggsave(plot = g, filename = paste0("Textures/Images/Margin/Margin_vs_non_margin/Texture/Clinical/Ly_in_margin_vs_nonmargin_", two1, "_", texture1, "_clin.png"), width = 5, height = 5, units = 'in', dpi = 300)
      
    }
    
  }
  
}


################################################### SURVIVAL ################################################################################


# Add OS time
df1 <- df[grep(pattern = "N:SAMP:", x = rownames(df)),]
df1 <- df1 %>% rownames_to_column()
df1 <- df1[apply(df1[,-1], 1, function(x) !all(x==0)),]
df1 <- df1[apply(df1[,-1], 1, function(x) !all(x==1)),]
# df1 <- df1 %>% filter_all(any_vars(!is.na(.)))
df1 <- df1 %>% dplyr::filter(!is.na(rowname))
df2 <- df1 %>% 
  pivot_longer(!rowname, names_to = "col1", values_to = "col2") %>% 
  pivot_wider(names_from = "rowname", values_from = "col2") %>%
  dplyr::rename(tcga_id = col1)
## Convert numeric data to binary data
colnames(df2) <- gsub("N:SAMP:", "", colnames(df2))
## Join
tcga_kirc1 <- tcga_kirc0 %>%
  dplyr::select(ID, Margin, NonMargin, Texture, FC, vital_status_Dead, `texture_normal_%`) %>%
  dplyr::left_join(df2 %>% dplyr::select(ID = tcga_id, OS_Time))


# Loop
for (texture1 in unique(tcga_kirc1$Texture) ) {
  print(texture1)
  
  # Reset data
  if (texture1 == "Normal") {
    
    tcga_kirc2 <- tcga_kirc1 %>%
      dplyr::filter(`texture_normal_%` > 1) %>%
      dplyr::filter(Texture %in% texture1)
    
  } else {
    tcga_kirc2 <- tcga_kirc1 %>%
      dplyr::filter(Texture %in% texture1)
  }
  
  # Categorize FC
  tcga_kirc2 <- tcga_kirc2 %>%
    dplyr::mutate(
      FC = ntile(FC, 3),
      Dead = ifelse(is.na(OS_Time), NA,
                                 ifelse(OS_Time > 5 * 365.25, 0, vital_status_Dead)),
      OS_Time = ifelse(is.na(OS_Time), NA,
                       ifelse(OS_Time > 5 * 365.25, 5 / 30.43, OS_Time/30.43))
    )
  
  fit1 <- survfit(Surv(tcga_kirc2$OS_Time, tcga_kirc2$Dead)~tcga_kirc2$FC, data = tcga_kirc2)
  
  
  #ggsurvplot with 3 groups
  g <- ggsurvplot(fit1,
             data = tcga_kirc2,
             # fun = "event", #cumulative incidence
             palette = "npg",      # Use JCO journal color palette
             size = 2,   #line thickness
             ggtheme = theme_minimal(), #theme
             font.main = c(20, "bold", "black"), #title font
             font.x = c(20, "bold", "black"), #x-axis font
             font.y = c(20, "bold", "black"), #y-axis font
             font.tickslab = c(20, "bold", "black"), #axis numbering font
             conf.int = FALSE, #confidence interval
             pval = TRUE, #p-value
             pval.size = 8, #p-value size
             # pval.coord = c(1500, 0.125), #p-value coordinate ## 1500 for MR4, 1250 for MMR and 500 for CCyR
             risk.table.pos = "out",
             tables.y.text = FALSE,
             tables.theme = theme_cleantable(),
             risk.table.fontsize = 6,
             risk.table = c("absolute"), #risk table ('absolute', 'percentage', 'abs_pct')
             risk.table.col = "strata", #risk table color
             risk.table.height = 0.25,
             fontsize = 4,
             xlab = "Time (months)",
             ylab = "Survival rate (%)",
             break.time.by = 6, # break time axis by 250, ##1000 for MR4, 500 for MMR, 250 for CCyR
             # legend = c(0.8, 0.175),
             censor = TRUE,
             legend.title = paste0("Margin/NonMargin FC"),
             font.legend = c(12, "bold", "black"),    #font voi olla esim. "bold" tai "plain"
             legend.labs = c("Low", "Intermediate", "High"),
             size = 1)  # change line size
  ggsave(plot = print(g), filename = paste0("Textures/Images/Margin/Margin_vs_non_margin/Texture/Clinical/KM_Ly_in_margin_vs_nonmargin_FC_", texture1, "_.png"), width = 6, height = 6, units = 'in', dpi = 300)
  
}





#################
## Balloonplot ##
################# 

# Reset data
tcga_kirc0 <- tcga_kirc

# Read statistics
if (exists("pvalue_df0")) {rm(pvalue_df0)}
for (texture1 in c("Blood", "Stroma", "Normal", "Other")) {
  pvalue_df <- read_xlsx(paste0("Textures/Images/Margin/Margin_vs_non_margin/Texture/Clinical/Ly_in_margin_vs_nonmargin_", texture1, "_clinical.xlsx"))
  pvalue_df$Texture = texture1
  if (exists("pvalue_df0")) {
    pvalue_df0 <- rbind(pvalue_df0, pvalue_df)
  } else {
    pvalue_df0 <- pvalue_df
  }
}


# Find variables with 2 values excluding NA
unique_values_0 = unique_values_2[40:85]
unique_values_1 = unique_values_2[1:39]


# Function to capiralize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


# Loop balloonplots
for (unique_value_name in 1:length(list(unique_values_0, unique_values_1))) {
  
  # Filter data
  unique_value <- list(unique_values_0, unique_values_1)[[unique_value_name]]
  
  # Reset data
  pvalue_df1 <- pvalue_df0 %>%
    # Filter
    dplyr::filter(genes %in% unique_value) %>%
    # Order
    dplyr::mutate(genes = factor(genes, levels = unique_value)) %>%
    dplyr::arrange(desc(genes)) %>%
    # Filter
    dplyr::filter(!genes %in% c("Eosinophils_1", "Neutrophils_1")) %>%
    dplyr::filter(!str_detect(genes, "__")) %>%
    dplyr::filter(!str_detect(genes, "Time")) %>%
    dplyr::filter(!str_detect(genes, "therapy_outcome")) %>%
    # Rename
    dplyr::mutate(genes = gsub("herpes", "Herpes", genes),
                  genes = gsub("stage_", "", genes),
                  genes = gsub("_", " ", genes))
  
  
  # Calculate FC and -log10 P value
  pvalue_df1 <- pvalue_df1 %>%
    mutate(FC = median1 / median2) %>%
    mutate(p_adjusted = p.adjust(p = pvalue, method = "BH")) %>%
    mutate(
      neg_log10_p = -log10(pvalue),
      neg_log10_p_adj = -log10(p_adjusted),
      # neg_log10_p_adj = ifelse(neg_log10_p_adj>10, 10, neg_log10_p_adj),
      p_adjusted_cat = ifelse(p_adjusted<0.001, 0.001,
                              ifelse(p_adjusted<0.01, 0.01,
                                     ifelse(p_adjusted<0.05, 0.05,
                                            ifelse(p_adjusted<0.1, 0.1, "ns")))),
      pvalue_cat = ifelse(pvalue<0.001, 0.001,
                          ifelse(pvalue<0.01, 0.01,
                                 ifelse(pvalue<0.05, 0.05,
                                        ifelse(pvalue<0.1, 0.1, "ns")))),
      FC_log = log(FC),
      Variables = gsub("_", " ", genes),
      Variables = firstup(Variables)
    )
  
  # Balloonplot
  p <- ggballoonplot(pvalue_df1, x = "Texture", y = "Variables",
                     fill = "FC_log",
                     size = "neg_log10_p",
                     # size.range = c(1, 10),
                     ggtheme = theme_bw()) +
    scale_y_discrete(limits = unique(pvalue_df1$Variables)) +
    # scale_size(breaks = c(0, 1, 2, 3), range = c(1, 16), limits = c(1, 16)) +
    # scale_size(breaks = c(exp(0.05), 2, 3), labels = c(0.05,0.01,0.001), range = c(1, 10), limits = c(exp(0.05), max(-log10(0.001), max(pvalue_df1$neg_log10_p_adj, na.rm = TRUE)))) +
    scale_size(breaks = c(exp(0.05), 2, 3), labels = c(0.05,0.01,0.001), range = c(1, 10), limits = c(exp(0.05), max(-log10(0.001), max(pvalue_df1$neg_log10_p, na.rm = TRUE)))) +
    # scale_size(breaks = c(0, 1, 2, 3), range = c(1, 10), limits = c(1, max(pvalue_df1$neg_log10_p_adj, na.rm = TRUE))) +
    # ceiling(max(pvalue_df1$neg_log10_p))
    # scale_fill_viridis_c(option = "C") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    # gradient_fill(c("blue", "white", "red")) +
    guides(size = guide_legend(title="P-value", nrow = 3, title.vjust = 0.5),
           fill = guide_colorbar(title="LOG10(FC)", title.vjust = 0.75)) +
    font("xy.text", size = 10, color = "black", face="plain") +
    theme(axis.title.y = element_text(size=12, colour="black", face="bold", angle = 90),
          axis.text.x = element_text(colour="black", face="bold", angle = 45, hjust = 1.0, vjust = 1.0),
          axis.title.x = element_text(size=12, colour="black", face="bold"),
          axis.text.y = element_text(colour="black"),
          legend.position = "right",
          # legend.position = c(1.3, 1.3),
          # legend.justification = ,
          legend.direction = "vertical",
          # legend.box="vertical",
          legend.margin=margin(),
          legend.text = element_text(size=12, colour="black", face="bold"),
          legend.title = element_text(size=12, colour="black", face="bold"))
  p
  
  ggsave(plot = p, filename = paste0("Textures/Images/Margin/Margin_balloonplot/Balloonplot_Texture_", unique_value_name, ".png"), width = 5, height = nrow(pvalue_df1)/25, units = 'in', dpi = 300)
  
}

