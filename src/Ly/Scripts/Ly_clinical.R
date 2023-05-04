rm(list=ls())

print("Start Ly/Scripts/Ly_clinical.R")

# Load libraries
library(readxl)
library(tidyverse)
library(data.table, warn.conflicts = TRUE)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(RColorBrewer)
library(fastDummies)
library(survival)
library(survminer)
library(sjPlot)

# Set theme for barplots
set_theme(theme_bw(), axis.title.color = "black", axis.title.size = 1.1, axis.textcolor = "black", axis.textsize = 1.05, geom.label.size = 3.5, legend.pos="bottom")


############################# LOAD DATA ##########################################################################################################


# Image data
tcga_kirc <- read_xlsx("../data/image_analysis_results_final.xlsx")

# Normalize by removing empty
tcga_kirc$`texture_cancer_%` <- 100*tcga_kirc$texture_cancer / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`texture_normal_%` <- 100*tcga_kirc$texture_normal / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)
tcga_kirc$`margin_texture_normal_%` <- 100*tcga_kirc$margin_texture_normal / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_cancer + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)
tcga_kirc$`margin_texture_blood_%` <- 100*tcga_kirc$margin_texture_blood / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_cancer + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)
tcga_kirc$`margin_texture_cancer_%` <- 100*tcga_kirc$margin_texture_cancer / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_cancer + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)
tcga_kirc$`margin_texture_stroma_%` <- 100*tcga_kirc$margin_texture_stroma / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_cancer + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)
tcga_kirc$`margin_texture_other_%` <- 100*tcga_kirc$margin_texture_other / (tcga_kirc$margin_texture_blood + tcga_kirc$margin_texture_cancer + tcga_kirc$margin_texture_normal + tcga_kirc$margin_texture_stroma + tcga_kirc$margin_texture_other)
tcga_kirc$`non_margin_texture_normal_%` <- 100*tcga_kirc$non_margin_texture_normal / (tcga_kirc$non_margin_texture_blood + tcga_kirc$non_margin_texture_cancer + tcga_kirc$non_margin_texture_normal + tcga_kirc$non_margin_texture_stroma + tcga_kirc$non_margin_texture_other)
tcga_kirc$`non_margin_texture_blood_%` <- 100*tcga_kirc$non_margin_texture_blood / (tcga_kirc$non_margin_texture_blood + tcga_kirc$non_margin_texture_cancer + tcga_kirc$non_margin_texture_normal + tcga_kirc$non_margin_texture_stroma + tcga_kirc$non_margin_texture_other)
tcga_kirc$`non_margin_texture_cancer_%` <- 100*tcga_kirc$non_margin_texture_cancer / (tcga_kirc$non_margin_texture_blood + tcga_kirc$non_margin_texture_cancer + tcga_kirc$non_margin_texture_normal + tcga_kirc$non_margin_texture_stroma + tcga_kirc$non_margin_texture_other)
tcga_kirc$`non_margin_texture_normal_%` <- 100*tcga_kirc$non_margin_texture_normal / (tcga_kirc$non_margin_texture_blood + tcga_kirc$non_margin_texture_cancer + tcga_kirc$non_margin_texture_normal + tcga_kirc$non_margin_texture_stroma + tcga_kirc$non_margin_texture_other)
tcga_kirc$`non_margin_texture_stroma_%` <- 100*tcga_kirc$non_margin_texture_stroma / (tcga_kirc$non_margin_texture_blood + tcga_kirc$non_margin_texture_cancer + tcga_kirc$non_margin_texture_normal + tcga_kirc$non_margin_texture_stroma + tcga_kirc$non_margin_texture_other)
tcga_kirc$`non_margin_texture_other_%` <- 100*tcga_kirc$non_margin_texture_other / (tcga_kirc$non_margin_texture_blood + tcga_kirc$non_margin_texture_cancer + tcga_kirc$non_margin_texture_normal + tcga_kirc$non_margin_texture_stroma + tcga_kirc$non_margin_texture_other)
tcga_kirc$inf_bin_lymphocytes_total <- 100*(tcga_kirc$bin_lymphocytes_blood + tcga_kirc$bin_lymphocytes_normal + tcga_kirc$bin_lymphocytes_cancer + tcga_kirc$bin_lymphocytes_stroma + tcga_kirc$bin_lymphocytes_other) / (tcga_kirc$texture_blood + tcga_kirc$texture_cancer + tcga_kirc$texture_normal + tcga_kirc$texture_stroma + tcga_kirc$texture_other)

# Filter samples with >5% cancer normal tissue
tcga_kirc <- tcga_kirc %>%
  dplyr::filter(`texture_cancer_%` > 5) %>%
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


# Save
tcga_kirc0 <- tcga_kirc


############################# MODIFY CLINICAL DATA ##########################################################################################################


# Rename
colnames(tcga_kirc)[grep(pattern = "^inf_bin_lymphocytes_[[:print:]]*", colnames(tcga_kirc))] <- c("Blood", "Cancer", "Normal", "Stroma", "Other", "Total")
# textures <- c("Blood", "Cancer", "Stroma", "Other", "Total")
textures <- c("Blood", "Cancer", "Normal", "Stroma", "Other", "Total")

# Petri's data
df <- readRDS("../data/clinical_transcriptome.rds") %>%
  # Filter patients
  dplyr::select(one_of(tcga_kirc$tcga_id))

# Clean clinical data
df1 <- df[grep(pattern = "CLIN:", x = rownames(df)),]
df1 <- df1 %>% rownames_to_column()
df1 <- df1[apply(df1[,-1], 1, function(x) !all(x==0)),]
df1 <- df1[apply(df1[,-1], 1, function(x) !all(x==1)),]
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
df1 <- df1 %>% dplyr::filter(!is.na(rowname))
df2 <- df1 %>% 
  pivot_longer(!rowname, names_to = "col1", values_to = "col2") %>% 
  pivot_wider(names_from = "rowname", values_from = "col2") %>%
  dplyr::rename(tcga_id = col1)
## Convert numeric data to binary data
# df2[2:ncol(df2)] <- sapply(df2[2:ncol(df2)], function(x) ifelse(x>median(x, na.rm=TRUE), 1, 0))
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
tcga_kirc$Lymph_node_metastasis <- ifelse(tcga_kirc$Lymph_node_metastasis == "NX", NA, tcga_kirc$Lymph_node_metastasis)

## CLIN:pathologic_M
cols1 = grep(pattern = "CLIN:pathologic_M", colnames(tcga_kirc))
mc = max.col(tcga_kirc[cols1], ties.method = "first")
tcga_kirc$Organ_metastasis <- gsub("B:CLIN:pathologic_M_", "", colnames(tcga_kirc[cols1])[mc])
tcga_kirc$Organ_metastasis <- ifelse(tcga_kirc$Organ_metastasis == "MX", NA, tcga_kirc$Organ_metastasis)

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
tcga_kirc <- tcga_kirc %>% dplyr::mutate(WBC = factor(WBC, levels = c("Low", "Normal", "Elevated", NA)))

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
tcga_kirc$therapy_outcome <- ifelse(tcga_kirc$therapy_outcome == "Stable Disease", NA, tcga_kirc$therapy_outcome)
tcga_kirc <- tcga_kirc %>% mutate(therapy_outcome = factor(therapy_outcome, levels = c("Complete Remission/Response", "Progressive Disease", NA)))

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
tcga_kirc <- tcga_kirc %>% dplyr::mutate(Karnofsky = factor(Karnofsky, levels = c("≤80", ">80")))

## CLIN:neoplasm_histologic_grade
cols1 = grep(pattern = "CLIN:neoplasm_histologic_grade", colnames(tcga_kirc))
mc = max.col(tcga_kirc[cols1], ties.method = "first")
tcga_kirc$grade <- gsub("B:CLIN:neoplasm_histologic_grade_", "", colnames(tcga_kirc[cols1])[mc])
tcga_kirc$grade <- ifelse(tcga_kirc$grade == "GX", NA, tcga_kirc$grade)
tcga_kirc <- tcga_kirc %>% dplyr::mutate(grade = factor(grade, levels = c("G1", "G2", "G3", "G4")))

## CLIN:tobacco_smoking_history
cols1 = grep(pattern = "CLIN:tobacco_smoking_history", colnames(tcga_kirc))
mc = max.col(tcga_kirc[cols1], ties.method = "first")
tcga_kirc$smoking_status <- gsub("B:CLIN:tobacco_smoking_history_", "", colnames(tcga_kirc[cols1])[mc])
tcga_kirc$smoking_status <- ifelse(tcga_kirc$smoking_status == "[Unknown]", NA, tcga_kirc$smoking_status)
tcga_kirc$smoking_status <- ifelse(str_detect(tcga_kirc$smoking_status, "eformed"), "Ex-smoker", tcga_kirc$smoking_status)
tcga_kirc$smoking_status <- gsub("Lifelong ", "", tcga_kirc$smoking_status)
tcga_kirc <- tcga_kirc %>% dplyr::mutate(smoking_status = factor(smoking_status, levels = c("Non-smoker", "Current smoker", "Ex-smoker", NA)))

## CLIN:lactate_dehydrogenase
cols1 = grep(pattern = "CLIN:lactate_dehydrogenase", colnames(tcga_kirc))
mc = max.col(tcga_kirc[cols1], ties.method = "first")
tcga_kirc$LDH <- gsub("B:CLIN:lactate_dehydrogenase_result_", "", colnames(tcga_kirc[cols1])[mc])
tcga_kirc$LDH <- ifelse(tcga_kirc$LDH == "[Unknown]", NA, tcga_kirc$LDH)
tcga_kirc <- tcga_kirc %>% dplyr::mutate(LDH = factor(LDH, levels = c("Normal", "Elevated", NA)))

## CLIN:hemoglobin_result
cols1 = grep(pattern = "CLIN:hemoglobin_result", colnames(tcga_kirc))
mc = max.col(tcga_kirc[cols1], ties.method = "first")
tcga_kirc$HB <- gsub("B:CLIN:hemoglobin_result_", "", colnames(tcga_kirc[cols1])[mc])
tcga_kirc$HB <- ifelse(tcga_kirc$HB == "[Unknown]", NA, tcga_kirc$HB)
tcga_kirc <- tcga_kirc %>% dplyr::mutate(HB = factor(HB, levels = c("Low", "Normal", "Elevated", NA)))

## CLIN:serum_calcium_result
cols1 = grep(pattern = "CLIN:serum_calcium_result", colnames(tcga_kirc))
mc = max.col(tcga_kirc[cols1], ties.method = "first")
tcga_kirc$calcium <- gsub("B:CLIN:serum_calcium_result_", "", colnames(tcga_kirc[cols1])[mc])
tcga_kirc$calcium <- ifelse(tcga_kirc$calcium == "[Unknown]", NA, tcga_kirc$calcium)
tcga_kirc <- tcga_kirc %>% dplyr::mutate(calcium = factor(calcium, levels = c("Low", "Normal", "Elevated", NA)))

## CLIN:platelet_qualitative
cols1 = grep(pattern = "CLIN:platelet_qualitative", colnames(tcga_kirc))
mc = max.col(tcga_kirc[cols1], ties.method = "first")
tcga_kirc$PLT <- gsub("B:CLIN:platelet_qualitative_result_", "", colnames(tcga_kirc[cols1])[mc])
tcga_kirc$PLT <- ifelse(tcga_kirc$PLT == "[Unknown]", NA, tcga_kirc$PLT)
tcga_kirc <- tcga_kirc %>% dplyr::mutate(PLT = factor(PLT, levels = c("Low", "Normal", "Elevated", NA)))


# Summary
sapply(tcga_kirc[, c(which( colnames(tcga_kirc)=="Tumor_size" ) : which( colnames(tcga_kirc)=="PLT" ) )], function(x) length(unique(x)))
sapply(tcga_kirc[, c(which( colnames(tcga_kirc)=="Tumor_size" ) : which( colnames(tcga_kirc)=="PLT" ) )], function(x) unique(x))
twos = c("gender", "LDH", "ECOG", "Karnofsky", "Age", "Lymph_node_metastasis", "Organ_metastasis",
         "Age_High", "gender_MALE", "therapy_outcome_Progressive Disease", "ECOG_1", "Lymph_node_metastasis_N1", "Organ_metastasis_M1", "Karnofsky_>80",
         colnames(df2)[2:ncol(df2)])
threes = c("WBC", "smoking_status", "HB", "calcium", "PLT")  # "Organ_metastasis_"
fours = c("Tumor_size", "stage", "grade") #"race"


# Convert variables to binary format
## Variables to hotencode
to_hotencode <- tcga_kirc[c("gender", "LDH", "ECOG", "Karnofsky", "Age", "Lymph_node_metastasis", "Organ_metastasis", "therapy_outcome", threes, fours)]
### k-cluster make one hot encoder
to_hotencode_var <- dummy_cols(tcga_kirc[colnames(to_hotencode)]) %>%
  dplyr::select(-ends_with("_NA"))
# to_hotencode_var <- to_hotencode_var %>% dplyr::select(-(1:ncol(to_hotencode)))
### Replace NA with 0
# for(i in 1:ncol(to_hotencode_var)){
#   to_hotencode_var[is.na(to_hotencode_var[,i]), i] = 0
# }
### Remove duplicate binary variables
to_hotencode_var <- to_hotencode_var %>%
  dplyr::select(-c(Age_Low, gender_FEMALE, `therapy_outcome_Progressive Disease`, Lymph_node_metastasis_N0, Organ_metastasis_M0, LDH_Normal, ECOG_0, `Karnofsky_≤80`))
to_hotencode_var <- to_hotencode_var %>%
  dplyr::select(-one_of(colnames(tcga_kirc)))
to_hotencode_var <- cbind(tcga_kirc %>%
                            dplyr::select(ID = tcga_id), to_hotencode_var)
### Join
b <- ncol(tcga_kirc)
# to_hotencode_var <- cbind(tcga_kirc %>% dplyr::select(tcga_id) %>% dplyr::rename(ID = tcga_id), to_hotencode_var)
tcga_kirc <- tcga_kirc %>%
  dplyr::rename(ID = tcga_id) %>%
  dplyr::left_join(to_hotencode_var)


############################# PLOT ##########################################################################################################


# To percent
tcga_kirc[textures] <- sapply(tcga_kirc[textures], function(x) 100*x)

# Arrange columns
tcga_kirc <- tcga_kirc %>%
  dplyr::relocate("neoantigen_num", .after = "TGF-beta_Response") %>%
  dplyr::relocate("indel_num", "immunogenic_indel_num", "numberOfImmunogenicMutation", "Genome_doublings", .after = "Indel_Neoantigens") %>%
  dplyr::rename(Neoantigens = neoantigen_num,
                Immunogenic_Mutations = numberOfImmunogenicMutation,
                Indels = indel_num,
                Immunogenic_Indels = immunogenic_indel_num
  )


# Loop over all textures
## Save data
tcga_kirc0 <- tcga_kirc


## Immunological and texture data (Wilcoxon test)
for (texture1 in textures) {
  # texture1 <- textures[2]
  print(texture1)
  
  # Reset data
  if (texture1 == "Normal") {
    
    tcga_kirc <- tcga_kirc0 %>%
      dplyr::filter(`texture_normal_%` > 1)
    
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
  
  # Find column number
  findcolnumber <- function(df, thecolumnname){
    which(colnames(df) == deparse(substitute(thecolumnname)))
  }
  
  
  # Identify immunological and genomic variables
  unique_values_2_imm = sapply(tcga_kirc[(findcolnumber(tcga_kirc, Leukocyte_Fraction)):(findcolnumber(tcga_kirc, Nrf2_Score))], function(x) length(unique(x[!is.na(x)])))
  unique_values_2_imm = names(unique_values_2_imm)
  ## Remove character variables
  nums <- unlist(lapply(tcga_kirc[unique_values_2_imm], is.numeric))
  unique_values_2_imm <- unique_values_2_imm[nums]
  ## Remove uninteresting variables
  unique_values_2_imm <- unique_values_2_imm[-which(unique_values_2_imm %in% c("Wound_Healing", "Lymphocyte_Infiltration_Signature_Score", "Macrophage_Regulation", "Lymphocytes", "Number_of_Segments", "Fraction_Altered",
                                                                               "BCR_Evenness", "BCR_Shannon", "TCR_Evenness", "TCR_Shannon", "CTA_Score", "Aneuploidy_Score",
                                                                               "SNV_Neoantigens", "Indel_Neoantigens", "Immunogenic_Indels", "Immunogenic_Mutations", "Silent_Mutation_Rate","Nonsilent_Mutation_Rate","Homologous_Recombination_Defects",
                                                                               "OS_Time", "PFI_Time", "Macrophages_M0", "Mast_Cells", "Mast_Cells_Activated", "Mast_Cells_Resting",
                                                                               "T_Cells_Follicular_Helper", "T_Cells_CD4_Memory_Resting", "B_Cells_Naive",
                                                                               "NK_Cells_Resting", "Dendritic_Cells_Resting", "Dendritic_Cells_Activated", "Neutrophils_1", "Eosinophils_1",
                                                                               "CMV", "EBV", "HBV", "HCV", "HHV", "HIV", "HPV", "HTLV", "MCV",
                                                                               "therapy_outcome_Complete Remission/Response", "therapy_outcome_Stable Disease", "therapy_outcome_Progressive Disease",
                                                                               "purity", "ploidy","Coverage_for_80_power", "Genome_doublings",
                                                                               "Cancer_DNA_fraction", "Subclonal_genome_fraction", "numberOfNonSynonymousSNP",
                                                                               "numberOfPeptideTested", "numberOfBindingPMHC", "numberOfBindingExpressedPMHC",
                                                                               "saimiriine_herpesvirus", "AS", "ASprime", "LOH_n_seg", "LOH_frac_altered", "Nrf2_Score"))]
  
  
  # Multiple comparison (two-sided, unpaired, Mann-Whitney U test)
  # multiple_t_tests_p_value <- lapply(tcga_kirc[(findcolnumber(tcga_kirc, Age)+1):(ncol(tcga_kirc)-1)],
  multiple_t_tests_p_value <- lapply(tcga_kirc[unique_values_2_imm],
                                     function(x) wilcox.test(x ~ tcga_kirc$FC, correct=FALSE, exact=FALSE, estimate=TRUE, na.rm=TRUE))
  ## P-values can be extracted from the result object
  pvalue <- data.frame(p.value = sapply(multiple_t_tests_p_value, getElement, name = "p.value"))
  ## Create a matrix and dataframe of the p-values
  pvalue_mat <- as.matrix(pvalue)
  pvalue_df <- data.frame(pvalue)
  ## Add the p values to a new dataframe
  pvalue_df$p_adjusted <- p.adjust(p = pvalue_mat, method = "BH")
  ## Rownames to column
  pvalue_df <- pvalue_df %>%
    rownames_to_column() %>%
    rename(genes = rowname,
           pvalue = p.value
    )
  rownames(pvalue_df) = NULL
  
  # Median
  tcga_kirc <- as.data.frame(tcga_kirc)
  pvalue_df$median1 <- NA
  pvalue_df$median2 <- NA
  for (i in 1:length(unique_values_2_imm)) {
    pvalue_df[i,]$median1 <- median(tcga_kirc[tcga_kirc$FC=="High",][unique_values_2_imm][[i]], na.rm=TRUE)
    pvalue_df[i,]$median2 <- median(tcga_kirc[tcga_kirc$FC=="Low",][unique_values_2_imm][[i]], na.rm=TRUE)
  }
  
  # Arrange
  pvalue_df <- pvalue_df %>%
    arrange(pvalue)
  
  # Export data
  dir.create(paste0("Ly/Images/Clinical/", texture1))
  writexl::write_xlsx(pvalue_df, paste0("Ly/Images/Clinical/", texture1, "/Ly_", texture1, "_immunological.xlsx"))
  
  # Filter variables with p_adjusted1 < 0.10
  vars <- pvalue_df %>%
    dplyr::filter(pvalue < 0.10) %>%
    dplyr::select(genes) %>%
    dplyr::mutate(genes = gsub("_[[:print:]]*", "", genes)) %>%
    dplyr::distinct()
  
  tmp = data.frame(twos = unique_values_2_imm)
  
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
  ## Remove NA
  a <- a %>% dplyr::filter(!is.na(twos))
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
      dplyr::mutate(FC = as.factor(FC))
    
    # Adjust p values
    pairwise.test = tcga_kirc1 %>%
      wilcox_test(two1~FC) %>%
      adjust_pvalue(method = 'BH') %>%
      # mutate(p.adj = round(p.adj, 2))
      mutate(p = ifelse(p < 0.001, "***",
                        ifelse(p < 0.01, "**",
                               ifelse(p < 0.05, "*", round(p, 2)))))
    
    
    # Plot
    g <- ggplot(tcga_kirc1, aes(x = FC, y = two1)) +
      geom_jitter(size=5, width = 0.2, aes(fill=FC), shape = 21, color = "black") +
      geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
      labs(x=paste0("Lymphocyte proportion in ", texture1), y=toupper(gsub("_mutation", "", two1))) +
      theme_bw() +
      theme(axis.text.x = element_text(size=12, colour = "black"),
            axis.text.y = element_text(size=12, colour = "black"),
            axis.title.y = element_text(size=14, face="bold", colour = "black"),
            axis.title.x = element_text(size=14, face="bold", colour = "black"),
            legend.position = "none") +
      scale_fill_brewer(palette = c("Set1")) +
      stat_pvalue_manual(
        pairwise.test,
        label = "p",
        bracket.size = 1.5,
        size = 5,
        # y.position = c(1.12*a, 1.24*a, 1.3*a, 1.0*a, 1.18*a, 1.06*a)
        y.position = 1.05*max(tcga_kirc1$two1, na.rm = TRUE)
      )
    ggsave(plot = g, filename = paste0("Ly/Images/Clinical/", texture1, "/Ly_", gsub(">", "", two1), "_", texture1, "_clin.png"), width = 5, height = 5, units = 'in', dpi = 300)
    
  }
  
}


## Clinical data (Chisq test)
for (texture1 in textures) {
  # texture1 <- textures[2]
  print(texture1)
  
  # Reset data
  if (texture1 == "Normal") {
    
    tcga_kirc <- tcga_kirc0 %>%
      dplyr::filter(`texture_normal_%` > 1)
    
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
  
  # Find column number
  findcolnumber <- function(df, thecolumnname){
    which(colnames(df) == deparse(substitute(thecolumnname)))
  }
  
  
  # Find variables with 2 values excluding NA
  unique_values_2_cli = sapply(tcga_kirc[(findcolnumber(tcga_kirc, Leukocyte_Fraction)):(ncol(tcga_kirc))], function(x) length(unique(x[!is.na(x)])))
  unique_values_2_cli = names(unique_values_2_cli[unique_values_2_cli==2])
  ## Remove character variables
  nums <- unlist(lapply(tcga_kirc[unique_values_2_cli], is.numeric))
  unique_values_2_cli <- unique_values_2_cli[nums]
  ## Remove uninteresting variables
  unique_values_2_cli <- unique_values_2_cli[-which(unique_values_2_cli %in% c("Wound_Healing", "Lymphocyte_Infiltration_Signature_Score", "Macrophage_Regulation", "Lymphocytes", "Number_of_Segments", "Fraction_Altered",
                                                                               "BCR_Evenness", "BCR_Shannon", "TCR_Evenness", "TCR_Shannon", "CTA_Score", "Aneuploidy_Score",
                                                                               "SNV_Neoantigens", "Indel_Neoantigens", "Immunogenic_Indels", "Immunogenic_Mutations", "Silent_Mutation_Rate","Nonsilent_Mutation_Rate","Homologous_Recombination_Defects",
                                                                               "OS_Time", "PFI_Time", "Macrophages_M0", "Mast_Cells", "Mast_Cells_Activated", "Mast_Cells_Resting",
                                                                               "T_Cells_Follicular_Helper", "T_Cells_CD4_Memory_Resting", "B_Cells_Naive",
                                                                               "NK_Cells_Resting", "Dendritic_Cells_Resting", "Dendritic_Cells_Activated", "Neutrophils_1", "Eosinophils_1",
                                                                               "CMV", "EBV", "HBV", "HCV", "HHV", "HIV", "HPV", "HTLV", "MCV",
                                                                               "therapy_outcome_Complete Remission/Response", "therapy_outcome_Stable Disease", "therapy_outcome_Progressive Disease",
                                                                               "purity", "ploidy","Coverage_for_80_power",
                                                                               "Cancer_DNA_fraction", "Subclonal_genome_fraction", "numberOfNonSynonymousSNP",
                                                                               "numberOfPeptideTested", "numberOfBindingPMHC", "numberOfBindingExpressedPMHC",
                                                                               "saimiriine_herpesvirus", "AS", "ASprime", "LOH_n_seg", "LOH_frac_altered", "Nrf2_Score"))]
  
  
  # Multiple comparison (chisq test)
  # multiple_t_tests_p_value <- lapply(tcga_kirc[(findcolnumber(tcga_kirc, Age)+1):(ncol(tcga_kirc)-1)],
  multiple_t_tests_p_value <- lapply(tcga_kirc[unique_values_2_cli],
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
  dir.create(paste0("Ly/Images/Clinical/", texture1))
  writexl::write_xlsx(pvalue_df, paste0("Ly/Images/Clinical/", texture1, "/Ly_", texture1, "_clinical.xlsx"))
  
  
  # Filter variables with p_adjusted1 < 0.10
  vars <- pvalue_df %>%
    dplyr::filter(pvalue < 0.10) %>%
    dplyr::select(genes) %>%
    dplyr::mutate(genes = gsub("_[[:print:]]*", "", genes)) %>%
    dplyr::distinct()
  
  tmp = data.frame(twos = unique_values_2_cli)
  
  if (nrow(vars) > 0) {
    
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
    ## Remove NA
    a <- a %>%
      dplyr::filter(!is.na(twos))
    b <- as.character(a$twos)
    
    
    # Plot
    for (two1 in b) {
      # two1 = b[1]
      
      # Reset data
      tcga_kirc1 <- tcga_kirc
      
      # Assign value
      tcga_kirc1$two1 <- ifelse(tcga_kirc1[[two1]] == 1, TRUE,
                                ifelse(tcga_kirc1[[two1]] == 0, FALSE, NA))
      
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
      
      ggsave(plot = g, filename = paste0("Ly/Images/Clinical/", texture1, "/Ly_", gsub(">", "", two1), "_", texture1, "_clin_barplot.png"), width = 4, height = 4.5, units = 'in', dpi = 300)
      
    }
    
  }
  
}


###################################### SURVIVAL #################################################################################################


# Add OS time
df1 <- df[grep(pattern = "N:SAMP:", x = rownames(df)),]
df1 <- df1 %>% rownames_to_column()
df1 <- df1[apply(df1[,-1], 1, function(x) !all(x==0)),]
df1 <- df1[apply(df1[,-1], 1, function(x) !all(x==1)),]
df1 <- df1 %>% dplyr::filter(!is.na(rowname))
df2 <- df1 %>% 
  pivot_longer(!rowname, names_to = "col1", values_to = "col2") %>% 
  pivot_wider(names_from = "rowname", values_from = "col2") %>%
  dplyr::rename(tcga_id = col1)
## Convert numeric data to binary data
colnames(df2) <- gsub("N:SAMP:", "", colnames(df2))
## Join
tcga_kirc1 <- tcga_kirc0 %>%
  dplyr::select(ID, vital_status, `texture_normal_%`, !!textures, tissue_source_site) %>%
  dplyr::left_join(df2 %>% dplyr::select(ID = tcga_id, OS_Time)) %>%
  dplyr::rename(vital_status_Dead = vital_status) %>%
  dplyr::mutate(vital_status_Dead = ifelse(vital_status_Dead == "Dead", 1, 0))


# Loop
for (texture1 in textures) {
  print(texture1)
  
  # Reset data
  if (texture1 == "Normal") {
    
    tcga_kirc2 <- tcga_kirc1 %>%
      dplyr::filter(`texture_normal_%` > 1)
    
  } else {
    
    tcga_kirc2 <- tcga_kirc1
    
  }
  
  # Divide patients into low and high by clinical centers
  tcga_kirc2$FC <- tcga_kirc2[[texture1]]
  tcga_kirc2 <- tcga_kirc2 %>%
    group_by(tissue_source_site) %>%
    mutate(FC = (FC > median(FC, na.rm = TRUE)) + 1) %>%
    mutate(FC = ifelse(FC == 2, "High", "Low")) %>%
    ungroup()
  
  
  # Categorize FC
  tcga_kirc2 <- tcga_kirc2 %>%
    dplyr::mutate(
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
                  palette = "Set1",      # Use JCO journal color palette
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
                  risk.table.fontsize = 5,
                  risk.table = c("absolute"), #risk table ('absolute', 'percentage', 'abs_pct')
                  risk.table.col = "strata", #risk table color
                  risk.table.height = 0.25,
                  fontsize = 5,
                  xlab = "Time (months)",
                  ylab = "Survival rate (%)",
                  break.time.by = 12, # break time axis by 250, ##1000 for MR4, 500 for MMR, 250 for CCyR
                  # legend = c(0.8, 0.175),
                  censor = TRUE,
                  legend.title = paste0("Lymphocyte % in ", tolower(texture1)),
                  font.legend = c(12, "bold", "black"),    #font voi olla esim. "bold" tai "plain"
                  legend.labs = c("High", "Low"),
                  size = 1)  # change line size
  ggsave(plot = print(g), filename = paste0("Ly/Images/Clinical/KM/KM_Ly_", texture1, ".png"), width = 6, height = 6, units = 'in', dpi = 300)
  
}



#################
## Balloonplot ##
################# 


# Reset data
tcga_kirc0 <- tcga_kirc

# Read statistics
if (exists("pvalue_df0")) {rm(pvalue_df0)}
for (texture1 in textures) {
  pvalue_df <- read_xlsx(paste0("Ly/Images/Clinical/", texture1, "/Ly_", texture1, "_immunological.xlsx"))
  pvalue_df$Texture = texture1
  if (exists("pvalue_df0")) {
    pvalue_df0 <- rbind(pvalue_df0, pvalue_df)
  } else {
    pvalue_df0 <- pvalue_df
  }
}

# Read statistics
for (texture1 in textures) {
  pvalue_df <- read_xlsx(paste0("Ly/Images/Clinical/", texture1, "/Ly_", texture1, "_clinical.xlsx"))
  pvalue_df$Texture = texture1
  if (exists("pvalue_df0")) {
    pvalue_df0 <- rbind(pvalue_df0, pvalue_df)
  } else {
    pvalue_df0 <- pvalue_df
  }
}

# Find variables with 2 values excluding NA
unique_values_0 = unique_values_2_cli
unique_values_1 = unique_values_2_imm


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
    dplyr::filter(!str_detect(genes, "_NX")) %>%
    dplyr::filter(!str_detect(genes, "_MX")) %>%
    dplyr::filter(!str_detect(genes, "calcium_Low")) %>%
    dplyr::filter(!str_detect(genes, "HB_Elevated")) %>%
    dplyr::filter(!str_detect(genes, "Time")) %>%
    dplyr::filter(!str_detect(genes, "therapy_outcome")) %>%
    # Rename
    dplyr::mutate(genes = gsub("herpes", "Herpes", genes),
                  genes = gsub("stage_", "", genes),
                  genes = gsub("_", " ", genes),
                  genes = gsub("gamma", "Gamma", genes),
                  genes = gsub("delta", "Delta", genes),
                  genes = gsub("Genome d", "Genome D", genes))
  
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
    scale_size(breaks = c(exp(0.05), 2, 3), labels = c(0.05,0.01,0.001), range = c(1, 10), limits = c(exp(0.05), max(-log10(0.001), max(pvalue_df1$neg_log10_p_adj, na.rm = TRUE)))) +
    # scale_size(breaks = c(exp(0.1), 2, 3), labels = c(0.1,0.01,0.001), range = c(1, 10), limits = c(exp(0.1), max(-log10(0.001), max(pvalue_df1$neg_log10_p, na.rm = TRUE)))) +
    # scale_size(breaks = c(exp(0.05), 2, 3), labels = c(0.05,0.01,0.001), range = c(1, 10), limits = c(0.1, max(-log10(0.001), max(pvalue_df1$neg_log10_p, na.rm = TRUE)))) +
    # scale_size(breaks = c(0, 1, 2, 3), range = c(1, 10), limits = c(1, max(pvalue_df1$neg_log10_p_adj, na.rm = TRUE))) +
    # ceiling(max(pvalue_df1$neg_log10_p))
    # scale_fill_viridis_c(option = "C") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    # gradient_fill(c("blue", "white", "red")) +
    guides(size = guide_legend(title="P-value", nrow = 3, title.vjust = 0.5),
           fill = guide_colorbar(title="LOG10(FC)", title.vjust = 0.75)) +
    font("xy.text", size = 10, color = "black", face="plain") +
    # theme(axis.title.y = element_text(size=12, colour="black", face="bold", angle = 90),
    #       axis.title.x = element_text(size=12, colour="black", face="bold"),
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(colour="black", face="bold", angle = 45, hjust = 1.0, vjust = 1.0),
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
  
  ggsave(plot = p, filename = paste0("Ly/Images/Clinical/Balloonplot/Balloonplot_Ly_", unique_value_name, ".png"), width = 5, height = nrow(pvalue_df1)/25, units = 'in', dpi = 300)
  
}

