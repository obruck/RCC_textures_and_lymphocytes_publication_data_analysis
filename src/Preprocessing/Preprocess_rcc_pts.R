rm(list=ls())

# Libraries
library(tidyverse)
library(readxl)
library(writexl)


# Change working directory
## Get current file location
getCurrentFileLocation <-  function()
{
  this_file <- commandArgs() %>% 
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
  if (length(this_file)==0)
  {
    this_file <- rstudioapi::getSourceEditorContext()$path
  }
  return(dirname(this_file))
}
setwd(getCurrentFileLocation())


# Load data
## Images with correct resolution
tcga_kirc <- read_xlsx("../data/image_analysis_results_before_sample_exclusion.xlsx")


## Mutation data
tcga_kirc1 <- read_xlsx("../data/mutations_before_sample_exclusion.xlsx") %>%
  dplyr::rename(tcga_id = bcr_patient_barcode) %>%
  dplyr::filter(!is.na(tcga_id)) %>%
  distinct(tcga_id, .keep_all=TRUE)  # Remove duplicate patients


## Patients to remove (≤5% cancer, wrong histology and poor quality)
to_remove <- read_xlsx("../data/samples_to_exclude.xlsx")


# Filter
## Raw data + MMC2 + to_remove
tcga_kirc <- tcga_kirc %>%
  dplyr::filter(!tcga_id %in% to_remove$TCGAid)
tcga_kirc_mut <- tcga_kirc %>%
  dplyr::inner_join(tcga_kirc1)


# Save
write_xlsx(tcga_kirc, "../data/image_analysis_results_final.xlsx")
write_xlsx(tcga_kirc_mut, "../data/mutations_final.xlsx")
