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

# Preprocess data if not done earlier
source("./Preprocessing/Preprocess_rcc_pts.R")


# Run all the rest
rscripts = list.files(path = ".", pattern = ".R$", recursive = TRUE, full.names = TRUE)
rscripts = rscripts[grep(x = rscripts, pattern = "Preprocess", invert = TRUE)]
for (rscript in rscripts) {
  source(rscript)
}
