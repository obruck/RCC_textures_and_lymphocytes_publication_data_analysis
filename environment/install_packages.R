# Load package tidyverse
if (!require("tidyverse", quietly = TRUE))
  install.packages("tidyverse")


# Get current file location
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




# Wrapper to install packages from CRAN
pkgLoad_cran <- function( packages ) {
  
  packagecheck <- match( packages, utils::installed.packages()[,1] )
  
  packagestoinstall <- packages[ is.na( packagecheck ) ]
  
  if( length( packagestoinstall ) > 0L ) {
    utils::install.packages( packagestoinstall )
  } else {
    print( "All requested packages already installed" )
  }
  
  for( package in packages ) {
    suppressPackageStartupMessages(
      library( package, character.only = TRUE, quietly = TRUE )
    )
  }
  
}

# Install CRAN packages
pkgLoad_cran(packages = as.vector(unlist(read.table(file.path(getCurrentFileLocation(), "requirements_cran.txt")))))





# Wrapper to install packages from Bioconductor
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
pkgLoad_bio <- function( packages ) {
  
  packagecheck <- match( packages, utils::installed.packages()[,1] )
  
  packagestoinstall <- packages[ is.na( packagecheck ) ]
  
  if( length( packagestoinstall ) > 0L ) {
    BiocManager::install( packagestoinstall )
  } else {
    print( "All requested Bioconductor packages already installed" )
  }
  
  for( package in packages ) {
    suppressPackageStartupMessages(
      library( package, character.only = TRUE, quietly = TRUE )
    )
  }
  
}


# Install Bioconductor packages
pkgLoad_bio(packages = c("ComplexHeatmap", "fgsea"))

