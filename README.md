# Welcome to the data analysis repository of the "Integrative Analysis of Tissue Textures and Lymphocyte Infiltration in Renal Cell Carcinoma using Deep Learning" by Brummer Otso et al.


## BACKGROUND
These codes will help you reproduce all plots and statistical analyses of the publication.  


## USEFUL LINKS
**CODE**$~~~$The [codes](https://version.helsinki.fi/hus_hematology/tcga-kirc-immunology) to reproduce the texture and lymphocyte data.  
**DATA**$~~~$The annotated image data are located in [Zenodo](https://zenodo.org/deposit/6384627).  
**VISUALIZATION**$~~~$The [TissUUmaps platform](http://hruh-20.it.helsinki.fi/rcc_texture_lymphocytes/) to visualize the texture and lymphocyte data.



## INSTRUCTIONS

### Instructions if you operate with RStudio
1. Install RStudio (the analyses here have been made with the version 3.5.1)
2. Install necessary R packages either by running the Rscript `./src/install_packages.R`
3. Run the analyses at once by running the Rscript `./src/main.R`. This will
- preprocess data files found in `./data/`
- produce tables and images in `./src/path/Images/`


### Instructions if you operate with the terminal
1. Install RStudio (the analyses here have been made with the version 3.5.1)
2. `cd ./RCC_textures_and_lymphocytes_publication_data_analysis`
3. `Rscript src/install_packages.R`
4. Run the analyses at once by running `Rscript src/main.R`. This will
- preprocess data files found in `./data/`
- produce tables and images in `./src/path/Images/`


A detailed description of  the R and library versions can be found in `session_info.yaml`.
