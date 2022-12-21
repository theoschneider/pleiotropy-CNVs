# Pleiotropic mechanisms of copy number variants

This script (process.R) uses both Mendelian randomization and CNV-GWAS files to determine the effects :
- of the CNV on the different traits,
- of the traits on every other trait.  
  
The selected_phenotypes.tsv file contains the phenotypes of interest for region 16p11.2 BP4-5 that were used for this study as an example, but the code is easily adaptable to other regions.  
  
This script is still a work in progress, as well as this README.

## Usage

This program requires the installation of a few packages, which are all listed in the script.  
It also requires a list of IDs as an input, as well as summary stats for the traits of interest.  
The output is a dataframe with rows being very pair of traits, and columns being the different effects, standard errors and p-values.  


## Acknowledgments

Thank you to all the members of the Statistical Genetics Group for their help and advice. 
