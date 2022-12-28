# Pleiotropic mechanisms of copy number variants

This script (pleiotropy-CNVs.R) uses both Mendelian randomization and CNV-GWAS files to determine the effects:
- of the CNV on the different traits,
- of the traits on every other trait.  
  

## Usage

This program requires the installation of a few packages, which are all listed in the script.  

It also requires two files as an input:
- a list of IDs with 2 columns : "ID" is a list of IDs and "trait" is the corresponding phenotype in the probes files. IDs must be coming from the [IEU GWAS database](https://gwas.mrcieu.ac.uk).
- a file containing the top probe for each trait in the region of interest. Should have at least a "trait" column, with traits corresponding to the ones in the IDs file above, and a "effect_allele" column with the direction.  

Please note that a smaller scale example of these files was provided, with only 3 traits, to try the script.  
The output is a dataframe with rows being every pair of traits, and columns being the different effects, standard errors and p-values.  
The script also asks for a treshold and a MTC, to create a dataframe with significant comparisons, and one with insignificant comparisons.
After setting all the variables, the entire script can be run at once. 


## Acknowledgments

Thank you to all the members of the Statistical Genetics Group for their help and advice. 
