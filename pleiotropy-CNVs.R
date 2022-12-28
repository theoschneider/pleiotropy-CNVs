#### 1. Packages ####

# remotes::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
# remotes::install_github("MRCIEU/MRInstruments")
library(MRInstruments)
library(data.table) # For fread and fwrite.
library(ggplot2) # To save MR plots.





#### 2. Variables ####

input_directory <- "~/"
output_directory <- "~/"

IDs_filename <- "ids.txt" # Name of the file containing 2 columns : "ID" and "trait".
# "ID" is a list of IDs and "trait" is the corresponding phenotype in the probes files (see below).
# IDs must be coming from the IEU GWAS database. Can be found at https://gwas.mrcieu.ac.uk. 
# File must be in the input directory.

probes_filename <- "probes.tsv" # Name of the file containing the top probe for each trait in the region of interest. 
# Should have a "trait" column, with traits corresponding to the ones in the IDs file above. 
# Should have a "effect_allele" column with the direction.
# For example, effect_allele can be "additional copy", "duplication (> 2 copies)", or "deletion (< 2 copies)".
# File must be in the input directory. 

treshold <- 0.05 # Can be assigned to define the treshold of the different p-values.
MTC <- 3 # Minimal number of significant tests for multiple testing correction, maximum is 5 (set to 0 or 1 for no MTC).
direction <- "additional copy" # Can be other directions, has to be the same as the column in probes.

save.MR.plots <- FALSE # Save MR plots? False by default. 
# The package can create 4 plots per MR, in a new folder (can take considerable disk space).





#### 3. Index ####

# Read file and save IDs
index <- fread(file = paste0(input_directory, "/", IDs_filename))
input_IDs <- index$ID

# Read available outcomes from TwoSampleMR package
ao <- available_outcomes()

# Add the corresponding db_trait to index
for(i in 1:nrow(index)){
    index$db_trait[i] <- ao$trait[ao$id == index$ID[i]]
}





#### 4. Combinations ####

# Create all possible pairs of traits, in both directions 
combinations <- expand.grid(exposure = input_IDs, outcome = input_IDs, stringsAsFactors = FALSE)

# Remove the pairs of identical traits
combinations <- combinations[combinations$exposure != combinations$outcome,]

# Recreate the indexes for clarity
rownames(combinations) <- 1:nrow(combinations)





#### 5. Two Sample MR ####

# Create empty dataframe
MR_results <- data.frame()

# Create folder for plots
if(save.MR.plots){
    dir.create(paste0(output_directory, "/Plots"))
}

# For every pair of traits (this loop might take a considerable amount of time to run)
for(i in 1:nrow(combinations)){
    
    # Retrieve the ID
    exposure_id = combinations$exposure[i]
    outcome_id = combinations$outcome[i]
    
    # Extract instruments and outcome data
    exp_dat <- extract_instruments(outcomes = exposure_id,
                                   clump = TRUE,
                                   p1 = 5e-8,
                                   r2 = 0.001)
    out_dat <- extract_outcome_data(snps = exp_dat$SNP,
                                    outcomes = outcome_id)
    
    # Harmonise data
    dat <- harmonise_data(exp_dat, out_dat)
    
    # Run MR and add it to the dataframe
    dat_mr <- mr(dat)
    MR_results <- rbind(MR_results, dat_mr)
    
    # Create and save plots
    if(save.MR.plots){
        
        p1 <- mr_scatter_plot(dat_mr, dat)
        ggsave(plot = p1[[1]], filename = paste0(exposure_id, "_against_", outcome_id, "_01_scatter_plot.png"), path = paste0(output_directory, "/Plots"))
        
        res_single <- mr_singlesnp(dat)
        p2 <- mr_forest_plot(res_single)
        ggsave(plot = p2[[1]], filename = paste0(exposure_id, "_against_", outcome_id, "_02_forest_plot.png"), path = paste0(output_directory, "/Plots"))
        
        res_loo <- mr_leaveoneout(dat)
        p3 <- mr_leaveoneout_plot(res_loo)
        ggsave(plot = p3[[1]], filename = paste0(exposure_id, "_against_", outcome_id, "_03_loo_plot.png"), path = paste0(output_directory, "/Plots"))
        
        res_single <- mr_singlesnp(dat)
        p4 <- mr_funnel_plot(res_single)
        ggsave(plot = p4[[1]], filename = paste0(exposure_id, "_against_", outcome_id, "_04_funnel_plot.png"), path = paste0(output_directory, "/Plots"))
    }
}

# Save MRs in a file
fwrite(MR_results, file = paste0(output_directory, "/MR_results.tsv"), quote = FALSE, sep = "\t")

# Create empty dataframe for methods
MR_methods <- data.frame()

# Create an unique ID for each pair of traits
MR_results$collapse <- paste0(MR_results$id.exposure, MR_results$id.outcome)

# For each unique pair
for(pair in unique(MR_results$collapse)){
    
    # Retrieve the five methods used on this pair
    five_methods <- MR_results[MR_results$collapse == pair,]
    
    # Keep only the method with the best p-value
    line <- five_methods[five_methods$pval == min(five_methods$pval),][1,]
    
    # Add a column with the number of significant methods (with corrected p-value)
    line$number.methods <- sum(five_methods$pval <= (0.05/nrow(combinations)))
    
    # Add the line to the methods dataframe
    MR_methods <- rbind(MR_methods, line)
}

# Save the methods in a file
fwrite(MR_methods, file = paste0(output_directory, "/MR_methods.tsv"), quote = FALSE, sep = "\t")





#### 6. Probes ####

# Import probes and keep 1 per trait
probes <- as.data.frame(fread(file = paste0(input_directory, "/", probes_filename), sep = "\t"))

# Keep only the direction of interest
probes <- probes[probes$effect_allele == direction,]

# Keep one probe per trait
probes <- probes[duplicated(probes$trait) == F,]

# Re-create the index for clarity
rownames(probes) <- 1:nrow(probes)





#### 7. Comparisons ####

# Create all possible pairs of traits, in both directions
comparisons <- expand.grid(T1.id = input_IDs, T2.id = input_IDs, stringsAsFactors = FALSE)

# Remove the pairs of identical traits
comparisons <- comparisons[comparisons$T1.id != comparisons$T2.id,]

# Recreate the indexes for clarity
rownames(comparisons) <- 1:nrow(comparisons)

# For every pair of traits
for(i in 1:nrow(comparisons)){
    
    # take CNV on T1 from probes
    # take CNV on T2 from probes
    # take T1 on T2 from MR_results -> IVW method
    
    # Extract the IDs
    T1 <- comparisons$T1.id[i]
    T2 <- comparisons$T2.id[i]
    
    # Extract the effect of the CNV region on trait 1 (from the probes)
    CNVonT1 <- probes$beta[probes$trait == index$trait[index$ID == T1]]
    
    # Extract the effect of the CNV region on trait 1 (from the probes)
    CNVonT2 <- probes$beta[probes$trait == index$trait[index$ID == T2]]
    
    # Extract the effect of trait 1 on trait 2 (from the MRs)
    T1onT2 <- MR_results$b[MR_results$id.exposure == T1 & MR_results$id.outcome == T2 & MR_results$method == "Inverse variance weighted"]
    
    # Do the same for standard errors
    seCNVonT1 <- probes$standard_error[probes$trait == index$probes_trait[index$db_id == T1]]
    seCNVonT2 <- probes$standard_error[probes$trait == index$probes_trait[index$db_id == T2]]
    seT1onT2 <- MR_results$se[MR_results$id.exposure == T1 & MR_results$id.outcome == T2 & MR_results$method == "Inverse variance weighted"]
    
    # Calculate the effect of the CNV region on trait 2 by passing through trait 1
    CNVonT1onT2 <- CNVonT1 * T1onT2
    
    # Calculate the corresponding standard error
    new.se <- seCNVonT1*seT1onT2 + seCNVonT1*T1onT2^2 + seT1onT2*CNVonT1^2
    
    # Calculate T stat and its p-value
    Tstat <- ((CNVonT1 * T1onT2) - CNVonT2)/(sqrt(new.se + seCNVonT2))
    pvalue <- 2 * pnorm(-abs(Tstat))
    
    # Calculate Z stat and its p-value
    indirect.Tstat <- CNVonT1onT2 / new.se
    indirectPvalue <- 2 * pnorm(-abs(indirect.Tstat))
    
    # Add every value to the dataframe
    comparisons$T1.trait[i] <- index$trait[index$ID == T1]
    comparisons$T2.trait[i] <- index$trait[index$ID == T2]
    comparisons$beta.CNV.on.T1[i] <- CNVonT1 
    comparisons$beta.CNV.on.T2[i] <- CNVonT2 
    comparisons$alpha.T1.on.T2[i] <- T1onT2
    comparisons$beta.CNV.on.T1.on.T2[i] <- CNVonT1onT2
    comparisons$se.CNV.on.T1[i] <- seCNVonT1
    comparisons$se.CNV.on.T2[i] <- seCNVonT2
    comparisons$se.T1.on.T2[i] <- seT1onT2 
    comparisons$se.CNV.on.T1.on.T2[i] <- new.se
    comparisons$T.stat[i] <- Tstat
    comparisons$p.value[i] <- pvalue
    comparisons$indirect.T.stat[i] <- indirect.Tstat
    comparisons$indirect.p.value[i] <- indirectPvalue
    comparisons$proportion <- comparisons$beta.CNV.on.T1.on.T2 / comparisons$beta.CNV.on.T2
    comparisons$signif.methods[i] <- MR_methods$number.methods[MR_methods$id.exposure == T1 & MR_methods$id.outcome == T2]
}

# Save the comparisons in a file
fwrite(comparisons, file = paste0(output_directory, "/comparisons.tsv"))

# Keep the significant comparisons (according to treshold and MTC set at the beginning)
signif.comp <- comparisons[comparisons$indirect.p.value <= treshold & comparisons$p.value <= treshold & comparisons$signif.methods >= MTC,]

# Keep the insignificant as well 
insignif.comp <- comparisons[comparisons$indirect.p.value > treshold | comparisons$p.value > treshold | comparisons$signif.methods < MTC,]








