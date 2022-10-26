# install.packages("remotes")
# library(remotes)
# remotes::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
# remotes::install_github("MRCIEU/MRInstruments")
library(MRInstruments)
# install.packages("R.utils")
library(R.utils) # for gunzip
# install.packages("data.table")
library(data.table) # for fread and fwrite
library(ggplot2)
# install.packages("Cairo")
library(Cairo)
library(dplyr)

setwd("~/THÉO/SummaryMR")

'''

# Unzip all the bgz files (unnecessary)

filenames <- list.files(path = "~/THÉO/SummaryMR", pattern = "*.bgz")
for(file in filenames){
    gunzip(filename = file, destname = gsub(".bgz", "", file), remove = FALSE, ext = "bgz")
}


# Format data to be the same as WHR.ma (unnecessary)

filenames <- list.files(path = "~/THÉO/SummaryMR", pattern = ".tsv")

for(fichier in filenames){
    data <- as.data.frame(fread(fichier, sep = "\t"))
    colnames(data)[which(names(data) == "variant")] <- "SNP"
    colnames(data)[which(names(data) == "minor_allele")] <- "A1"
    colnames(data)[which(names(data) == "minor_AF")] <- "Freq"
    colnames(data)[which(names(data) == "beta")] <- "b"
    colnames(data)[which(names(data) == "pval")] <- "p"
    colnames(data)[which(names(data) == "n_complete_samples")] <- "N"
    data$A2 <- sapply(strsplit(data$SNP, ":"), `[`, 3)
    fwrite(data, file = gsub(".tsv", ".ma", fichier), quote = FALSE, sep = "\t")
}


# Create all possible combinations (unnecessary)

db_trait <- c("Body mass index (BMI)", "Whole body fat mass", "Weight", "Standing height", "Hand grip strength (right)", "Heel bone mineral density (BMD)", "Forced vital capacity (FVC)", "Platelet count", "Neutrophill count", "C-reactive protein (quantile)", "Creatinine (quantile)", "Cystatin C (quantile)", "Vitamin D", "Albumin (quantile)", "Alkaline phosphatase (quantile)", "Alanine aminotransferase (quantile)", "Aspartate aminotransferase (quantile)", "Gamma glutamyltransferase (quantile)", "Triglycerides (quantile)", "Glycated haemoglobin (quantile)", "IGF-1 (quantile)", "SHBG (quantile)", "Fluid intelligence score", "Age when periods started (menarche)", "Relative age of first facial hair", "WHRadjBMI")
combinations <- expand.grid(exposure = db_trait, outcome = db_trait, stringsAsFactors = FALSE)
combinations <- combinations[combinations$exposure != combinations$outcome,]
rownames(combinations) <- 1:nrow(combinations)

'''


# Add corresponding IDs (Neale Lab round 2, IEU if not available) (exception : WHR is adj for BMI, and pop is european instead of UK)

db_id <- c("ukb-b-3768", "ukb-d-30620_irnt", "ukb-d-30600_irnt", "ukb-d-30610_irnt", "ukb-d-30650_irnt", "ukb-b-19953", "ukb-d-30710_irnt", "ukb-d-30700_irnt", "ukb-d-30720_irnt", "ukb-b-5238", "ukb-b-7953", "ukb-d-30730_irnt", "ukb-d-30750_irnt", "ukb-b-10215", "ukb-b-8875", "ukb-d-30770_irnt", "ukb-d-30140_irnt", "ukb-d-30080_irnt", "ukb-b-5945", "ukb-d-30830_irnt", "ukb-b-10787", "ukb-d-30870_irnt", "ukb-d-30890_irnt", "ukb-b-11842", "ukb-b-19393", "ieu-a-79")
combinations <- expand.grid(exposure = db_id, outcome = db_id, stringsAsFactors = FALSE)
combinations <- combinations[combinations$exposure != combinations$outcome,]
rownames(combinations) <- 1:nrow(combinations)


# Running the MRs

knitr::opts_chunk$set(dpi = 500, fig.width = 7)

results <- data.frame()

for(i in 77:77){
    exposure_id = combinations$exposure[i]
    outcome_id = combinations$outcome[i]
    exp_dat <- extract_instruments(outcomes = exposure_id,
                                   clump = TRUE,
                                   p1 = 5e-8,
                                   r2 = 0.001)
    out_dat <- extract_outcome_data(snps = exp_dat$SNP,
                                    outcomes = outcome_id)
    dat <- harmonise_data(exp_dat, out_dat)
    dat_mr <- mr(dat)
    results <- rbind(results, dat_mr)
    mr_report(dat, output_path = "~/THÉO/SummaryMR/1Results")
}

fwrite(results, file = "1Results/results.tsv", append = TRUE, quote = FALSE, sep = "\t")

# 107 seconds for 5 MR = 22 seconds per MR, total of 650*22 = 14300 seconds = 4 hours



# Import results as dataframe, remove duplicates and save in a file

results_df <- as.data.frame(fread("1Results/results.tsv", sep = "\t"))

new_results_df <- distinct(results_df, id.exposure, id.outcome, method, .keep_all = TRUE)

rownames(new_results_df) <- 1:nrow(new_results_df)

fwrite(new_results_df, file = "1Results/new_results.tsv", quote = FALSE, sep = "\t")





















