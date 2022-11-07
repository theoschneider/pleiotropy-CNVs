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
library(sys)
#install.packages("svglite")
library(svglite)

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

'''

# Create the index 

db_trait <- c("Body mass index (BMI)", "Whole body fat mass", "Weight", "Standing height", "Hand grip strength (right)", "Heel bone mineral density (BMD)", "Forced vital capacity (FVC)", "Platelet count", "Neutrophill count", "C-reactive protein (quantile)", "Creatinine (quantile)", "Cystatin C (quantile)", "Vitamin D", "Albumin (quantile)", "Alkaline phosphatase (quantile)", "Alanine aminotransferase (quantile)", "Aspartate aminotransferase (quantile)", "Gamma glutamyltransferase (quantile)", "Triglycerides (quantile)", "Glycated haemoglobin (quantile)", "IGF-1 (quantile)", "SHBG (quantile)", "Fluid intelligence score", "Age when periods started (menarche)", "Relative age of first facial hair", "WHRadjBMI")
db_id <- c("ukb-b-3768", "ukb-d-30620_irnt", "ukb-d-30600_irnt", "ukb-d-30610_irnt", "ukb-d-30650_irnt", "ukb-b-19953", "ukb-d-30710_irnt", "ukb-d-30700_irnt", "ukb-d-30720_irnt", "ukb-b-5238", "ukb-b-7953", "ukb-d-30730_irnt", "ukb-d-30750_irnt", "ukb-b-10215", "ukb-b-8875", "ukb-d-30770_irnt", "ukb-d-30140_irnt", "ukb-d-30080_irnt", "ukb-b-5945", "ukb-d-30830_irnt", "ukb-b-10787", "ukb-d-30870_irnt", "ukb-d-30890_irnt", "ukb-b-11842", "ukb-b-19393", "ieu-a-79")
index <- data.frame(db_trait = sort(db_trait), db_id = db_id)


# Create all combinations of IDs (Neale Lab round 2, IEU if not available) (exception : WHR is adj for BMI, and pop is european instead of UK)

combinations <- expand.grid(exposure = db_id, outcome = db_id, stringsAsFactors = FALSE)
combinations <- combinations[combinations$exposure != combinations$outcome,]
rownames(combinations) <- 1:nrow(combinations)


# Running the MRs

knitr::opts_chunk$set(dpi = 500, fig.width = 7)

results <- data.frame()

for(i in 501:650){
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





exposure_id = combinations$exposure[77]
outcome_id = combinations$outcome[77]
exp_dat <- extract_instruments(outcomes = exposure_id,
                               clump = TRUE,
                               p1 = 5e-8,
                               r2 = 0.001)
out_dat <- extract_outcome_data(snps = exp_dat$SNP,
                                outcomes = outcome_id)

which(table(out_dat$SNP) != 1)
exp_dat[exp_dat$SNP == "rs377381568",]
out_dat[out_dat$SNP == "rs377381568",]

dat <- harmonise_data(exp_dat, out_dat[-c(18, 172),])
dat <- harmonise_data(exp_dat, out_dat)

res <- mr(dat)
p1 <- mr_scatter_plot(res, dat)
p1[[1]]

res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]

res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]

res_single <- mr_singlesnp(dat)
p4 <- mr_funnel_plot(res_single)
p4[[1]]


mr_report(dat, output_path = "~/THÉO/SummaryMR/2Results/")

exp_dat_new$pval.exposure[exp_dat_new$pval.exposure < 1e-14] <- 1e-14

out_dat_new <- out_dat
out_dat_new$pval.outcome[out_dat_new$pval.exposure < 1e-14] <- 1e-14


dat_subset <- harmonise_data(exp_dat[1:200,], out_dat[1:200,])
mr_report(dat_subset, output_path = "~/THÉO/SummaryMR/2Results/")

# problème ligne 104 : pval 10^-98
# 10^-63 is also too low
# 10^-39 is fine
# 10^-56 is fine
# 10^57 is too low

# works with 106 numbers, not with 107 or higher


# Keep only the results where p ≥ 0.05/650 and sort them

results <- as.data.frame(fread("1Results/new_results.tsv", sep = "\t"))

bon_correction <- results[results$pval <= (0.05/650),]
bon_correction <- bon_correction[order(bon_correction$pval),]

fwrite(bon_correction, file = "1Results/results_bonferroni_correction.tsv", quote = FALSE, sep = "\t")


# Create a collapse column, containing the IDs pasted

bon_correction$collapse <- paste0(bon_correction$id.exposure, bon_correction$id.outcome)


# Create another dataframe with only the pairs that have 2 or more significant methods

significant <- data.frame()

for(pair in unique(bon_correction$collapse)){
    if(sum(bon_correction$collapse == pair) >= 2){
        significant <- rbind(significant, bon_correction[bon_correction$collapse == pair, 1:9])
    }
}

significant <- significant[order(significant$pval),]

fwrite(significant, file = "1Results/significant_MR.tsv", quote = FALSE, sep = "\t")


# Create file with number of significant methods and best p-value for each pair

results$collapse <- paste0(results$id.exposure, results$id.outcome)
eachpair <- data.frame()
for(pair in unique(results$collapse)){
    fiveMethods <- results[results$collapse == pair,]
    line <- fiveMethods[fiveMethods$pval == min(fiveMethods$pval),][1,]
    line$number.methods <- sum(fiveMethods$pval <= (0.05/650))
    eachpair <- rbind(eachpair, line)
}

colnames(eachpair)[5] <- "best.method"

sum(eachpair$number.methods >= 3)

eachpair <- eachpair[order(eachpair$number.methods, decreasing = T),]
fwrite(eachpair, file = "1Results/eachpair.tsv", quote = FALSE, sep = "\t")



# Try some multivariable MR

exposure1_id = "ukb-b-11842"    # poids
exposure2_id = "ukb-b-10787"    # taille
outcome_id = "ukb-b-19953"      # BMI
 

exp_dat <- mv_extract_exposures(id_exposure = db_id[-1], 
                                clump_r2 = 0.001,
                                pval_threshold = 5e-8)

out_dat <- extract_outcome_data(snps = exp_dat$SNP,
                                outcomes = db_id[1])

dat <- mv_harmonise_data(exp_dat, out_dat)

mv_example <- mv_basic(dat)

mv_example$result      # pvalues of each trait against the outcome
mv_example$plots[[1]]  # plot of SNP on outcome VS SNP on exp 1
mv_example$plots[[2]]  # plot of SNP on outcome VS SNP on exp 2


# Automatic MVMR

results <- as.data.frame(fread("1Results/significant_MR.tsv", sep = "\t"))

mvmr_df <- data.frame()


for(outcome in unique(results$id.outcome)){
    exposures = unique(results$id.exposure[results$id.outcome == outcome])
    exp_dat <- mv_extract_exposures(id_exposure = exposures,
                                    clump_r2 = 0.001,
                                    pval_threshold = 5e-8)
    
    out_dat <- extract_outcome_data(snps = exp_dat$SNP,
                                    outcomes = outcome)
    
    dat <- mv_harmonise_data(exp_dat, out_dat)
    
    mvmr <- mv_basic(dat)
    mvmr_df <- rbind(mvmr_df, mvmr$result)
    for(i in 1:length(mvmr$plots)){
        ggsave(filename = paste0("3MVMR/", mvmr$result$exposure[i], "_on_", mvmr$result$outcome[i], ".svg"), plot = mvmr$plots[[i]], width = 10, height = 8)
    }
}

fwrite(mvmr_df, file = "1Results/MVMR_results.tsv", quote = FALSE, sep = "\t")
























