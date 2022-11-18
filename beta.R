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
# install.packages("ggplot2")
library(ggplot2)
# install.packages("Cairo")
library(Cairo)
library(dplyr)
library(sys)
#install.packages("svglite")
library(svglite)
# install.packages("pheatmap")
library(pheatmap)
# install.packages("ggrepel")
library(ggrepel)
# BiocManager::install("biomaRt")
library(biomaRt)
# install.packages("scales")
library(scales)
# install.packages("gridExtra")
library(gridExtra)

setwd("~/THÉO/SummaryMR")


# Import probes and keep 1 per trait

probes <- as.data.frame(fread("probes.tsv", sep = "\t"))

probes_mirror <- probes[probes$effect_allele == "additional copy",]
probes_mirror <- probes_mirror[duplicated(probes$trait) == F,]
rownames(probes_mirror) <- 1:nrow(probes_mirror)


# Import results from MRs

MRs <- as.data.frame(fread("1Results/new_results.tsv", sep = "\t"))
MRs.methods <- as.data.frame(fread("1Results/eachpair.tsv", sep = "\t"))


# Create index, column for database and probes

db_trait <- c("Body mass index (BMI)", "Whole body fat mass", "Weight", "Standing height", "Hand grip strength (right)", "Heel bone mineral density (BMD)", "Forced vital capacity (FVC)", "Platelet count", "Neutrophill count", "C-reactive protein (quantile)", "Creatinine (quantile)", "Cystatin C (quantile)", "Vitamin D", "Albumin (quantile)", "Alkaline phosphatase (quantile)", "Alanine aminotransferase (quantile)", "Aspartate aminotransferase (quantile)", "Gamma glutamyltransferase (quantile)", "Triglycerides (quantile)", "Glycated haemoglobin (quantile)", "IGF-1 (quantile)", "SHBG (quantile)", "Fluid intelligence score", "Age when periods started (menarche)", "Relative age of first facial hair", "WHRadjBMI")
db_id <- c("ukb-b-3768", "ukb-d-30620_irnt", "ukb-d-30600_irnt", "ukb-d-30610_irnt", "ukb-d-30650_irnt", "ukb-b-19953", "ukb-d-30710_irnt", "ukb-d-30700_irnt", "ukb-d-30720_irnt", "ukb-b-5238", "ukb-b-7953", "ukb-d-30730_irnt", "ukb-d-30750_irnt", "ukb-b-10215", "ukb-b-8875", "ukb-d-30770_irnt", "ukb-d-30140_irnt", "ukb-d-30080_irnt", "ukb-b-5945", "ukb-d-30830_irnt", "ukb-b-10787", "ukb-d-30870_irnt", "ukb-d-30890_irnt", "ukb-b-11842", "ukb-b-19393", "ieu-a-79")
index <- data.frame(db_trait = sort(db_trait), db_id = db_id)

index$probes_trait <- c("age at menarche", "serum alanine aminotransferase measurement", "serum albumin measurement", "alkaline phosphatase measurement", "aspartate aminotransferase measurement", "body mass index", "C-reactive protein measurement", "creatinine measurement", "cystatin C measurement", "intelligence", "vital capacity", "serum gamma-glutamyl transferase measurement", "HbA1c measurement", "grip strength measurement", "heel bone mineral density", "insulin like growth factor measurement", "neutrophil count", "platelet count", "age at first facial hair", "sex hormone-binding globulin measurement", "body height", "triglyceride measurement", "vitamin D measurement", "body weight", "fat body mass", "waist-hip ratio")


# Create outcome dataframe

comparisons <- expand.grid(T1.id = db_id, T2.id = db_id, stringsAsFactors = FALSE)
comparisons <- comparisons[comparisons$T1 != comparisons$T2,]
rownames(comparisons) <- 1:nrow(comparisons)


for(i in 1:nrow(comparisons)){
    
    # take CNV on T1 from probes
    # take CNV on T2 from probes
    # take T1 on T2 from new_results -> IVW method
    
    T1 <- comparisons$T1.id[i]
    T2 <- comparisons$T2.id[i]
    
    CNVonT1 <- probes_mirror$beta[probes_mirror$trait == index$probes_trait[index$db_id == T1]]
    CNVonT2 <- probes_mirror$beta[probes_mirror$trait == index$probes_trait[index$db_id == T2]]
    T1onT2 <- MRs$b[MRs$id.exposure == T1 & MRs$id.outcome == T2 & MRs$method == "Inverse variance weighted"]
    CNVonT1onT2 <- CNVonT1 * T1onT2
    
    # Do the same for variances
    
    seCNVonT1 <- probes_mirror$standard_error[probes_mirror$trait == index$probes_trait[index$db_id == T1]]
    seCNVonT2 <- probes_mirror$standard_error[probes_mirror$trait == index$probes_trait[index$db_id == T2]]
    seT1onT2 <- MRs$se[MRs$id.exposure == T1 & MRs$id.outcome == T2 & MRs$method == "Inverse variance weighted"]
    
    # Perform the new se
    
    new.se <- seCNVonT1*seT1onT2 + seCNVonT1*T1onT2^2 + seT1onT2*CNVonT1^2

    # Calculate T stat and Z stat, and corresponding p-values
    
    Tstat <- ((CNVonT1 * T1onT2) - CNVonT2)/(sqrt(new.se + seCNVonT2))
    pvalue <- 2 * pnorm(-abs(Tstat))
    
    Zstat <- CNVonT1onT2 / new.se
    indirectPvalue <- 2 * pnorm(-abs(Zstat))
    
    
    # Add every value to the dataframe
    
    comparisons$T1.trait[i] <- index$probes_trait[index$db_id == T1]
    comparisons$T2.trait[i] <- index$probes_trait[index$db_id == T2]
    
    comparisons$beta.CNV.on.T1[i] <- CNVonT1 
    comparisons$beta.CNV.on.T2[i] <- CNVonT2 
    comparisons$beta.T1.on.T2[i] <- T1onT2
    comparisons$beta.CNV.on.T1.on.T2[i] <- CNVonT1onT2

    comparisons$se.CNV.on.T1[i] <- seCNVonT1
    comparisons$se.CNV.on.T2[i] <- seCNVonT2
    comparisons$se.T1.on.T2[i] <- seT1onT2 
    comparisons$se.CNV.on.T1.on.T2[i] <- new.se
    
    comparisons$T.stat[i] <- Tstat
    comparisons$p.value[i] <- pvalue
    
    comparisons$Z.stat[i] <- Zstat
    comparisons$indirect.p.value[i] <- indirectPvalue
    
    comparisons$signif.methods[i] <- MRs.methods$number.methods[MRs.methods$id.exposure == T1 & MRs.methods$id.outcome == T2]
    
}

signif.comp <- comparisons[comparisons$indirect.p.value <= 0.05 & comparisons$signif.methods >= 3,]
insignif.comp <- subset(comparisons, !(comparisons %in% signif.comp))


matrice <- matrix(nrow = 26, ncol = 26)
colnames(matrice) <- index$probes_trait
rownames(matrice) <- index$probes_trait

for(i in 1:nrow(signif.comp)){
    matrice[signif.comp$T2.trait[i], signif.comp$T1.trait[i]] <- signif.comp$beta.CNV.on.T2[i] - signif.comp$beta.CNV.on.T1.on.T2[i]
}

matrice_labels <- matrix(nrow = 26, ncol = 26)
colnames(matrice_labels) <- index$probes_trait
rownames(matrice_labels) <- index$probes_trait

for(i in 1:nrow(comparisons)){
    matrice_labels[comparisons$T2.trait[i], comparisons$T1.trait[i]] <- ifelse(comparisons$p.value[i] <= 0.05 & comparisons$indirect.p.value[i] <= 0.05, "*", "")
}

matrice[is.na(matrice)] <- 0
matrice_labels[is.na(matrice_labels)] <- ""

heatmap <- pheatmap(mat = matrice, cluster_rows = T, cluster_cols = F, # cannot do clustering when NAs are present
                    silent = T)

pheatmap(mat = matrice[,heatmap$tree_row$order], cluster_rows = T, cluster_cols = F, # cannot do clustering when NAs are present
         display_numbers = matrice_labels[,heatmap$tree_row$order], 
         color = hcl.colors(100, "Blue-Red 3"),
         number_color = "black",
         fontsize_number = 20,
         angle_col = 315,
         main = "Difference between direct and indirect effect",
         na_col = "black")


ggplot(data = signif.comp, aes(x = signif.comp$T2.trait)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5), plot.title = element_text(hjust = 0.5)) +
    ylim(-1, 1) +
    geom_point(y = signif.comp$beta.CNV.on.T1.on.T2, aes(color = signif.comp$T1.trait)) + # black dots for every trait 1
    geom_point(y = signif.comp$beta.CNV.on.T2, col = "red", shape = 18, size = 2) + # red diamonds for direct effect on trait 2
    labs(x = "Trait 2", y = "Effect of CNV on trait 2 passing by trait 1", title = "Direct and indirect effects of the CNV region on 26 associated traits") + 
    geom_text_repel(aes(label = signif.comp$T1.trait, x = signif.comp$T2.trait, y = signif.comp$beta.CNV.on.T1.on.T2), nudge_y = 0.03, size = 3)


# Locus zoom plot

mart.obj <- useMart(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl')

locus_genes <- getBM(c("entrezgene_id","external_gene_name","chromosome_name","start_position","end_position", "description", "phenotype_description"), 
                     filters = c("chromosome_name","start","end","with_entrezgene"), 
                     values = list(16, 29400000, 30400000, TRUE), 
                     mart = mart.obj)
unique_locus_genes <- distinct(locus_genes, external_gene_name, start_position, end_position, .keep_all = TRUE)
for(i in 1:nrow(unique_locus_genes)){
    unique_locus_genes$description[i] <- strsplit(unique_locus_genes$description[i], "\\[")[[1]]
}

probes_plot <- ggplot(data = probes_mirror, aes(x = base_pair_location)) +
    theme_bw() +
    geom_point(y = rescale(-log10(probes_mirror$p_value))) +
    ylim(0,1) +
    xlim(29400000,30400000) +
    labs(y = "-log10 of p-value (scaled)", title = "P-value of probes associated with the locus") +
    theme(plot.margin = margin(l = 0.5, unit = "cm"))

jitter <- rep(c(1,2), 25)

genes_plot <- ggplot(data = unique_locus_genes) +
    geom_segment(aes(x = start_position, xend = end_position, y = jitter, yend = jitter), linewidth = 3) +
    geom_text_repel(aes(label = description, x = start_position + ((end_position - start_position) / 2), y = jitter), size = 3, max.overlaps = 30) +
    theme_bw() +
    ylim(0,3) +
    xlim(29400000,30400000) +
    guides(y = "none") +
    labs(x = "Base pair location of the genes and the probes", y = "Genes") +
    theme(plot.margin = margin(l = 1.4, unit = "cm"))

grid.arrange(probes_plot, genes_plot, nrow = 2, heights = c(4,1))


# Pearson correlation





































