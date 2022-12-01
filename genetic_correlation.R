library(data.table)
library(dplyr)
library(ggplot2)

setwd("~/THÉO/SummaryMR")


db_trait <- c("Body mass index (BMI)", "Whole body fat mass", "Weight", "Standing height", "Hand grip strength (right)", "Heel bone mineral density (BMD)", "Forced vital capacity (FVC)", "Platelet count", "Neutrophill count", "C-reactive protein (quantile)", "Creatinine (quantile)", "Cystatin C (quantile)", "Vitamin D", "Albumin (quantile)", "Alkaline phosphatase (quantile)", "Alanine aminotransferase (quantile)", "Aspartate aminotransferase (quantile)", "Gamma glutamyltransferase (quantile)", "Triglycerides (quantile)", "Glycated haemoglobin (quantile)", "IGF-1 (quantile)", "SHBG (quantile)", "Fluid intelligence score", "Age when periods started (menarche)", "Relative age of first facial hair", "WHRadjBMI")
db_id <- c("ukb-b-3768", "ukb-d-30620_irnt", "ukb-d-30600_irnt", "ukb-d-30610_irnt", "ukb-d-30650_irnt", "ukb-b-19953", "ukb-d-30710_irnt", "ukb-d-30700_irnt", "ukb-d-30720_irnt", "ukb-b-5238", "ukb-b-7953", "ukb-d-30730_irnt", "ukb-d-30750_irnt", "ukb-b-10215", "ukb-b-8875", "ukb-d-30770_irnt", "ukb-d-30140_irnt", "ukb-d-30080_irnt", "ukb-b-5945", "ukb-d-30830_irnt", "ukb-b-10787", "ukb-d-30870_irnt", "ukb-d-30890_irnt", "ukb-b-11842", "ukb-b-19393", "ieu-a-79")
index <- data.frame(db_trait = sort(db_trait), db_id = db_id)

index$probes_trait <- c("age at menarche", "serum alanine aminotransferase measurement", "serum albumin measurement", "alkaline phosphatase measurement", "aspartate aminotransferase measurement", "body mass index", "C-reactive protein measurement", "creatinine measurement", "cystatin C measurement", "intelligence", "vital capacity", "serum gamma-glutamyl transferase measurement", "HbA1c measurement", "grip strength measurement", "heel bone mineral density", "insulin like growth factor measurement", "neutrophil count", "platelet count", "age at first facial hair", "sex hormone-binding globulin measurement", "body height", "triglyceride measurement", "vitamin D measurement", "body weight", "fat body mass", "waist-hip ratio")
index$db_dim <- c("Age at menarche", "ALT", "ALB", "ALP", "AST", "BMI", "CRP", "Creatinin", "Cystatin C", "Intelligence", "FVC", "GGT", "HbA1c", "Grip strength", "Heel bone mineral density", "IGF-1", "Neutrophils", "Platelets", "Age at first facial hair", "SHBG", "Height", "Triglycerides", "VitD", "Weight", "Body fat mass", "WHR")

index$neale_id <- c("2714", "30620_irnt", "30600_irnt", "30610_irnt", "30650_irnt", "21001_irnt", "30710_irnt", "30700_irnt", "30720_irnt", "20016_irnt", "3062_irnt", "30730_irnt", "30750_irnt", "47_irnt", "3148_irnt", "30770_irnt", "30140_irnt", "30080_irnt", "2375", "30830_irnt", "50_irnt", "30870_irnt", "100021_irnt", "21002_irnt", "23100_irnt", NA)

gen_cor <- fread("geno_correlation_sig.txt")
gen_cor_int <- gen_cor[p1 %in% index$neale_id & p2 %in% index$neale_id, c(3, 14, 15)] 
# In theory, there are 300 values of correlation (25*24 / 2)
# However there are only 66

# ggplot(data = gen_cor_int, aes(x = description_p1, y = description_p2, fill = rg)) + 
#     geom_tile() +
#     scale_fill_gradient2(low = "red", high = "blue", mid = "white", midpoint = 0, na.value = "grey50") +
#     theme_classic() +
#     theme(axis.text.x = element_text(angle = -45, vjust = 0, hjust = 0))

matrice <- matrix(nrow = 26, ncol = 26)
colnames(matrice) <- index$db_trait
rownames(matrice) <- index$db_trait

for(i in 1:nrow(gen_cor_int)){
    matrice[gen_cor_int$description_p1[i], gen_cor_int$description_p2[i]] <- gen_cor_int$rg[i]
    matrice[gen_cor_int$description_p2[i], gen_cor_int$description_p1[i]] <- gen_cor_int$rg[i]
}

matrice[upper.tri(matrice)] <- 0

bk1 <- c(seq(-0.4, -0.1, by = 0.1), -0.001)
bk2 <- c(0.001, seq(0.1, 1, by = 0.1))
bk <- c(bk1, bk2)  #combine the break limits for purpose of graphing

my_palette <- c(colorRampPalette(colors = c("darkblue", "white"))(n = length(bk1) - 1),
                "white", "white",
                c(colorRampPalette(colors = c("white", "darkred"))(n = length(bk2) - 1)))


p <- pheatmap(mat = matrice, cluster_rows = F, cluster_cols = F, # cannot do clustering when NAs are present
         color = my_palette,
         breaks = bk,
         fontsize_number = 20,
         angle_col = 315,
         main = "Genetic correlation between traits",
         na_col = "grey80")

ggsave(filename = "test cor.png", plot = p, width = 11, height = 10, dpi = 500)

# Photoshop pour mettre l'axe y sur la gauche












