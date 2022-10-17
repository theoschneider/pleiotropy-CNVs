setwd("~/Desktop/Uni/First step project/GWAS-Neale")
# install.packages("stringr")
library(stringr)

manifest <- read.csv(file = "UKBB-GWAS-Manifest.tsv", sep = "\t", header = T)
traits <- c("Body mass index (BMI)", "Whole body fat mass", "Weight", "Standing height", "Hand grip strength (right)", "Heel bone mineral density (BMD)", "Forced vital capacity (FVC)", "Platelet count", "Neutrophill count", "C-reactive protein (quantile)", "Creatinine (quantile)", "Cystatin C (quantile)", "Vitamin D", "Albumin (quantile)", "Alkaline phosphatase (quantile)", "Alanine aminotransferase (quantile)", "Aspartate aminotransferase (quantile)", "Gamma glutamyltransferase (quantile)", "Triglycerides (quantile)", "Glycated haemoglobin (quantile)", "IGF-1 (quantile)", "SHBG (quantile)", "Fluid intelligence score", "Age when periods started (menarche)", "Relative age of first facial hair")
# only 25 traits because WHR is a separate file (already dl)

selected_phenotypes <- manifest[grepl(pattern = "irnt", x = manifest$Phenotype.Code) == T & manifest$Sex == "both_sexes" & manifest$Phenotype.Description %in% traits == T,]
# irnt = inversed rank-normal transformed

# weight and BMI appear 2 times :
selected_phenotypes <- selected_phenotypes[-c(5, 7),]


# menarche and facial hair are  missing because there is no irnt

selected_phenotypes <- rbind(selected_phenotypes, na.omit(manifest[manifest$Phenotype.Description == "Age when periods started (menarche)" & manifest$Sex == "both_sexes",]))
selected_phenotypes <- rbind(selected_phenotypes, na.omit(manifest[manifest$Phenotype.Description == "Relative age of first facial hair" & manifest$Sex == "both_sexes",]))


# check the phenotypes :

summary(traits %in% selected_phenotypes$Phenotype.Description)


# export df as tsv

write.table(selected_phenotypes, "selected_phenotypes.tsv", sep = "\t", row.names = F, col.names = T)


# create wget commands with right name and directory

commands = ""
for(link in selected_phenotypes$AWS.File){
    commands = paste(commands, "wget -O ~/THÉO/SummaryMR/", str_replace_all(selected_phenotypes[selected_phenotypes$AWS.File == link, 2], c(" " = "-", "\\(" = "", "\\)" = "")), ".tsv.bgz ", link, "\n", sep = "")
}


# export commands in a sh file

write.table(commands, "commands.sh", row.names = F, col.names = F, quote = F)







