setwd("~/Desktop/Uni/First step project")

folder <- list.files(path = "GWAS")
studies <- read.csv("studies.tsv", sep = "\t", header = T)


### Fichier différent pour chaque trait

for(i in 1:length(folder)){
    file <- read.csv(paste("GWAS/", folder[i], sep = ""), sep = "\t", header = T)
    CNVregion <- file[file$chromosome == "16" & file$base_pair_location >= 29400000 & file$base_pair_location <= 30400000,]
    write.table(CNVregion, paste("GWAS-CNV/", folder[i], sep = ""), sep = "\t", row.names = TRUE, col.names = TRUE)
}



### Tout dans un même fichier

AllTraits <- data.frame()

for(i in 1:length(folder)){
    file <- read.csv(paste("GWAS/", folder[i], sep = ""), sep = "\t", header = T)
    CNVregion <- file[file$chromosome == "16" & file$base_pair_location >= 29400000 & file$base_pair_location <= 30400000,]
    CNVregion$trait = studies$MAPPED_TRAIT[studies$STUDY.ACCESSION == sub("_.*", "", folder[i])]
    AllTraits <- rbind(AllTraits, CNVregion)
}

write.table(AllTraits, "GWAS-CNV/AllTraits.tsv", sep = "\t", row.names = F, col.names = T)



file <- read.csv(paste("GWAS/", folder[1], sep = ""), sep = "\t", header = T)
CNVregion <- file[file$chromosome == "16" & file$base_pair_location >= 29400000 & file$base_pair_location <= 30400000,]
CNVregion$trait = studies$MAPPED_TRAIT[studies$STUDY.ACCESSION == sub("_.*", "", folder[1])]







