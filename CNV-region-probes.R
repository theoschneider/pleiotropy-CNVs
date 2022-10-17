setwd("~/Desktop/Uni/First step project")

folder <- list.files(path = "GWAS")
studies <- read.csv("studies.tsv", sep = "\t", header = T)


### Tout dans un même fichier

# Isolate the right region, then find the smallest p-value in this region, then add the trait

AllTraits <- data.frame()

for(i in 1:length(folder)){
    file <- read.csv(paste("GWAS/", folder[i], sep = ""), sep = "\t", header = T)
    CNVregion <- file[file$chromosome == "16" 
                      & file$base_pair_location >= 29400000 
                      & file$base_pair_location <= 30400000,]
    CNVregion <- CNVregion[CNVregion$p_value == min(CNVregion$p_value),]
    CNVregion$trait = studies$MAPPED_TRAIT[studies$STUDY.ACCESSION == sub("_.*", "", folder[i])]
    AllTraits <- rbind(AllTraits, CNVregion)
}

write.table(AllTraits, "GWAS-CNV/AllTraits.tsv", sep = "\t", row.names = F, col.names = T)

length(unique(AllTraits$trait)) # check that we have 26 traits














