################################################
## DNA divers - organize data for basic plots ##
################################################
#######
## Read in raw MOTU tables
#######
# Disclaimer: some file names have 'ASVs' but the sequences have been clustered into MOTUs
## Phase 1 - all samples with tele02 marker
p1 <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/divers_methods_experiment/sintax_taxonomy/divmeth1_ASVs.csv")
## Phase 2 
## elas02 marker
p2e_aq <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/divers_methods_experiment_part2/sintax_output_MainTank/e_st_divmeth2_ASVs.csv") #shark tank aquarium
p2e_uk <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/divers_methods_experiment_part2/sintax_output_UKreferencedatabase/e_divmeth2_ASVs.csv") #UK samples
p2e_sk <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/St.Kilda_andothersamples/Elas02_run/other_sintax/elas_other_ASVs.csv") #UK samples accidently not sequenced/added on to a different project
## tele02 marker
p2t_aq <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/divers_methods_experiment_part2/sintax_output_CoralCavedatabase/t_cc_divmeth2_ASVs.csv") #coral cave aquarium
p2t_uk <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/divers_methods_experiment_part2/sintax_output_UKreferencedatabase/t_divmeth2_ASVs.csv") #UK samples

#######
## Remove contamination from MOTU tables
#######
#####
library(tidyverse)
library(janitor)
library(decontam)
library(phyloseq)
library(ggplot2)
#####
## 1st Sequencing Run
#####
# Convert to phyloseq objects 
# partition MOTU table to only include counts
count_tab <- p1 %>% select(c(id, sample.10A:sample.9H))
# list of seq_ids
p1_names <- as.data.frame(colnames(count_tab[-c(1)])) # 79 samples
colnames(p1_names)[1] <- "seq_id" # can use list to compare to sample data to make sure they're in the same order
# make MOTU ids row names
count_tab <- count_tab %>%  column_to_rownames(var = "id")
# read in metadata
p1_meta <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/divers_methods_experiment/p1_meta.csv")
colnames(p1_meta)[1] <- "seq_id"
sample_info_tab <- merge(p1_names,p1_meta, by="seq_id", all.x = T) %>% ##merge because 3 controls dropped out for having no reads assigned
                   column_to_rownames(var = "seq_id")
# partition MOTU table to only include taxa
tax_tab <- p1 %>% select(c(id, kingdom, phylum, class, order, family, genus, species)) %>% 
                  replace_na(list(kingdom = "unassigned", phylum = "unassigned", class = "unassigned", order = "unassigned", family = "unassigned", genus = "unassigned", species = "unassigned")) %>%
                  column_to_rownames(var = "id")
# create phyloseq object
OTU = otu_table(as.matrix(count_tab), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(tax_tab))
SAM = sample_data(sample_info_tab)
dataset <- merge_phyloseq(phyloseq(OTU, TAX), SAM)

# code taken from tutorial: https://benjjneb.github.io/decontam/vignettes/decontam_intro.html#identify-contaminants---frequency
# inspect library size
df <- as.data.frame(sample_data(dataset)) # Put sample_data into a ggplot-friendly data.frame
df$SampleReads <- sample_sums(dataset)
df <- df[order(df$SampleReads),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=SampleReads, color=sampletype2)) + geom_point() 
# identify contaminants through prevalence 
sample_data(dataset)$is.neg <- sample_data(dataset)$sampletype2 == "negative"
# Default prevalence threshold of 0.1
contamdf.prev <- isContaminant(dataset, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))
# Prevalence threshold of 0.5
contamdf.prev05 <- isContaminant(dataset, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)
# extract rows which are possible contaminants and merge with taxonomy to inspect
# 0.1
contamdf.prevT <- contamdf.prev[contamdf.prev$contaminant == TRUE,]
contamdf.prevT <- contamdf.prevT %>% rownames_to_column(var = "id")
contamdf.prevT <- merge(p1, contamdf.prevT, by="id", all.x = F)
# 0.5
contamdf.prev05T <- contamdf.prev05[contamdf.prev05$contaminant == TRUE,]
contamdf.prev05T <- contamdf.prev05T %>% rownames_to_column(var = "id")
contamdf.prev05T <- merge(p1, contamdf.prev05T, by="id", all.x = F)
## use 0.5 prevalence threshold
p1_prev <- p1[!(p1$id %in% contamdf.prev05T$id),]
write.csv(p1_prev, "C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/divers_methods_experiment/sintax_taxonomy/divmeth1_MOTU_decontam.csv")



