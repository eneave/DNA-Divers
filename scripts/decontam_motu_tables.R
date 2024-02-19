################################################
## DNA divers - organize data for basic plots ##
################################################
#######
## Read in MOTU tables
#######
# Disclaimer: some file names have 'ASVs' but the sequences have been clustered into MOTUs
## Phase 1 - all samples with tele02 marker
p1 <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/motu98_final_p1.csv")
p1 <- p1[c(2:90)]
## Phase 2 
## elas02 marker
p2e <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/motu98_final_p2e.csv")
p2e <- p2e[c(2:49)]
p2e_aq <- subset(p2e, select = -c(sample.6AeORK_eDNA_FBblank,sample.6BeORK_eDNAA_bottle1,
                                  sample.6CeORK_eDNAB_bottle2, sample.6DeORK_eDNAC_bottle3,
                                  sample.6EeORK_eDNAD_bottle4, sample.6FeORK_MPEtA_kurt,
                                  sample.6GeORK_MPEtB_mike, sample.6HeORK_MPEtC_lisa,
                                  sample.8AeORK_MPbB_mike, sample.8BeORK_MPbC_lisa))
p2e_uk <- subset(p2e, select = c(id, final_name, final_genus, final_order, final_class, pid, final_rank, 
                                  method_assign,total_reads, sequence,
                                  sample.11Be_EBMay_13extblank, sample.11Ce_EBJun_12extblank,  
                                  sample.11De_EBJun_19extblank,   
                                  sample.11Fe_negativePCRcontrol, sample.6BeORK_eDNAA_bottle1, 
                                  sample.6CeORK_eDNAB_bottle2, sample.6DeORK_eDNAC_bottle3,
                                  sample.6EeORK_eDNAD_bottle4, sample.6FeORK_MPEtA_kurt,
                                  sample.6GeORK_MPEtB_mike, sample.6HeORK_MPEtC_lisa,
                                  sample.8AeORK_MPbB_mike, sample.8BeORK_MPbC_lisa))

## tele02 marker
p2t <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/motu98_final_p2t.csv")
p2t <- p2t[c(2:50)]
p2t_aq <- subset(p2t, select = c(id, final_name,final_genus, final_order, final_class, 
                                 pid, final_rank, method_assign,                 
                                 total_reads, sequence, sample.5Et10BLUE_MPEtA_,
                                 sample.5Ft10BLUE_MPEtB_, sample.5Gt30BLUE_MPEtA_,        
                                 sample.5Ht30BLUE_MPEtB_, sample.6At60BLUE_MPEtA_,
                                 sample.6Bt60BLUE_MPEtB_, sample.6Ct120BLUE_MPEtA,
                                 sample.6Dt120BLUE_MPEtB, sample.6Et240BLUE_MPEtA,
                                 sample.6Ft240BLUE_MPEtB, sample.6GtBLUE_FBEtblank,
                                 sample.8Dt_EBMay_13extblank, sample.8Et_EBJun_12extblank,
                                 sample.8Ft_EBJun_13extblank))
p2t_uk <- subset(p2t, select = -c(sample.5Et10BLUE_MPEtA_, sample.5Ft10BLUE_MPEtB_, sample.5Gt30BLUE_MPEtA_,        
                                 sample.5Ht30BLUE_MPEtB_, sample.6At60BLUE_MPEtA_,
                                 sample.6Bt60BLUE_MPEtB_, sample.6Ct120BLUE_MPEtA,
                                 sample.6Dt120BLUE_MPEtB, sample.6Et240BLUE_MPEtA,
                                 sample.6Ft240BLUE_MPEtB, sample.6GtBLUE_FBEtblank,
                                 sample.4Ht_positivePCRcontrol))  
## Phase 3 - all samples with elas02 marker
p3 <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/motu98_final_p3.csv")
p3 <- p3[c(2:54)]

#######
## Remove contamination from MOTU tables using decontam prevalence method
#######
#####
library(tidyverse)
library(janitor)
library(decontam)
library(phyloseq)
library(ggplot2)
#####
## 1st Sequencing Run - All UK Samples
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
## LOOK at supp table1 and decide if I should be subsetting metadata from that table
p1_meta <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/metadata/p1_meta.csv")
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
## save a list of the identified contaminants
#write.csv(contamdf.prev05T, "C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/divers_methods_experiment/sintax_taxonomy/divmeth1_contamdf_prev05T.csv")
## remove generic named data, since same names will be used for different MOTU tables
rm(count_tab,sample_info_tab,tax_tab,OTU,TAX,SAM,dataset,df)
rm(contamdf.prev,contamdf.prev05,contamdf.prev05T,contamdf.prevT)

#####
#####
## 2nd Sequencing Run - Aquarium Samples - Shark tank
#####
# Convert to phyloseq objects 
# partition MOTU table to only include counts
# don't select extraction blank from June 19 because it's not usable
View(p2e_aq)
p2e_aq_v1 <- p2e_aq %>% select(c(id, sample.10Ae60BLUE_MPEtB:sample.11Ce_EBJun_12extblank,sample.11EeBLUE_FBEt_fieldb,sample.11Fe_negativePCRcontrol,sample.8CeBLUE_eDNA_FBblank:sample.9He60BLUE_MPEtA_)) %>%
                        column_to_rownames(var = "id") # make MOTU ids row names
count_tab <- p2e_aq_v1[rowSums(p2e_aq_v1[])>0,] # remove rows that sum to zero (they are there due to removal of non-aquarium samples)
p2e_aq_v2 <- count_tab %>% rownames_to_column(var = "id") # have a version of count_tab with a MOTU ID column
p2e_aq_v3 <- merge(p2e_aq_v2,p2e_aq[c(2:18)], by="id") # add taxonomy back 
# list of seq_ids
p2e_aq_names <- as.data.frame(colnames(p2e_aq_v2[-c(1)])) # 27 samples
colnames(p2e_aq_names)[1] <- "seq_id" # can use list to compare to sample data to make sure they're in the same order
# read in metadata
p2e_aq_meta <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/divers_methods_experiment_part2/p2e_aq_meta.csv")
colnames(p2e_aq_meta)[1] <- "seq_id"
sample_info_tab <- merge(p2e_aq_names,p2e_aq_meta, by="seq_id", all.x = T) %>% ##merge to make sure we have metadata for all of your samples
  column_to_rownames(var = "seq_id")
# partition MOTU table to only include taxa
tax_tab <- p2e_aq_v3 %>% select(c(id, kingdom, phylum, class, order, family, genus, species)) %>% 
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
contamdf.prevT <- merge(p2e_aq_v3, contamdf.prevT, by="id", all.x = F) 
# 0.5
contamdf.prev05T <- contamdf.prev05[contamdf.prev05$contaminant == TRUE,]
contamdf.prev05T <- contamdf.prev05T %>% rownames_to_column(var = "id")
contamdf.prev05T <- merge(p2e_aq_v3, contamdf.prev05T, by="id", all.x = F) 
## use 0.5 prevalence threshold
p2e_aq_prev <- p2e_aq_v3[!(p2e_aq_v3$id %in% contamdf.prev05T$id),] ## zero contaminants
write.csv(p2e_aq_prev, "C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/divers_methods_experiment_part2/sintax_output_MainTank/divmeth2_sharktank_MOTU_decontam.csv")
## remove generic named data, since same names will be used for different MOTU tables
rm(count_tab,sample_info_tab,tax_tab,OTU,TAX,SAM,dataset,df) 
rm(contamdf.prev,contamdf.prev05,contamdf.prev05T,contamdf.prevT)
#####
#####
## 2nd Sequencing Run - Aquarium Samples - Coral cave tank
#####
# Convert to phyloseq objects 
# partition MOTU table to only include counts (of samples and negative controls)
View(p2t_aq)
count_tab <- p2t_aq %>% select(c(id, sample.5Et10BLUE_MPEtA_:sample.8Dt_EBMay_13extblank)) %>%
  column_to_rownames(var = "id") # make MOTU ids row names
#count_tab <- p2e_aq_v1[rowSums(p2e_aq_v1[])>0,] # remove rows that sum to zero; no need to do this for coral cave; already done manually 
p2t_aq_v1 <- p2t_aq %>% select(c(id, sample.5Et10BLUE_MPEtA_:sample.8Dt_EBMay_13extblank))# have a version of count_tab with a MOTU ID column
#p2e_aq_v3 <- merge(p2e_aq_v2,p2e_aq[c(2:18)], by="id") # add taxonomy back; no need since all taxa have been kept
# list of seq_ids
p2t_aq_names <- as.data.frame(colnames(p2t_aq_v1[-c(1)])) # 12 samples + negative controls
colnames(p2t_aq_names)[1] <- "seq_id" # can use list to compare to sample data to make sure they're in the same order
# read in metadata
p2t_aq_meta <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/divers_methods_experiment_part2/p2t_aq_meta.csv")
colnames(p2t_aq_meta)[1] <- "seq_id" 
sample_info_tab <- merge(p2t_aq_names,p2t_aq_meta, by="seq_id", all.x = T) %>% ##merge to make sure we have metadata for all of your samples
  column_to_rownames(var = "seq_id")
# partition MOTU table to only include taxa
tax_tab <- p2t_aq %>% select(c(id, kingdom, phylum, class, order, family, genus, species)) %>% 
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
#contamdf.prevT <- contamdf.prev[contamdf.prev$contaminant == TRUE,]
#contamdf.prevT <- contamdf.prevT %>% rownames_to_column(var = "id")
#contamdf.prevT <- merge(p1, contamdf.prevT, by="id", all.x = F)
# 0.5
contamdf.prev05T <- contamdf.prev05[contamdf.prev05$contaminant == TRUE,]
contamdf.prev05T <- contamdf.prev05T %>% rownames_to_column(var = "id")
contamdf.prev05T <- merge(p2t_aq, contamdf.prev05T, by="id", all.x = F)
View(contamdf.prev05T)
## use 0.5 prevalence threshold
p2t_aq_prev <- p2t_aq[!(p2t_aq$id %in% contamdf.prev05T$id),] ## 11 contaminants
write.csv(p2t_aq_prev, "C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/divers_methods_experiment_part2/sintax_output_CoralCavedatabase/divmeth2_coralcave_MOTU_decontam.csv")
## remove generic named data, since same names will be used for different MOTU tables
rm(count_tab,sample_info_tab,tax_tab,OTU,TAX,SAM,dataset,df)
rm(contamdf.prev,contamdf.prev05,contamdf.prev05T,contamdf.prevT)

#####
## 2nd Sequencing Run - UK samples - Tele02 primer 
#####
# Convert to phyloseq objects 
# partition MOTU table to only include counts
View(p2t_uk)
p2t_uk_v1 <- p2t_uk %>% select(c(id, sample.4AtLIV_eDNA_FBblank:sample.5DtLIV_MPbA_Wendy, sample.6HtNEW_MPEt1_unis:sample.8CtNEW_MPb6_unis, sample.8Et_EBJun_12extblank:sample.8HtLIV_FBEtblank)) %>%
  column_to_rownames(var = "id") # make MOTU ids row names
count_tab <- p2t_uk_v1[rowSums(p2t_uk_v1[])>0,] # remove rows that sum to zero (they are there due to removal of non-aquarium samples)
p2t_uk_v2 <- count_tab %>% rownames_to_column(var = "id") # have a version of count_tab with a MOTU ID column
p2t_uk_v3 <- merge(p2t_uk_v2,p2t_uk[c(2:18)], by="id") # add taxonomy back 
# list of seq_ids
p2t_uk_names <- as.data.frame(colnames(p2t_uk_v2[-c(1)])) # 27 samples
colnames(p2t_uk_names)[1] <- "seq_id" # can use list to compare to sample data to make sure they're in the same order
# read in metadata
p2t_uk_meta <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/divers_methods_experiment_part2/p2t_uk_meta.csv")
colnames(p2t_uk_meta)[1] <- "seq_id"
sample_info_tab <- merge(p2t_uk_names,p2t_uk_meta, by="seq_id", all.x = T) %>% ##merge to make sure we have metadata for all of your samples
  column_to_rownames(var = "seq_id")
# partition MOTU table to only include taxa
tax_tab <- p2t_uk_v3 %>% select(c(id, kingdom, phylum, class, order, family, genus, species)) %>% 
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
contamdf.prevT <- merge(p2t_uk_v3, contamdf.prevT, by="id", all.x = F)
# 0.5
contamdf.prev05T <- contamdf.prev05[contamdf.prev05$contaminant == TRUE,]
contamdf.prev05T <- contamdf.prev05T %>% rownames_to_column(var = "id")
contamdf.prev05T <- merge(p2t_uk_v3, contamdf.prev05T, by="id", all.x = F)
## use 0.5 prevalence threshold
p2t_uk_prev <- p2t_uk_v3[!(p2t_uk_v3$id %in% contamdf.prev05T$id),] ## zero contaminants
write.csv(p2t_uk_prev, "C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/divers_methods_experiment_part2/sintax_output_UKreferencedatabase/t_divmeth2_UK_MOTU_decontam.csv")
## remove generic named data, since same names will be used for different MOTU tables
rm(count_tab,sample_info_tab,tax_tab,OTU,TAX,SAM,dataset,df) 
rm(contamdf.prev,contamdf.prev05,contamdf.prev05T,contamdf.prevT)

#####
## 2nd Sequencing Run - UK samples - Elas02 primer  
#####
##
# Convert to phyloseq objects 
# partition MOTU table to only include counts
View(p2e_uk)
p2e_uk2 <- merge(p2e_uk, p2e_sk, all = TRUE)
p2e_uk2[c(20:65)][is.na(p2e_uk2[c(20:65)])] <- 0

#p2e_uk_v1 <- p2e_uk %>% select(c(id, sample.11De_EBJun_19extblank, sample.11Fe_negativePCRcontrol:sample.8BeORK_MPbC_lisa)) 
#p2e_sk_v1 <- p2e_sk %>% select(c(id,sample.12A_ORKMP_boatEtA:sample.12G_ORKMPbA_kurt)) 

p2e_uk_v1 <- p2e_uk2 %>% select(c(id, sample.11De_EBJun_19extblank, 
                                  sample.11Fe_negativePCRcontrol:sample.8BeORK_MPbC_lisa,
                                  sample.12A_ORKMP_boatEtA:sample.12G_ORKMPbA_kurt)) %>%
                              column_to_rownames(var = "id")

count_tab <- p2e_uk_v1[rowSums(p2e_uk_v1[])>0,] # remove rows that sum to zero (they are there due to removal of non-aquarium samples)

p2e_uk_v2 <- count_tab %>% rownames_to_column(var = "id") # have a version of count_tab with a MOTU ID column
p2e_uk_v3 <- merge(p2e_uk_v2,p2e_uk2[c(1:18)], by="id") # add taxonomy back 
# list of seq_ids
p2e_uk_names <- as.data.frame(colnames(p2e_uk_v2[-c(1)])) # 19 samples
colnames(p2e_uk_names)[1] <- "seq_id" # can use list to compare to sample data to make sure they're in the same order
# read in metadata
p2e_uk_meta <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/divers_methods_experiment_part2/p2e_uk_meta.csv") 
colnames(p2e_uk_meta)[1] <- "seq_id"
sample_info_tab <- merge(p2e_uk_names,p2e_uk_meta, by="seq_id", all.x = T) %>% ##merge to make sure we have metadata for all of your samples
  column_to_rownames(var = "seq_id")
# partition MOTU table to only include taxa
tax_tab <- p2e_uk_v3 %>% select(c(id, kingdom, phylum, class, order, family, genus, species)) %>% 
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
contamdf.prevT <- merge(p2e_uk_v3, contamdf.prevT, by="id", all.x = F)
# 0.5
contamdf.prev05T <- contamdf.prev05[contamdf.prev05$contaminant == TRUE,]
contamdf.prev05T <- contamdf.prev05T %>% rownames_to_column(var = "id")
contamdf.prev05T <- merge(p2e_uk_v3, contamdf.prev05T, by="id", all.x = F)
## use 0.5 prevalence threshold
p2e_uk_prev <- p2e_uk_v3[!(p2e_uk_v3$id %in% contamdf.prev05T$id),] ## zero contaminants
write.csv(p2e_uk_prev, "C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/divers_methods_experiment_part2/sintax_output_UKreferencedatabase/e_divmeth2_UK_MOTU_decontam.csv")
## remove generic named data, since same names will be used for different MOTU tables
rm(count_tab,sample_info_tab,tax_tab,OTU,TAX,SAM,dataset,df) 
rm(contamdf.prev,contamdf.prev05,contamdf.prev05T,contamdf.prevT)




