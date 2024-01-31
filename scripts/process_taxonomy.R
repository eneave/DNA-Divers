############################################################
## Process taxonomic output from ecotag, sintax, & blastn ##
############################################################
# load packages
library(splitstackshape)
library(tidyverse)

#####
## Load data
#####
# ecotag
# open ecotag output fasta file notepad, copy to excel and sort by motu ids
# keep motu ids and taxonomic info, remove sequences
# save as a csv file
eco <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/divmeth1_SWARM1_nohuman_ecotag_noseq.csv") # sequence run 1
#eco <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/e_divemeth2_SWARM1_ecotag_noseq.csv") # sequence run 2, elas02 library
#eco <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/t_divemeth2_SWARM1_ecotag_noseq.csv") # sequence run 2, tele02 library


# sintax
# Read .tsv file generated from using the sintax option in vsearch
tax <- read.csv(file = "C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/SWARM1_sintax_v2_26_1_v258.tsv", sep = ';', header = FALSE) # sequence run 1
#tax <- read.csv(file = "C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/e_SWARM1_sintax_output_ALL.tsv", sep = ';', header = FALSE) # sequence run 2, elas02 library
#tax <- read.csv(file = "C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/NEEDTODO.tsv", sep = ';', header = FALSE) # sequence run 2, tele02 library

# blastn
# Read .tsv file generated from using the blastn option in blastn with a word size of 7
bn <- read.csv(file = "C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/NEEDTODO.tsv", sep = "", header = TRUE) # sequence run 1
#bn <- read.csv(file = "C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/elas02_divers2_blast_v3.tsv", sep = "", header = TRUE) # sequence run 2, elas02 library

# read counts
# Read count data csv file generated from obitab into R
abund <- read.delim("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/divmeth1_SWARM1_output.counts.csv", sep=";", header=T) # sequence run 1
#abund <- read.delim("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/e_divmeth2_SWARM1_output.counts.csv", sep=";", header=T) # sequence run 2, elas02 library
#abund <- read.delim("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/e_divmeth2_SWARM1_output.counts.csv", sep=";", header=T) # sequence run 2, elas02 library
#abund <- read.delim("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/e_other_SWARM1_output.counts.csv", sep=";", header=T) # sequence run 3


#####
## Process ecotag output
#####
# separate by semi-colons into columns
eco_v1 <- eco %>% 
  unnest(id) %>%
  separate(id, c("id", "id_status", "family", "species_name", "best_match",
                 "taxid_by_db", "rank_by_db", "scientific name",
                 "match_count", "rank", "taxid", "species",
                 "order_name", "best_identity", "scientific_name_by_db",
                 "spcies_list", "genus_name", "family_name", "genus", "order"), sep = "\\;",
           convert = TRUE, extra = "merge")
eco_v1 <- eco_v1 %>% separate(id, c("id", "count"), sep = " ")

# separate the best rank scientific name
eco_v2 <- cSplit(eco_v1, c("scientific name"), sep = c("="))
eco_v2 <- cSplit(eco_v2, c("best_identity"), sep = c(":"))
eco_v2 <- cSplit(eco_v2, c("best_identity_2"), sep = c("}"))
eco_v2 <- cSplit(eco_v2, c("rank"), sep = c("="))
eco_v2 <- cSplit(eco_v2, c("species_name"), sep = c("="))
eco_v2 <- cSplit(eco_v2, c("genus_name"), sep = c("="))
eco_v2 <- cSplit(eco_v2, c("family_name"), sep = c("="))
eco_v2 <- cSplit(eco_v2, c("order_name"), sep = c("="))

# fix motu ids
eco_v2$id <-gsub(">","",as.character(eco_v2$id))

# remove unnecessary columns
eco_v3 <- subset(eco_v2, select= -c(best_identity_1,
                                            rank_1, species_name_1, genus_name_1,
                                            family_name_1, order_name_1, count, id_status,
                                            family, taxid_by_db, rank_by_db, match_count,
                                            taxid, species, scientific_name_by_db, genus, order))
eco_v3 <- subset(eco_v3, select= -c(`scientific name_1`))

# rename columns to meaningful names
eco_clean <- eco_v3 %>% 
  rename(scientific_name = `scientific name_2`,
         best_identity = best_identity_2_1,
         rank = rank_2,
         species_name = species_name_2,
         genus_name = genus_name_2, 
         family_name = family_name_2,
         order_name = order_name_2)

#####
## Process sintax output
#####
tax2 <- tax %>%
  unnest(V1) %>%
  separate(V1, into = c("all_assign", "qc_assign_70"), 
           sep = "\\+", convert = TRUE, extra = "merge")

tax3 <- tax2 %>% separate(all_assign, c("id", "assign"), sep = 14)

tax4 <- tax3 %>%
  select(id, assign, qc_assign_70) %>%
  unnest(assign) %>%
  separate(assign, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), 
           sep = ",", convert = TRUE, extra = "merge")

# separate percent identity from taxonomic information
tax5 <- cSplit(tax4, c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = c("("))

# remove non-descriptive characters
tax5$kingdom_1<-gsub("k:","",as.character(tax5$kingdom_1))
tax5$kingdom_2<-gsub(")","",as.character(tax5$kingdom_2))
tax5$phylum_1<-gsub("p:","",as.character(tax5$phylum_1))
tax5$phylum_2<-gsub(")","",as.character(tax5$phylum_2))
tax5$class_1<-gsub("c:","",as.character(tax5$class_1))
tax5$class_2<-gsub(")","",as.character(tax5$class_2))
tax5$order_1<-gsub("o:","",as.character(tax5$order_1))
tax5$order_2<-gsub(")","",as.character(tax5$order_2))
tax5$family_1<-gsub("f:","",as.character(tax5$family_1))
tax5$family_2<-gsub(")","",as.character(tax5$family_2))
tax5$genus_1<-gsub("g:","",as.character(tax5$genus_1))
tax5$genus_2<-gsub(")","",as.character(tax5$genus_2))
tax5$species_1<-gsub("s:","",as.character(tax5$species_1))
tax5$species_2<-gsub(")","",as.character(tax5$species_2))
tax5$species_1<-gsub("_"," ",as.character(tax5$species_1))

# rename columns to meaningful names
tax_clean <- tax5 %>% 
  rename(
    kingdom = kingdom_1,
    k.id = kingdom_2,
    phylum = phylum_1,
    p.id = phylum_2,
    class = class_1,
    c.id = class_2, 
    order = order_1,
    o.id = order_2,
    family = family_1,
    f.id = family_2,
    genus = genus_1,
    g.id = genus_2,
    species = species_1,
    s.id = species_2 
  )

#####
## Process blastn output
#####
bn1 <- bn %>%
  select(id, blastDbid, blastEvalue, blastLength, blastPident, blastNident, blastScore, blastBitscore) %>%
  unnest(blastDbid) %>%
  separate(blastDbid, into = c("accession","taxa"), 
           sep = ";", convert = TRUE, extra = "merge")

bn2 <- bn1 %>%
  select(id, accession, taxa, blastEvalue, blastLength, blastPident, blastNident, blastScore, blastBitscore) %>%
  unnest(taxa) %>%
  separate(taxa, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), 
           sep = ",", convert = TRUE, extra = "merge")

# remove columns not needed
bn3 <- subset(bn2, select= -c(kingdom, phylum, class, order, family))

# remove unecessary characters and columns
bn3$genusBlast<-gsub("g:","",as.character(bn3$genus))
bn3$species1<-gsub("s:","",as.character(bn3$species))
bn3$speciesBlast<-gsub("_"," ",as.character(bn3$species1))

bnAll <- subset(bn3, select= -c(genus, species, species1))

#####
## Filter blastn output
#####
# keep maximum percent identities for each motu or all percent identities which are 100
bnident <- bnAll %>%
  group_by(id) %>%
  filter(blastPident == max(blastPident) | blastPident == 100) 

# keep all unique motu and species combinations
bntopscore <-bnident %>%
  group_by(id, speciesBlast) %>%
  slice_max(blastBitscore, n = 1, with_ties = FALSE)

# if there is only one species match, keep that, if multiple, assign to genus
bngroup <- bntopscore %>% group_by(id) %>% count(id)
bntaxa <- merge(bngroup, bntopscore, by.x = "id")
bntaxa$taxaBlast <- ifelse(bntaxa$n==1, bntaxa$speciesBlast, bntaxa$genusBlast)

# Remove duplicates based on id
bn_clean <- bntaxa[!duplicated(bntaxa$id), ]
# clean up columns that no longer have context
bn_clean <- subset(bn_clean, select = -c(accession, genusBlast, speciesBlast))

#####
## Merge all taxonomy and counts
#####
# Remove columns with unnecessary variables
abund <- subset(abund, select = -c(definition, ali_length, count, cut, direction, experiment, forward_match,                 
                                        forward_primer, forward_score, forward_tag, goodali, mode,
                                        reverse_match, reverse_primer, reverse_score, reverse_tag, score,                         
                                        score_norm, seq_a_deletion, seq_a_insertion, seq_a_mismatch,                
                                        seq_a_single, seq_ab_match, seq_b_deletion, seq_b_insertion,               
                                        seq_b_mismatch, seq_b_single, seq_length, seq_length_ori, seq_rank, 
                                        status, subsequence, cluster_weight))

# Join the dataframes to make raw motu table 
# put all data frames into list
df_list <- list(eco_clean, tax_clean, bn_clean, abund)
# merge all data frames in list
motu_all <- df_list %>% reduce(full_join, by="id")

# save as a csv file
#write.csv(motu_all, "C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/motu_all_CHANGE.csv")
write.csv(motu_all, "C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/motu_all_p2e.csv") # sequence run 2, elas02 library

#####
## Filter to assign final taxonomy
#####
motu_all_ident <- subset(motu_all, motu_all$best_identity>0.70)

# filter human
motu_all_ident$final_name1 <- ifelse((motu_all_ident$taxaBlast=="Homo sapiens" & motu_all_ident$order_name=="Primates") |
                                (motu_all_ident$taxaBlast=="Homo sapiens" & motu_all_ident$blastPident==100), "Homo sapiens", NA)
motu_all_ident$pid1 <- ifelse((motu_all_ident$taxaBlast=="Homo sapiens" & motu_all_ident$order_name=="Primates") |
                                (motu_all_ident$taxaBlast=="Homo sapiens" & motu_all_ident$blastPident==100), motu_all_ident$blastPident/100, NA)
motu_all_ident$rank1 <- ifelse((motu_all_ident$taxaBlast=="Homo sapiens" & motu_all_ident$order_name=="Primates") |
                                 (motu_all_ident$taxaBlast=="Homo sapiens" & motu_all_ident$blastPident==100), "species", NA)

# filter three-way consensus
# species level
motu_all_ident$final_name2 <- ifelse(motu_all_ident$scientific_name== motu_all_ident$species & motu_all_ident$species== motu_all_ident$taxaBlast,
                                     motu_all_ident$species, NA)
motu_all_ident$pid2 <- ifelse(motu_all_ident$scientific_name== motu_all_ident$species & motu_all_ident$species== motu_all_ident$taxaBlast,
                              motu_all_ident$s.id, NA)
motu_all_ident$rank2 <- ifelse(motu_all_ident$scientific_name== motu_all_ident$species & motu_all_ident$species== motu_all_ident$taxaBlast,
                              "species", NA)
# genus level
motu_all_ident$final_name3 <- ifelse(motu_all_ident$scientific_name== motu_all_ident$genus & motu_all_ident$genus== motu_all_ident$taxaBlast,
                                     motu_all_ident$genus, NA)
motu_all_ident$pid3 <- ifelse(motu_all_ident$scientific_name== motu_all_ident$genus & motu_all_ident$genus== motu_all_ident$taxaBlast,
                              motu_all_ident$best_identity, NA)
motu_all_ident$rank3 <- ifelse(motu_all_ident$scientific_name== motu_all_ident$genus & motu_all_ident$genus== motu_all_ident$taxaBlast,
                               "genus", NA)

# filter two-way consensus
# sintax and blast - species level
motu_all_ident$final_name4 <- ifelse(motu_all_ident$s.id>0.95 & motu_all_ident$species==motu_all_ident$taxaBlast,
                                     motu_all_ident$species, NA)
motu_all_ident$pid4 <- ifelse(motu_all_ident$s.id>0.95 & motu_all_ident$species==motu_all_ident$taxaBlast,
                                     motu_all_ident$s.id, NA)
motu_all_ident$rank4 <- ifelse(motu_all_ident$s.id>0.95 & motu_all_ident$species==motu_all_ident$taxaBlast,
                                     "species", NA)

# Compile filters into master taxonomy
# collapse filters above
motu_all_ident$final_name <- ifelse(is.na(motu_all_ident$final_name1)==FALSE, motu_all_ident$final_name1, 
                                    ifelse(is.na(motu_all_ident$final_name2)==FALSE, motu_all_ident$final_name2, 
                                           ifelse(is.na(motu_all_ident$final_name3)==FALSE, motu_all_ident$final_name3, 
                                                  motu_all_ident$final_name4)))
motu_all_ident$pid <- ifelse(is.na(motu_all_ident$pid1)==FALSE, motu_all_ident$pid1, 
                                    ifelse(is.na(motu_all_ident$pid2)==FALSE, motu_all_ident$pid2, 
                                           ifelse(is.na(motu_all_ident$pid3)==FALSE, motu_all_ident$pid3, 
                                                  motu_all_ident$pid4)))
motu_all_ident$final_rank <- ifelse(is.na(motu_all_ident$rank1)==FALSE, motu_all_ident$rank1, 
                             ifelse(is.na(motu_all_ident$rank2)==FALSE, motu_all_ident$rank2, 
                                    ifelse(is.na(motu_all_ident$rank3)==FALSE, motu_all_ident$rank3, 
                                           motu_all_ident$rank4)))

# any remaining NAs are filled with results from ecotag
motu_all_ident$final_name <- ifelse(is.na(motu_all_ident$final_name)==TRUE, motu_all_ident$scientific_name, motu_all_ident$final_name)
motu_all_ident$pid <- ifelse(is.na(motu_all_ident$pid)==TRUE, motu_all_ident$best_identity, motu_all_ident$pid)
motu_all_ident$final_rank <- ifelse(is.na(motu_all_ident$final_rank)==TRUE, motu_all_ident$rank, motu_all_ident$final_rank)

# save as a csv file
#write.csv(motu_all_ident, "C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/motu_all_CHANGE.csv")
write.csv(motu_all_ident, "C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/motu70_all_p2e.csv") # sequence run 2, elas02 library

#####
## Final motu table (not collapsed or decontaminated)
#####
# keep assignments 98% or higher
motu98 <- subset(motu_all_ident, motu_all_ident$pid>=0.98)
# remove unnecessary columns
motu98 <- subset(motu98, select = -c(best_match, spcies_list, scientific_name,
                                             best_identity, rank, species_name, genus_name, 
                                             family_name, order_name, qc_assign_70, kingdom,
                                             k.id, phylum, p.id, class, c.id, order, o.id, 
                                             family, f.id, genus, g.id, species, s.id, n,
                                             blastEvalue, blastLength, blastPident, blastNident, 
                                             blastScore, blastBitscore, taxaBlast, final_name1,
                                             final_name2, final_name3, final_name4, pid1, pid2,
                                             pid3, pid4, rank1, rank2, rank3, rank4))
# reorder columns
motu98 <- motu98 %>% relocate(c(final_name, pid, final_rank, total_reads, sequence), .after = id)
# save as a csv file
#write.csv(motu98, "C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/motu98_final_CHANGE.csv")
write.csv(motu98, "C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/motu98_final_p2e.csv") # sequence run 2, elas02 library

#####
# Calculate taxonomy & assignment statistics
#####

# taxonomy stats
# names of statistics calculated
stat_all <- c("total_motus", "lib_reads", NA, NA)
stat_70 <- c("motus70", "reads70", "percentmotus70", "percentreads70")
stat_98 <- c("motus98", "reads98", "percentmotus98", "percentreads98")

# value calculations
# totals
total_motus = nrow(motu_all)
lib_reads = sum(motu_all$total_reads)
# motus and data after removing everything <70%
motus70 = nrow(motu_all_ident)
reads70 = sum(motu_all_ident$total_reads)
percentmotus70 = (motus70/total_motus)*100
percentreads70 = (reads70/lib_reads)*100
# motus and data after reassigning taxa and keeping >=98%
motus98 = nrow(motu98)
reads98 = sum(motu98$total_reads)
percentmotus98 = (motus98/total_motus)*100
percentreads98 = (reads98/lib_reads)*100

# compile calculated values
value_all <- c(total_motus, lib_reads, NA, NA)
value_70 <- c(motus70, reads70, percentmotus70, percentreads70)
value_98 <- c(motus98, reads98, percentmotus98, percentreads98)

# make table
df <- data.frame(stat_all, value_all, stat_70, value_70, stat_98, value_98)

print(df)

# compare assignments
stat <- c("motus70","species_3way", "percentmotuspecies_3way", "genus_3way", "percentmotugenus_3way",
          "species_2way", "percentmotuspecies_2way", "human", "percentmotu_human",
          "ecotag", "percentmotu_ecotag")

# value calculations
# 3 way consensus
species_3way = sum(!is.na(motu_all_ident$final_name2))
percentmotuspecies_3way = (species_3way/motus70) * 100
genus_3way = sum(!is.na(motu_all_ident$final_name3))
percentmotugenus_3way = (genus_3way/motus70) * 100
# 2 way consensus
species_2way = sum(!is.na(motu_all_ident$final_name4))
percentmotuspecies_2way = (species_2way/motus70) * 100
# human
human = sum(!is.na(motu_all_ident$final_name1))
percentmotu_human = (human/motus70) * 100
# default to ecotag
ecotag = motus70 - (species_3way + genus_3way + species_2way + human)
percentmotu_ecotag = (ecotag/motus70) * 100

value <- c(motus70, species_3way, percentmotuspecies_3way, genus_3way, percentmotugenus_3way,
           species_2way, percentmotuspecies_2way, human, percentmotu_human,
           ecotag, percentmotu_ecotag)

df2 <- data.frame(stat, value)

print(df2)

#####
# clean up (optional)
#####
rm(eco, eco_v1, eco_v2, eco_v3) #ecotag processing
rm(tax, tax2, tax3, tax4, tax5) #sintax processing
rm(bn, bn1, bn2, bn3) #blastn processing
# keep bntaxa because it can be used to see blast hits that have the same score
rm(bnAll, bngroup, bnident, bntopscore) #more blastn processing
rm(abund, bn_clean, bntaxa, df_list, eco_clean, motu_all, motu_all_ident, tax_clean) #all except final table
