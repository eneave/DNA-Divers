############################################################
## Process taxonomic output from ecotag, sintax, & blastn ##
############################################################
# load packages
library(splitstackshape)
library(tidyverse)

####################################
## Sequence run 2, elas02 library ##
####################################
#####
## Process ecotag output
#####
# open ecotag output fasta file notepad, copy to excel and sort by motu ids
# keep motu ids and taxonomic info, remove sequences
# save as a csv file
eco_p2e <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/e_divemeth2_SWARM1_ecotag_noseq.csv")

# separate by semi-colons into columns
eco_p2e_v1 <- eco_p2e %>% 
  unnest(id) %>%
  separate(id, c("id", "id_status", "family", "species_name", "best_match",
                 "taxid_by_db", "rank_by_db", "scientific name",
                 "match_count", "rank", "taxid", "species",
                 "order_name", "best_identity", "scientific_name_by_db",
                 "spcies_list", "genus_name", "family_name", "genus", "order"), sep = "\\;",
           convert = TRUE, extra = "merge")
eco_p2e_v1 <- eco_p2e_v1 %>% separate(id, c("id", "count"), sep = " ")

# separate the best rank scientific name
eco_p2e_v2 <- cSplit(eco_p2e_v1, c("scientific name"), sep = c("="))
eco_p2e_v2 <- cSplit(eco_p2e_v2, c("best_identity"), sep = c(":"))
eco_p2e_v2 <- cSplit(eco_p2e_v2, c("best_identity_2"), sep = c("}"))
eco_p2e_v2 <- cSplit(eco_p2e_v2, c("rank"), sep = c("="))
eco_p2e_v2 <- cSplit(eco_p2e_v2, c("species_name"), sep = c("="))
eco_p2e_v2 <- cSplit(eco_p2e_v2, c("genus_name"), sep = c("="))
eco_p2e_v2 <- cSplit(eco_p2e_v2, c("family_name"), sep = c("="))
eco_p2e_v2 <- cSplit(eco_p2e_v2, c("order_name"), sep = c("="))

# fix motu ids
eco_p2e_v2$id <-gsub(">","",as.character(eco_p2e_v2$id))

# remove unnecessary columns
eco_p2e_v3 <- subset(eco_p2e_v2, select= -c(best_identity_1,
                                            rank_1, species_name_1, genus_name_1,
                                            family_name_1, order_name_1, count, id_status,
                                            family, taxid_by_db, rank_by_db, match_count,
                                            taxid, species, scientific_name_by_db, genus, order))
eco_p2e_v3 <- subset(eco_p2e_v3, select= -c(4))

# rename columns to meaningful names
ecop2e <- eco_p2e_v3 %>% 
  rename(scientific_name = 4,
         best_identity = best_identity_2_1,
         rank = rank_2,
         species_name = species_name_2,
         genus_name = genus_name_2, 
         family_name = family_name_2,
         order_name = order_name_2)

# remove unnecessary data
rm(eco_p2e, eco_p2e_v1, eco_p2e_v2, eco_p2e_v3)

#####
## Process sintax output
#####
# Read .tsv file generated from using the sintax option in vsearch
tax <- read.csv(file = "C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/e_SWARM1_sintax_output_ALL.tsv", sep = ';', header = FALSE)

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
taxp2e <- tax5 %>% 
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
# clean up
rm(tax, tax2, tax3, tax4, tax5)

#####
## Process blastn output
#####
bn <- read.csv(file = "C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/elas02_divers2_blast_v3.tsv", sep = "", header = TRUE)

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

# clean up
rm(bn, bn1, bn2, bn3)

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
bnp2e <- bntaxa[!duplicated(bntaxa$id), ]
# clean up columns that no longer have context
bnp2e <- subset(bnp2e, select = -c(accession, genusBlast, speciesBlast))

# clean up
rm(bnAll, bngroup, bnident, bntopscore)

#####
## Merge all taxonomy and counts
#####

# Read tab file generated from obitab into R
abundp2e <- read.delim("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/e_divmeth2_SWARM1_output.counts.csv", sep=";", header=T)

# Remove columns with unnecessary variables
abundp2e <- subset(abundp2e, select = -c(definition, ali_length, count, cut, direction, experiment, forward_match,                 
                                        forward_primer, forward_score, forward_tag, goodali, mode,
                                        reverse_match, reverse_primer, reverse_score, reverse_tag, score,                         
                                        score_norm, seq_a_deletion, seq_a_insertion, seq_a_mismatch,                
                                        seq_a_single, seq_ab_match, seq_b_deletion, seq_b_insertion,               
                                        seq_b_mismatch, seq_b_single, seq_length, seq_length_ori, seq_rank, 
                                        status, subsequence, cluster_weight))

# Join the dataframes to make raw motu table 
# put all data frames into list
df_list <- list(ecop2e, taxp2e, bnp2e, abundp2e)
# merge all data frames in list
motu_p2e <- df_list %>% reduce(full_join, by="id")
# save as a csv file
write.csv(motu_p2e, "C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/motu_p2e.csv")

#####
## Filter to assign final taxonomy
#####

motu_p2e_ident <- subset(motu_p2e, motu_p2e$best_identity>0.70)

# filter human
motu_p2e_ident$final_name1 <- ifelse((motu_p2e_ident$taxaBlast=="Homo sapiens" & motu_p2e_ident$order_name=="Primates") |
                                (motu_p2e_ident$taxaBlast=="Homo sapiens" & motu_p2e_ident$blastPident==100), "Homo sapiens", NA)
motu_p2e_ident$pid1 <- ifelse((motu_p2e_ident$taxaBlast=="Homo sapiens" & motu_p2e_ident$order_name=="Primates") |
                                (motu_p2e_ident$taxaBlast=="Homo sapiens" & motu_p2e_ident$blastPident==100), motu_p2e_ident$blastPident/100, NA)
motu_p2e_ident$rank1 <- ifelse((motu_p2e_ident$taxaBlast=="Homo sapiens" & motu_p2e_ident$order_name=="Primates") |
                                 (motu_p2e_ident$taxaBlast=="Homo sapiens" & motu_p2e_ident$blastPident==100), "species", NA)

# filter three-way consensus
# species level
motu_p2e_ident$final_name2 <- ifelse(motu_p2e_ident$scientific_name== motu_p2e_ident$species & motu_p2e_ident$species== motu_p2e_ident$taxaBlast,
                                     motu_p2e_ident$species, NA)
motu_p2e_ident$pid2 <- ifelse(motu_p2e_ident$scientific_name== motu_p2e_ident$species & motu_p2e_ident$species== motu_p2e_ident$taxaBlast,
                              motu_p2e_ident$s.id, NA)
motu_p2e_ident$rank2 <- ifelse(motu_p2e_ident$scientific_name== motu_p2e_ident$species & motu_p2e_ident$species== motu_p2e_ident$taxaBlast,
                              "species", NA)
# genus level
motu_p2e_ident$final_name3 <- ifelse(motu_p2e_ident$scientific_name== motu_p2e_ident$genus & motu_p2e_ident$genus== motu_p2e_ident$taxaBlast,
                                     motu_p2e_ident$genus, NA)
motu_p2e_ident$pid3 <- ifelse(motu_p2e_ident$scientific_name== motu_p2e_ident$genus & motu_p2e_ident$genus== motu_p2e_ident$taxaBlast,
                              motu_p2e_ident$best_identity, NA)
motu_p2e_ident$rank3 <- ifelse(motu_p2e_ident$scientific_name== motu_p2e_ident$genus & motu_p2e_ident$genus== motu_p2e_ident$taxaBlast,
                               "genus", NA)

# filter two-way consensus
# sintax and blast - species level
motu_p2e_ident$final_name4 <- ifelse(motu_p2e_ident$s.id>0.95 & motu_p2e_ident$species==motu_p2e_ident$taxaBlast,
                                     motu_p2e_ident$species, NA)
motu_p2e_ident$pid4 <- ifelse(motu_p2e_ident$s.id>0.95 & motu_p2e_ident$species==motu_p2e_ident$taxaBlast,
                                     motu_p2e_ident$s.id, NA)
motu_p2e_ident$rank4 <- ifelse(motu_p2e_ident$s.id>0.95 & motu_p2e_ident$species==motu_p2e_ident$taxaBlast,
                                     "species", NA)

# Compile filters into master taxonomy
# collapse filters above
motu_p2e_ident$final_name <- ifelse(is.na(motu_p2e_ident$final_name1)==FALSE, motu_p2e_ident$final_name1, 
                                    ifelse(is.na(motu_p2e_ident$final_name2)==FALSE, motu_p2e_ident$final_name2, 
                                           ifelse(is.na(motu_p2e_ident$final_name3)==FALSE, motu_p2e_ident$final_name3, 
                                                  motu_p2e_ident$final_name4)))
motu_p2e_ident$pid <- ifelse(is.na(motu_p2e_ident$pid1)==FALSE, motu_p2e_ident$pid1, 
                                    ifelse(is.na(motu_p2e_ident$pid2)==FALSE, motu_p2e_ident$pid2, 
                                           ifelse(is.na(motu_p2e_ident$pid3)==FALSE, motu_p2e_ident$pid3, 
                                                  motu_p2e_ident$pid4)))
motu_p2e_ident$final_rank <- ifelse(is.na(motu_p2e_ident$rank1)==FALSE, motu_p2e_ident$rank1, 
                             ifelse(is.na(motu_p2e_ident$rank2)==FALSE, motu_p2e_ident$rank2, 
                                    ifelse(is.na(motu_p2e_ident$rank3)==FALSE, motu_p2e_ident$rank3, 
                                           motu_p2e_ident$rank4)))

# any remaining NAs are filled with results from ecotag
motu_p2e_ident$final_name <- ifelse(is.na(motu_p2e_ident$final_name)==TRUE, motu_p2e_ident$scientific_name, motu_p2e_ident$final_name)
motu_p2e_ident$pid <- ifelse(is.na(motu_p2e_ident$pid)==TRUE, motu_p2e_ident$best_identity, motu_p2e_ident$pid)
motu_p2e_ident$final_rank <- ifelse(is.na(motu_p2e_ident$final_rank)==TRUE, motu_p2e_ident$rank, motu_p2e_ident$final_rank)

# save as a csv file
#write.csv(motu_p2e_ident, "C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/motu_p2e_ident.csv")

#####
## Final motu table (not collapsed or decontaminated)
#####

# keep assignments 98% or higher
p2e_motu98 <- subset(motu_p2e_ident, motu_p2e_ident$pid>=0.98)
# remove unnecessary columns
p2e_motu98 <- subset(p2e_motu98, select = -c(best_match, spcies_list, scientific_name,
                                             best_identity, rank, species_name, genus_name, 
                                             family_name, order_name, qc_assign_70, kingdom,
                                             k.id, phylum, p.id, class, c.id, order, o.id, 
                                             family, f.id, genus, g.id, species, s.id, n,
                                             blastEvalue, blastLength, blastPident, blastNident, 
                                             blastScore, blastBitscore, taxaBlast, final_name1,
                                             final_name2, final_name3, final_name4, pid1, pid2,
                                             pid3, pid4, rank1, rank2, rank3, rank4))
# reorder columns
p2e_motu98 <- p2e_motu98 %>% relocate(c(final_name, pid, final_rank, total_reads, sequence), .after = id)
# save as a csv file
write.csv(p2e_motu98, "C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/p2e_motu98.csv")
# clean up
rm(abundp2e, bnp2e, bntaxa, df_list, ecop2e, motu_p2e, motu_p2e_ident, taxp2e)






