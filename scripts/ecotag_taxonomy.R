###############################################
## Ecotag output with custom-ecopcr database ##
###############################################

library(splitstackshape)
library(tidyverse)

# open ecotag output fasta file in excel and sort by motu ids
# keep motu ids and taxonomic info, remove sequences
# save as a csv file

eco_p2e <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/e_divemeth2_SWARM1_nohuman.ecotag_noseq.csv")

#eco_p2e_v1 <- as.data.frame(eco_p2e)

# separate by semi-colon
eco_p2e_v1 <- eco_p2e %>% 
                  unnest(id) %>%
                  separate(id, c("id", "size", "count", "id_status",
                                  "family", "species_name", "best_match",
                                  "taxid_by_db", "rank_by_db", "scientific name",
                                  "match_count", "rank", "taxid", "species",
                                  "order_name", "best_identity", "scientific_name_by_db",
                                  "spcies_list", "genus_name", "family_name", "genus", "order"), sep = "\\;",
                                  convert = TRUE, extra = "merge")

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

# remove unnecessary rows
eco_p2e_v3 <- subset(eco_p2e_v2, select= -c(best_identity_1,
                            rank_1, species_name_1, genus_name_1,
                            family_name_1, order_name_1))
eco_p2e_v3 <- subset(eco_p2e_v3, select= -c(3,5,11,14,15,16))

# rename columns to meaningful names
ecop2e <- eco_p2e_v3 %>% 
  rename(scientific_name = 11,
    best_identity = best_identity_2_1,
    rank = rank_2,
    species_name = species_name_2,
    genus_name = genus_name_2, 
    family_name = family_name_2,
    order_name = order_name_2)

# remove unnecessary data
rm(eco_p2e, eco_p2e_v1, eco_p2e_v2, eco_p2e_v3)

# combine ecotag and sintax taxanomy
# BE CAREFUL, THIS SINTX FILE IS NOT UPDATED WITH THE NEW v258 UK reference database/ MADE WITH SEED FILE
# ALSO UK REFERENCE DATABASE USED ON SAMPLES THAT ARE ALSO FROM THE AQUARIUM
p2e_uk <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/e_divmeth2_ASVs.csv") 

# Join the dataframes to make a master motu table
p2e_merge <- merge(x=ecop2e, y=p2e_uk ,by="id", all.x=TRUE)

write.csv(p2e_merge, "C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/p2_e_masterMOTU.csv")
## NEED TO THINK OF HOW TO MERGE WITH THE DIFFERENT SUBSETS OF REFERENCE DATABASES FOR AQUARIUM SAMPLES

