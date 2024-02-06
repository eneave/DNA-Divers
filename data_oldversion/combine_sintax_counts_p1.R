################################################################
## Combine sintax taxonomy w/ counts per sample (from obitab) ##
################################################################

## set working directory
setwd("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers")
#####
## wrangling taxonomy file generated from sintax
#####

## Read .tsv file generated from using the sintax option in vsearch into R

tax <- read.csv(file = 'C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/SWARM1_sintax_output_v258.tsv', sep = ';', header = FALSE)

## separate characters and values into columns
library(tidyverse)

tax2 <- tax %>%
          select(V1, V2, V3) %>%
          unnest(V3) %>%
          separate(V3, into = c("all_assign", "qc_assign_70"), 
                 sep = "\\+", convert = TRUE, extra = "merge")

tax3 <- tax2 %>%
  select(V1, V2, all_assign, qc_assign_70) %>%
  unnest(all_assign) %>%
  separate(all_assign, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), 
           sep = ",", convert = TRUE, extra = "merge")

## separate percent identity from taxonomic information
library(splitstackshape)

tax4 <- cSplit(tax3, c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = c("("))

## remove non-descriptive characters
tax4$kingdom_1<-gsub("k:","",as.character(tax4$kingdom_1))
tax4$kingdom_2<-gsub(")","",as.character(tax4$kingdom_2))
tax4$phylum_1<-gsub("p:","",as.character(tax4$phylum_1))
tax4$phylum_2<-gsub(")","",as.character(tax4$phylum_2))
tax4$class_1<-gsub("c:","",as.character(tax4$class_1))
tax4$class_2<-gsub(")","",as.character(tax4$class_2))
tax4$order_1<-gsub("o:","",as.character(tax4$order_1))
tax4$order_2<-gsub(")","",as.character(tax4$order_2))
tax4$family_1<-gsub("f:","",as.character(tax4$family_1))
tax4$family_2<-gsub(")","",as.character(tax4$family_2))
tax4$genus_1<-gsub("g:","",as.character(tax4$genus_1))
tax4$genus_2<-gsub(")","",as.character(tax4$genus_2))
tax4$species_1<-gsub("s:","",as.character(tax4$species_1))
tax4$species_2<-gsub(")","",as.character(tax4$species_2))

tax4$species_1<-gsub("_"," ",as.character(tax4$species_1))

## rename columns to meaningful names
taxafinal <- tax4 %>% 
                rename(
                        id = V1,
                        size = V2,
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
## clean up
rm(tax, tax2, tax3, tax4)
#####
## Combine abundance file from obitab/owi_recount_swarm script 
## with wrangled taxonomy information
## remove unnecessary information
#####

## Read tab file generated from obitab into R
abund <- read.delim("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/divmeth1_SWARM1_output.counts.csv", sep=";", header=T)

## Join the dataframes to make table of ASVs
asvs <- merge(x=taxafinal,y=abund,by="id",all.x=TRUE)

## Remove rows with NA vales in total_reads column
asvs2 <- asvs %>% drop_na(total_reads)

## Remove columns with unnecessary variables
asvs3 <- asvs2[,-c(18:28,108:129)]

## Write new csv
write.csv(asvs3, file="C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/divers_methods_experiment/sintax_taxonomy/divmeth1_ASVs.csv")
write.csv(asvs3, file="C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/divmeth1_ASVs.csv")

## Remove unnecessary data
rm(asvs, asvs2,taxafinal)

