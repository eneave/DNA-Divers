################################################################
## Combine sintax taxonomy w/ counts per sample (from obitab) ##
################################################################

## set working directory
setwd("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/divers_methods_experiment/sintax_taxonomy")

#####
## wrangling taxonomy file generated from sintax
#####

## Read .tsv file generated from using the sintax option in vsearch into R

tax <- read.csv(file = 'SWARM1_sintax_v2_26_1_v258.tsv', sep = ';', header = FALSE)

# previous file
taxo <- read.csv(file = 'SWARM1_sintax_output.tsv', sep = ';', header = FALSE)


## separate characters and values into columns
library(tidyverse)
library(splitstackshape)

tax2 <- tax %>%
          #select(V1, V2, V3) %>%
          unnest(V1) %>%
          separate(V1, into = c("all_assign", "qc_assign_70"), 
                 sep = "\\+", convert = TRUE, extra = "merge")

tax3 <- tax2 %>% separate(all_assign, c("id", "assign"), sep = 14)


tax4 <- tax3 %>%
  select(id, assign, qc_assign_70) %>%
  unnest(assign) %>%
  separate(assign, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), 
           sep = ",", convert = TRUE, extra = "merge")

## separate percent identity from taxonomic information

tax5 <- cSplit(tax4, c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = c("("))

## remove non-descriptive characters
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

## rename columns to meaningful names
taxafinal <- tax5 %>% 
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
## clean up
rm(tax, tax2, tax3, tax4, tax5)
## Note: taxafinal has so many sequences because it's the seeds file of the 
## SWARM output, not the counts file; seeds file was correct
## SHOULD RUN SINTAX AGAIN AND SEE WHETHER THIS MAKES A DIFFERENCE

#####
## Combine abundance file from obitab/owi_recount_swarm script 
## with wrangled taxonomy information
## remove unnecessary information
#####

## Read tab file generated from obitab into R
abund <- read.delim("divmeth1_SWARM1_output.counts.csv", sep=";", header=T)

## Join the dataframes to make table of ASVs
asvs <- merge(x=taxafinal,y=abund,by="id",all.x=TRUE)

## Remove columns with unnecessary variables
asvs2 <- asvs[,-c(17:27,107:128)]

## Write new csv
write.csv(asvs2, file="C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/divers_methods_experiment/sintax_taxonomy/divmeth1_ASVs_FINAL.csv")

## Remove unnecessary data
rm(asvs, asvs2,taxafinal)

