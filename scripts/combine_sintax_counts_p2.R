################################################################
## Combine sintax taxonomy w/ counts per sample (from obitab) ##
################################################################

## set working directory
setwd("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers")
#####
## wrangling taxonomy file generated from sintax
#####

## Read .tsv file generated from using the sintax option in vsearch into R

taxe <- read.csv(file = 'C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/e_SWARM1_sintax_output_miyav258.tsv', sep = ';', header = FALSE)
taxt <- read.csv(file = 'C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/t_SWARM1_sintax_output_taberletv258.tsv', sep = ';', header = FALSE)


## separate characters and values into columns
library(tidyverse)

taxe2 <- taxe %>%
          select(V1, V2, V3) %>%
          unnest(V3) %>%
          separate(V3, into = c("all_assign", "qc_assign_70"), 
                 sep = "\\+", convert = TRUE, extra = "merge")

taxe3 <- taxe2 %>%
  select(V1, V2, all_assign, qc_assign_70) %>%
  unnest(all_assign) %>%
  separate(all_assign, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), 
           sep = ",", convert = TRUE, extra = "merge")

taxt2 <- taxt %>%
  select(V1, V2, V3) %>%
  unnest(V3) %>%
  separate(V3, into = c("all_assign", "qc_assign_70"), 
           sep = "\\+", convert = TRUE, extra = "merge")

taxt3 <- taxt2 %>%
  select(V1, V2, all_assign, qc_assign_70) %>%
  unnest(all_assign) %>%
  separate(all_assign, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), 
           sep = ",", convert = TRUE, extra = "merge")

## separate percent identity from taxonomic information
library(splitstackshape)

taxe4 <- cSplit(taxe3, c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = c("("))

taxt4 <- cSplit(taxt3, c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = c("("))


## remove non-descriptive characters
taxe4$kingdom_1<-gsub("k:","",as.character(taxe4$kingdom_1))
taxe4$kingdom_2<-gsub(")","",as.character(taxe4$kingdom_2))
taxe4$phylum_1<-gsub("p:","",as.character(taxe4$phylum_1))
taxe4$phylum_2<-gsub(")","",as.character(taxe4$phylum_2))
taxe4$class_1<-gsub("c:","",as.character(taxe4$class_1))
taxe4$class_2<-gsub(")","",as.character(taxe4$class_2))
taxe4$order_1<-gsub("o:","",as.character(taxe4$order_1))
taxe4$order_2<-gsub(")","",as.character(taxe4$order_2))
taxe4$family_1<-gsub("f:","",as.character(taxe4$family_1))
taxe4$family_2<-gsub(")","",as.character(taxe4$family_2))
taxe4$genus_1<-gsub("g:","",as.character(taxe4$genus_1))
taxe4$genus_2<-gsub(")","",as.character(taxe4$genus_2))
taxe4$species_1<-gsub("s:","",as.character(taxe4$species_1))
taxe4$species_2<-gsub(")","",as.character(taxe4$species_2))
taxe4$species_1<-gsub("_"," ",as.character(taxe4$species_1))

taxt4$kingdom_1<-gsub("k:","",as.character(taxt4$kingdom_1))
taxt4$kingdom_2<-gsub(")","",as.character(taxt4$kingdom_2))
taxt4$phylum_1<-gsub("p:","",as.character(taxt4$phylum_1))
taxt4$phylum_2<-gsub(")","",as.character(taxt4$phylum_2))
taxt4$class_1<-gsub("c:","",as.character(taxt4$class_1))
taxt4$class_2<-gsub(")","",as.character(taxt4$class_2))
taxt4$order_1<-gsub("o:","",as.character(taxt4$order_1))
taxt4$order_2<-gsub(")","",as.character(taxt4$order_2))
taxt4$family_1<-gsub("f:","",as.character(taxt4$family_1))
taxt4$family_2<-gsub(")","",as.character(taxt4$family_2))
taxt4$genus_1<-gsub("g:","",as.character(taxt4$genus_1))
taxt4$genus_2<-gsub(")","",as.character(taxt4$genus_2))
taxt4$species_1<-gsub("s:","",as.character(taxt4$species_1))
taxt4$species_2<-gsub(")","",as.character(taxt4$species_2))
taxt4$species_1<-gsub("_"," ",as.character(taxt4$species_1))

## rename columns to meaningful names
taxaefinal <- taxe4 %>% 
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

taxatfinal <- taxt4 %>% 
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
#rm(tax, tax2, tax3, tax4)
## Note: taxafinal has so many sequences because it's the seeds file of the 
## SWARM output, not the counts file; seeds file was correct
## SHOULD RUN SINTAX AGAIN AND SEE WHETHER THIS MAKES A DIFFERENCE

#####
## Combine abundance file from obitab/owi_recount_swarm script 
## with wrangled taxonomy information
## remove unnecessary information
#####

## Read tab file generated from obitab into R
abunde <- read.delim("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/e_divmeth2_SWARM1_output.counts.csv", sep=";", header=T)

abundt <- read.delim("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/t_divmeth2_SWARM1_output.counts.csv", sep=";", header=T)

## Join the dataframes to make table of ASVs
asvse <- merge(x=taxaefinal,y=abunde,by="id",all.x=TRUE)

asvst <- merge(x=taxatfinal,y=abundt,by="id",all.x=TRUE)

## Remove rows with NA vales in total_reads column
asvs2e <- asvse %>% drop_na(total_reads)

asvs2t <- asvst %>% drop_na(total_reads)

## Remove columns with unnecessary variables
asvs3e <- asvs2e[,-c(18:28,67:88)]

asvs3t <- asvs2t[,-c(18:28,75:96)]

## Write new csv
write.csv(asvs3e, file="C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/divers_methods_experiment_part2/sintax_output_UKreferencedatabase/e_divmeth2_ASVs.csv")
write.csv(asvs3e, file="C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/e_divmeth2_ASVs.csv")


write.csv(asvs3t, file="C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/divers_methods_experiment_part2/sintax_output_UKreferencedatabase/t_divmeth2_ASVs.csv")
write.csv(asvs3t, file="C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/t_divmeth2_ASVs.csv")


## Remove unnecessary data
rm(asvs, asvs2,taxafinal)

