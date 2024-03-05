#############################################################
## Spotting and displaying species of conservation concern ##
#############################################################

library(tidyverse)

#####
# load data  
#####
# data that has 0.001% read counts removed
# tele02 sequence run 1
p1_prev <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/decontam/p1_decontam_001.csv")
# elas02 ocean display sequence run 2
p2e_aq_prev <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/decontam/p2e_aq_decontam_001.csv")
# tele02 UK sequence run 2
p2t_uk_prev <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/decontam/p2t_uk_decontam_001.csv")
# elas02 Orkney sequence run 2
p2e_uk_prev <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/decontam/p2e_uk_decontam_001.csv")

# sample metadata
meta <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/metadata/supp_table_1.csv")

#####
# load and clean data that has not has 0.001% reads removed
#####
# need to remove 0.001% read counts
# tele02 coral cave sequence run 2
p2t_aq_prev <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/decontam/p2t_aq_decontam.csv")
all_reads <- p2t_aq_prev %>%
  mutate(sum = rowSums(across(c(sample.5Et10BLUE_MPEtA_:sample.8Ft_EBJun_13extblank)))) 
total <- sum(all_reads$sum)
# calculate 0.001% of reads
remove <- total*0.00001
rm(all_reads)
# replace 0.001% of reads with zero
var <- subset(meta$seq_id, meta$sequence.run==2 & meta$primer=="tele02" & meta$location_abbreviation=="BLUE")
p2t_aq_prev[,var][p2t_aq_prev[,var] <= 4] <- 0
# total reads per MOTU
p2t_aq_prev <- p2t_aq_prev %>% mutate(total_reads = rowSums(.[12:25]))
p2t_aq_prev <- subset(p2t_aq_prev, p2t_aq_prev$total_reads>0)
# save
write.csv(p2t_aq_prev, "C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/decontam/p2t_aq_decontam_001.csv")


#  elas02 sequence run 3
p3_prev <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/decontam/p3_decontam.csv")
all_reads <- p3_prev %>%
  mutate(sum = rowSums(across(c(sample.10A_SAtra1001:sample.9H_SA_SL7F2A_c))))
total <- sum(all_reads$sum) - 525000 #human assignments estimate
# calculate 0.001% of reads
remove <- total*0.00001
rm(all_reads)
# replace 0.001% of reads with zero
var <- subset(meta$seq_id, meta$sequence.run==3)
p3_prev[,var][p3_prev[,var] <= 2] <- 0
# total reads per MOTU
p3_prev <- p3_prev %>% mutate(total_reads = rowSums(.[12:54]))
p3_prev <- subset(p3_prev, p3_prev$total_reads>0)
# save
write.csv(p3_prev, "C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/decontam/p3_decontam_001.csv")

# tele02 sequence run 4
p4_prev <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/decontam/p4_decontam.csv")
# add samples back
p4_prev["sample.2B_PLYM_Fi_Jul22"] <- 0 # sample has no reads
p4_prev["sample.3A_B_CLASNOR_Jul22"] <- 0 # sample has no reads
p4_prev["sample.4A_B_CLASNOR_Jul22"] <- 0 # sample (bead preserved) has no reads
p4_prev["sample.5A_B_CLAS_NOR_Jul22"] <- 0 # sample has no reads
p4_prev["sample.6A_B_CLAS_NOR_Jul22"] <- 0 # sample has no reads
p4_prev["sample.7C_PLYM_Fi_Jul22"] <- 0 # sample has no reads
p4_prev["sample.8A_PLYM_Fi_Jul22"] <- 0 # sample has no reads
all_reads <- p4_prev %>%
  mutate(sum = rowSums(across(c(sample.1A_B_CLASNOR_Jul22:sample.pcr_positive_2_Jul22))))
total <- sum(all_reads$sum) - 680000 #add an estimate of human reads
# calculate 0.001% of reads
remove <- total*0.00001
rm(all_reads)
# replace 0.001% of reads with zero
var <- subset(meta$seq_id, meta$sequence.run==4)
p4_prev[,var][p4_prev[,var] <= 2] <- 0
# total reads per MOTU
p4_prev <- p4_prev %>% mutate(total_reads = rowSums(.[12:52]))
p4_prev <- subset(p4_prev, p4_prev$total_reads>0)
# save
write.csv(p4_prev, "C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/decontam/p4_decontam_001.csv")

