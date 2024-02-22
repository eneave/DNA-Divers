################################################
## Add Taxonomic metadata to aquarium samples ##
################################################

library(ggplot2)
library(tidyverse)
library(janitor)

#####
## Shark tank
#####
p2e_aq_prev <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/decontam/p2e_aq_decontam.csv")
#p2t_aq_prev <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/decontam/p2t_aq_decontam.csv")
# total reads per MOTU
p2e_aq_prev <- p2e_aq_prev %>% mutate(total_reads = rowSums(.[3:30]))

#####
# subset elasmobranchs
#####
# fix known mislabeling in taxonomy
p2e_aq_prev$final_class <- ifelse(p2e_aq_prev$final_order=="Myliobatiformes",
                                  "Elasmobranchii",
                                  p2e_aq_prev$final_class)
# subset Elasmobranchs
p2e_aq_elas <- subset(p2e_aq_prev, final_class=="Elasmobranchii")
# Remove genus-level assignments for species detected (low amount of reads)
p2e_aq_elas <- subset(p2e_aq_elas, final_name!="Dipturus" & final_name!="Hypanus")
# Remove contamination from aquarium (tide pool exhibit) or other samples
p2e_aq_elas <- subset(p2e_aq_elas, final_name!="Scyliorhinus canicula")
# convert to long dataframe
long_aq_elas <- p2e_aq_elas %>%
  pivot_longer(cols = sample.10Ae60BLUE_MPEtB:sample.9He60BLUE_MPEtA_,
               names_to = "seq_id",
               values_to = "reads")
# add metadata to long dataframe
meta <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/metadata/supp_table_1.csv")
long_aq_elas <- merge(long_aq_elas, meta, by="seq_id")

# add proportional read counts
long_aq_elas$prc <- long_aq_elas$reads/long_aq_elas$total_reads
# specify sample type more clearly
long_aq_elas$type2 <-ifelse(long_aq_elas$type=="eDNA", "Syringe Filter", 
                              ifelse(long_aq_elas$type=="MP" & (long_aq_elas$time==65|long_aq_elas$time==50), "Diver MP", "Soak MP"))


# Fix taxonomy based on inventory knowledge 
long_aq_elas$manual_name <- ifelse(long_aq_elas$final_name=="Chiloscyllium griseum", "Chiloscyllium sp.",
                                    ifelse(long_aq_elas$final_name=="Heterodontus", "Heterodontus sp.",
                                           ifelse(long_aq_elas$final_name=="Orectolobus japonicus", "Orectolobus sp.",
                                                  long_aq_elas$final_name)))

# quick stacked bar plot of sharks
windows()
ggplot(long_aq_elas , aes(x = seq_id, y = prc, fill = manual_name)) + 
  geom_bar(stat = "identity", color = "black", position = "fill") +
  facet_grid(. ~ type2, scales = "free") +
  labs(title="Proportional Read counts of Elasmobranchs \ndetected in the Ocean display by different methods",
       x ="Sample", y = "Proportional Read Counts (PRC)") +
  theme(axis.text.x = element_text(angle = 90))

#####

#####
## Subset and prepare other detections to add to the stacked bar charts
#####

p2e_aq_notelas <- subset(p2e_aq_prev, final_class!="Elasmobranchii")
# Remove domestic species (low amount of reads)
p2e_aq_notelas <- subset(p2e_aq_notelas, final_order!="Artiodactyla" & final_order!="Galliformes"
                      & final_order!="Carnivora" & final_order!="Pelecaniformes")

# NEED TO FIGURE OUT WHICH FISH ARE FOOD AND WHICH ARE TANK RESIDENTS

# Remove contamination from aquarium (tide pool exhibit) or other samples
p2e_aq_elas <- subset(p2e_aq_elas, final_name!="Scyliorhinus canicula")
# convert to long dataframe
long_aq_elas <- p2e_aq_elas %>%
  pivot_longer(cols = sample.10Ae60BLUE_MPEtB:sample.9He60BLUE_MPEtA_,
               names_to = "seq_id",
               values_to = "reads")
# add metadata to long dataframe
meta <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/metadata/supp_table_1.csv")
long_aq_elas <- merge(long_aq_elas, meta, by="seq_id")

# add proportional read counts
long_aq_elas$prc <- long_aq_elas$reads/long_aq_elas$total_reads
# specify sample type more clearly
long_aq_elas$type2 <-ifelse(long_aq_elas$type=="eDNA", "Syringe Filter", 
                            ifelse(long_aq_elas$type=="MP" & (long_aq_elas$time==65|long_aq_elas$time==50), "Diver MP", "Soak MP"))

# Fix taxonomy based on inventory knowledge 
long_aq_elas$manual_name <- ifelse(long_aq_elas$final_name=="Chiloscyllium griseum", "Chiloscyllium sp.",
                                   ifelse(long_aq_elas$final_name=="Heterodontus", "Heterodontus sp.",
                                          ifelse(long_aq_elas$final_name=="Orectolobus japonicus", "Orectolobus sp.",
                                                 long_aq_elas$final_name)))



#####
# OLD CODE BELOW HERE
#####
# explore data including all detections
#####

# subset species identity >0.95
p2e_aq_95 <- subset(p2e_aq_prev, s.id>0.95)

# subset species identity to 1, see if it differs much from >0.95
p2e_aq_100 <- subset(p2e_aq_prev, s.id==1)

##### 
#95% identity to species-level
#####
# assign taxonomy for presentation for all taxa
p2e_aq_95$manual_taxo1 <- ifelse(p2e_aq_95$species=="Chiloscyllium griseum", "Chiloscyllium sp.",
                          ifelse(p2e_aq_95$genus=="Heterodontus", "Heterodontus sp.",
                          ifelse(p2e_aq_95$genus=="Orectolobus", "Orectolobus sp.",
                          ifelse(p2e_aq_95$genus=="Scomber", "Food",
                          ifelse(p2e_aq_95$genus=="Salmo", "Food",
                          ifelse(p2e_aq_95$genus=="Sprattus", "Food",
                          ifelse(p2e_aq_95$genus=="Hippoglossus", "Food",
                          ifelse(p2e_aq_95$genus=="Pseudocaranx", "Food",
                          ifelse(p2e_aq_95$genus=="Clupea", "Food",
                          ifelse(p2e_aq_95$genus=="Gadus", "Food",
                          ifelse(p2e_aq_95$genus=="Melanogrammus", "Food",
                          ifelse(p2e_aq_95$genus=="Euthynnus", "Food",
                          ifelse(p2e_aq_95$genus=="Sardina", "Food",
                          ifelse(p2e_aq_95$species=="Carcharhinus melanopterus", "Carcharhinus melanopterus",
                          ifelse(p2e_aq_95$species=="Carcharias taurus", "Carcharias taurus",
                          ifelse(p2e_aq_95$species=="Chiloscyllium punctatum", "Chiloscyllium punctatum",
                          ifelse(p2e_aq_95$species=="Ginglymostoma cirratum", "Ginglymostoma cirratum",
                          ifelse(p2e_aq_95$species=="Glaucostegus cemiculus", "Glaucostegus cemiculus",
                          ifelse(p2e_aq_95$species=="Hypanus americanus", "Hypanus americanus",
                          ifelse(p2e_aq_95$species=="Stegostoma tigrinum", "Stegostoma tigrinum",
                                 "Teleostei"))))))))))))))))))))

# convert to long dataframe
long_aq_95 <- p2e_aq_95 %>%
  pivot_longer(cols = sample.10Ae60BLUE_MPEtB:sample.9He60BLUE_MPEtA_,
               names_to = "long_id",
               names_prefix = "sample.",
               values_to = "reads")
# add metadata to long dataframe
long_aq_95 <- merge(long_aq_95, p2e_aq_meta, by="long_id")
# remove controls 
all95samples <- subset(long_aq_95, sampletype1=="sample")
# add proportional read counts
all95samples$prc <- all95samples$reads/all95samples$total_reads
# specify sample type more clearly
all95samples$type2 <-ifelse(all95samples$type=="eDNA", "Syringe Filter", 
                             ifelse(all95samples$type=="MP" & (all95samples$time==65|all95samples$time==50), "Diver MP", "Soak MP"))
#####

##### 
#100% identity to species-level
#####
# assign taxonomy for presentation for all taxa
p2e_aq_100$manual_taxo1 <- ifelse(p2e_aq_100$species=="Chiloscyllium griseum", "Chiloscyllium sp.",
                                 ifelse(p2e_aq_100$genus=="Heterodontus", "Heterodontus sp.",
                                        ifelse(p2e_aq_100$genus=="Orectolobus", "Orectolobus sp.",
                                               ifelse(p2e_aq_100$genus=="Scomber", "Food",
                                                      ifelse(p2e_aq_100$genus=="Salmo", "Food",
                                                             ifelse(p2e_aq_100$genus=="Sprattus", "Food",
                                                                    ifelse(p2e_aq_100$genus=="Hippoglossus", "Food",
                                                                           ifelse(p2e_aq_100$genus=="Pseudocaranx", "Food",
                                                                                  ifelse(p2e_aq_100$genus=="Clupea", "Food",
                                                                                         ifelse(p2e_aq_100$genus=="Gadus", "Food",
                                                                                                ifelse(p2e_aq_100$genus=="Melanogrammus", "Food",
                                                                                                       ifelse(p2e_aq_100$genus=="Euthynnus", "Food",
                                                                                                              ifelse(p2e_aq_100$genus=="Sardina", "Food",
                                                                                                                     ifelse(p2e_aq_100$species=="Carcharhinus melanopterus", "Carcharhinus melanopterus",
                                                                                                                            ifelse(p2e_aq_100$species=="Carcharias taurus", "Carcharias taurus",
                                                                                                                                   ifelse(p2e_aq_100$species=="Chiloscyllium punctatum", "Chiloscyllium punctatum",
                                                                                                                                          ifelse(p2e_aq_100$species=="Ginglymostoma cirratum", "Ginglymostoma cirratum",
                                                                                                                                                 ifelse(p2e_aq_100$species=="Glaucostegus cemiculus", "Glaucostegus cemiculus",
                                                                                                                                                        ifelse(p2e_aq_100$species=="Hypanus americanus", "Hypanus americanus",
                                                                                                                                                               ifelse(p2e_aq_100$species=="Stegostoma tigrinum", "Stegostoma tigrinum",
                                                                                                                                                                      "Teleostei"))))))))))))))))))))

# convert to long dataframe
long_aq_100 <- p2e_aq_100 %>%
  pivot_longer(cols = sample.10Ae60BLUE_MPEtB:sample.9He60BLUE_MPEtA_,
               names_to = "long_id",
               names_prefix = "sample.",
               values_to = "reads")
# add metadata to long dataframe
long_aq_100 <- merge(long_aq_100, p2e_aq_meta, by="long_id")
# remove controls 
all100samples <- subset(long_aq_100, sampletype1=="sample")
# add proportional read counts
all100samples$prc <- all100samples$reads/all100samples$total_reads
# specify sample type more clearly
all100samples$type2 <-ifelse(all100samples$type=="eDNA", "Syringe Filter", 
                            ifelse(all100samples$type=="MP" & (all100samples$time==65|all100samples$time==50), "Diver MP", "Soak MP"))
#####

## quick plot of everything
ggplot(all100samples , aes(x = long_id, y = prc, fill = manual_taxo1)) + 
  geom_bar(stat = "identity", color = "black", position = "fill") +
  facet_grid(. ~ type2, scales = "free") +
  labs(title="Proportional Read counts of Elasmobranchs & Teleosts \ndetected in the Main Tank by different methods \n 100% s.id",
       x ="Sample", y = "Proportional Read Counts (PRC)") +
  theme(axis.text.x = element_text(angle = 90))

#####
# Funky heatmap for shark tank
#####

library(funkyheatmap)
library(dplyr, warn.conflicts = FALSE)
library(tibble, warn.conflicts = FALSE)
# example data
#data(mtcars)
#funky_heatmap(mtcars)
#funky_heatmap(data, column_info = column_info, expand = list(xmax = 4))

## Need to collapse dataframe by manual taxonomy
  hm_e_aq <- p2e_aq_100 %>% 
  group_by(manual_taxo1) %>%
  summarize(sample.8DeBLUE_eDNAA_bottle1 = sum(sample.8DeBLUE_eDNAA_bottle1),
            sample.8EeBLUE_eDNAB_bottle2 = sum(sample.8EeBLUE_eDNAB_bottle2),
            sample.8FeBLUE_eDNAC_bottle3 = sum(sample.8FeBLUE_eDNAC_bottle3),
            sample.8GeBLUE_eDNAD_bottle4 = sum(sample.8GeBLUE_eDNAD_bottle4),
            sample.8HeBLUE_MPEtA_RA = sum(sample.8HeBLUE_MPEtA_RA),
            sample.9AeBLUE_MPEtB_RB = sum(sample.9AeBLUE_MPEtB_RB),
            sample.9BeBLUE_MPEtC_DA = sum(sample.9BeBLUE_MPEtC_DA),
            sample.9CeBLUE_MPEtD_DB = sum(sample.9CeBLUE_MPEtD_DB),
            sample.10FeBLUE_MPEtA_RA = sum(sample.10FeBLUE_MPEtA_RA),
            sample.10GeBLUE_MPEtB_RB = sum(sample.10GeBLUE_MPEtB_RB),
            sample.10HeBLUE_MPEtC_DA = sum(sample.10HeBLUE_MPEtC_DA),
            sample.11AeBLUE_MPEtD_DB = sum(sample.11AeBLUE_MPEtD_DB),
            sample.9De10BLUE_MPEtA_ = sum(sample.9De10BLUE_MPEtA_),
            sample.9Ee10BLUE_MPEtB_ = sum(sample.9Ee10BLUE_MPEtB_),
            sample.9Fe30BLUE_MPEtA_ = sum(sample.9Fe30BLUE_MPEtA_),
            sample.9Ge30BLUE_MPEtB_ = sum(sample.9Ge30BLUE_MPEtB_),
            sample.9He60BLUE_MPEtA_ = sum(sample.9He60BLUE_MPEtA_),
            sample.10Ae60BLUE_MPEtB = sum(sample.10Ae60BLUE_MPEtB),
            sample.10Be120BLUE_MPEtA = sum(sample.10Be120BLUE_MPEtA),
            sample.10Ce120BLUE_MPEtB = sum(sample.10Ce120BLUE_MPEtB),
            sample.10Ee240BLUE_MPEtA = sum(sample.10Ee240BLUE_MPEtA),
            sample.10De240BLUE_MPEtB = sum(sample.10De240BLUE_MPEtB))
  
## transform dataframe for funkyheatmap
hm_e_aqt <- t(hm_e_aq)
hm_e_aqt <- hm_e_aqt %>% row_to_names(row_number = 1) %>%
            as.data.frame()
hm_e_aqt2 <- sapply(hm_e_aqt, as.numeric)  
hm_e_aqt2 <- as.data.frame(hm_e_aqt2)
hm_e_aqt2$id <- rownames(hm_e_aqt)

## Add more data for heatmap
# reads per sample
hm_e_aqt2$reads_per_sample <- rowSums(hm_e_aqt2[c(1:12)])
# total motus 
#elas_motus <- p2e_aq_100 %>% 
#  count(manual_taxo1) 

# motus per samples
## STILL NEED TO CALCULATE THIS

# proportional reads counts
## NEED TO CALCULATE THIS

######
## If columns are samples
#####
column_info <- tribble(
  ~id,     ~group,         ~name,                      ~geom,        ~palette,    ~options,
  "id",    "",             "",                         "text",       NA,          list(hjust = 0, width = 6),
  "reads_per_sample",   "overall",      "Reads",           "bar",        "palette1",  list(width = 4, legend = FALSE),
  #"cyl",   "overall",      "Number of cylinders",      "bar",        "palette2",  list(width = 4, legend = FALSE),
  "sample.8DeBLUE_eDNAA_bottle1",  "eDNA Bottle",       "1",    "funkyrect",  "palette2",  lst(),
  "sample.8EeBLUE_eDNAB_bottle2",    "eDNA Bottle",       "2",         "funkyrect",  "palette2",  lst(),
  "sample.8FeBLUE_eDNAC_bottle3",  "eDNA Bottle",       "3",          "funkyrect",  "palette2",  lst(),
  "sample.8GeBLUE_eDNAD_bottle4",    "eDNA Bottle",       "4",        "funkyrect",  "palette2",  lst(),
  "sample.8HeBLUE_MPEtA_RA",  "Dive 1",       "1a",            "funkyrect",     "palette3",  lst(),
  "sample.9AeBLUE_MPEtB_RB",    "Dive 1",       "1b",                   "funkyrect",     "palette3",  lst(),
  "sample.9BeBLUE_MPEtC_DA",    "Dive 1",       "2a",             "funkyrect",     "palette3",  lst(),
  "sample.9CeBLUE_MPEtD_DB",  "Dive 1",       "2b",          "funkyrect",     "palette3",  lst(),
  "sample.10FeBLUE_MPEtA_RA",  "Dive 2",       "1a",            "funkyrect",     "palette3",  lst(),
  "sample.10GeBLUE_MPEtB_RB",  "Dive 2",       "1b",            "funkyrect",     "palette3",  lst(),
  "sample.10HeBLUE_MPEtC_DA",  "Dive 2",       "2a",            "funkyrect",     "palette3",  lst(),
  "sample.11AeBLUE_MPEtD_DB",  "Dive 2",       "2b",            "funkyrect",     "palette3",  lst(),
  "sample.9De10BLUE_MPEtA_",  "10",       "1",    "funkyrect",  "palette2",  lst(),
  "sample.9Ee10BLUE_MPEtB_",    "10",       "2",         "funkyrect",  "palette4",  lst(),
  "sample.9Fe30BLUE_MPEtA_",  "30",       "1",          "funkyrect",  "palette4",  lst(),
  "sample.9Ge30BLUE_MPEtB_",    "30",       "2",        "funkyrect",  "palette4",  lst(),
  "sample.9He60BLUE_MPEtA_",  "60",       "1",    "funkyrect",  "palette4",  lst(),
  "sample.10Ae60BLUE_MPEtB",    "60",       "2",         "funkyrect",  "palette4",  lst(),
  "sample.10Be120BLUE_MPEtA",  "120",       "1",          "funkyrect",  "palette4",  lst(),
  "sample.10Ce120BLUE_MPEtB",    "120",       "2",        "funkyrect",  "palette4",  lst(),
  "sample.10Ee240BLUE_MPEtA",  "240",       "1",    "funkyrect",  "palette4",  lst(),
  "sample.10De240BLUE_MPEtB",    "240",       "2",         "funkyrect",  "palette4",  lst()
)
#####

######
## If columns are sharks
#####
column_info <- tribble(
  ~id,     ~group,         ~name,                      ~geom,        ~palette,    ~options,
  "id",    "",             "",                         "text",       NA,          list(hjust = 0, width = 6),
  "reads_per_sample",   "overall",      "Reads",           "bar",        "palette1",  list(width = 4, legend = FALSE),
  #"cyl",   "overall",      "Number of cylinders",      "bar",        "palette2",  list(width = 4, legend = FALSE),
  "Carcharhinus melanopterus",  "Group1",       "Carcharhinus melanopterus",    "funkyrect",  "palette2",  lst(),
  "Carcharias taurus",  "Group1",       "Carcharias taurus",    "funkyrect",  "palette2",  lst(),
  "Chiloscyllium punctatum",  "Group1",       "Chiloscyllium punctatum",    "funkyrect",  "palette2",  lst(),
  "Chiloscyllium sp.",  "Group1",       "Chiloscyllium sp.",    "funkyrect",  "palette2",  lst(),
  "Ginglymostoma cirratum",  "Group1",       "Ginglymostoma cirratum",    "funkyrect",  "palette2",  lst(),
  "Glaucostegus cemiculus",  "Group1",       "Glaucostegus cemiculus",    "funkyrect",  "palette2",  lst(),
  "Heterodontus sp.",  "Group1",       "Heterodontus sp.",    "funkyrect",  "palette2",  lst(),
  "Hypanus americanus",  "Group1",       "Hypanus americanus",    "funkyrect",  "palette2",  lst(),
  "Orectolobus sp.",  "Group1",       "Orectolobus sp.",    "funkyrect",  "palette2",  lst(),
  "Stegostoma tigrinum",  "Group1",       "Stegostoma tigrinum",    "funkyrect",  "palette2",  lst()
)

#####
# Convert shark reads to percentages
#####
## including Teleosts
hm_e_aqt2$Cm_p <- (hm_e_aqt2[c(1)]/hm_e_aqt2$reads_per_sample)*100
hm_e_aqt2$Ct_p <- (hm_e_aqt2[c(2)]/hm_e_aqt2$reads_per_sample)*100
hm_e_aqt2$Cp_p <- (hm_e_aqt2[c(3)]/hm_e_aqt2$reads_per_sample)*100
hm_e_aqt2$Cg_p <- (hm_e_aqt2[c(4)]/hm_e_aqt2$reads_per_sample)*100
hm_e_aqt2$food_p <- (hm_e_aqt2[c(5)]/hm_e_aqt2$reads_per_sample)*100
hm_e_aqt2$Gc_p <- (hm_e_aqt2[c(6)]/hm_e_aqt2$reads_per_sample)*100
hm_e_aqt2$Glc_p <- (hm_e_aqt2[c(7)]/hm_e_aqt2$reads_per_sample)*100
hm_e_aqt2$Hg_p <- (hm_e_aqt2[c(8)]/hm_e_aqt2$reads_per_sample)*100
hm_e_aqt2$Ha_p <- (hm_e_aqt2[c(9)]/hm_e_aqt2$reads_per_sample)*100
hm_e_aqt2$Og_p <- (hm_e_aqt2[c(10)]/hm_e_aqt2$reads_per_sample)*100
hm_e_aqt2$St_p <- (hm_e_aqt2[c(11)]/hm_e_aqt2$reads_per_sample)*100
hm_e_aqt2$T_p <- (hm_e_aqt2[c(12)]/hm_e_aqt2$reads_per_sample)*100

column_info <- tribble(
  ~id,     ~group,         ~name,                      ~geom,        ~palette,    ~options,
  "id",    "",             "",                         "text",       NA,          list(hjust = 0, width = 6),
  "reads_per_sample",   "overall",      "Reads",           "bar",        "palette1",  list(width = 4, legend = FALSE),
  #"cyl",   "overall",      "Number of cylinders",      "bar",        "palette2",  list(width = 4, legend = FALSE),
  "Cm_p",  "Group1",       "Carcharhinus melanopterus",    "funkyrect",  "Reds",  lst(),
  "Ct_p",  "Group1",       "Carcharias taurus",    "funkyrect",  "Reds",  lst(),
  "Cp_p",  "Group1",       "Chiloscyllium punctatum",    "funkyrect",  "Reds",  lst(),
  "Cg_p",  "Group1",       "Chiloscyllium sp.",    "funkyrect",  "Reds",  lst(),
  "Gc_p",  "Group1",       "Ginglymostoma cirratum",    "funkyrect",  "Reds",  lst(),
  "Glc_p",  "Group1",       "Glaucostegus cemiculus",    "funkyrect",  "Reds",  lst(),
  "Hg_p",  "Group1",       "Heterodontus sp.",    "funkyrect",  "Reds",  lst(),
  "Ha_p",  "Group1",       "Hypanus americanus",    "funkyrect",  "Reds",  lst(),
  "Og_p",  "Group1",       "Orectolobus sp.",    "funkyrect",  "Reds",  lst(),
  "St_p",  "Group1",       "Stegostoma tigrinum",    "funkyrect",  "Reds",  lst(),
  "T_p",  "Group1",       "Teleostei",    "funkyrect",  "Reds",  lst(),
  "food_p",  "Group1",       "Food",    "funkyrect",  "Reds",  lst()
)

## plot
windows()
funky_heatmap(hm_e_aqt2, column_info = column_info, expand = list(xmax = 4))


# calculate proportional read counts for just sharks
hm_e_aqt2$elasreads_per_sample <- rowSums(hm_e_aqt2[c(1:4,6:11)])

hm_e_aqt2$Cm_pe <- (hm_e_aqt2[c(1)]/hm_e_aqt2$elasreads_per_sample)*100
hm_e_aqt2$Ct_pe <- (hm_e_aqt2[c(2)]/hm_e_aqt2$elasreads_per_sample)*100
hm_e_aqt2$Cp_pe <- (hm_e_aqt2[c(3)]/hm_e_aqt2$elasreads_per_sample)*100
hm_e_aqt2$Cg_pe <- (hm_e_aqt2[c(4)]/hm_e_aqt2$elasreads_per_sample)*100
hm_e_aqt2$Gc_pe <- (hm_e_aqt2[c(6)]/hm_e_aqt2$elasreads_per_sample)*100
hm_e_aqt2$Glc_pe <- (hm_e_aqt2[c(7)]/hm_e_aqt2$elasreads_per_sample)*100
hm_e_aqt2$Hg_pe <- (hm_e_aqt2[c(8)]/hm_e_aqt2$elasreads_per_sample)*100
hm_e_aqt2$Ha_pe <- (hm_e_aqt2[c(9)]/hm_e_aqt2$elasreads_per_sample)*100
hm_e_aqt2$Og_pe <- (hm_e_aqt2[c(10)]/hm_e_aqt2$elasreads_per_sample)*100
hm_e_aqt2$St_pe <- (hm_e_aqt2[c(11)]/hm_e_aqt2$elasreads_per_sample)*100

hm_e_aqt3 <- hm_e_aqt2[c(13,28:37)]


#funkyheatmap

column_info2 <- tribble(
  ~id,     ~group,         ~name,                      ~geom,        ~palette,    ~options,
  "id",    "",             "",                         "text",       NA,          list(hjust = 0, width = 6),
  #"elasreads_per_sample",   "overall",      "Reads",           "bar",        "palette1",  list(width = 4, legend = FALSE),
  #"cyl",   "overall",      "Number of cylinders",      "bar",        "palette2",  list(width = 4, legend = FALSE),
  "Cm_pe",  "Group1",       "Carcharhinus melanopterus",    "funkyrect",  "stability",  lst(),
  "Ct_pe",  "Group1",       "Carcharias taurus",    "funkyrect",  "stability",  lst(),
  "Cp_pe",  "Group1",       "Chiloscyllium punctatum",    "funkyrect",  "stability",  lst(),
  "Cg_pe",  "Group1",       "Chiloscyllium sp.",    "funkyrect",  "stability",  lst(),
  "Gc_pe",  "Group1",       "Ginglymostoma cirratum",    "funkyrect",  "stability",  lst(),
  "Glc_pe",  "Group1",       "Glaucostegus cemiculus",    "funkyrect",  "stability",  lst(),
  "Hg_pe",  "Group1",       "Heterodontus sp.",    "funkyrect",  "stability",  lst(),
  "Ha_pe",  "Group1",       "Hypanus americanus",    "funkyrect",  "stability",  lst(),
  "Og_pe",  "Group1",       "Orectolobus sp.",    "funkyrect",  "stability",  lst(),
  "St_pe",  "Group1",       "Stegostoma tigrinum",    "funkyrect",  "stability",  lst(),
)

# find suitable color palettes
data("dynbenchmark_data")
palettes <- dynbenchmark_data$palettes

windows()
funky_heatmap(hm_e_aqt3, column_info = column_info2, expand = list(xmax = 4))

#####
## Bubble plot since funkyheatmaps are proving tricky
####
# bubbleplot

ggplot(elas100samples, aes(x = long_id, y = manual_taxo, size = prc)) +
  geom_point(pch = 21, fill = "black") +
  scale_size_continuous(name = "Proportional \nread counts (%)",
                        range = c(1, 6),
                        limits = c(0,1),
                        breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                        labels = c("< 1%","1-20%","21-40%","41-60%","61%-80%", "> 80%")) +
  facet_grid(. ~ type2, scales = "free") +
  theme(axis.text.x = element_text(angle = 90))

# bubbleplot, add colour

ggplot(elas100samples, aes(x = long_id, y = manual_taxo, size = prc)) +
  geom_point(pch = 21, aes(fill = prc)) +
  scale_fill_gradientn(colours = c("#fde725", "#7ad151", "#22a884", "#2a788e", "#414487", "#440154"),
                       limits = c(0,1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
                       #guide = guide_colourbar(reverse = TRUE, barwidth = 1, barheight = 20), 
                       name = "Proportional \nread counts (%)") +
  scale_size_continuous(name = "Proportional \nread counts (%)",
                        range = c(1, 6),
                        limits = c(0,1),
                        breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                        labels = c("< 1%","1-20%","21-40%","41-60%","61%-80%", "> 80%")) +
  facet_grid(. ~ type2, scales = "free") +
  theme(axis.text.x = element_text(angle = 90))

# convert to long dataframe
hm_e_aqt3_l <- hm_e_aqt3 %>%
  pivot_longer(c(2:11), names_to = "elas", values_to = "prc") 
hm_e_aqt3_l[is.na(hm_e_aqt3_l)] <- 0
hm_e_aqt3_l$prc <- rowSums(hm_e_aqt3_l[c(3:11),])

windows()
ggplot(hm_e_aqt3_l, aes(x = id, y = elas, size = prc)) +
  geom_point(pch = 21, fill = "black") +
  scale_size_continuous(name = "Proportional \nread counts (%)",
                        range = c(1, 6),
                        limits = c(0,1),
                        breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                        labels = c("< 1%","1-20%","21-40%","41-60%","61%-80%", "> 80%"))  


#####
## Coral cave tank
#####
# FUNKY HEATMAPS TUTORIAL FOR LATER 
# https://www.youtube.com/watch?v=9XbxL-Is22k

# total reads per MOTU
p2t_aq_prev <- p2t_aq_prev %>% mutate(total_reads = rowSums(.[18:30]))


# subset species identity 100%
p2t_aq_100 <- subset(p2t_aq_prev, s.id==1)









