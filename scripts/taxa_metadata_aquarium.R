################################################
## Add Taxonomic metadata to aquarium samples ##
################################################

library(ggplot2)
library(tidyverse)
library(janitor)

#####
## Shark tank
#####
# total reads per MOTU
p2e_aq_prev <- p2e_aq_prev %>% mutate(total_reads = rowSums(.[2:28]))

#####
# subset elasmobranchs
#####
p2e_aq_elas <- subset(p2e_aq_prev, class=="Elasmobranchii")

##### 
#95% identity to species-level
#####
# subset species identity >0.95
p2e_aq_elas95 <- subset(p2e_aq_elas, s.id>0.95)

# assign best-match taxonomy for Elasmobranchs
p2e_aq_elas95$manual_taxo <- ifelse(p2e_aq_elas95$species=="Chiloscyllium griseum", "Chiloscyllium sp.",
                                    ifelse(p2e_aq_elas95$genus=="Heterodontus", "Heterodontus sp.",
                                           ifelse(p2e_aq_elas95$genus=="Orectolobus", "Orectolobus sp.",
                                                  p2e_aq_elas95$species)))
p2e_aq_elas95$manual_taxo.id <- ifelse(p2e_aq_elas95$species=="Chiloscyllium griseum", p2e_aq_elas95$g.id,
                                       ifelse(p2e_aq_elas95$genus=="Heterodontus", p2e_aq_elas95$g.id,
                                              ifelse(p2e_aq_elas95$genus=="Orectolobus", p2e_aq_elas95$g.id,
                                                     p2e_aq_elas95$s.id)))

# convert to long dataframe
long_aq_elas95 <- p2e_aq_elas95 %>%
  pivot_longer(cols = sample.10Ae60BLUE_MPEtB:sample.9He60BLUE_MPEtA_,
               names_to = "long_id",
               names_prefix = "sample.",
               values_to = "reads")
# add metadata to long dataframe
long_aq_elas95 <- merge(long_aq_elas95, p2e_aq_meta, by="long_id")
# remove controls 
elas95samples <- subset(long_aq_elas95, sampletype1=="sample")
# add proportional read counts
elas95samples$prc <- elas95samples$reads/elas95samples$total_reads
# specify sample type more clearly
elas95samples$type2 <-ifelse(elas95samples$type=="eDNA", "Syringe Filter", 
                              ifelse(elas95samples$type=="MP" & (elas95samples$time==65|elas95samples$time==50), "Diver MP", "Soak MP"))
#####

##### 
#100% identity to species level
#####
# subset species identity = 1
p2e_aq_elas100 <- subset(p2e_aq_elas, s.id==1)

# assign best-match taxonomy for Elasmobranchs
p2e_aq_elas100$manual_taxo <- ifelse(p2e_aq_elas100$species=="Chiloscyllium griseum", "Chiloscyllium sp.",
                                    ifelse(p2e_aq_elas100$genus=="Heterodontus", "Heterodontus sp.",
                                           ifelse(p2e_aq_elas100$genus=="Orectolobus", "Orectolobus sp.",
                                                  p2e_aq_elas100$species)))
p2e_aq_elas100$manual_taxo.id <- ifelse(p2e_aq_elas100$species=="Chiloscyllium griseum", p2e_aq_elas100$g.id,
                                       ifelse(p2e_aq_elas100$genus=="Heterodontus", p2e_aq_elas100$g.id,
                                              ifelse(p2e_aq_elas100$genus=="Orectolobus", p2e_aq_elas100$g.id,
                                                     p2e_aq_elas100$s.id)))

# convert to long dataframe
long_aq_elas100 <- p2e_aq_elas100 %>%
  pivot_longer(cols = sample.10Ae60BLUE_MPEtB:sample.9He60BLUE_MPEtA_,
               names_to = "long_id",
               names_prefix = "sample.",
               values_to = "reads")
# add metadata to long dataframe
long_aq_elas100 <- merge(long_aq_elas100, p2e_aq_meta, by="long_id")
# remove controls 
elas100samples <- subset(long_aq_elas100, sampletype1=="sample")
# add proportional read counts
elas100samples$prc <- elas100samples$reads/elas100samples$total_reads
# specify sample type more clearly
elas100samples$type2 <-ifelse(elas100samples$type=="eDNA", "Syringe Filter", 
                             ifelse(elas100samples$type=="MP" & (elas100samples$time==65|elas100samples$time==50), "Diver MP", "Soak MP"))
#####

## quick stacked bar plot of sharks
ggplot(elas100samples , aes(x = long_id, y = prc, fill = manual_taxo)) + 
  geom_bar(stat = "identity", color = "black", position = "fill") +
  facet_grid(. ~ type2, scales = "free") +
  labs(title="Proportional Read counts of Elasmobranchs \ndetected in the Main Tank by different methods \n100% s.id",
       x ="Sample", y = "Proportional Read Counts (PRC)") +
  theme(axis.text.x = element_text(angle = 90))

#####
# explore data including all detentions
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


##
windows()
funky_heatmap(hm_e_aqt2, column_info = column_info, expand = list(xmax = 4))


#####
## Coral cave tank
#####
# FUNKY HEATMAPS TUTORIAL FOR LATER 
# https://www.youtube.com/watch?v=9XbxL-Is22k

# total reads per MOTU
p2t_aq_prev <- p2t_aq_prev %>% mutate(total_reads = rowSums(.[18:30]))


# subset species identity 100%
p2t_aq_100 <- subset(p2t_aq_prev, s.id==1)









