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
# load and clean data that has not has 0.001% reads removed; in prep for this plot or other supplement figures
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


#####
# Prep data for CCA and presence-absence IUCN vulnerable species
#####

# make motu tables long
# sequence run 1 tele02
l1 <- p1_prev %>% pivot_longer(cols = sample.10A:sample.3A,
                                 names_to = "sample",
                                 names_prefix = "sample.",
                                 values_to = "reads")
l1 <- subset(l1, select = -c(X.1,X))
# UK elas02
l2 <- p2e_uk_prev %>% pivot_longer(cols = sample.11Be_EBMay_13extblank:sample.8BeORK_MPbC_lisa,
                                 names_to = "sample",
                                 names_prefix = "sample.",
                                 values_to = "reads")
l2 <- subset(l2, select = -c(X.1,X))
# UK tele02
l3 <- p2t_uk_prev %>% pivot_longer(cols = sample.4AtLIV_eDNA_FBblank:sample.8HtLIV_FBEtblank,
                                 names_to = "sample",
                                 names_prefix = "sample.",
                                 values_to = "reads")
l3 <- subset(l3, select = -c(X.1,X))
# sequence run 3 elas02
l4 <- p3_prev %>% pivot_longer(cols = sample.10A_SAtra1001:sample.9H_SA_SL7F2A_c,
                                 names_to = "sample",
                                 names_prefix = "sample.",
                                 values_to = "reads")
l4 <- subset(l4, select = -X)
# sequence run 4 tele02
l5 <- p4_prev %>% pivot_longer(cols = sample.1A_B_CLASNOR_Jul22:sample.8A_PLYM_Fi_Jul22,
                               names_to = "sample",
                               names_prefix = "sample.",
                               values_to = "reads")
l5 <- subset(l5, select = -X)

# combine long datafame
l6 <- rbind(l1,l2,l3,l4,l5)
rm(l1,l2,l3,l4,l5)

# wide motu table
w1 <- l6 %>% pivot_wider(names_from = "sample",
                         names_prefix = "sample.",
                         values_from = "reads",
                         values_fill = 0)


# Refine list of samples of relevant for these data visualizations by sample type
list <- subset(meta$seq_id, (meta$via=="snorkel" | meta$via=="SCUBA") & meta$location_abbreviation!="BLUE")

# Master MOTU table for real-world MP samples from divers/snorkelers
all_motu <- w1 %>% select("id","final_name","final_genus","final_order","final_class", 
                                "pid","final_rank","method_assign","total_reads","sequence",
                                all_of(list)) 

all_motu <- all_motu %>%
  group_by(final_rank,final_class,final_order,final_genus,final_name) %>%
  summarise(across(c(total_reads, sample.1A:sample.6B_B_CLAS_NOR_Jul22), sum))

all_motu <- all_motu %>% filter(total_reads > 0)

#####
# Fix MOTU table taxonomy
#####

# Remove contaminanats
# Dama dama; not possible in Orkney
all_motu <- subset(all_motu, all_motu$final_name!="Dama dama")
# Bos mutus; not possible in Orkney
all_motu <- subset(all_motu, all_motu$final_name!="Bos mutus")
# Equus asinus; not possible in Liverpool
all_motu <- subset(all_motu, all_motu$final_name!="Equus asinus")
# Carcharias taurus; not possible in Newcastle
all_motu <- subset(all_motu, all_motu$final_name!="Carcharias taurus")
# Galeus melastomus; unlikely for Dorset seagrass beds
all_motu <- subset(all_motu, all_motu$final_name!="Galeus melastomus")
# Carcharhinus melanopterus; not in Orkney
all_motu <- subset(all_motu, all_motu$final_name!="Carcharhinus melanopterus")
# Columba palumbus and Columba livia and other birds; needs to be changed to Aves
all_motu$final_class <- ifelse(all_motu$final_name=="Columba palumbus" | all_motu$final_name=="Columba livia" |
                           all_motu$final_name=="Corvus" | all_motu$final_name=="Gavia stellata", "Aves", all_motu$final_class)
# Rhynchobatus djiddensis
all_motu$final_class <- ifelse(all_motu$final_name=="Rhynchobatus djiddensis", "Elasmobranchii", all_motu$final_class) 
# Caretta caretta - change class, currently with fish
all_motu$final_class <- ifelse(all_motu$final_name=="Caretta caretta", "Reptilia", all_motu$final_class) 
# fix mammals
all_motu$final_class <- ifelse(all_motu$final_name=="Equus", "Mammalia", all_motu$final_class)
# Chimaera monstrosa; unlikely
all_motu <- subset(all_motu, all_motu$final_name!="Chimaera monstrosa")
# fix synonym problem Taeniura meyeni
all_motu$final_name <- ifelse(all_motu$final_name=="Taeniura meyeni", "Taeniurops meyeni", all_motu$final_name)
all_motu$final_genus <- ifelse(all_motu$final_genus=="Taeniura", "Taeniurops", all_motu$final_genus)
# Pangasianodon hypophthalmus; positive control
all_motu <- subset(all_motu, all_motu$final_name!="Pangasianodon hypophthalmus")
all_motu <- subset(all_motu, all_motu$final_name!="Hypanus americanus")
# Gymnothorax kidako should not be present in the Newcastle sample
all_motu[all_motu$final_name=="Gymnothorax kidako", "sample.7BtNEW_MPEt2_unis"] <- 0


# Subset genus and species level assignments
final_motu <- subset(all_motu, all_motu$final_rank=="genus" | all_motu$final_rank=="species")
# Might want to collapse further by species, genus, and class...
# excluding order since different notation is causing genus detections to not group
final_motu <- final_motu %>%
  group_by(final_rank,final_class,final_genus,final_name) %>%
  summarise(across(c(total_reads, sample.1A:sample.6B_B_CLAS_NOR_Jul22), sum))
# recalculate total reads column, to be sure...
library(plyr)
final_motu <- final_motu %>% mutate(total_reads = rowSums(.[6:151]))
final_motu <- subset(final_motu, final_motu$total_reads>0)


## Make different plots for different purposes
# remove domestic species
final_motu <- subset(final_motu, final_motu$final_name!="Equus")
final_motu <- subset(final_motu, final_motu$final_name!="Canis")
final_motu <- subset(final_motu, final_motu$final_name!="Ovis")
final_motu <- subset(final_motu, final_motu$final_name!="Felis catus")
final_motu <- subset(final_motu, final_motu$final_name!="Sus scrofa")

# remove human
final_motu <- subset(final_motu, final_motu$final_name!="Homo sapiens")
final_motu <- subset(final_motu, final_motu$final_name!="Pan troglodytes")


write.csv(final_motu, "C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/final_motu_mp_nature_cca.csv")

#####
## Prepare data for plot
#####

# iucn information
stat <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/metadata/supp_table_3.csv")

# subset species
sp_motu <- subset(final_motu, final_motu$final_rank=="species")
sp_motu <- merge(sp_motu, stat, by="final_name")

# subset vulnerable species
vsp_motu <- subset(sp_motu, sp_motu$iucn!="LC")
vsp_motu <- subset(vsp_motu, vsp_motu$iucn!="NE")

# make long dataframe
l_vsp <- vsp_motu %>% pivot_longer(cols = sample.1A:sample.6B_B_CLAS_NOR_Jul22,
                               names_to = "seq_id",
                               values_to = "reads")
l_vsp <- merge(l_vsp, meta, by.x="seq_id")

# collapse by country
lc_vsp <- l_vsp %>%
  group_by(final_class, final_name, country, iucn) %>%
  summarise(sum_read = sum(reads))


#####
## Bubble plot 
#####
lc_fig <- lc_vsp

# make adjustments for the figure

# remove 6 reads of Atlantic salmon in CA, most likely food
lc_fig[lc_fig$final_name=="Salmo salar" & lc_fig$country=="USA", "sum_read"] <- 0

# remove Atlantic horse mackeral because it's not regionally threatened
lc_fig <- subset(lc_fig, lc_fig$final_name!="Trachurus trachurus")

# remove Jordan
lc_fig <- subset(lc_fig, lc_fig$country!="Jordan")

# transform read counts
lc_fig$log_reads <- log(lc_fig$sum_read)
lc_fig$sqrt_reads <- sqrt(lc_fig$sum_read)

# NT "#ADFF2F"
# new facet label names for classes
class.labs <- c("Teleosts", "Elasmobranchs", "M", "R")
names(class.labs) <- c("Actinopterygii", "Elasmobranchii", "Mammalia", "Reptilia")

bubble <-
ggplot(lc_fig) + 
  geom_point(aes(x = country, y = final_name,fill = iucn, size = sqrt_reads), pch = 21) +
  scale_fill_manual(values = c(CR = "#FF0000", DD = "#808080", 
                               EN = "#FFA500", NT = "#CDFF3B", VU = "#FFFF00"),
                    breaks = c("CR", "EN", "VU", "NT", "DD"),
                    name = "IUCN") +
  scale_size_continuous(name = expression(sqrt("reads")),
                        range = c(3, 9),
                        limits = c(1,200),
                        breaks = c(1, 60, 120, 180, 200),
                        labels = c("1","60","120","180", "200")) + 
  guides(fill = guide_legend(override.aes = list(size = 4))) +
  labs(x = "", y ="") + 
  scale_y_discrete(limits = rev) +
  #facet_grid(. ~ final_class, scales = "free", space = "free", labeller = labeller(bioregion_nmds = region.labs)) +
  facet_grid(final_class ~., scales = "free", space = "free", labeller = labeller(final_class = class.labs)) +
  theme_light() +
  theme(#panel.grid.major = element_line(colour = "white"),
        #panel.grid.minor = element_line(colour = "white"), 
        axis.ticks = element_line(colour = "#000000"),
        strip.text.y = element_text(colour = "#000000"),
        axis.text.x = element_text(colour = "#000000", angle = 35, hjust = 1),
        axis.text.y = element_text(colour = "#000000", face = "italic"),
        text = element_text(size = 14),
        legend.direction = "vertical", 
        legend.box = "vertical") 


ggsave(filename=c("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_writing/figures/bubble.jpg"), 
       plot = bubble, width = 6, height = 7.5, units = "in")


