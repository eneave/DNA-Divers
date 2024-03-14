##############################################
## Compare eDNA bottles to diver metaprobes ##
##############################################
library(plyr)
library(tidyverse)

#####
# load data
#####
# tele02 Orkney sequence run 1
p1_prev <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/decontam/p1_decontam.csv")
# elas02 shark tank sequence run 2
p2e_aq_prev <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/decontam/p2e_aq_decontam_001.csv")
# tele02 Liverpool sequence run 2
p2t_uk_prev <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/decontam/p2t_uk_decontam.csv")
# elas02 Orkney sequence run 2
p2e_uk_prev <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/decontam/p2e_uk_decontam.csv")
# sample metadata
meta <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/metadata/supp_table_1.csv")

#####
# extra decontamination step for data 0.001% of reads, read removal
#####
# calculate total reads and number of reads to remove

# tele02 run 1
# add samples back 
p1_prev["sample.8C"] <- 0 # sample (bead preserved) has no reads
p1_prev["sample.3A"] <- 0 # sample has no reads
all_reads <- p1_prev %>%
  mutate(sum = rowSums(across(c(sample.10A:sample.8C))))
total <- sum(all_reads$sum)
# calculate 0.001% of reads
remove <- total*0.00001
rm(all_reads)
# replace 0.001% of reads with zero
var <- subset(meta$seq_id, meta$sequence.run==1)
p1_prev[,var][p1_prev[,var] <= 6] <- 0
# total reads per MOTU
p1_prev <- p1_prev %>% mutate(total_reads = rowSums(.[12:92]))
p1_prev <- subset(p1_prev, p1_prev$total_reads>0)
# save
write.csv(p1_prev, "C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/decontam/p1_decontam_001.csv")


# elas02 shark tank sequence run 2
# cleaned in taxa_metadata_aquarium.R


# tele02 sequence run 2
all_reads <- p2t_uk_prev %>%
  mutate(sum = rowSums(across(c(sample.4AtLIV_eDNA_FBblank:sample.8HtLIV_FBEtblank))))
total <- sum(all_reads$sum)
# calculate 0.001% of reads
remove <- total*0.00001
rm(all_reads)
# replace 0.001% of reads with zero
var <- subset(meta$seq_id, meta$sequence.run==2 & meta$primer=="tele02" & meta$location_abbreviation!="BLUE")
p2t_uk_prev[,var][p2t_uk_prev[,var] <= 8] <- 0
# total reads per MOTU
p2t_uk_prev <- p2t_uk_prev %>% mutate(total_reads = rowSums(.[3:29]))
p2t_uk_prev <- subset(p2t_uk_prev, p2t_uk_prev$total_reads>0)
# save
write.csv(p2t_uk_prev, "C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/decontam/p2t_uk_decontam_001.csv")


# elas02 run 2
all_reads <- p2e_uk_prev %>%
  mutate(sum = rowSums(across(c(sample.11Be_EBMay_13extblank:sample.8BeORK_MPbC_lisa))))
total <- sum(all_reads$sum)
# calculate 0.001% of reads
remove <- total*0.00001
rm(all_reads)
# replace 0.001% of reads with zero
var <- subset(meta$seq_id, meta$sequence.run==2 & meta$primer=="elas02" & meta$location_abbreviation!="BLUE")
p2e_uk_prev[,var][p2e_uk_prev[,var] <= 2] <- 0
# total reads per MOTU
p2e_uk_prev <- p2e_uk_prev %>% mutate(total_reads = rowSums(.[3:15]))
p2e_uk_prev <- subset(p2e_uk_prev, p2e_uk_prev$total_reads>0)
# save
write.csv(p2e_uk_prev, "C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/decontam/p2e_uk_decontam_001.csv")


# elas02 run 3
p3_prev <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/decontam/p3_decontam.csv")
# replace 0.001% of reads with zero
var <- subset(meta$seq_id, meta$sequence.run==3)
p3_prev[,var][p3_prev[,var] <= 2] <- 0  ## match elas02 run 2 p2e_uk_prev because only using one sample
# total reads per MOTU
p3_prev <- p3_prev %>% mutate(total_reads = rowSums(.[12:54]))
p3_prev <- subset(p3_prev, p3_prev$total_reads>0)

# subset sample that is needed from run 3
samp <- p3_prev[c("id", "final_name","final_genus","final_order","final_class", 
                  "pid","final_rank","method_assign","total_reads","sequence",
                  "sample.12G_ORKMPbA_kurt")]


# extract data for this visualization
# Refine list of samples of relevant for these data visualizations by sample type
orkt <- subset(meta$seq_id, meta$site.name=="SMS Bayern, Scapa Flow, Orkney" & (meta$type=="eDNA" | (meta$extraction=="QBT" & meta$perservation=="Et")))
blue <- subset(meta$seq_id, meta$site.name=="Ocean Exhibit, Blue Planet, Ellesmere Port" & (meta$via=="SCUBA" | meta$type=="eDNA"))
livt <- subset(meta$seq_id, meta$site.name=="Dukes Dock, Liverpool"& (meta$type=="eDNA" | (meta$extraction=="QBT" & meta$perservation=="Et")))
orke <- subset(meta$seq_id, meta$site.name=="SMS Brummer, Scapa Flow, Orkney" &  (meta$type=="eDNA" | (meta$extraction=="QBT" & meta$perservation=="Et")))

# tele02 Orkney sequence run 1

orkt_motu <- p1_prev %>% select("id","final_name","final_genus","final_order","final_class", 
                                "pid","final_rank","method_assign","total_reads","sequence",
                                all_of(orkt)) 
# elas02 shark tank sequence run 2
blue_motu <- p2e_aq_prev %>% select("id","final_name","final_genus","final_order","final_class", 
                                    "pid","final_rank","method_assign","total_reads","sequence",
                                    all_of(blue)) 
# tele02 Liverpool sequence run 2
livt_motu <- p2t_uk_prev %>% select("id","final_name","final_genus","final_order","final_class", 
                                    "pid","final_rank","method_assign","total_reads","sequence",    
                                    all_of(livt)) 

# elas02 Orkney sequence run 2
# adding sample that was left off of run 2, so added to run 3
p2e_uk_prev["sample.12G_ORKMPbA_kurt"] <- 0
p2e_uk_prev_1 <- rbind.fill(p2e_uk_prev, samp)
p2e_uk_prev_1 <- p2e_uk_prev_1 %>% replace(is.na(.), 0)
p2e_uk_prev_1 <- p2e_uk_prev_1 %>% mutate(total_reads = rowSums(.[3:15]))
orke_motu <- p2e_uk_prev_1 %>% select("id","final_name","final_genus","final_order","final_class", 
                                    "pid","final_rank","method_assign","total_reads","sequence",  
                                    all_of(orke)) 
rm(p1_prev,p2e_aq_prev,p2t_uk_prev,p2e_uk_prev,p2e_uk_prev_1, samp, p3_prev)


#####
# Preparing data in long format and collapsing for nmds
#####

# make motu tables long
l1 <- orkt_motu %>% pivot_longer(cols = sample.7A:sample.9E,
                           names_to = "sample",
                           names_prefix = "sample.",
                           values_to = "reads")
# Blue planet
l2 <- blue_motu %>% pivot_longer(cols = sample.10FeBLUE_MPEtA_RA:sample.9CeBLUE_MPEtD_DB,
                                 names_to = "sample",
                                 names_prefix = "sample.",
                                 values_to = "reads")
# manually fix elasmobranch taxonomy
l2$final_name <- ifelse(l2$final_name=="Chiloscyllium griseum", "Chiloscyllium",
          ifelse(l2$final_name=="Heterodontus", "Heterodontus",
              ifelse(l2$final_name=="Orectolobus japonicus", "Orectolobus",
                  l2$final_name)))
l2 <- subset(l2, l2$final_name!="Dicentrarchus labrax" |
                           l2$final_name!="Hypophthalmichthys" |
                           l2$final_name!="Molva molva" |
                           l2$final_name!="Salmo salar" |
                           l2$final_name!="Trisopterus minutus")

l3 <- livt_motu %>% pivot_longer(cols = sample.4BtLIV_eDNA_seA_bottle1:sample.5AtLIV_MPEtC_Cath,
                                 names_to = "sample",
                                 names_prefix = "sample.",
                                 values_to = "reads")
l4 <- orke_motu %>% pivot_longer(cols = sample.6BeORK_eDNAA_bottle1:sample.6HeORK_MPEtC_lisa,
                                 names_to = "sample",
                                 names_prefix = "sample.",
                                 values_to = "reads")
l5 <- rbind(l1,l2,l3,l4)

# master wide mou table
w1 <- l5 %>% pivot_wider(names_from = "sample",
                   names_prefix = "sample.",
                   values_from = "reads",
                   values_fill = 0)

w1$total_reads <- rowSums(w1[c(11:46)])

w2 <- w1 %>% filter(total_reads > 0)

w3 <- w2 %>%
  group_by(final_rank,final_class,final_order,final_genus,final_name) %>%
  summarise(across(c(total_reads, sample.7A:sample.6HeORK_MPEtC_lisa), sum))

# group again just by final name because some have different family names
# due to taxonomic assignment from different references so not collapsing properly
# Pan troglodytes; misassigned human
w3$final_name <- ifelse(w3$final_name=="Pan troglodytes", "Homo sapiens", w3$final_name)
w4 <- w3 %>%
  group_by(final_rank,final_class,final_name) %>%
  summarise(across(c(total_reads, sample.7A:sample.6HeORK_MPEtC_lisa), sum))
w4 <- subset(w4, w4$final_rank=="genus" | w4$final_rank=="species")

  
# Inspect for "contamination" or errors in taxonomy 
# decide what to keep or remove for the NMDS etc.

# Dama dama; not possible in Orkney
w4 <- subset(w4, w4$final_name!="Dama dama")
# Bos mutus; not possible in Orkney
w4 <- subset(w4, w4$final_name!="Bos mutus")
# Equus asinus; not possible in Liverpool
w4 <- subset(w4, w4$final_name!="Equus asinus")
# Columba palumbus and Columba livia; needs to be changed to Aves
w4$final_class <- ifelse(w4$final_name=="Columba palumbus" | w4$final_name=="Columba livia" |
                           w4$final_name=="Corvus", "Aves", w4$final_class)
# fix mammals
w4$final_class <- ifelse(w4$final_name=="Equus", "Mammalia", w4$final_class)
# Chimaera monstrosa; unlikely
w4 <- subset(w4, w4$final_name!="Chimaera monstrosa")
# Pangasianodon hypophthalmus; positive control
w4 <- subset(w4, w4$final_name!="Pangasianodon hypophthalmus")
# Hypanus sp.; genus level assignment was removed during other analysis; low reads
w4 <- subset(w4, w4$final_name!="Hypanus")

## Make different plots for different purposes
# remove domestic species
w4 <- subset(w4, w4$final_name!="Equus")
w4 <- subset(w4, w4$final_name!="Canis")
w4 <- subset(w4, w4$final_name!="Ovis")
w4 <- subset(w4, w4$final_name!="Felis catus")
w4 <- subset(w4, w4$final_name!="Sus scrofa")

# remove human
w5 <- subset(w4, w4$final_name!="Homo sapiens")


# Interesting to note: 	
#Osmerus eperlanus; European smelt

#####
## Prepare data for NMDS
#####
library(vegan)
library(sjmisc)
library(janitor)

w6 <- w5[c(3,5:40)]
w7 <- w6 %>%
  rotate_df(w6) %>%
  row_to_names(row_number = 1) 
colnames(w7)[1] <- "seq_id"

# extract abundance data and prepare for MDS
dat <- w7[,2:ncol(w7)]
# convert from character to numeric dataframe
dat <- as.data.frame(sapply(dat, as.numeric)) 

# presence-absence dataframe
dat_pa <- as.data.frame(ifelse(dat > 0, 1, 0))
# Hellinger's standardization
dat_hell <- decostand(dat, method = "hellinger")

# Calculate distances
# Bray-Curtis dissimilarity matrix
bray_dat <- vegdist(dat_hell, method = "bray")
# Calculate jaccard dissimilariy matrix
jac_dat <- vegdist(dat_pa, method = "jaccard", binary =  TRUE)

#####
## NMDS & extract data for plots
#####

set.seed(1995)
ord_bray <- metaMDS(bray_dat, distance = "bray", trymax = 1000)
plot(ord_bray)

set.seed(1996)
ord_jac <- metaMDS(jac_dat, distance = "jaccard", trymax = 1000)
plot(ord_jac)

## extract site scores from NMDS and add metadata
cem_bray <- as.data.frame(scores(ord_bray))  
cem_bray$seq_id <- w7$seq_id
cem_bray <- merge(cem_bray, meta, by.x = "seq_id")

cem_jac <- as.data.frame(scores(ord_jac))  
cem_jac$seq_id <- w7$seq_id
cem_jac <- merge(cem_jac, meta, by.x = "seq_id")

#####
## NMDS plots
#####
library(ggplot2)
library(ggthemes)

# basic
ggplot(cem_bray) +
  geom_point(aes(x = NMDS1, y = NMDS2, color = site.name, shape = type))

ggplot(cem_jac) +
  geom_point(aes(x = NMDS1, y = NMDS2, color = site.name, shape = type))

#prepare data for more informative plot
cem_bray$type <- factor(cem_bray$type, levels = c("eDNA", "MP"))
cem_bray$primer <- factor(cem_bray$primer, levels = c("elas02", "tele02"))
cem_bray$site.name <- factor(cem_bray$site.name, levels = c("Ocean Exhibit, Blue Planet, Ellesmere Port",
                                                            "Dukes Dock, Liverpool",                     
                                                            "SMS Brummer, Scapa Flow, Orkney",           
                                                            "SMS Bayern, Scapa Flow, Orkney"))
#prepare data for more informative plot
cem_jac$type <- factor(cem_jac$type, levels = c("eDNA", "MP"))
cem_jac$primer <- factor(cem_jac$primer, levels = c("elas02", "tele02"))
cem_jac$site.name <- factor(cem_jac$site.name, levels = c("Ocean Exhibit, Blue Planet, Ellesmere Port",
                                                            "Dukes Dock, Liverpool",                     
                                                            "SMS Brummer, Scapa Flow, Orkney",           
                                                            "SMS Bayern, Scapa Flow, Orkney"))

# Bray curtis groups so strong that all the points are on top of each other
# sample type
cem_bray_nmds <-
ggplot(cem_bray, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(fill = site.name, shape = type), size = 5, alpha = 0.5) +
  scale_shape_manual(values = c(21, 24)) +
  scale_color_colorblind()
ggsave(filename=c("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_writing/figures/nmds_bray_cem.jpg"), 
       plot = cem_bray_nmds, width = 8, height = 6.5, units = "in")

# primer type
ggplot(cem_bray,
       aes(x = NMDS1, y = NMDS2)
) +
  geom_point(aes(fill = location, shape = primer), size = 5, alpha = 0.8) +
  scale_shape_manual(values = c(21, 24)) +
  scale_color_colorblind()


# prepare plot for paper

cem_jac_nmds <-
ggplot(cem_jac,  aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(shape = type, fill = site.name, colour = primer), alpha= 0.7, size = 6, stroke = 1.5) + 
  scale_shape_manual(values = c(24,21),labels = c("Syringe Filter", "Diver MP")) +
  scale_fill_manual(values = c("#1E88E5","#D81B60","#FFC107", "#004D40"),
                    labels = c("Ocean Display, Blue Planet", "Dukes Dock, Liverpool", 
                               "SMS Brummer, Orkney","SMS Bayern, Orkney")) +
  scale_colour_manual(values = c("black", "white"),labels = c("Elas02", "Tele02")) +
  labs(x = "NMDS1", colour = "Primer", y = "NMDS2") +
  guides(fill = guide_legend("Site", override.aes = list(shape = 21, colour = "darkgrey")),
         shape = guide_legend("Sample Type", override.aes = list(fill = "#4477AA", colour = "darkgrey"))) +
  theme(axis.text.y = element_text(colour = "black", size = 16, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 16), 
        legend.text = element_text(size = 16, colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 18), 
        axis.title.x = element_text(face = "bold", size = 18, colour = "black"), 
        legend.title = element_text(size = 18, colour = "black", face = "bold")) 

ggsave(filename=c("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_writing/figures/nmds_jac_cem.jpg"), 
       plot = cem_jac_nmds, width = 8, height = 6.5, units = "in")

#####
## Prepare data for richness boxplots
#####

## NEED TO EDIT THIS SECTION FOR OLD CODE

cem_pa <- cbind(w7[,1], dat_pa)
colnames(cem_pa)[1] <- "seq_id"

cem_pa$richness <- rowSums(cem_pa[,2:ncol(cem_pa)])

cem_pa <- merge(meta, cem_pa, by = "seq_id")

loc.labs <- c("Blue Planet", "Liverpool", "Orkney")
names(loc.labs) <- c("Blue Planet Aquarium", "Liverpool", "Orkney")

# Might want to change this later, figure out how to
type.labs <- c("Bottle", "Diver MP")
names(type.labs) <- c("eDNA", "Diver MP")

# make boxplots
cem_boxplots <-
ggplot(cem_pa, aes(x = type, y = richness)) + 
  geom_boxplot(aes(fill = site.name), alpha = 0.7, outlier.shape = NA) +
  geom_jitter(shape = 19, position = position_jitter(0.3)) +
  facet_wrap(~ location, labeller = labeller(location = loc.labs)) +
  scale_fill_manual(values = c( "#D81B60", "#1E88E5", "#004D40", "#FFC107"),
                    labels = c("Dukes Dock, Liverpool", "Ocean Exhibit, Blue Planet",
                               "SMS Bayern, Orkney", "SMS Brummer, Orkney")) +
  labs(x = "Sample Type", fill = "Site", y = "Species Richness") +
  theme(axis.text.y = element_text(colour = "black", size = 16, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 16), 
        legend.text = element_text(size = 16, colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 18), 
        axis.title.x = element_text(face = "bold", size = 18, colour = "black"), 
        legend.title = element_text(size = 18, colour = "black", face = "bold"),
        strip.text = element_text(colour = "black", size = 16, face = "bold"))

ggsave(filename=c("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_writing/figures/boxplots_cem.jpg"), 
       plot = cem_boxplots, width = 8, height = 6.5, units = "in")


#####
## Prepare data for Venn diagrams
#####

vcem <- w5

# add together columns by groups
vcem$bayern_bottle <- vcem$sample.7A + vcem$sample.7B +
                   vcem$sample.7C + vcem$sample.7D 
vcem$brummer_bottle <- vcem$sample.6BeORK_eDNAA_bottle1 + vcem$sample.6CeORK_eDNAB_bottle2 +
                   vcem$sample.6DeORK_eDNAC_bottle3 + vcem$sample.6EeORK_eDNAD_bottle4

vcem$bayern_diver <- vcem$sample.7E + vcem$sample.7G +
                  vcem$sample.8D + vcem$sample.8F +
                  vcem$sample.9C + vcem$sample.9E 
vcem$brummer_diver <- vcem$sample.6FeORK_MPEtA_kurt + vcem$sample.6GeORK_MPEtB_mike +   
                  vcem$sample.6HeORK_MPEtC_lisa    

vcem$blue_bottle <- vcem$sample.8DeBLUE_eDNAA_bottle1 + vcem$sample.8EeBLUE_eDNAB_bottle2 +
                    vcem$sample.8FeBLUE_eDNAC_bottle3 + vcem$sample.8GeBLUE_eDNAD_bottle4

vcem$blue_diver <- vcem$sample.10FeBLUE_MPEtA_RA + vcem$sample.10GeBLUE_MPEtB_RB +
                   vcem$sample.10HeBLUE_MPEtC_DA + vcem$sample.11AeBLUE_MPEtD_DB +
                   vcem$sample.8HeBLUE_MPEtA_RA + vcem$sample.9AeBLUE_MPEtB_RB +
                   vcem$sample.9BeBLUE_MPEtC_DA + vcem$sample.9CeBLUE_MPEtD_DB

vcem$liv_bottle <- vcem$sample.4BtLIV_eDNA_seA_bottle1 + vcem$sample.4CtLIV_eDNA_seB_bottle2 +
                   vcem$sample.4DtLIV_eDNA_DeC_bottle3 + vcem$sample.4EtLIV_eDNA_DeD_bottle4

vcem$liv_diver <- vcem$sample.4FtLIV_MPEtA_Wendy + vcem$sample.4GtLIV_MPEtB_Rosie +
                  vcem$sample.5AtLIV_MPEtC_Cath

# group; remove taxa that both == 0
vcem_bayern <- vcem[c("final_name","bayern_bottle", "bayern_diver")]
vcem_bayern <- vcem_bayern[rowSums(vcem_bayern[c(2,3)])>0,] 

vcem_brummer <- vcem[c("final_name","brummer_bottle", "brummer_diver")]
vcem_brummer <- vcem_brummer[rowSums(vcem_brummer[c(2,3)])>0,] 

vcem_blue <- vcem[c("final_name","blue_bottle", "blue_diver")]
vcem_blue <- vcem_blue[rowSums(vcem_blue[c(2,3)])>0,] 

vcem_liv <- vcem[c("final_name","liv_bottle", "liv_diver")]
vcem_liv <- vcem_liv[rowSums(vcem_liv[c(2,3)])>0,] 

## convert to true or false
vcem_bayern$bayern_bottle <- ifelse(vcem_bayern$bayern_bottle==0, FALSE, TRUE)
vcem_bayern$bayern_diver <- ifelse(vcem_bayern$bayern_diver==0, FALSE, TRUE)

vcem_brummer$brummer_bottle <- ifelse(vcem_brummer$brummer_bottle==0, FALSE, TRUE)
vcem_brummer$brummer_diver <- ifelse(vcem_brummer$brummer_diver==0, FALSE, TRUE)

vcem_blue$blue_bottle <- ifelse(vcem_blue$blue_bottle==0, FALSE, TRUE)
vcem_blue$blue_diver <- ifelse(vcem_blue$blue_diver==0, FALSE, TRUE)

vcem_liv$liv_bottle <- ifelse(vcem_liv$liv_bottle==0, FALSE, TRUE)
vcem_liv$liv_diver <- ifelse(vcem_liv$liv_diver==0, FALSE, TRUE)


#####
## Venn diagrams
#####

library(ggvenn)
library(cowplot)

names(vcem_bayern) <- c("Taxa","Filter (4)","MP (6)")
v4 <- ggvenn(vcem_bayern, c(A = "Filter (4)", B = "MP (6)"),
       set_name_size = 5, text_size = 5, 
       fill_color = c("#004D40", "#99FFEE")) 
ggsave(filename=c("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_writing/figures/venn_4.jpg"), 
       plot = v4, width = 6, height = 5, units = "in")

names(vcem_brummer) <- c("Taxa","Filter (4)","MP (3)")
v3 <- ggvenn(vcem_brummer, c(A = "Filter (4)", B = "MP (3)"),
       set_name_size = 5, text_size = 5, 
       fill_color = c("#FFC107", "#FFE699")) 

names(vcem_blue) <- c("Taxa","Filter (4)","MP (8)")
v1 <- ggvenn(vcem_blue, c(A = "Filter (4)", B = "MP (8)"),
       set_name_size = 5, text_size = 5,
       fill_color = c("#1E88E5", "#A4CFF4")) 

names(vcem_liv) <- c("Taxa","Filter (4)","MP (3)")
v2 <- ggvenn(vcem_liv, c(A = "Filter (4)", B = "MP (3)"),
       set_name_size = 5, text_size = 5,
       fill_color = c("#D81B60", "#EE77A2")) 

av <- plot_grid(v1,v2,v3,v4, labels = "AUTO", ncol = 2)


ggsave(filename=c("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_writing/figures/venn_cem.jpg"), 
       plot = av, width = 10, height = 9.5, units = "in")

avn <- plot_grid(av, cem_jac_nmds, labels = c("", "E"), ncol = 1)

ggsave(filename=c("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_writing/figures/figure2.jpg"), 
       plot = avn, width = 10.5, height = 14, units = "in")
