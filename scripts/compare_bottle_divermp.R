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
all_reads <- p3_prev %>%
  mutate(sum = rowSums(across(c(sample.10A_SAtra1001:sample.9H_SA_SL7F2A_c))))
total <- sum(all_reads$sum)
# calculate 0.001% of reads
remove <- total*0.00001
rm(all_reads)
# replace 0.001% of reads with zero
var <- subset(meta$seq_id, meta$sequence.run==3)
p3_prev[,var][p3_prev[,var] <= 7] <- 0
# total reads per MOTU
p3_prev <- p3_prev %>% mutate(total_reads = rowSums(.[12:54]))
p3_prev <- subset(p3_prev, p3_prev$total_reads>0)
# save
write.csv(p3_prev, "C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/decontam/p3_decontam_001.csv")

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
l2 <- blue_motu %>% pivot_longer(cols = sample.10FeBLUE_MPEtA_RA:sample.9CeBLUE_MPEtD_DB,
                                 names_to = "sample",
                                 names_prefix = "sample.",
                                 values_to = "reads")
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
w4 <- w3 %>%
  group_by(final_rank,final_class,final_name) %>%
  summarise(across(c(total_reads, sample.7A:sample.6HeORK_MPEtC_lisa), sum))
  
## NEED TO INSPECT FOR "CONTAMINATION"
## decide what to keep or remove for the NMDS etc.

#####
## Prepare data for NMDS
#####
library(vegan)
library(sjmisc)
library(janitor)

w4 <- w3[c(4:40)]
w5 <- w4 %>%
  rotate_df(w4) %>%
  row_to_names(row_number = 1) 
colnames(w5)[1] <- "seq_id"

# extract abundance data and prepare for MDS
dat <- w5[,2:ncol(w5)]
# convert from character to numeric dataframe
dat <- as.data.frame(sapply(dat, as.numeric)) 

# presence-absence dataframe
dat_pa <- as.data.frame(ifelse(dat > 0, 1, 0))
# Hellinger's standardization
dat_hell <- decostand(dat, method = "hellinger")

# Calculate distances
# Bray-Curtis dissimilarity matrix
bray_dat <- vegdist(dat_hell, method = "bray")
## Calculate jaccard dissimilariy matrix
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
cem_bray$seq_id <- w5$seq_id
cem_bray <- merge(cem_bray, cem_meta, by.x = "seq_id")

cem_jac <- as.data.frame(scores(ord_jac))  
cem_jac$seq_id <- w5$seq_id
cem_jac <- merge(cem_jac, cem_meta, by.x = "seq_id")

#####
## NMDS plots
#####
library(ggplot2)
library(ggthemes)

# basic
ggplot(cem_bray) +
  geom_point(aes(x = NMDS1, y = NMDS2, color = location, shape = type))

ggplot(cem_jac) +
  geom_point(aes(x = NMDS1, y = NMDS2, color = location, shape = type))

#prepare data for more informative plot
cem_bray$type <- factor(cem_bray$type, levels = c("eDNA", "MP"))
cem_bray$primer <- factor(cem_bray$primer, levels = c("elas02", "tele02"))
cem_bray$location <- factor(cem_bray$location, levels = c("Blue Planet", "Liverpool", "Orkney"))


# make more informative

# sample type
ggplot(cem_bray,
       aes(x = NMDS1, y = NMDS2)
       ) +
  geom_point(aes(fill = location, shape = type), size = 5, alpha = 0.8) +
  scale_shape_manual(values = c(21, 24)) +
  scale_color_colorblind()
# primer type
ggplot(cem_bray,
       aes(x = NMDS1, y = NMDS2)
) +
  geom_point(aes(fill = location, shape = primer), size = 5, alpha = 0.8) +
  scale_shape_manual(values = c(21, 24)) +
  scale_color_colorblind()

#####
# have a go at a spider plot
#library(ggordiplots)
#gg_ordiplot(ord_bray, groups = cem_bray$location, 
#            spiders=TRUE, ellipse=FALSE, plot=TRUE) # groups aren't labeled properly

#gg_ordiplot(ord_jac, groups = cem_bray$location, 
#            spiders=TRUE, ellipse=FALSE, plot=TRUE) # groups aren't labeled properly
#####

# make plot for paper

# primer is shape, sample type is colour
#ggplot(cem_bray,  aes(x = NMDS1, y = NMDS2)) + 
#  geom_point(aes(shape = primer, fill = location, colour = type), alpha= 0.7, size = 4, stroke = 1.5) + 
#  scale_shape_manual(values = c(24,21),labels = c("Elas02", "Tele02")) +
#  scale_fill_manual(values = c("#4477AA","#EE6677","#228833"),labels = c("Blue Planet", "Liverpool", "Orkney")) +
#  scale_colour_manual(values = c("black", "darkgrey"),labels = c("eDNA", "MP")) +
#  labs(x = "NMDS1", colour = "Sample Type", y = "NMDS2", shape = "Primer") +
#  guides(fill = guide_legend("Location", override.aes = list(shape = 21, colour = "darkgrey")),
#         shape = guide_legend("Primer", override.aes = list(fill = "#4477AA", colour = "darkgrey")))

cem_bray_nmds <-
ggplot(cem_bray,  aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(shape = type, fill = location, colour = primer), alpha= 0.7, size = 4, stroke = 1.5) + 
  scale_shape_manual(values = c(24,21),labels = c("eDNA", "MP")) +
  scale_fill_manual(values = c("#4477AA","#EE6677","#228833"),labels = c("Blue Planet", "Liverpool", "Orkney")) +
  scale_colour_manual(values = c("black", "darkgrey"),labels = c("Elas02", "Tele02")) +
  labs(x = "NMDS1", colour = "Primer", y = "NMDS2") +
  guides(fill = guide_legend("Location", override.aes = list(shape = 21, colour = "darkgrey")),
         shape = guide_legend("Sample Type", override.aes = list(fill = "#4477AA", colour = "darkgrey"))) +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold")) 

ggsave(filename=c("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_writing/figures/nmds_bray_cem.jpg"), 
       plot = cem_bray_nmds, width = 8, height = 6.5, units = "in")

#####
## Prepare data for richness boxplots
#####

cem_pa <- cbind(w5[,1], dat_pa)
colnames(cem_pa)[1] <- "seq_id"

cem_pa$richness <- rowSums(cem_pa[,2:ncol(cem_pa)])

cem_pa <- merge(cem_meta, cem_pa, by = "seq_id")

# make boxplots
cem_boxplots <-
ggplot(cem_pa, aes(x = type, y = richness)) + 
  geom_boxplot(aes(fill = location), alpha = 0.7, outlier.shape = NA) +
  geom_jitter(shape = 19, position = position_jitter(0.3)) +
  facet_wrap(~ location) +
  scale_fill_manual(values = c("#4477AA","#EE6677","#228833"),labels = c("Blue Planet", "Liverpool", "Orkney")) +
  labs(x = "Sample Type", fill = "Location", y = "Species Richness") +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"),
        strip.text = element_text(colour = "black", size = 12, face = "bold"))

ggsave(filename=c("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_writing/figures/boxplots_cem.jpg"), 
       plot = cem_boxplots, width = 8, height = 6.5, units = "in")


#####
## Prepare data for Venn diagrams
#####

vcem <- w3

# add together columns by groups
vcem$ork_bottle <- vcem$sample.7A + vcem$sample.7B +
                   vcem$sample.7C + vcem$sample.7D +
                   vcem$sample.6BeORK_eDNAA_bottle1 + vcem$sample.6CeORK_eDNAB_bottle2 +
                   vcem$sample.6DeORK_eDNAC_bottle3 + vcem$sample.6EeORK_eDNAD_bottle4

vcem$ork_diver <- vcem$sample.7E + vcem$sample.7G +
                  vcem$sample.8D + vcem$sample.8F +
                  vcem$sample.9C + vcem$sample.9E +
                  vcem$sample.6FeORK_MPEtA_kurt + vcem$sample.6GeORK_MPEtB_mike +   
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

## group
vcem_ork <- vcem[c(4,42,43)]
vcem_ork <- vcem_ork[rowSums(vcem_ork[c(2,3)])>0,] 

vcem_blue <- vcem[c(4,44,45)]
vcem_blue <- vcem_blue[rowSums(vcem_blue[c(2,3)])>0,] 

vcem_liv <- vcem[c(4,46,47)]
vcem_liv <- vcem_liv[rowSums(vcem_liv[c(2,3)])>0,] 

## convert to true or false
vcem_ork$ork_bottle <- ifelse(vcem_ork$ork_bottle==0, FALSE, TRUE)
vcem_ork$ork_diver <- ifelse(vcem_ork$ork_diver==0, FALSE, TRUE)

vcem_blue$blue_bottle <- ifelse(vcem_blue$blue_bottle==0, FALSE, TRUE)
vcem_blue$blue_diver <- ifelse(vcem_blue$blue_diver==0, FALSE, TRUE)

vcem_liv$liv_bottle <- ifelse(vcem_liv$liv_bottle==0, FALSE, TRUE)
vcem_liv$liv_diver <- ifelse(vcem_liv$liv_diver==0, FALSE, TRUE)


#####
## Venn diagrams
#####

library(ggvenn)

ggvenn(vcem_ork, c(A = "ork_bottle", B = "ork_diver"),
       set_name_size = 10, text_size = 5) 

ggvenn(vcem_blue, c(A = "blue_bottle", B = "blue_diver"),
       set_name_size = 10, text_size = 5) 

ggvenn(vcem_liv, c(A = "liv_bottle", B = "liv_diver"),
       set_name_size = 10, text_size = 5) 

