##############################################
## Compare eDNA bottles to diver metaprobes ##
##############################################
library(tidyverse)

#load data
# tele02 Orkney sequence run 1
p1_prev <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/divers_methods_experiment/sintax_taxonomy/divmeth1_MOTU_decontam.csv")
# elas02 shark tank sequence run 2
p2e_aq_prev <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/divers_methods_experiment_part2/sintax_output_MainTank/divmeth2_sharktank_MOTU_decontam.csv")
# tele02 Liverpool sequence run 2
p2t_uk_prev <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/divers_methods_experiment_part2/sintax_output_UKreferencedatabase/t_divmeth2_UK_MOTU_decontam.csv")
# elas02 Orkney sequence run 2
p2e_uk_prev <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/divers_methods_experiment_part2/sintax_output_UKreferencedatabase/e_divmeth2_UK_MOTU_decontam.csv")


# select metadata of samples for comparison of eDNA and MP
## might want to add some lines of code to get to samp_meta
cem_meta <- samp_meta[c("49":"53","55","60","62","67","69","89":"105","114":"120","144":"151"),]
cem_meta <- cem_meta %>% drop_na(seq_id)

# list of sample ids
orkt <- subset(cem_meta$seq_id, cem_meta$site=="Bayern")
blue <- subset(cem_meta$seq_id, cem_meta$location=="Blue Planet")
livt <- subset(cem_meta$seq_id, cem_meta$location=="Liverpool")
orke <- subset(cem_meta$seq_id, cem_meta$site=="SMS Brummer")

# extract data for this visualization
# tele02 Orkney sequence run 1
orkt_motu <- p1_prev %>% select("id","class","c.id","order","o.id","family","f.id","genus","g.id","species","s.id",
                        "sample.7A", "sample.7B", "sample.7C", "sample.7D", "sample.7E", "sample.7G",
                        "sample.8D", "sample.8F", "sample.9C", "sample.9E") %>%
                        filter(s.id==1)

# elas02 shark tank sequence run 2
blue_motu <- p2e_aq_prev %>% select("id","class","c.id","order","o.id","family","f.id","genus","g.id","species","s.id",
                                    "sample.10FeBLUE_MPEtA_RA",     "sample.10GeBLUE_MPEtB_RB",    
                                    "sample.10HeBLUE_MPEtC_DA",    "sample.11AeBLUE_MPEtD_DB",    
                                    "sample.8DeBLUE_eDNAA_bottle1", "sample.8EeBLUE_eDNAB_bottle2",
                                    "sample.8FeBLUE_eDNAC_bottle3", "sample.8GeBLUE_eDNAD_bottle4",
                                    "sample.8HeBLUE_MPEtA_RA",      "sample.9AeBLUE_MPEtB_RB",     
                                    "sample.9BeBLUE_MPEtC_DA",      "sample.9CeBLUE_MPEtD_DB") %>%
                                    filter(s.id==1)

# tele02 Liverpool sequence run 2
livt_motu <- p2t_uk_prev %>% select("id","class","c.id","order","o.id","family","f.id","genus","g.id","species","s.id",
                                    "sample.4BtLIV_eDNA_seA_bottle1", "sample.4CtLIV_eDNA_seB_bottle2",
                                    "sample.4DtLIV_eDNA_DeC_bottle3", "sample.4EtLIV_eDNA_DeD_bottle4",
                                    "sample.4FtLIV_MPEtA_Wendy",      "sample.4GtLIV_MPEtB_Rosie",     
                                    "sample.5AtLIV_MPEtC_Cath") %>%
                                    filter(s.id==1)

# elas02 Orkney sequence run 2
orke_motu <- p2e_uk_prev %>% select("id","class","c.id","order","o.id","family","f.id","genus","g.id","species","s.id",
                                    "sample.6BeORK_eDNAA_bottle1", "sample.6CeORK_eDNAB_bottle2",
                                    "sample.6DeORK_eDNAC_bottle3", "sample.6EeORK_eDNAD_bottle4",
                                    "sample.6FeORK_MPEtA_kurt",    "sample.6GeORK_MPEtB_mike",   
                                    "sample.6HeORK_MPEtC_lisa") %>%
                                    filter(s.id==1)
rm(p1_prev,p2e_aq_prev,p2t_uk_prev,p2e_uk_prev)


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

w1$total_reads <- rowSums(w1[c(12:47)])

w2 <- w1 %>% filter(total_reads > 0)

w3 <- w2 %>%
  group_by(order,family,genus,species) %>%
  summarise(across(c(sample.7A:total_reads), sum))


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


# have a go at a spider plot
library(ggordiplots)

gg_ordiplot(ord_bray, groups = cem_bray$location, 
            spiders=TRUE, ellipse=FALSE, plot=TRUE) # groups aren't labeled properly

gg_ordiplot(ord_jac, groups = cem_bray$location, 
            spiders=TRUE, ellipse=FALSE, plot=TRUE) # groups aren't labeled properly


