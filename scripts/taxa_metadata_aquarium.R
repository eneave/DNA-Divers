################################################
## Add Taxonomic metadata to aquarium samples ##
################################################

library(ggplot2)
library(tidyverse)
library(janitor)
library(ggthemes)

#####
## Shark tank
#####
p2e_aq_prev <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/decontam/p2e_aq_decontam.csv")
#p2t_aq_prev <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/decontam/p2t_aq_decontam.csv")

# calculate total reads and number of reads to remove
all_reads <- p2e_aq_prev %>%
        mutate(sum = rowSums(across(c(sample.10Ae60BLUE_MPEtB:sample.9He60BLUE_MPEtA_))))
total <- sum(all_reads$sum)
# calculate 0.001% of reads
remove <- total*0.00001
rm(all_reads)

# replace 0.001% of reads with zero
var <- c("sample.10Ae60BLUE_MPEtB", "sample.10Be120BLUE_MPEtA",      
        "sample.10Ce120BLUE_MPEtB", "sample.10De240BLUE_MPEtB",      
         "sample.10Ee240BLUE_MPEtA", "sample.10FeBLUE_MPEtA_RA",      
          "sample.10GeBLUE_MPEtB_RB", "sample.10HeBLUE_MPEtC_DA",      
          "sample.11AeBLUE_MPEtD_DB", "sample.11Be_EBMay_13extblank",  
          "sample.11Ce_EBJun_12extblank", "sample.11De_EBJun_19extblank",  
          "sample.11EeBLUE_FBEt_fieldb", "sample.11Fe_negativePCRcontrol",
          "sample.8CeBLUE_eDNA_FBblank", "sample.8DeBLUE_eDNAA_bottle1",  
          "sample.8EeBLUE_eDNAB_bottle2", "sample.8FeBLUE_eDNAC_bottle3",  
          "sample.8GeBLUE_eDNAD_bottle4", "sample.8HeBLUE_MPEtA_RA",       
           "sample.9AeBLUE_MPEtB_RB", "sample.9BeBLUE_MPEtC_DA",       
          "sample.9CeBLUE_MPEtD_DB", "sample.9De10BLUE_MPEtA_",       
         "sample.9Ee10BLUE_MPEtB_", "sample.9Fe30BLUE_MPEtA_",       
          "sample.9Ge30BLUE_MPEtB_", "sample.9He60BLUE_MPEtA_")
p2e_aq_prev[,var][p2e_aq_prev[,var] <= 5] <- 0

# total reads per MOTU
p2e_aq_prev <- p2e_aq_prev %>% mutate(total_reads = rowSums(.[3:30]))
p2e_aq_prev <- subset(p2e_aq_prev, p2e_aq_prev$total_reads>0)

write.csv(p2e_aq_prev, "C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/decontam/p2e_aq_decontam_001.csv")

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
# Calculate total sample reads
# CHECK TO SEE IF THIS IS CORRECT
long_aq_elas %>% 
  group_by(seq_id) %>% 
  summarise(total_sample_reads = sum(reads))
# add proportional read counts NOT CALCULATED CORRECTLY - NEEDS TOTAL SAMPLE READS
#long_aq_elas$prc <- long_aq_elas$reads/long_aq_elas$total_reads
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
ggplot(long_aq_elas , aes(x = seq_id, y = reads, fill = manual_name)) + 
#ggplot(long_aq_elas , aes(x = seq_id, y = prc, fill = manual_name)) + 
  geom_bar(stat = "identity", color = "black", position = "fill") +
  facet_grid(. ~ type2, scales = "free") +
  labs(title="Proportional Read counts of Elasmobranchs \ndetected in the Ocean display by different methods",
       x ="Sample", y = "Proportional Read Counts (PRC)") +
  theme(axis.text.x = element_text(angle = 90))

#####
## Subset and prepare other detections to add to the stacked bar charts
#####

p2e_aq_notelas <- subset(p2e_aq_prev, final_class!="Elasmobranchii")
# Remove domestic species (low amount of reads)
p2e_aq_notelas <- subset(p2e_aq_notelas, final_order!="Artiodactyla" & final_order!="Galliformes"
                      & final_order!="Carnivora" & final_order!="Pelecaniformes")
# Remove contamination from aquarium (tide pool exhibit) or other samples or could be food but we don't know
p2e_aq_notelas <- subset(p2e_aq_notelas, final_name!="Dicentrarchus labrax" |
                           final_name!="Hypophthalmichthys" |
                           final_name!="Molva molva" |
                           final_name!="Salmo salar" |
                           final_name!="Trisopterus minutus")

# convert to long dataframe
long_aq_notelas <- p2e_aq_notelas %>%
  pivot_longer(cols = sample.10Ae60BLUE_MPEtB:sample.9He60BLUE_MPEtA_,
               names_to = "seq_id",
               values_to = "reads")
# add metadata to long dataframe
meta2 <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/metadata/supp_table_2.csv")
long_aq_notelas <- merge(long_aq_notelas, meta, by="seq_id")
long_aq_notelas <- merge(long_aq_notelas, meta2, by="final_name")


# add proportional read counts
#long_aq_notelas$prc <- long_aq_notelas$reads/long_aq_notelas$total_reads NOT CALCULATED CORRECTLY - NEEDS TOTAL SAMPLE READS
# specify sample type more clearly
long_aq_notelas$type2 <-ifelse(long_aq_notelas$type=="eDNA", "Syringe Filter", 
                            ifelse(long_aq_notelas$type=="MP" & (long_aq_notelas$time==65|long_aq_notelas$time==50), "Diver MP", "Soak MP"))


# quick stacked bar plot of fish
windows()
ggplot(long_aq_notelas , aes(x = seq_id, y = reads, fill = final_name)) + 
  geom_bar(stat = "identity", color = "black", position = "fill") +
  facet_grid(. ~ type2, scales = "free") +
  labs(title="Proportional Read counts of Teleosts \ndetected in the Ocean display by different methods",
       x ="Sample", y = "Proportional Read Counts (PRC)") +
  theme(axis.text.x = element_text(angle = 90))

# quick stacked bar plot of category fish
windows()
ggplot(long_aq_notelas , aes(x = seq_id, y = reads, fill = category2)) + 
  geom_bar(stat = "identity", color = "black", position = "fill") +
  facet_grid(. ~ type2, scales = "free") +
  labs(title="Proportional Read counts of Teleosts \ndetected in the Ocean display by different methods",
       x ="Sample", y = "Proportional Read Counts (PRC)") +
  theme(axis.text.x = element_text(angle = 90))


#####
## Make a combined long dataframe of Elasmobranchs, Teleosts and Human
#####
long_aq <- bind_rows(long_aq_notelas, long_aq_elas)
# make manual taxonomy for teleosts category 2
long_aq$manual_name <- if_else(is.na(long_aq$manual_name), long_aq$category2, long_aq$manual_name)
# fix proportional read counts
long_aq$prc <- long_aq$reads/long_aq$total_reads

# quick stacked bar plot of ocean display
windows()
ggplot(long_aq , aes(x = seq_id, y = reads, fill = manual_name)) + 
  geom_bar(stat = "identity", color = "black", position = "fill") +
  facet_grid(. ~ type2, scales = "free") +
  labs(x ="Sample", y = "Proportional Read Counts (PRC)") +
  theme(axis.text.x = element_text(angle = 90))

long_aq$manual_name2 <- ifelse(long_aq$manual_name=="human", "Human",
                               ifelse(long_aq$manual_name=="food", "Food",
                                      ifelse(long_aq$manual_name=="inventory", "Teleost",
                                            "Elasmobranch")))
 
# quick stacked bar plot of ocean display with manual name two
windows()
ggplot(long_aq , aes(x = seq_id, y = reads, fill = manual_name2)) + 
  geom_bar(stat = "identity", color = "black", position = "fill") +
  facet_grid(. ~ type2, scales = "free") +
  labs(x ="Sample", y = "Proportional Read Counts (PRC)") +
  theme(axis.text.x = element_text(angle = 90))


#####
## Make a combined long dataframe of Elasmobranchs, and previous combined long dataframe
#####
# make a column that you will facet_grid with
long_aq$df_from <- "Ocean Display"
long_aq_elas$df_from <- "Elasmobranchs Only"
# master long df
master_aq <- bind_rows(long_aq, long_aq_elas)
# update manual name 2
master_aq$manual_name2 <- if_else(is.na(master_aq$manual_name2), master_aq$manual_name, master_aq$manual_name2)
  
# set up factors
levels(master_aq$df_from) <- c("Elasmobranchs", "Ocean Display")
master_aq$manual_name2 <- factor(master_aq$manual_name2, levels = c("Carcharhinus melanopterus", "Carcharias taurus",
                                                                    "Chiloscyllium punctatum", "Chiloscyllium sp.",
                                                                    "Ginglymostoma cirratum", "Glaucostegus cemiculus",
                                                                    "Heterodontus sp.", "Hypanus americanus", "Orectolobus sp.",
                                                                    "Stegostoma tigrinum", "Elasmobranch", "Teleost", 
                                                                    "Food", "Human"))
master_aq$type2 <- factor(master_aq$type2, levels = c("Syringe Filter", "Diver MP", "Soak MP"))
# set colours
# colourblind safe colour palettes
# https://davidmathlogic.com/colorblind/#%23332288-%23117733-%2344AA99-%2388CCEE-%23DDCC77-%23CC6677-%23AA4499-%23882255-%234c0268-%23D99B0B

cols <- c("Carcharhinus melanopterus" = "#332288", 
          "Carcharias taurus" = "#117733",
          "Chiloscyllium punctatum" = "#44AA99", 
          "Chiloscyllium sp." = "#88CCEE",
          "Ginglymostoma cirratum" = "#D99B0B", 
          "Glaucostegus cemiculus" = "#DDCC77",
          "Heterodontus sp." = "#CC6677", 
          "Hypanus americanus" = "#AA4499", 
          "Orectolobus sp." = "#882255",
          "Stegostoma tigrinum" = "#4C0268", 
          "Elasmobranch" = "#000000", 
          "Teleost" = "#404040", 
          "Food" = "#808080", 
          "Human" = "#ffffff")


# version 1
#windows()
ggplot(master_aq, aes(x = reorder(seq_id, time), y = reads, fill = manual_name2)) + 
  geom_bar(stat = "identity", position = "fill") +
  #geom_bar(stat = "identity", color = "black", position = "fill") +
  scale_fill_manual(values = c(cols)) +
  facet_grid(df_from ~ type2, scales = "free") +
  labs(x ="Sample", y = "Proportional Read Counts (PRC)") +
  theme(axis.text.x = element_text(angle = 90)) 


master_aq <- mutate(master_aq,
                    seq_id2 = case_when(
                      replicate =="1A_RA" ~ "Dive 1A", 
                      replicate =="1B_RB" ~ "Dive 1B", 
                      replicate =="1C_DA" ~ "Dive 1C", 
                      replicate =="1D_DB" ~ "Dive 1D",
                      replicate =="2A_RA" ~ "Dive 2A", 
                      replicate =="2B_RB" ~ "Dive 2B",
                      replicate =="2C_DA" ~ "Dive 2C", 
                      replicate =="2D_DB" ~ "Dive 2D", 
                      time == 10 & replicate=="A" ~ "10 A",
                      time == 10 & replicate=="B" ~ "10 B",
                      time == 30 & replicate=="A" ~ "30 A",
                      time == 30 & replicate=="B" ~ "30 B",
                      time == 60 & replicate=="A" ~ "60 A",
                      time == 60 & replicate=="B" ~ "60 B",
                      time == 120 & replicate=="A" ~ "120 A",
                      time == 120 & replicate=="B" ~ "120 B",
                      time == 240 & replicate=="A" ~ "240 A",
                      time == 240 & replicate=="B" ~ "240 B",
                      type == "eDNA" & replicate== "A" ~ "Bottle A",
                      type == "eDNA" & replicate== "B" ~ "Bottle B",
                      type == "eDNA" & replicate== "C" ~ "Bottle C",
                      type == "eDNA" & replicate== "D" ~ "Bottle D",
                      TRUE ~ NA # This is for all other values 
                    ))                # not covered by the above.



# version 2
aq_barplot2 <-
#aq_barplot2a <-
ggplot(master_aq, aes(x = reorder(seq_id2, time), y = reads, fill = manual_name2)) + 
  #geom_bar(stat = "identity", position = "fill") + #2a
  geom_bar(stat = "identity", color = "black", position = "fill") + #2
  scale_fill_manual(values = c(cols)) +
  facet_grid(df_from ~ type2, scales = "free") +
  labs(x ="Sample", y = "Proportional Read Counts (PRC)") +
  theme(axis.text.x = element_text(angle = 90)) 
ggsave(filename = c("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_writing/figures/aq_barplot2.jpg"),
       plot = aq_barplot2,  width = 7, height = 7, units = "in")
#ggsave(filename = c("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_writing/figures/aq_barplot2a.jpg"),
#       plot = aq_barplot2a,  width = 7, height = 7, units = "in")


master_aq <- mutate(master_aq,
                  seq_id3 = case_when(
                    replicate =="1A_RA" | replicate =="1B_RB"  ~ "Dive 1A", 
                    replicate =="1C_DA" | replicate =="1D_DB"  ~ "Dive 1B", 
                    replicate =="2A_RA" | replicate =="2B_RB"  ~ "Dive 2A", 
                    replicate =="2C_DA" | replicate =="2D_DB"  ~ "Dive 2B", 
                    time == 10 ~ "10",
                    time == 30 ~ "30",
                    time == 60 ~ "60",
                    time == 120 ~ "120",
                    time == 240 ~ "240",
                    type == "eDNA" & replicate== "A" ~ "Bottle A",
                    type == "eDNA" & replicate== "B" ~ "Bottle B",
                    type == "eDNA" & replicate== "C" ~ "Bottle C",
                    type == "eDNA" & replicate== "D" ~ "Bottle D",
                    TRUE ~ NA # This is for all other values 
                  ))                # not covered by the above.


# version 3
aq_barplot3 <-
#aq_barplot3a <-  
ggplot(master_aq, aes(x = reorder(seq_id3, time), y = reads, fill = manual_name2)) + 
  #geom_bar(stat = "identity", position = "fill") + #3a
  geom_bar(stat = "identity", color = "black", size=0.075, position = "fill") + #3
  scale_fill_manual(values = c(cols)) +
  facet_grid(df_from ~ type2, scales = "free") +
  labs(x ="Sample", y = "Proportional Read Counts per MOTU") +
  theme(axis.text.x = element_text(angle = 90)) +
  theme_bw() +
  guides(fill=guide_legend(title="Taxa"))
ggsave(filename = c("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_writing/figures/aq_barplot3_BEST.jpg"),
       plot = aq_barplot3,  width = 7, height = 7, units = "in")
#ggsave(filename = c("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_writing/figures/aq_barplot3a.jpg"),
#       plot = aq_barplot3a,  width = 7, height = 7, units = "in")


#####
## Plot Motu Richness, Species richness, & reads per sample
#####
# need to add same grouping variable names to long_aq
long_aq2 <- mutate(long_aq,
                    seq_id2 = case_when(
                      replicate =="1A_RA" ~ "Dive 1A", 
                      replicate =="1B_RB" ~ "Dive 1B", 
                      replicate =="1C_DA" ~ "Dive 1C", 
                      replicate =="1D_DB" ~ "Dive 1D",
                      replicate =="2A_RA" ~ "Dive 2A", 
                      replicate =="2B_RB" ~ "Dive 2B",
                      replicate =="2C_DA" ~ "Dive 2C", 
                      replicate =="2D_DB" ~ "Dive 2D", 
                      time == 10 & replicate=="A" ~ "10 A",
                      time == 10 & replicate=="B" ~ "10 B",
                      time == 30 & replicate=="A" ~ "30 A",
                      time == 30 & replicate=="B" ~ "30 B",
                      time == 60 & replicate=="A" ~ "60 A",
                      time == 60 & replicate=="B" ~ "60 B",
                      time == 120 & replicate=="A" ~ "120 A",
                      time == 120 & replicate=="B" ~ "120 B",
                      time == 240 & replicate=="A" ~ "240 A",
                      time == 240 & replicate=="B" ~ "240 B",
                      type == "eDNA" & replicate== "A" ~ "Bottle A",
                      type == "eDNA" & replicate== "B" ~ "Bottle B",
                      type == "eDNA" & replicate== "C" ~ "Bottle C",
                      type == "eDNA" & replicate== "D" ~ "Bottle D",
                      TRUE ~ NA # This is for all other values 
                    )) 
long_aq2 <- mutate(long_aq,
                    seq_id3 = case_when(
                      replicate =="1A_RA" | replicate =="1B_RB"  ~ "Dive 1A", 
                      replicate =="1C_DA" | replicate =="1D_DB"  ~ "Dive 1B", 
                      replicate =="2A_RA" | replicate =="2B_RB"  ~ "Dive 2A", 
                      replicate =="2C_DA" | replicate =="2D_DB"  ~ "Dive 2B", 
                      time == 10 ~ "10",
                      time == 30 ~ "30",
                      time == 60 ~ "60",
                      time == 120 ~ "120",
                      time == 240 ~ "240",
                      type == "eDNA" & replicate== "A" ~ "Bottle A",
                      type == "eDNA" & replicate== "B" ~ "Bottle B",
                      type == "eDNA" & replicate== "C" ~ "Bottle C",
                      type == "eDNA" & replicate== "D" ~ "Bottle D",
                      TRUE ~ NA # This is for all other values 
                    )) 
# remove rows where reads = 0
long_aq2 <- subset(long_aq, long_aq$reads>0)

# Calculate unique Motus
# calculate the number of motus per version 2 grouping
motus_2 <- long_aq2 %>% count(seq_id2, id, sort = TRUE)
motus_2 <- motus_2 %>% count(seq_id2, sort = TRUE)
# calculate the number of motus per version 3 grouping
motus_3 <- long_aq2 %>% count(seq_id3, id, sort = TRUE)
motus_3 <- motus_3 %>% count(seq_id3, sort = TRUE)

# Calculate unique taxa
# calculate the number of taxa per version 2 grouping
name_2 <- long_aq2 %>% count(seq_id2, final_name, sort = TRUE)
name_2 <- name_2 %>% count(seq_id2, sort = TRUE)
# calculate the number of taxa per version 3 grouping
name_3 <- long_aq2 %>% count(seq_id3, final_name, sort = TRUE)
name_3 <- name_3 %>% count(seq_id3, sort = TRUE)

# Calculate reads per sample or per combined replicate
## MAYBE COME BACK TO THIS LATER
#samp_total <- as.data.frame(colSums(p2e_aq_prev[,c(3:30)]))

# Calculate the number of Elamobrachs
# subset the elasmobranchs
long_aq3 <- subset(long_aq2, long_aq2$manual_name2=="Elasmobranch")
# calculate the number of elas per version 2 grouping
elas_2 <- long_aq3 %>% count(seq_id2, manual_name, sort = TRUE)
elas_2 <- elas_2 %>% count(seq_id2, sort = TRUE)
# calculate the number of elas per version 3 grouping
elas_3 <- long_aq3 %>% count(seq_id3, manual_name, sort = TRUE)
elas_3 <- elas_3 %>% count(seq_id3, sort = TRUE)


# Combine motus and unique taxa to make alpha diversity plots
ad2 <- merge(motus_2, name_2, by="seq_id2")
ad2 <- merge(ad2, elas_2, by="seq_id2")
colnames(ad2)[2] <- "MOTUs2"
colnames(ad2)[3] <- "Richness2"
colnames(ad2)[4] <- "Elas2"
ad3 <- merge(motus_3, name_3, by="seq_id3")
ad3 <- merge(ad3, elas_3, by="seq_id3")
colnames(ad3)[2] <- "MOTUs3"
colnames(ad3)[3] <- "Richness3"
colnames(ad3)[4] <- "Elas3"
final_aq <- merge(master_aq, ad2, by="seq_id2")
final_aq <- merge(final_aq, ad3, by="seq_id3")

# make plot for version 2 sample names, with both MOTU and Reads
#library(patchwork)
#m2 <-
#ggplot(final_aq, aes(x=seq_id2, y=MOTUs)) +
#  geom_segment(aes(x=seq_id2, xend=seq_id2, y=0, yend=MOTUs), color="grey") +
#  geom_point(aes(x=seq_id2, y=MOTUs), color="orange", size=4) +
#  facet_grid( ~ type2, scales = "free") +
#  labs(x ="Sample", y = "MOTUs")

#m3 <-
#ggplot(final_aq, aes(x=seq_id2, y=Richness)) +
#  geom_segment( aes(x=seq_id2, xend=seq_id2, y=0, yend=Richness), color="grey") +
#  geom_point( color="green", size=4) +
#  facet_grid( ~ type2, scales = "free") +
#  labs(x ="Sample", y = "Unique Taxa")

# did not work how I wanted it to
#m2 + m3

# set factor levels
final_aq$seq_id2 <- factor(final_aq$seq_id2, levels=c("Bottle A", "Bottle B","Bottle C", "Bottle D",
                                                      "Dive 1A",  "Dive 1B",  "Dive 1C",  "Dive 1D",
                                                      "Dive 2A",  "Dive 2B",  "Dive 2C",  "Dive 2D",
                                                      "10 A", "10 B", "30 A", "30 B", "60 A", "60 B",     
                                                      "120 A", "120 B", "240 A", "240 B"))

# lollipop plot 2
aq_dotplot2 <-
ggplot(final_aq) +
  geom_segment(aes(x=seq_id2, xend=seq_id2, y=0, yend=MOTUs2), color="grey") +
  geom_point(aes(x=seq_id2, y=MOTUs2, color="MOTUs"), size=4) +
  geom_point(aes(x=seq_id2, y=Richness2, color="Taxa"), size=4) +
  geom_point(aes(x=seq_id2, y=Elas2, color="Elasmobranch"), size=4) +
  scale_colour_manual(values = c("#FFC20A", "#0C7BDC", "#000000"),
                      limits = c("MOTUs", "Taxa", "Elasmobranch"),
                      name = "") +
  facet_grid( ~ type2, scales = "free") +
  labs(x ="Sample", y = "") +
  theme(axis.text.x = element_text(angle = 90)) 
ggsave(filename = c("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_writing/figures/aq_dotplot2.jpg"),
       plot = aq_dotplot2,  width = 7, height = 3.5, units = "in")

# set factor levels
final_aq$seq_id3 <- factor(final_aq$seq_id3, levels=c("Bottle A", "Bottle B","Bottle C", "Bottle D",
                                                      "Dive 1A",  "Dive 1B",  "Dive 2A",  "Dive 2B",
                                                      "10", "30", "60", "120", "240"))



# lollipop plot 3
aq_dotplot3 <-
  ggplot(final_aq) +
  geom_segment(aes(x=seq_id3, xend=seq_id3, y=0, yend=MOTUs3), color="grey") +
  geom_point(aes(x=seq_id3, y=MOTUs3, color="MOTUs"), size=4) +
  geom_point(aes(x=seq_id3, y=Richness3, color="Taxa"), size=4) +
  geom_point(aes(x=seq_id3, y=Elas3, color="Elasmobranch"), size=4) +
  scale_colour_manual(values = c("#FFC20A", "#0C7BDC", "#000000"),
                      limits = c("MOTUs", "Taxa", "Elasmobranch"),
                      name = "") +
  facet_grid( ~ type2, scales = "free") +
  labs(x ="Sample", y = "") +
  theme(axis.text.x = element_text(angle = 90)) +
  theme_bw() 
ggsave(filename = c("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_writing/figures/aq_dotplot3.jpg"),
       plot = aq_dotplot3,  width = 7, height = 3.5, units = "in")




#####
## OLD CODE BELOW

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
#####
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









