###############################################################
## Repeat plots from taxa_metadata_aquarium.R for Coral cave ##
###############################################################

library(ggplot2)
library(tidyverse)
library(janitor)
library(ggthemes)

p2t_aq_prev <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/decontam/p2t_aq_decontam_001.csv")

# Separate out teleosts
p2t_aq_all <- subset(p2t_aq_prev, final_class=="Actinopterygii" | final_name=="Homo sapiens") 
p2t_aq_all <- subset(p2t_aq_all, select = -c(X.1,X)) #get rid of meaningless row name columns


# remove contamination
p2t_aq_all <- subset(p2t_aq_all, final_name!="Pangasianodon hypophthalmus") #positive control
p2t_aq_all <- subset(p2t_aq_all, final_name!="Leporidae") #hare, mislabeled at Actinopterygii


# add inventory information
meta2 <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/metadata/supp_table_3.csv")
p2t_aq_all <- merge(p2t_aq_all, meta2, by="final_name")

# make long dataframe
long_aqt_all <- p2t_aq_all %>%
  pivot_longer(cols = sample.5Et10BLUE_MPEtA_:sample.6Ft240BLUE_MPEtB, names_to = "seq_id", values_to = "reads")

# figure out what the top 10 fish are
# collapse by country
top10<- long_aqt_all %>%
  group_by(final_name, category2) %>%
  summarise(sum_read = sum(reads))
top10list <- data.frame(
  final_name = c("Gymnothorax kidako", "Abudefduf", "Platax", "Chaetodon", "Diagramma",
                  "Chrysiptera cyanea", "Arothron", "Myripristis", "Siganus", "Chaetodon auriga"),
  stringsAsFactors = FALSE)
rm(top10)

# add metadata to long dataframe
meta <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/metadata/supp_table_1.csv")
long_aqt_all <- merge(long_aqt_all, meta, by.x="seq_id")


# make a subset of the long dataframe that is just the top 10 fish
long_aqt_10 <- merge(long_aqt_all, top10list, by.x ="final_name")

# quick plots to have a look
ggplot(long_aqt_10 , aes(x = seq_id, y = reads, fill = final_name)) + 
  geom_bar(stat = "identity", color = "black", position = "fill") +
  facet_grid(. ~ via, scales = "free") +
  labs(x ="Sample", y = "Proportional Read Counts (PRC)") +
  theme(axis.text.x = element_text(angle = 90))

ggplot(long_aqt_all , aes(x = seq_id, y = reads, fill = category2)) + 
  geom_bar(stat = "identity", color = "black", position = "fill") +
  facet_grid(. ~ via, scales = "free") +
  labs(x ="Sample", y = "Proportional Read Counts (PRC)") +
  theme(axis.text.x = element_text(angle = 90))

#####
## Make a combined long dataframe of Top 10 fish, all fish and Human
#####
# make a column that you will facet_grid with
long_aqt_all$df_from <- "Coral Cave"
long_aqt_10$df_from <- "Top 10 Teleosts"
# master long df
master_aq <- bind_rows(long_aqt_all, long_aqt_10)

# set colours
cols <- c("Gymnothorax kidako" = "#332288", 
          "Abudefduf" = "#117733",
          "Platax" = "#44AA99", 
          "Chaetodon" = "#88CCEE",
          "Diagramma" = "#D99B0B", 
          "Chrysiptera cyanea" = "#DDCC77",
          "Arothron" = "#CC6677", 
          "Myripristis" = "#AA4499", 
          "Siganus" = "#882255",
          "Chaetodon auriga" = "#4C0268", 
          "Teleost" = "#000000", 
          "Food" = "#808080", 
          "Human" = "#ffffff")

# assign appropriate names
master_aq$manual_name <- ifelse(master_aq$category2=="human", "Human",
                                ifelse(master_aq$category2=="food", "Food",
                                       ifelse(master_aq$df_from=="Coral Cave", "Inventory", NA)))
# update manual name 2 EDIT THE TOP 10 FISHES NAME
master_aq$manual_name2 <- ifelse(is.na(master_aq$manual_name), master_aq$final_name, master_aq$manual_name)
master_aq$manual_name2 <- ifelse(master_aq$df_from=="Coral Cave" & (master_aq$final_name=="Gymnothorax kidako" |
                                 master_aq$final_name=="Gymnothorax kidako" |master_aq$final_name=="Gymnothorax kidako" |
                                 master_aq$final_name=="Gymnothorax kidako" |master_aq$final_name=="Gymnothorax kidako" |
                                 master_aq$final_name=="Gymnothorax kidako" |master_aq$final_name=="Gymnothorax kidako" |
                                 master_aq$final_name=="Gymnothorax kidako" |master_aq$final_name=="Gymnothorax kidako" |
                                master_aq$final_name=="Gymnothorax kidako"), "Top 10", master_aq$manual_name2))


ggplot(master_aq, aes(x = seq_id, y = reads, fill = manual_name2)) + 
  #geom_bar(stat = "identity", position = "fill") + #3a
  geom_bar(stat = "identity", color = "black", size=0.075, position = "fill") + #3
  scale_fill_manual(values = c(cols)) +
  facet_grid(df_from ~ type2, scales = "free") +
  labs(x ="Sample", y = "Proportional Read Counts per MOTU") +
  theme(axis.text.x = element_text(angle = 90)) +
  theme_bw() +
  guides(fill=guide_legend(title="Taxa"))



ggplot(master_aq, aes(x = reorder(seq_id3, time), y = reads, fill = manual_name2)) + 
  geom_bar(stat = "identity", color = "black", size=0.075, position = "fill") + #3
  scale_fill_manual(values = c(cols), 
                    labels = c("Elasmobranch", "Teleost", "Food", "Human", 
                               expression(italic("Carcharhinus melanopterus")), 
                               expression(italic("Carcharias taurus")),
                               expression(italic("Chiloscyllium punctatum")), 
                               expression(italic("Chiloscyllium")~plain("sp.")),
                               expression(italic("Ginglymostoma cirratum")), 
                               expression(italic("Glaucostegus cemiculus")),
                               expression(italic("Heterodontus")~plain("sp.")), 
                               expression(italic("Hypanus americanus")), 
                               expression(italic("Orectolobus")~plain("sp.")),
                               expression(italic("Stegostoma tigrinum")))) +
  facet_grid(df_from ~ type2, scales = "free") +
  labs(x ="Sample", y = "Proportional Read Counts") +
  theme_bw() +
  theme(legend.text.align = 0,
        strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 12, colour ="black", face = "bold"), 
        legend.text = element_text(size = 12, colour ="black"), 
        axis.text.x = element_text(colour = "black", angle = 35, hjust = 1, size = 12, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
        axis.title.y = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"),
        axis.ticks = element_line(colour ="black", linewidth = 2)) +
  guides(fill=guide_legend(title="Taxa"))





