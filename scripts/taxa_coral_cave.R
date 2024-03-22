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
          "Top 10" = "#000000", 
          "Food" = "#808080", 
          "Human" = "#ffffff")

# assign appropriate names
master_aq$manual_name <- ifelse(master_aq$category2=="human", "Human",
                                ifelse(master_aq$category2=="food", "Food",
                                       ifelse(master_aq$df_from=="Coral Cave", "Inventory", NA)))
# update manual name 2 EDIT THE TOP 10 FISHES NAME
master_aq$manual_name2 <- ifelse(is.na(master_aq$manual_name), master_aq$final_name, master_aq$manual_name)
master_aq$manual_name2 <- ifelse(master_aq$df_from=="Coral Cave" & (master_aq$final_name=="Gymnothorax kidako" |
                                 master_aq$final_name=="Abudefduf" |master_aq$final_name=="Platax" |
                                 master_aq$final_name=="Chaetodon" |master_aq$final_name=="Diagramma" |
                                 master_aq$final_name=="Chrysiptera cyanea" |master_aq$final_name=="Arothron" |
                                 master_aq$final_name=="Myripristis" |master_aq$final_name=="Siganus" |
                                master_aq$final_name=="Chaetodon auriga"), "Top 10", master_aq$manual_name2)

supp_fig_coralcave <-
ggplot(master_aq, aes(x = seq_id, y = reads, fill = manual_name2)) + 
  geom_bar(stat = "identity", color = "black", size=0.075, position = "fill") + #3
  scale_fill_manual(values = c(cols), breaks = c("Top 10", "Food", "Human", "Gymnothorax kidako",
                                                 "Abudefduf", "Platax", "Chaetodon", "Diagramma",
                                                 "Chrysiptera cyanea", "Arothron", "Myripristis",
                                                 "Siganus", "Chaetodon auriga"),
                    labels = c("Top 10 Teleosts", "Food", "Human", 
                               expression(italic("Gymnothorax kidako")), 
                               expression(italic("Abudefduf")~plain("sp.")),
                               expression(italic("Platax")~plain("sp.")), 
                               expression(italic("Chaetodon")~plain("sp.")),
                               expression(italic("Diagramma")~plain("sp.")), 
                               expression(italic("Chrysiptera cyanea")),
                               expression(italic("Arothron")~plain("sp.")), 
                               expression(italic("Myripristis")~plain("sp.")), 
                               expression(italic("Siganus")~plain("sp.")),
                               expression(italic("Chaetodon auriga")))) +
  scale_x_discrete(labels=c("sample.5Et10BLUE_MPEtA_" = "10",
                            "sample.5Ft10BLUE_MPEtB_" = "10",
                            "sample.5Gt30BLUE_MPEtA_" = "30",
                            "sample.5Ht30BLUE_MPEtB_" = "30",
                            "sample.6At60BLUE_MPEtA_" = "60",
                            "sample.6Bt60BLUE_MPEtB_" = "60",
                            "sample.6Ct120BLUE_MPEtA" = "120",
                            "sample.6Dt120BLUE_MPEtB" = "120",
                            "sample.6Et240BLUE_MPEtA" = "240",
                            "sample.6Ft240BLUE_MPEtB" = "240")) +
  facet_grid(df_from ~ via, scales = "free") +
  labs(x ="Time", y = "Proportional Read Counts per MOTU") +
  theme_bw() +
  theme(legend.text.align = 0) +
  guides(fill=guide_legend(title="Taxa"))

ggsave(filename = c("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_writing/figures/sfig_coralcavebar.jpg"),
       plot = supp_fig_coralcave,  width = 7, height = 7, units = "in")


#####
## Make lollipop plot for coral cave
#####

# remove rows where reads = 0
long_aqt_2 <- subset(long_aqt_all, long_aqt_all$reads>0)

# Calculate unique Motus
# calculate the number of motus per replicate
motus_t <- long_aqt_2 %>% count(seq_id, id, sort = TRUE)
motus_t <- motus_t %>% count(seq_id, sort = TRUE)

# Calculate unique taxa
# calculate the number of taxa per replicate
name_t <- long_aqt_2 %>% count(seq_id, final_name, sort = TRUE)
name_t <- name_t %>% count(seq_id, sort = TRUE)

# make dataframe
# Combine motus and unique taxa to make alpha diversity plots
ct <- merge(motus_t, name_t, by="seq_id")
ct <- mutate(ct, seq_id3 = case_when(
                     seq_id =="sample.5Et10BLUE_MPEtA_" | seq_id =="sample.5Ft10BLUE_MPEtB_"  ~ 10, 
                     seq_id =="sample.5Gt30BLUE_MPEtA_" | seq_id =="sample.5Ht30BLUE_MPEtB_"  ~ 30, 
                     seq_id =="sample.6At60BLUE_MPEtA_" | seq_id =="sample.6Bt60BLUE_MPEtB_"  ~ 60,  
                     seq_id =="sample.6Ct120BLUE_MPEtA" | seq_id =="sample.6Dt120BLUE_MPEtB"  ~ 120, 
                     seq_id =="sample.6Et240BLUE_MPEtA" | seq_id =="sample.6Ft240BLUE_MPEtB"  ~ 240, 
                     TRUE ~ NA # This is for all other values 
                   )) 
colnames(ct)[2] <- "MOTUs2"
colnames(ct)[3] <- "Richness2"


# combine dataframe from shark tank data
head(ad2)

ad2 <- mutate(ad2, seq_id3 = case_when(
                                  seq_id2=="10 A" ~ 10,
                                  seq_id2=="10 B" ~ 10,
                                  seq_id2=="30 A" ~ 30,
                                  seq_id2=="30 B" ~ 30,
                                  seq_id2=="60 A" ~ 60,
                                  seq_id2=="60 B" ~ 60,
                                  seq_id2=="120 A" ~ 120,
                                  seq_id2=="120 B" ~ 120,
                                  seq_id2=="240 A" ~ 240,
                                  seq_id2=="240 B" ~ 240,
                                  seq_id2=="Dive 1A" ~ 50,
                                  seq_id2=="Dive 1B" ~ 50,
                                  seq_id2=="Dive 1C" ~ 50,
                                  seq_id2=="Dive 1D" ~ 50,
                                  seq_id2=="Dive 2A" ~ 65,
                                  seq_id2=="Dive 2B" ~ 65,
                                  seq_id2=="Dive 2C" ~ 65,
                                  seq_id2=="Dive 2D" ~ 65,
                                  TRUE ~ NA # This is for all other values 
                                  )) 

time_ex <- rbind(ct[c(2:4)], ad2[c(2,3,5)])
time_ex <- 
time_ex %>% drop_na()
time_ex$location <- rep(c("Ocean Exhibit","Coral Cave","Ocean Exhibit"),times=c(10,10,8))

time_ex2 <- time_ex[c(1:20),]
#time_ex2$location <- rep(c("Ocean Exhibit","Coral Cave"),each=10)

# not including dives
time_lm <-
ggplot() +
  geom_point(data = time_ex2, aes(x = seq_id3, y = MOTUs2, color="MOTUs", shape=location), size = 4, alpha=0.7) +
  geom_point(data = time_ex2, aes(x = seq_id3, y = Richness2, color="Taxa", shape=location), size = 4, alpha=0.7) +
  scale_colour_manual(values = c("#FFC20A", "#0C7BDC"),
                      limits = c("MOTUs", "Taxa"),
                      name = "") +
  geom_smooth(data = time_ex2, aes(x = seq_id3, y = MOTUs2, color="MOTUs"),method="lm") +
  geom_smooth(data = time_ex2, aes(x = seq_id3, y = Richness2, color="Taxa"),method="lm") +
  labs(x="Time", y="") +
  theme_bw() +
  theme(axis.text.y = element_text(colour = "black", size = 16, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 16), 
        legend.text = element_text(size = 16, colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 18), 
        axis.title.x = element_text(face = "bold", size = 18, colour = "black"), 
        legend.title = element_text(size = 18, colour = "black", face = "bold"))

ggsave(filename=c("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_writing/figures/time_lm.jpg"), 
       plot = time_lm, width = 8, height = 6.5, units = "in")



lm_tm = lm(MOTUs2 ~ seq_id3, time_ex2)
summary(lm_tm)

lm_tr = lm(Richness2 ~ seq_id3, time_ex2)
summary(lm_tr)

lm_tmr = lm(MOTUs2 + Richness2 ~ seq_id3, time_ex2)
summary(lm_tmr)


time_lm_dv <-
  ggplot() +
  geom_point(data = time_ex, aes(x = seq_id3, y = MOTUs2, color="MOTUs", shape=location), size = 4, alpha=0.7) +
  geom_point(data = time_ex, aes(x = seq_id3, y = Richness2, color="Taxa", shape=location), size = 4, alpha=0.7) +
  scale_colour_manual(values = c("#FFC20A", "#0C7BDC"),
                      limits = c("MOTUs", "Taxa"),
                      name = "") +
  geom_smooth(data = time_ex, aes(x = seq_id3, y = MOTUs2, color="MOTUs"),method="lm") +
  geom_smooth(data = time_ex, aes(x = seq_id3, y = Richness2, color="Taxa"),method="lm") +
  geom_vline(xintercept=50, linetype="dashed", color = "red", size=1) + #indicating data from dives
  geom_vline(xintercept=65, linetype="dashed", color = "red", size=1) +
  labs(x="Time", y="") +
  theme_bw() +
  theme(axis.text.y = element_text(colour = "black", size = 16, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 16), 
        legend.text = element_text(size = 16, colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 18), 
        axis.title.x = element_text(face = "bold", size = 18, colour = "black"), 
        legend.title = element_text(size = 18, colour = "black", face = "bold"))

ggsave(filename=c("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_writing/figures/time_lm_dv.jpg"), 
       plot = time_lm_dv, width = 8, height = 6.5, units = "in")



lm_tmd = lm(MOTUs2 ~ seq_id3, time_ex)
summary(lm_tmd)

lm_trd = lm(Richness2 ~ seq_id3, time_ex)
summary(lm_trd)

