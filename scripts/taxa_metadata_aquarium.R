################################################
## Add Taxonomic metadata to aquarium samples ##
################################################
#####
## Shark tank
#####
# total reads per MOTU
p2e_aq_prev <- p2e_aq_prev %>% mutate(total_reads = rowSums(.[2:28]))

#####
# subset elasmobranchs
#####
p2e_aq_elas <- subset(p2e_aq_prev, class=="Elasmobranchii")

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

## quick stacked bar plot of sharks
ggplot(elas95samples , aes(x = long_id, y = prc, fill = manual_taxo)) + 
  geom_bar(stat = "identity", color = "black", position = "fill") +
  facet_grid(. ~ type2, scales = "free") +
  labs(title="Proportional Read counts of Elasmobranchs \ndetected in the Main Tank by different methods",
       x ="Sample", y = "Proportional Read Counts (PRC)") +
  theme(axis.text.x = element_text(angle = 90))

#####
# explore data including all detentions
#####

# subset species identity >0.95
p2e_aq_95 <- subset(p2e_aq_prev, s.id>0.95)

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

## quick plot of everything
ggplot(all95samples , aes(x = long_id, y = prc, fill = manual_taxo1)) + 
  geom_bar(stat = "identity", color = "black", position = "fill") +
  facet_grid(. ~ type2, scales = "free") +
  labs(title="Proportional Read counts of Elasmobranchs & Teleosts \ndetected in the Main Tank by different methods",
       x ="Sample", y = "Proportional Read Counts (PRC)") +
  theme(axis.text.x = element_text(angle = 90))

#####
## Coral cave tank
#####
# FUNKY HEATMAPS TUTORIAL FOR LATER 
# https://www.youtube.com/watch?v=9XbxL-Is22k











