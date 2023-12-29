################################################
## Add Taxonomic metadata to aquarium samples ##
################################################
## Shark tank
# total reads per MOTU
p2e_aq_v3 <- p2e_aq_v3 %>% mutate(total_reads = rowSums(.[2:28]))

# subset elasmobranchs
p2e_aq_elas <- subset(p2e_aq_v3, class=="Elasmobranchii")

# subset species identity >0.95
p2e_aq_elas95 <- subset(p2e_aq_elas, s.id>0.95)

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

## quick plot of sharks
windows()
ggplot(elas95samples , aes(x = long_id, y = prc, fill = species)) + 
  geom_bar(stat = "identity", color = "black", position = "fill") +
  theme(axis.text.x = element_text(angle = 90))
  

# subset species identity >0.95
p2e_aq_95 <- subset(p2e_aq_v3, s.id>0.95)

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

## quick plot of everything
windows()
ggplot(all95samples , aes(x = long_id, y = reads, fill = species)) + 
  geom_bar(stat = "identity", color = "black", position = "fill") +
  theme(axis.text.x = element_text(angle = 90))






