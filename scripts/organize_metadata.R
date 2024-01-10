##########################################################
## DNA divers - organize metadata for tables & analysis ##
##########################################################
library(tidyverse)

# load data
# sequence run 1 - all UK tele02 primer
p1_meta <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/divers_methods_experiment/p1_meta.csv")
colnames(p1_meta)[1] <- "seq_id"
p1_meta$run <- "1"
#sequence run 2
# shark tank elas02 primer
p2e_aq_meta <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/divers_methods_experiment_part2/p2e_aq_meta.csv")
colnames(p2e_aq_meta)[1] <- "seq_id"
p2e_aq_meta$run <- "2"
# coral cave tele02 primer
p2t_aq_meta <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/divers_methods_experiment_part2/p2t_aq_meta.csv")
colnames(p2t_aq_meta)[1] <- "seq_id" 
p2t_aq_meta$run <- "2"
# UK tele02 primer
p2t_uk_meta <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/divers_methods_experiment_part2/p2t_uk_meta.csv")
colnames(p2t_uk_meta)[1] <- "seq_id"
p2t_uk_meta$run <- "2"
# UK elas02 primer
p2e_uk_meta <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/divers_methods_experiment_part2/p2e_uk_meta.csv") 
colnames(p2e_uk_meta)[1] <- "seq_id"
p2e_uk_meta$run <- "2"
# sequence run 3 - all elas02 primer shared with another project
p3_meta <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/St.Kilda_andothersamples/Elas02_run/p3_meta.csv") 
colnames(p3_meta)[1] <- "seq_id"
p3_meta$run <- "3"


# combine
all_meta <- rbind(p1_meta, p2e_aq_meta, p2e_uk_meta, p2t_aq_meta, p2t_uk_meta, p3_meta)

# just samples
samp_meta <- subset(all_meta, all_meta$sampletype1=="sample")
write.csv(samp_meta, "C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/sample_metadata.csv")

rm(p1_meta, p2e_aq_meta, p2e_uk_meta, p2t_aq_meta, p2t_uk_meta, p3_meta)

