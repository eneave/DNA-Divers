############################
## Descriptive statistics ##
############################

# Read in decontaminated data
# note that human reads and controls are still included
p1_prev <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/decontam/p1_decontam.csv")
p2e_aq_prev <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/decontam/p2e_aq_decontam.csv")
p2t_aq_prev <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/decontam/p2t_aq_decontam.csv")
p2t_uk_prev <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/decontam/p2t_uk_decontam.csv")
p2e_uk_prev <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/decontam/p2e_uk_decontam.csv")
p3_prev <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/decontam/p3_decontam.csv")
p4_prev <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/decontam/p4_decontam.csv")




l1 <- orkt_motu %>% pivot_longer(cols = sample.7A:sample.9E,
                                 names_to = "sample",
                                 names_prefix = "sample.",
                                 values_to = "reads")