################################################
## DNA divers - organize data for basic plots ##
################################################

## Read in raw MOTU tables
# Disclaimer: some file names have 'ASVs' but the sequences have been clustered into MOTUs
## Phase 1 - all samples with tele02 marker
p1 <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/divers_methods_experiment/sintax_taxonomy/divmeth1_ASVs.csv")
## Phase 2 
## elas02 marker
p2e_aq <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/divers_methods_experiment_part2/sintax_output_MainTank/e_st_divmeth2_ASVs.csv") #shark tank aquarium
p2e_uk <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/divers_methods_experiment_part2/sintax_output_UKreferencedatabase/e_divmeth2_ASVs.csv") #UK samples
p2e_sk <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/St.Kilda_andothersamples/Elas02_run/other_sintax/elas_other_ASVs.csv") #UK samples accidently not sequenced/added on to a different project
## tele02 marker
p2t_aq <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/divers_methods_experiment_part2/sintax_output_CoralCavedatabase/t_cc_divmeth2_ASVs.csv") #coral cave aquarium
p2t_uk <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/DNAdivers/divers_methods_experiment_part2/sintax_output_UKreferencedatabase/t_divmeth2_ASVs.csv") #UK samples
