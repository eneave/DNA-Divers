#########
## CCA ##
#########

# load data
# sample metadata
meta <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/metadata/supp_table_1.csv")
# MOTU table of diver MP from nature
mps <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/final_motu_mp_nature_cca.csv")
mps <- subset(mps, select = -X)

# remove 6 reads salmon from California sample
mps[mps$final_name=="Salmo salar", "sample.7G_CA_diveA_a"] <- 0

library(vegan)
library(tidyverse)
library(sjmisc)
library(janitor)

# remove samples that have no data (once human reads and contam removed)
mps2 <- mps[, colSums(mps != 0) > 0]

mps2 <- mps2[c(4,6:ncol(mps2))]
mps3 <- mps2 %>%
  rotate_df(mps2) %>%
  row_to_names(row_number = 1) 
colnames(mps3)[1] <- "seq_id"

# extract abundance data and prepare for MDS
dat <- mps3[,2:ncol(mps3)]
# convert from character to numeric dataframe
dat <- as.data.frame(sapply(dat, as.numeric)) 

# Hellinger's standardization
dat_hell <- decostand(dat, method = "hellinger")

# squaroot transformation
dat_sq <- sqrt(dat)

# cca
mps.ca.hell = cca(dat_hell)

plot(mps.ca.hell)


mps.ca = cca(dat_sq)

print(mps.ca, digits = 3)
eigenvals(mps.ca)
scores(mps.ca, display = "sites")

plot(mps.ca)

# This looks like a much better tutorial
# https://medium.com/@saurav12das/cca-plot-using-ggplot2-125159f13bbd
