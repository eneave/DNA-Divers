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

# rotate so that species are rows
mps3 <- mps2 %>%
  rotate_df(mps2) %>%
  row_to_names(row_number = 1) 
colnames(mps3)[1] <- "seq_id"

# prepare the environmental variables
mps4 <- merge(mps3, meta, by.x="seq_id")
env <- mps4[c("latitude", "ocean_basin")]

# extract abundance data and prepare for MDS
dat <- mps4[c(2:213)]
# convert from character to numeric dataframe
dat <- as.data.frame(sapply(dat, as.numeric)) 

# Hellinger's standardization
dat_hell <- decostand(dat, method = "hellinger")

# squaroot transformation
dat_sq <- sqrt(dat)

# presence-absence data
dat_pa <- as.data.frame(ifelse(dat > 0, 1, 0))

# CCA
#v1
mps.ca.hell <- cca(dat_hell,env)

plot(mps.ca.hell)

#v2
mps.ca <- cca(dat_sq, env)

print(mps.ca, digits = 3)
eigenvals(mps.ca)
scores(mps.ca, display = "sites")

plot(mps.ca)

#v3
mps.ca.pa <- cca(dat_pa,env)

plot(mps.ca.pa)

# This tutorial shows how to extract species scores
# https://medium.com/@saurav12das/cca-plot-using-ggplot2-125159f13bbd

library(ggplot2)
library(ggrepel)

#extracting the data as data frame; env data
veg_1 = as.data.frame(mps.ca$CCA$biplot)
veg_1["env"] = row.names(veg_1)

#extracting the data; taxa
veg_2 = as.data.frame(mps.ca$CCA$v)
veg_2["genus"] = row.names(veg_2)

#extracting the data; sitea
veg_3 = as.data.frame(mps.ca$CCA$u)
veg_3["seq_id"] = mps4[c(1)]
veg_3 <- merge(veg_3, meta, by.x="seq_id")


#rm(veg_1, veg_2)

ggplot() +
  geom_point(data = veg_2, aes(x = CCA1, y = CCA2), color = "red") +
  geom_point(data = veg_1, aes(x = CCA1, y = CCA2), color = "blue") +
  geom_text_repel(data = veg_2,
                  aes(x = CCA1, y = CCA2, label = veg_2$genus),
                  nudge_y = -0.05) +
  theme_bw() +
  geom_segment(
    data = veg_1,
    aes(
      x = 0,
      y = 0,
      xend = CCA1,
      yend = CCA2
    ),
    arrow = arrow(length = unit(0.25, "cm"))
  ) +
  geom_text_repel(
    data = veg_1,
    aes(x = CCA1, y = CCA2, label = veg_1$env),
    nudge_y = -0.05,
    color = "blue",
    size = 5
  ) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))

# see what happens when you extract site scores

ggplot() +
  geom_point(data = veg_3, aes(x = CCA1, y = CCA2), color = "red") +
  geom_point(data = veg_1, aes(x = CCA1, y = CCA2), color = "blue") +
  geom_text_repel(data = veg_3,
                  aes(x = CCA1, y = CCA2, label = site),
                  nudge_y = -0.05,
                  max.overlaps = 100) +
  theme_bw() +
  geom_segment(
    data = veg_1,
    aes(
      x = 0,
      y = 0,
      xend = CCA1,
      yend = CCA2
    ),
    arrow = arrow(length = unit(0.25, "cm"))
  ) +
  geom_text_repel(
    data = veg_1,
    aes(x = CCA1, y = CCA2, label = veg_1$env),
    nudge_y = -0.05,
    color = "blue",
    size = 5
  ) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))


# change aestetics to see if the plot becomes more meaningful

jitter <- position_jitter(width = 0.2, height = 0.2)

cca_rough <-
ggplot() +
  geom_point(data = veg_3, aes(x = CCA1, y = CCA2, color = country),
             size = 5, alpha=0.2, position = jitter) +
  geom_point(data = veg_2, aes(x = CCA1, y = CCA2), color = "black",
             size = 1) +
  geom_text_repel(data = veg_2,
                  aes(x = CCA1, y = CCA2, label = veg_2$genus),
                  nudge_y = -0.05,
                  max.overlaps = 40) +
  #geom_text_repel(data = veg_3,
  #                aes(x = CCA1, y = CCA2, label = site),
  #                nudge_y = -0.05,
  #                max.overlaps = 75) +
  theme_bw() +
  geom_segment(
    data = veg_1,
    aes(
      x = 0,
      y = 0,
      xend = CCA1,
      yend = CCA2
    ),
    arrow = arrow(length = unit(0.25, "cm"))
  ) +
  geom_text_repel(
    data = veg_1,
    aes(x = CCA1, y = CCA2, label = veg_1$env),
    nudge_y = -0.05,
    color = "blue",
    size = 5
  ) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))

ggsave(filename=c("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_writing/figures/cca_rough.jpg"), 
       plot = cca_rough, width = 8, height = 6.5, units = "in")

#####
# repeat CCA for just the North Atlantic samples
#####

# prepare the environmental variables
dat_na <- subset(mps4, mps4$latitude>50)
#env_na <- dat_na[c("latitude", "ices_ecor", "ices_area")]
env_na <- dat_na[c("latitude", "ices_ecor")]


# extract abundance data and prepare for MDS
na_samples <- as.data.frame(dat_na$seq_id)
dat_na <- dat_na[c(2:213)]
# convert from character to numeric dataframe
dat_na <- as.data.frame(sapply(dat_na, as.numeric)) 
# remove columns which add to zero
dat_naf <- dat_na[colSums(dat_na == 0) != nrow(dat_na)]
# keep only species-level assignments
dat_nas <- dat_naf[c(10:86)]

# check rowsums, remove any samples with no reads
rowSums(dat_nas) #need to remove row 19
dat_nas2 <- dat_nas[c(-19),]
na_samples2 <- na_samples[c(-19),]
env_na2 <- env_na[c(-19),]
#rm(env_na)

#rownames(dat_nas2) <- c(1:93) 
#rownames(env_na2) <- c(1:93) 

# squaroot transformation
dat_nas_sq <- sqrt(dat_nas2)

# cca
mps.nas.ca <- cca(dat_nas_sq, env_na2)
plot(mps.nas.ca)

#extracting the data as data frame; env data
veg_1na = as.data.frame(mps.nas.ca$CCA$biplot)
veg_1na["env"] = row.names(veg_1na)
sqrt((veg_1na$CCA1)**2+(veg_1na$CCA2)**2)
veg_1naf = veg_1na[c(1,4,5),]

#extracting the data; taxa
veg_2na = as.data.frame(mps.nas.ca$CCA$v)
veg_2na["genus"] = row.names(veg_2na)

#extracting the data; sitea
veg_3na = as.data.frame(mps.nas.ca$CCA$u)
veg_3na["seq_id"] = na_samples2
veg_3na <- merge(veg_3na, meta, by.x="seq_id")

# plot with ggplot2

ggplot() +
  geom_point(data = veg_3na, aes(x = CCA1, y = CCA2, color = ices_area),
             size = 5, alpha=0.2, position = jitter) +
  geom_point(data = veg_2na, aes(x = CCA1, y = CCA2), color = "black",
             size = 1) +
  geom_text_repel(data = veg_2na,
                  aes(x = CCA1, y = CCA2, label = veg_2na$genus),
                  nudge_y = -0.05,
                  max.overlaps = 5) +
  theme_bw() +
  geom_segment(data = veg_1na, aes(
      x = 0,
      y = 0,
      xend = CCA1,
      yend = CCA2
    ),
    arrow = arrow(length = unit(0.25, "cm"))
  ) +
  geom_text_repel(
    data = veg_1na,
    aes(x = CCA1, y = CCA2, label = veg_1na$env),
    nudge_y = -0.05,
    color = "blue",
    size = 5
  ) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))



##########
## NMDS ## #could not reach a solution due to stress
##########

# just have a looksies

# Calculate distances
# Bray-Curtis dissimilarity matrix
bray_dat <- vegdist(dat_hell, method = "bray")
# Calculate jaccard dissimilariy matrix
jac_dat <- vegdist(dat_pa, method = "jaccard", binary =  TRUE)

# could not reach solution in 1000 tries
set.seed(432)
ord_bray <- metaMDS(bray_dat, distance = "bray", trymax = 1000)
plot(ord_bray)
# could not reach solution in 1000 tries
set.seed(92)
ord_jac <- metaMDS(jac_dat, distance = "jaccard", trymax = 1000)
plot(ord_jac)

# try for just North Atlantic
dat_hell_na <- decostand(dat_nas2, method = "hellinger")

# Bray-Curtis dissimilarity matrix
bray_dat_na <- vegdist(dat_hell_na, method = "bray")

set.seed(32)
ord_bray_na <- metaMDS(bray_dat_na, distance = "bray", trymax = 1000)
plot(ord_bray_na)

#########
## PCA ## # Not a good way to visualize community data
#########

pca_all <- rda(dat_hell)

summary(pca_all)
screeplot(pca_all)
plot(pca_all)

# extract data
pca_site <- as.data.frame(scores(pca_all, choices=1:2, "sites"))
pca_site["seq_id"] <- mps4[c(1)]
pca_site <- merge(pca_site, meta, by.x="seq_id")

pca_spec <- as.data.frame(scores(pca_all, choices=1:2, "species"))
pca_spec$final_name <- row.names(pca_spec)


# plot with ggplot2

ggplot() +
  geom_point(data = pca_site, aes(x = PC1, y = PC2, color = country),
             size = 5, alpha=0.2) +
  geom_point(data = pca_spec, aes(x = PC1, y = PC2), color = "green",
            size = 1) +
  geom_text(data = pca_spec, aes(x = PC1, y = PC2, label = final_name))


##########
## PCoA ##
##########

# Calculate distances
# Bray-Curtis dissimilarity matrix
bray_dat <- vegdist(dat_hell, method = "bray")
# Calculate jaccard dissimilariy matrix
jac_dat <- vegdist(dat_pa, method = "jaccard", binary =  TRUE)

cmdscale(bray_dat) %>% head()

cmdscale(bray_dat) %>% plot()
cmdscale(jac_dat) %>% plot()

#https://www.youtube.com/watch?v=G5Qckqq5Erw
