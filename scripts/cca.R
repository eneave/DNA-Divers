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
#env <- mps4[c("latitude", "ocean_basin")]

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

#####
## cca with model building
#####

# prepare all environmental variables for model building
envm <- mps4[c(215:239)]
# remove incomplete/redundant variables
envm <- subset(envm, select = -c(type,replicate,PCR_primer_well,weight_g,
                time,Diver,Partnership,Main.Contact,ices_ecor,ices_area, ices_ea,
                location, site.name, collection_date, country, longitude))
glimpse(envm)

# this will be used as the scope in the add1() command when model building
mpam = cca(dat_sq ~ ., data = envm)

# make a blank model (containing an intercept only) to be used at the base for the model building
mpam.cca = cca(dat_sq ~ 1, data = envm)

# use add1 command to see if any variables are likely to significantly affect the model
add1(mpam.cca, scope = formula(mpam), test = "permutation")

# add ocean
mpam.cca = cca(dat_sq ~ ocean_basin, data = envm)

# repeat
add1(mpam.cca, scope = formula(mpam), test = "permutation")

# add latitude
mpam.cca = cca(dat_sq ~ ocean_basin + latitude, data = envm)

# repeat
add1(mpam.cca, scope = formula(mpam), test = "permutation")

# add via
mpam.cca = cca(dat_sq ~ ocean_basin + latitude + via, data = envm)
# no more variables significant

anova(mpam.cca, by = "terms")
#Model: cca(formula = dat_sq ~ ocean_basin + latitude + via, data = envm)
#Df ChiSquare      F Pr(>F)    
#ocean_basin   2    1.9760 5.0657  0.001 ***
#  latitude      1    0.8184 4.1962  0.001 ***
#  via           1    0.3591 1.8414  0.002 ** 
#  Residual    107   20.8691  

# BUT, snorkel samples were only from one location, so remove from model
mpam.cca = cca(dat_sq ~ ocean_basin + latitude, data = envm)

anova(mpam.cca, by = "terms")
#Model: cca(formula = dat_sq ~ ocean_basin + latitude, data = envm)
#Df ChiSquare      F Pr(>F)    
#ocean_basin   2    1.9760 5.0265  0.001 ***
#  latitude      1    0.8184 4.1638  0.001 ***
#  Residual    108   21.2283          

plot(mpam.cca)

# extract and plot data
library(ggplot2)
library(ggrepel)

#extracting the data as data frame; env data
mpam_1 = as.data.frame(mpam.cca$CCA$biplot)
mpam_1["env"] = row.names(mpam_1)

#extracting the data; taxa
mpam_2 = as.data.frame(mpam.cca$CCA$v)
mpam_2["taxa"] = row.names(mpam_2)

#extracting the data; sites
mpam_3 = as.data.frame(mpam.cca$CCA$u)
mpam_3["seq_id"] = mps4[c(1)]
mpam_3 <- merge(mpam_3, meta, by.x="seq_id")


# plot with ggplot2

# use jitter to make the points more visable
jitter <- position_jitter(width = 0.2, height = 0.2)
jitter2 <- position_jitter(width = 0.1, height = 0.1)


cca_all <-
  ggplot() +
  geom_point(data = mpam_3, aes(x = CCA1, y = CCA2, color = country),
             size = 5, alpha=0.4, position = jitter) +
  scale_colour_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                                 "#0072B2", "#D55E00", "#CC79A7")) +
  geom_point(data = mpam_2, aes(x = CCA1, y = CCA2), color = "black",
             size = 1, alpha=0.4, position = jitter2) +
  geom_text_repel(data = mpam_2, aes(x = CCA1, y = CCA2, label = mpam_2$taxa),
                  nudge_y = -0.05, max.overlaps = 30, size = 3) +
  theme_bw() +
  geom_segment(data = mpam_1, aes(x = 0,y = 0,xend = CCA1,yend = CCA2),
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_text_repel(data = mpam_1,aes(x = CCA1, y = CCA2, label = c("N. Atlantic", "N. Pacific", "latitude")),
    nudge_y = -0.05, color = "black", size = 4.5) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12))

ggsave(filename=c("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_writing/figures/cca_all.jpg"), 
       plot = cca_all, width = 8, height = 6.5, units = "in")

#####
## just N. Atlantic sites
#####

# prepare the environmental variables
dat_na <- subset(mps4, mps4$latitude>50)
dat_na_n <- subset(mps4, mps4$latitude>50) # keep for names

# prepare all environmental variables for model building
envam <- dat_na[c(215:239)]
# remove incomplete/redundant variables
envam <- subset(envam, select = -c(type,replicate,PCR_primer_well,weight_g, weight_range_mg,
                                 time,Diver,Partnership,Main.Contact,ices_ea, ocean_basin,
                                 location, site.name, collection_date, country, longitude))
glimpse(envam)


# extract abundance data and prepare for MDS; remove Amphipiron genus & Gnathanodon speciosus (3 reads)
dat_na <- dat_na[c(2:4,6:88,90:213)]
# convert from character to numeric dataframe
dat_na <- as.data.frame(sapply(dat_na, as.numeric)) 
# remove columns which add to zero
dat_naf <- dat_na[colSums(dat_na == 0) != nrow(dat_na)]
# keep only species-level assignments
#dat_nas <- dat_naf[c(9:86)]

# check rowsums, remove any samples with no reads
#rowSums(dat_nas) #need to remove row 19
#dat_nas2 <- dat_nas[c(-19),]
rowSums(dat_naf)
#envam <- envam[c(-19),]

# squaroot transformation
#dat_nas_sq <- sqrt(dat_nas2)
#dat_na_sq <- sqrt(dat_na2)
dat_naf_sq <- sqrt(dat_naf)

# build the model
# this will be used as the scope in the add1() command when model building
#mpna = cca(dat_nas_sq ~ ., data = envam)
mpna = cca(dat_naf_sq ~ ., data = envam)

# make a blank model (containing an intercept only) to be used at the base for the model building
#mpna.cca = cca(dat_nas_sq ~ 1, data = envam)

mpna.cca = cca(dat_naf_sq ~ 1, data = envam)

# use add1 command to see if any variables are likely to significantly affect the model
add1(mpna.cca, scope = formula(mpna), test = "permutation")

# add ices area
#mpna.cca = cca(dat_nas_sq ~ ices_area, data = envam)
mpna.cca = cca(dat_naf_sq ~ ices_area, data = envam)

# repeat
add1(mpna.cca, scope = formula(mpna), test = "permutation")

# model fits better when latitude is also added
#mpna.cca = cca(dat_nas_sq ~ ices_area + latitude, data = envam)
mpna.cca = cca(dat_naf_sq ~ ices_area + latitude, data = envam)

# run stats on the variables in the model
anova(mpna.cca, by = "terms")
#Model: cca(formula = dat_naf_sq ~ ices_area + latitude, data = envam)
#Df ChiSquare      F Pr(>F)    
#ices_area  5    3.0978 5.1515  0.001 ***
#  latitude   1    0.1656 1.3772  0.190    
#Residual  85   10.2229          


#extracting the data as data frame; env data
mpna_1 = as.data.frame(mpna.cca$CCA$biplot)
mpna_1["env"] = row.names(mpna_1)

#extracting the data; taxa
mpna_2 = as.data.frame(mpna.cca$CCA$v)
mpna_2["taxa"] = row.names(mpna_2)

#extracting the data; sites
mpna_3 = as.data.frame(mpna.cca$CCA$u)
mpna_3["seq_id"] = dat_na_n[c(1)]
mpna_3 <- merge(mpna_3, meta, by.x="seq_id")

# fix assignment of Platichthys 
mpna_2$taxa <- ifelse(mpna_2$taxa=="Platichthys stellatus", "Platichthys flesus", mpna_2$taxa)
# fix assignment of Chirolophis
mpna_2$taxa <- ifelse(mpna_2$taxa=="Chirolophis japonicus", "Chirolophis ascanii", mpna_2$taxa)
# fix assignment of Limanda - we can deduce that it is common dab
mpna_2$taxa <- ifelse(mpna_2$taxa=="Limanda", "Limanda limanda", mpna_2$taxa)

# plot with ggplot2

# use jitter to make the points more visable
jitter <- position_jitter(width = 0.2, height = 0.2)
jitter2 <- position_jitter(width = 0.05, height = 0.05)

# arrow labels
al1 <- mpna_1[c(1,2,3),]
al2 <- mpna_1[6,]

cca_na <-
  ggplot() +
  geom_point(data = mpna_3, aes(x = CCA1, y = CCA2, color = location),
             size = 5, alpha=0.6, position = jitter) +
  scale_colour_manual(values = c("#332288", "#117733", "#44AA99", "#88CCEE", 
                                 "#DDCC77", "#CC6677", "#AA4499", "#882255")) +
  geom_point(data = mpna_2, aes(x = CCA1, y = CCA2), color = "black",
             size = 1, alpha=0.4) +
  geom_text_repel(data = mpna_2, aes(x = CCA1, y = CCA2, label = mpna_2$taxa),
                  nudge_y = -0.05, max.overlaps = 4, size = 3) +
  theme_bw() +
  geom_segment(data = mpna_1, aes(x = 0,y = 0,xend = CCA1,yend = CCA2),
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_text(data = al1,aes(x = CCA1 + 0.2, y = CCA2 -0.2, 
            label = c("4.a","4.b", "7.a,7.d,7.e"))) +
  geom_text(data = al2,aes(x = CCA1 + 0.2, y = CCA2 +0.2, 
                          label = c("latitude"))) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12))

ggsave(filename=c("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_writing/figures/cca_na.jpg"), 
       plot = cca_na, width = 8, height = 6.5, units = "in")




#####
## old cca all code below
#####
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
env_na <- dat_na[c("latitude", "ices_ecor", "ices_area")]
#env_na <- dat_na[c("latitude", "ices_ecor")]
#env_na <- dat_na[c("latitude", "ices_area")]

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
mps.nas.ca <- cca(dat_nas_sq ~ latitude + ices_area, data = env_na2)
plot(mps.nas.ca)


# run stats on the variables in the model
anova(mps.nas.ca, by = "terms")
#Model: cca(formula = dat_nas_sq ~ latitude + ices_area, data = env_na2)
#Df ChiSquare      F Pr(>F)    
#latitude   1    0.7450 6.1955  0.001 ***
#  ices_area  5    2.5047 4.1661  0.001 ***
#  Residual  84   10.1005  

#extracting the data as data frame; env data
veg_1na = as.data.frame(mps.nas.ca$CCA$biplot)
veg_1na["env"] = row.names(veg_1na)

#extracting the data; taxa
veg_2na = as.data.frame(mps.nas.ca$CCA$v)
veg_2na["genus"] = row.names(veg_2na)

#extracting the data; sitea
veg_3na = as.data.frame(mps.nas.ca$CCA$u)
veg_3na["seq_id"] = na_samples2
veg_3na <- merge(veg_3na, meta, by.x="seq_id")

# plot with ggplot2

ggplot() +
  geom_point(data = veg_3na, aes(x = CCA1, y = CCA2, color = country),
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
