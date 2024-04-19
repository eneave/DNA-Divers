############################################
## Understand factors effecting diversity ##
############################################
library(tidyverse)
library(vegan)
#####
# load data
#####
# samples from nature
nat_tab <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/final_motu_mp_nature_cca.csv")
# make minor edits to csv
nat_tab  <- subset(nat_tab, select = -X)
# remove 6 reads salmon from California sample
nat_tab[nat_tab$final_name=="Salmo salar", "sample.7G_CA_diveA_a"] <- 0


# elas02 shark tank sequence run 2
p2e_tab <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/decontam/p2e_aq_decontam_001.csv")
# Separate out teleosts and elasmobranchs
p2e_tab <- subset(p2e_tab, final_class!="Mammalia") 
p2e_tab <- subset(p2e_tab, select = -c(X.1,X,sample.11Be_EBMay_13extblank,
                                       sample.11Ce_EBJun_12extblank,sample.11De_EBJun_19extblank,
                                       sample.11EeBLUE_FBEt_fieldb,sample.11Fe_negativePCRcontrol,
                                       sample.8CeBLUE_eDNA_FBblank,sample.8DeBLUE_eDNAA_bottle1,
                                       sample.8EeBLUE_eDNAB_bottle2,sample.8FeBLUE_eDNAC_bottle3,
                                       sample.8GeBLUE_eDNAD_bottle4)) #get rid of meaningless row name columns

# tele02 coral cave sequencr run 2
p2t_tab <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/decontam/p2t_aq_decontam_001.csv")
# Separate out teleosts
p2t_tab <- subset(p2t_tab, final_class=="Actinopterygii") 
p2t_tab <- subset(p2t_tab, select = -c(X.1,X,sample.6GtBLUE_FBEtblank,
                                       sample.8Dt_EBMay_13extblank, sample.8Et_EBJun_12extblank,
                                       sample.8Ft_EBJun_13extblank)) #get rid of meaningless row name columns
# remove contamination
p2t_tab <- subset(p2t_tab, final_name!="Pangasianodon hypophthalmus") #positive control
p2t_tab <- subset(p2t_tab, final_name!="Leporidae") #hare, mislabeled at Actinopterygii

# sample metadata
meta <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/metadata/supp_table_1.csv")


#####
# Combine all metaprobe data together
#####

ntl <- nat_tab %>% pivot_longer(cols = sample.1A:sample.6B_B_CLAS_NOR_Jul22,
                               names_to = "seq_id",
                               #names_prefix = "sample.",
                               values_to = "reads")

pel <- p2e_tab %>% pivot_longer(cols = sample.10Ae60BLUE_MPEtB:sample.9He60BLUE_MPEtA_,
                                names_to = "seq_id",
                                #names_prefix = "sample.",
                                values_to = "reads")

ptl <- p2t_tab %>% pivot_longer(cols = sample.5Et10BLUE_MPEtA_:sample.6Ft240BLUE_MPEtB,
                                names_to = "seq_id",
                                #names_prefix = "sample.",
                                values_to = "reads")

# remove unnecessary columns
pel <- subset(pel, select = -c(id,final_order,pid,method_assign,sequence))
ptl <- subset(ptl, select = -c(id,final_order,pid,method_assign,sequence))


# combine long datafame
mp <- rbind(ntl,pel,ptl)
rm(ntl,pel,ptl)

# wide motu table
mp1 <- mp %>% group_by(seq_id, final_name) %>%
  summarise(reads = sum(reads)) %>%
  tidyr::pivot_wider(names_from = final_name, values_from = reads, 
                     values_fill = 0)


# species table 
mp0 <- mp1[c(2:276)]

# species data and metadata combined
mp2 <- merge(mp1, meta, by="seq_id")
#mp2$richness <- meta0$richness

write.csv(mp2, "C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/final_motu_mp_all.csv")


# new metadata 
meta0 <- mp2[c(1,277:304)]

#####
# Alpha diversity indices
#####

meta0$richness <- specnumber(mp0)
meta0$shan <- diversity(mp0, "shannon")
meta0$simp <- diversity(mp0, "simpson")
meta0$unsimp <- simpson.unb(mp0) #unbiased simpson

# calculate number of reads per sample for target taxa (vertebrates; non-human contaminants removed)
meta0$targetreads <- rowSums(mp0)

#####
# Sub-model
#####

## Only data for which exact weight of input gauze was measured
g1 <- subset(meta0, weight_g>0)

#####
model1a <- glm(richness ~ country/site.name/mp_id + via + preservation*extraction + weight_g + gauze,  data = g1)
summary.glm(model1a)
anova(model1a)

# check assumptions of model
# 1. normally distributed residuals
sresid <- (model1a$residuals - mean(model1a$residuals))/sd(model1a$residuals)

hist(sresid, freq=F)
lines(density(sresid, adjust=1))

qqnorm(sresid, cex=1.8, pch=20) #yes

# 2. the variances of residuals are homogenous

plot(sresid ~ model1a$fitted.values, pch=20, cex=2, cex.lab=1.5)

# 3. No collinearity of independant variables

#plot(g1[c("country","site.name","mp_id","via","preservation","extraction","weight_g","gauze")], panel = panel.smooth)

library(car)
durbinWatsonTest(model1a) #glm
#####

library(lme4)
#GLMM
baseline = lmer(richness ~ 1 + (1|mp_id), data=g1)
summary(baseline)

baseline2 = lmer(richness ~ 1 + (1|country/site.name/mp_id), data=g1)
summary(baseline2)

# fixed gauze,via
# random country/site.name/mp_id, input weight, extraction, preservation

fixed_gv = lmer(richness ~ 1 + gauze + via + (1 |country/site.name/mp_id), data=g1)
summary(fixed_gv)


# fixed extraction and preservation, input weight, extraction, preservation
# random input weight, extraction, preservation

fixed_gv = lmer(richness ~ 1 + gauze + via + weight_g + extraction + preservation + (1 + weight_g + extraction + preservation |country/site.name/mp_id), data=g1)
summary(fixed_gv)

#####

## Only data for which exact weight of input gauze was measured; from one sample
g2 <- subset(meta0, location=="Dorset" & extraction=="MUT")

## add in data for exact amount of lysis buffer pushed through the filter
g2$lysis <- c(750, 600, 500, 500, 700, 300,
       2000, 3000, 1500, 2000, 1750, 1750,
       3500, 2500, 3200, 5000, 3000, 3000)


baseline = lmer(richness ~ 1 + (1|mp_id), data=g2)
summary(baseline)


# fixed time
# random input weight

#fixed_tw = lmer(richness ~ weight_g + time + (weight_g |mp_id), data=g2)
#summary(fixed_tw)

# FINAL MODEL
fixed_tw2 = lmer(richness ~ weight_g + lysis + time + (weight_g + lysis |mp_id), data=g2)
summary(fixed_tw2)

#fixed_tw3 = lmer(richness ~ weight_g + lysis + targetreads + time + (weight_g + lysis + targetreads |mp_id), data=g2)
#summary(fixed_tw3)

weight <-
ggplot(g2,aes(x=mp_id,y=richness,col=mp_id)) + 
  geom_jitter() +
  geom_boxplot(alpha=0.2) + 
  facet_wrap(~ weight_class)
ggsave(filename=c("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_writing/figures/weight_exp.jpg"), 
       plot = weight, width = 4, height = 4, units = "in")

#####
# NOT SURE THIS IS DONE CORRECTLY
## Data to compare the effects of extraction and preservation
g3 <- subset(meta0, mp_id=="PLYM_1"| mp_id=="PLYM_2"| mp_id=="PLYM_3"| 
               mp_id=="DOR_1"| mp_id=="DOR_2"| mp_id=="LIV_1"| mp_id=="LIV_2"| 
               mp_id=="LIV_3"| mp_id=="NEW_1"| mp_id=="NEW_2"| mp_id=="NEW_3"| 
               mp_id=="NEW_4"| mp_id=="NEW_5"| mp_id=="NEW_6" | mp_id=="ORK_1"| 
               mp_id=="ORK_2"| mp_id=="ORK_3"| 
               mp_id=="ORK_4"| mp_id=="ORK_5"| mp_id=="ORK_6")
g3 <- subset(g3, weight_class=="standard")


baseline = lmer(richness ~ 1 + (1|country/site.name/mp_id), data=g3)
summary(baseline)


# fixed time, preservation
# random sequence run, extraction

fixed_pe = lmer(richness ~ 1 + time + preservation + extraction + sequence.run + (1 + extraction + preservation + sequence.run |country/site.name/mp_id), data=g3)
summary(fixed_pe)

fixed_pe2 = lmer(richness ~ time + preservation + extraction + (extraction + preservation |country/site.name/sequence.run/mp_id), data=g3)
summary(fixed_pe2)


# Would like to generate p-values for model output
library(lmerTest)


# Preservation-extraction model
coef(summary(as(baseline,"merModLmerTest")))
coef(summary(as(fixed_pe,"merModLmerTest")))

ggplot(g3,aes(x=mp_id,y=richness,col=site.name)) + 
  geom_jitter() + 
  geom_boxplot(alpha=0.2) + 
  facet_wrap(extraction ~ preservation)

#####
# test how extractions effect richness, only samples where extraction varies within

g4 <- subset(meta0, mp_id=="PLYM_2"| mp_id=="PLYM_3"| 
               mp_id=="DOR_1"| mp_id=="DOR_2"| mp_id=="ORK_1"| 
               mp_id=="ORK_2"| mp_id=="ORK_3")
g4 <- subset(g4, weight_class=="standard" & preservation=="Et" & sequence.run=="1")



baseline = lmer(richness ~ 1 + (1|country/site.name/mp_id), data=g4)
summary(baseline)


# fixed time
# random sequence run, extraction


# interesting; when you add the ices_area the intercept is no longer significant
#fixed_e = lmer(richness ~ extraction + (extraction |country/site.name/mp_id), data=g4)
#summary(fixed_e)

# FINAL MODEL
fixed_e1 = lmer(richness ~ extraction + site.name + (extraction |mp_id), data=g4)
summary(fixed_e1)

#fixed_e2 = lmer(shan ~ 1 + extraction + (1 + extraction |country/site.name/mp_id), data=g4)
#summary(fixed_e2)

#fixed_e2a = lmer(shan ~ extraction + site.name + (extraction |mp_id), data=g4)
#summary(fixed_e2a)

#fixed_e3 = lmer(targetreads ~ 1 + extraction + (1 + extraction |country/site.name/mp_id), data=g4)
#summary(fixed_e3)

#coef(summary(as(fixed_e,"merModLmerTest")))

extraction <-
ggplot(g4,aes(x=mp_id,y=richness,col=site.name)) + 
  geom_jitter() + 
  geom_boxplot(alpha=0.2) + 
  facet_wrap(~ extraction)
ggsave(filename=c("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_writing/figures/extraction_exp.jpg"), 
       plot = extraction, width = 4, height = 4, units = "in")


#ggplot(g4,aes(x=mp_id,y=shan,col=site.name)) + 
#  geom_jitter() + 
#  geom_boxplot(alpha=0.2) + 
#  facet_wrap(~ extraction)

#ggplot(g4,aes(x=mp_id,y=targetreads,col=site.name)) + 
#  geom_jitter() + 
#  geom_boxplot(alpha=0.2) + 
#  facet_wrap(~ extraction)

#####
# test how preservation effects richness, only samples where preservation varies within

g5 <- subset(meta0, mp_id=="LIV_1"| mp_id=="LIV_2"| 
               mp_id=="LIV_3"| mp_id=="ORK_1"| 
               mp_id=="ORK_2"| mp_id=="ORK_3" | mp_id=="ORK_4"| 
               mp_id=="ORK_5"| mp_id=="ORK_6" | mp_id=="NEW_1"| 
               mp_id=="NEW_2"| mp_id=="NEW_3" | mp_id=="NEW_4"| 
               mp_id=="NEW_5"| mp_id=="NEW_6")
g5 <- subset(g5, weight_class=="standard" & extraction=="QBT")


#fixed_p = lmer(richness ~ preservation + (preservation|country/site.name/mp_id), data=g5)
#summary(fixed_p)

#FINAL MODEL
fixed_p1 = lmer(richness ~ preservation + site.name + sequence.run + (preservation|mp_id), data=g5)
summary(fixed_p1)


#ggplot(g5,aes(x=mp_id,y=richness,col=site.name)) + 
#  geom_jitter() + 
# geom_boxplot(alpha=0.2) + 
#  facet_wrap(~preservation)

