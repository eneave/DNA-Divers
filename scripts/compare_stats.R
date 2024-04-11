####################################################
## Statistics for bottle and metaprobe comparison ##
####################################################

#library(betapart)
library(tidyverse)
library(vegan)

###############
## PERMANOVA ##
###############

## Calculate beta-dispersion - does variation in data vary between different sampling groups

#location
cem_bd1 <- betadisper(jac_dat, cem_jac$site.name)
anova(cem_bd1)
permutest(cem_bd1)

#type
cem_bd2 <- betadisper(jac_dat, cem_jac$type)
anova(cem_bd2)
permutest(cem_bd2)

#primer
#cem_bd3 <- betadisper(jac_dat, cem_jac$primer)
#anova(cem_bd3)

#location*type
cem_jac2 <- cem_jac %>%
                mutate(site_type=paste0(site.name,type))
cem_bd4 <- betadisper(jac_dat, cem_jac2$site_type)
anova(cem_bd4)
permutest(cem_bd4, pairwise=TRUE)

## Permanova

str(jac_dat) #check that this is a distance matrix

# compare sample type, nested by site.name
test1 <- adonis2(jac_dat ~ type, data = cem_jac, permutations = 999, strata = cem_jac$site.name)

# compare sample type and primer, nested by site.name
#test2 <- adonis2(jac_dat ~ type*primer, data = cem_jac, permutations = 999, strata = cem_jac$site.name)

# compare sample site.name, nested by type
test3 <- adonis2(jac_dat ~ site.name, data = cem_jac, permutations = 999, strata = cem_jac$type)






