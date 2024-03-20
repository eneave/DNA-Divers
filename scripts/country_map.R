#######################################
## DNA Divers map to add to CCA plot ##
#######################################

# load data
loc <- read.csv("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_dnadivers/DNA-Divers/data/metadata/supp_table_1.csv")

# necessary packages
install.packages(c("cowplot", "googleway", "ggplot2", "ggrepel", 
                   "ggspatial", "libwgeom", "sf", "rnaturalearth", "rnaturalearthdata"))
# load packages
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(ggspatial)
library(ggrepel)
library(maptools)
library(tidyverse)

# download shape files for map
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

## Map of North Atlantic sites sampled

loc_na <- subset(loc, loc$ocean_basin=="North Atlantic" & loc$location!="Ireland")
loc_na2 <- subset(loc_na, loc_na$country!="Cape Verde")
loc_na3 <- loc_na2 %>%
              group_by(location, site.name, country, latitude, longitude, type, via) %>%
              summarise(sequence.run = mean(sequence.run))
loc_na4 <- loc_na3 %>%
  group_by(location, latitude, longitude) %>%
  summarise(sequence.run = mean(sequence.run))

#windows()
na <-
ggplot(data = world) +
  geom_sf(color = "black", fill = "lightgrey") +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.1, "in"), pad_y = unit(0.1, "in"),
                         style = north_arrow_fancy_orienteering) +
  geom_point(data = loc_na4, aes(x=longitude, y=latitude, color = location), size = 5, alpha=0.6) +
  scale_colour_manual(values = c("#332288", "#117733", "darkgrey", "#44AA99", "#88CCEE", 
                                 "#DDCC77", "#CC6677", "#AA4499", "#882255")) +
  coord_sf(xlim = c(-10, 12), ylim = c(46, 65), expand = FALSE) +
  #labs(colour= "Dive Locations") + 
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "azure"),
        axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, colour ="black"), 
        legend.position = "right", 
        legend.title = element_text(size = 14, colour = "black", face = "bold"))  

library(cowplot)
figure3bm <- plot_grid(na, labels = "B", ncol = 1)

ggsave(filename=c("C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/chapter3_writing/figures/na_map_B.jpg"), 
       plot = figure3bm, width = 6, height = 8, units = "in")


## make world map

data(wrld_simpl)
myCountries = wrld_simpl@data$NAME %in% c("United Kingdom", "Norway", "United States", "South Africa", "Cape Verde", "Jordan")

windows()
plot(wrld_simpl, col = c(gray(.80), "red")[myCountries+1])

## combine maps

library(cowplot)

#windows()
#plot_grid(sam, wsam, align = c("v"), axis = c("l"))

ggsave(filename=c("C:/Users/BESENEAV/OneDrive - Liverpool John Moores University/PhD/DNAdivers/R/map/samples_processed.jpeg"), 
       plot = sam, width = 6, height = 6, dpi = 1000, units = "in")

ggsave(filename=c("C:/Users/BESENEAV/OneDrive - Liverpool John Moores University/PhD/DNAdivers/R/map/countries_sampled.jpeg"), 
       plot = wsam, width = 7, height = 3, dpi = 1000, units = "in")

