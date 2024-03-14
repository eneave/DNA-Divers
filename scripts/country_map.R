#######################################
## DNA Divers map to add to CCA plot ##
#######################################

# load data
loc <- read.csv("C:/Users/BESENEAV/OneDrive - Liverpool John Moores University/PhD/DNAdivers/R/map/divesite_metadata.csv")

# necessary packages
install.packages(c("cowplot", "googleway", "ggplot2", "ggrepel", 
                   "ggspatial", "libwgeom", "sf", "rnaturalearth", "rnaturalearthdata"))
# load packages
library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library(ggplot2)
library(ggspatial)
library(ggrepel)
library(maptools)

# download shape files for map
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

## Map of sites sampled

#windows()
sam <-
ggplot(data = world) +
  geom_sf(color = "black", fill = "lightgrey") +
  annotation_scale(location = "br", width_hint = 0.75) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.75, "in"),
                         style = north_arrow_fancy_orienteering) +
  geom_point(data = loc, aes(x=londd, y=latdd, color = map_name), size = 5) +
  geom_label_repel(data = loc, 
                   aes(x=londd, y=latdd, label = sampletot),
                   fontface = 'bold', color = 'black',
                   box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.5, "lines"),
                   segment.color = 'black'
  ) +
  coord_sf(xlim = c(-10, 12), ylim = c(46, 65), expand = FALSE) +
  labs(colour= "Dive Locations") + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        #axis.line.x = element_line(colour=c("#CC79A7")),
        #axis.line.y = element_line(colour=c("#CC79A7")),
        #legend.direction = "vertical", legend.box = "horizontal",
        legend.title = element_text(size = 12, colour = "black", face = "bold"), 
        legend.text = element_text(size = 11, face ="bold", colour ="black"))  

## make world map

data(wrld_simpl)
myCountries = wrld_simpl@data$NAME %in% c("United Kingdom", "Norway", "United States", "South Africa")

wsam <- plot(wrld_simpl, col = c(gray(.80), "red")[myCountries+1])

## combine maps

library(cowplot)

#windows()
#plot_grid(sam, wsam, align = c("v"), axis = c("l"))

ggsave(filename=c("C:/Users/BESENEAV/OneDrive - Liverpool John Moores University/PhD/DNAdivers/R/map/samples_processed.jpeg"), 
       plot = sam, width = 6, height = 6, dpi = 1000, units = "in")

ggsave(filename=c("C:/Users/BESENEAV/OneDrive - Liverpool John Moores University/PhD/DNAdivers/R/map/countries_sampled.jpeg"), 
       plot = wsam, width = 7, height = 3, dpi = 1000, units = "in")

