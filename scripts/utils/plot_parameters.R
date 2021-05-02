library(ggplot2)
library(RColorBrewer)
library(ggthemes)
library(patchwork)
library(grid)
library(ggpubr)
library(viridis)

# For plotting purposes ---------------------------------------------------

phases <- list(W1=c("Danian", "Selandian-Thanetian"),
               H = c("Ypresian"),
               W2 = c("Lutetian", "Bartonian", "Priabonian"),
               C = c("Rupelian", "Chattian", "Lower Miocene", 
                     "Middle Miocene", "Upper Miocene", "Pliocene"),
               I = c("Pleistocene", "Holocene")
)

##colour palette ####
pal <- c("t"=rgb(241,88,0, max=255), "nt"=rgb(0,120,190, max=255))
clim_clr <- c("#d73232", "#a93d32","#fc7e3c", "#8aaccf", "#321e72")
names(clim_clr) <- names(phases)
ldg_clr <- c("black", "#2A58A1", "#2E8B55")
gp <- gpar(fontsize = 16, fontface = "bold") #labels

#for labelling purposes ####
group_names <- c("all"="All", 
                 "D" ="Diatoms", 
                 "F" = "Foraminifera", 
                 "N" = "Nannofossils", 
                 "R" = "Radiolaria")

scenario_names <- c("A"="Scenario A", "B"="Scenario B")

per_names <- c("W1" = "Warmhouse I", "H"="Hothouse", "W2" = "Warmhouse II", "C" = "Coolhouse", "I"="Icehouse")

env_names <- c("t"="Tropical", "nt"="Extratropical")
#names(clr) <- names(group_names)

#theme
theme_hce <- ggthemes::theme_hc(base_size = 14) + theme(panel.grid.major.y = element_blank(), 
                                                        axis.title = element_text(face="bold"), 
                                                        strip.background = element_rect(fill=alpha("lightgrey", 0.2)), 
                                                        legend.title = element_text(face="bold"),
                                                        legend.position = "top")

#### geogscale
data(stages, package="divDyn")
syscol <- stages %>% group_by(system, sys, systemCol) %>%
  summarise(bottom=max(bottom), top=min(top), mid=(bottom+top)/2) %>%
  arrange(bottom) 

syscol[1,]$system <- ""

ints <- stages
ints$per <- stages$per <- categorize(stages$stage, phases) 
ints <- ints[!is.na(ints$per),]

ints <- ints %>% 
  mutate(per = factor(per, levels=names(phases))) %>% 
  group_by(per) %>% 
  summarise(bottom=max(bottom), top=min(top), mid=(bottom+top)/2)

add_scale <- function(p, y=-15) {
  p +
    annotate("rect", xmin=ints$bottom, xmax=ints$top, ymin=0, ymax=y, 
             fill=clim_clr[c(3,2,3,4,5)], col="white") +
    annotate("text", x=ints$mid, y=y/2, 
             label=names(phases),  size=3,
             hjust=0.5, fontface="bold", col="white") +
    annotate("rect", xmin=syscol$bottom, xmax=syscol$top, 
             ymin=y, ymax=y*2, col="white", fill="grey95") +
  annotate("text", x=syscol$mid, y=y+y/2, label=syscol$system, 
           hjust=0.5, fontface="bold", col="grey20") 
  }


  
