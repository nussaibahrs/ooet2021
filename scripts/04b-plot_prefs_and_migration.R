# This script generates the plots of origination, extinction and 
# diversitfication preferences, and dispersal proportions across 
# climate zones over the last 66 million years, divided into climate 
# phases, as per the analyses in `04a-preferences_migration.R`.

# load libraries ---------------------------------------------------------
library(tidyverse)
library(divDyn)
library(magrittr)
library(reshape2)

source(file.path("scripts","utils", "functions.R"))
source(file.path("scripts","utils", "plot_parameters.R"))

# Load results ---------------------------------------------------------------
load(file.path("output", "prefs_sp.RData"))
prefs2$per <- factor(prefs2$per, levels=c("all", "W1","H", "W2", "C", "I"))

load(file.path("output", "mig_overall_stage.RData"))
mig.tot$per <- factor(mig.tot$per, levels=c("all", "W1","H", "W2", "C", "I"))
# Plot --------------------------------------------------------------------

# * a. Origination preferences --------------------------------------------
prefs_all <- prefs2 %>%  filter(group=="all" & method=="sqs" & env=="env1")

prefs_ori <- prefs_all %>% 
  dplyr::select(nont, not, ont, ot, per) %>% 
  group_by(per) %>% 
  summarise_if(is.numeric, mean, na.rm=TRUE)

prefs_ori <-  prefs_ori %>% 
  group_by(per) %>% 
  do(oddsratio(.$ot, .$not, .$ont, .$nont))

p1 <- ggplot(prefs_ori, aes(x=per, y=OR, ymin=LowerCI, ymax=UpperCI, col=per)) +
  geom_hline(yintercept = 0, col="darkgrey", linetype="dashed")+
  geom_point(size=3) +
  geom_errorbar(width=0.05, size=1) +
  scale_x_discrete(labels=c("all"="All", per_names))+
  scale_color_manual(values=c("black", clim_clr), guide=FALSE) +
  labs(y="Origination Preferences") +
  coord_cartesian(xlim=c(1,6), ylim=c(0.5, -0.8), clip="off")+
  theme_hce +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1),
        plot.margin = unit(c(1,1,1,3), "lines"),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

lims <- layer_scales(p1)$y$range$range

w <- 0.6
y1 <- c(lims[1] + w, lims[2] - w)
x1 <- 0.3
p1 <- p1 + annotate("segment", x=-x1, xend=-x1, y=-y1[1], yend=-lims[1],  
              arrow = arrow(length = unit(0.03, "npc")), 
              col="darkgrey") +
  annotate("segment", x=-x1, xend=-x1, y=-y1[2], yend=-lims[2],  
           arrow = arrow(length = unit(0.03, "npc")), 
           col="darkgrey") +
  annotate("text", x=-(x1+0.05), y=-y1[1], label="Tropics", angle=90, 
           vjust=-0.2, hjust=1,
           col="darkgrey") +
  annotate("text", x=-(x1+0.05), y=-y1[2], label="Extratropics", angle=90, 
           vjust=-0.2,hjust=0,
           col="darkgrey") 
  

# * b. Extinction preferences ---------------------------------------------
prefs_ext <- prefs_all %>% 
  dplyr::select(nent, net, ent, et, per) %>% 
  group_by(per) %>% 
  summarise_if(is.numeric, mean, na.rm=TRUE) 

prefs_ext[1,-1] <- t(apply(prefs_ext[-1,-1],2,sum)) #calculate for all

prefs_ext <- prefs_ext%>% 
  group_by(per) %>% 
  do(oddsratio(.$et, .$net, .$ent, .$nent))


p2 <- ggplot(prefs_ext, aes(x=per, y=OR, ymin=LowerCI, ymax=UpperCI, col=per)) +
  geom_hline(yintercept = 0, col="darkgrey", linetype="dashed")+
  geom_point(size=3) +
  geom_errorbar(width=0.05, size=1) +
  scale_x_discrete(labels=c("all"="All", per_names))+
  scale_color_manual(values=c("black", clim_clr), guide=FALSE) +
  scale_y_continuous(trans="reverse") +
  labs(y="Extinction Preferences") +
  coord_cartesian(xlim=c(1,6), clip="off")+
  theme_hce +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1),
        plot.margin = unit(c(1,1,1,3), "lines"),
        axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)))

lims <- layer_scales(p2)$y$range$range
w <- 0.6
y1 <- c(lims[1] + w, lims[2] - w)
x1 <- 0.2

p2 <- p2 + annotate("segment", x=-x1, xend=-x1, y=-y1[1], yend=-lims[1],  
                 arrow = arrow(length = unit(0.03, "npc")), 
                 col="darkgrey") +
  annotate("segment", x=-x1, xend=-x1, y=-y1[2], yend=-lims[2],  
           arrow = arrow(length = unit(0.03, "npc")), 
           col="darkgrey") +
  annotate("text", x=-(x1+0.05), y=-y1[1], label="Tropics", angle=90, 
           vjust=-0.2, hjust=1,
           col="darkgrey") +
  annotate("text", x=-(x1+0.05), y=-y1[2], label="Extratropics", angle=90, 
           vjust=-0.2,hjust=0,
           col="darkgrey") 



# * Migration patterns ----------------------------------------------------
mig.tot2 <- mig.tot %>% 
  filter(method=="sqs" & env=="env1" & group == "all") %>% 
  dplyr::select(-method) %>% 
  reshape2::melt(id.vars=c("env", "group", "per")) %>%  
  mutate(value=as.numeric(value)) %>% 
  group_by(group, env, per, variable) %>% 
  do(boot.summary2(.$value)) %>% 
  ungroup() 

mig.tot2 <- mig.tot2 %>% 
  pivot_wider(
    names_from = variable,
    names_glue = "{variable}_{.value}",
    values_from = c(mean, lwr, upr)
  )

p3 <- ggplot(mig.tot2, 
       aes(
         #Proportion imported into the extratropics, i.e. exported from the tropics
         x=nt_mean, 
         xmin=nt_lwr,
         xmax=nt_upr,
         #Proportion imported into the tropics, i.e. exported from the extratropics
           y=t_mean,
       ymin=t_lwr,
       ymax=t_upr,
       
       col=per
       
       )
       )+ 
  geom_point(size=3) +
  geom_errorbarh(height=0, size=1)+
  geom_errorbar(width=0, size=1) +
  geom_abline(slope=1, intercept = 0, linetype="dashed", col="darkgrey") +
  coord_cartesian(xlim=c(0, 0.4), ylim=c(0,0.4)) +
  scale_color_manual(values=c("black", clim_clr), 
                     labels=c("all"="All", per_names), 
                     guide = guide_legend(ncol=2)) +
  labs(x="Proportion exported from the tropics", y="Proportion exported from the extratropics",
       col="Climate Phase") +
  theme_hce


# * Diversification -------------------------------------------------------
prefs_div <- prefs_all %>% 
  dplyr::select(nont, not, ont, ot, per) %>% 
  group_by(per) %>% 
  summarise_if(is.numeric, mean, na.rm=TRUE) %>% 
  left_join(
    prefs_all %>% 
      dplyr::select(nent, net, ent, et, per) %>% 
      group_by(per) %>% 
      summarise_if(is.numeric, mean, na.rm=TRUE)
  ) 

prefs_div[1,-1] <- t(apply(prefs_div[-1,-1],2,sum)) #calculate for all

prefs_div <- prefs_div%>% 
  group_by(per) %>% 
  do(oddsratio(.$ot, .$et, .$ont, .$ent))

p4 <- ggplot(prefs_div, aes(x=per, y=OR, ymin=UpperCI, ymax=LowerCI, col=per)) +
  geom_hline(yintercept = 0, col="darkgrey", linetype="dashed")+
  geom_point(size=3) +
  geom_errorbar(width=0.05, size=1) +
  scale_x_discrete(labels=c("all"="All", per_names)) +
  scale_color_manual(values=c("black", clim_clr), guide=FALSE) +
  labs(y="Diversification Preferences") +
  coord_cartesian(xlim=c(1,6), ylim=c(1.2, -1.2), clip="off")+
  theme_hce +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1),
        plot.margin = unit(c(1,1,1,3), "lines"),
        axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)))

lims <- layer_scales(p4)$y$range$range
w <- 0.5
y1 <- c(lims[1] + w, lims[2] - w)
x1 <- 0.2

p4 <- p4 + annotate("segment", x=-x1, xend=-x1, y=-y1[1], yend=-lims[1],  
                    arrow = arrow(length = unit(0.03, "npc")), 
                    col="darkgrey") +
  annotate("segment", x=-x1, xend=-x1, y=-y1[2], yend=-lims[2],  
           arrow = arrow(length = unit(0.03, "npc")), 
           col="darkgrey") +
  annotate("text", x=-(x1+0.05), y=-y1[1], label="Tropics", angle=90, 
           vjust=-0.2, hjust=1,
           col="darkgrey") +
  annotate("text", x=-(x1+0.05), y=-y1[2], label="Extratropics", angle=90, 
           vjust=-0.2,hjust=0,
           col="darkgrey") 

# Save plot ---------------------------------------------------------------

svg(file.path("figs", "Fig_03_preferences_migration_overall.svg"), 
    w=12, h=12)
p1+p2 + p4 + p3 +
  plot_annotation(tag_prefix = "(", tag_levels = "a", tag_suffix = ")") +
  plot_layout(nrow=2, guides = "collect") &
  theme(legend.position='bottom')
dev.off()
