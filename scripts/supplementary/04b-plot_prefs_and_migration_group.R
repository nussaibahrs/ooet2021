# Plots origination & extinction preferences for each plankton group

# load libraries ---------------------------------------------------------
library(tidyverse)
library(divDyn)
library(magrittr)
library(reshape2)

source(file.path("code", "functions.R"))
source(file.path("code", "plot_parameters.R"))

# Load results ---------------------------------------------------------------
load(file.path("output", "prefs_sp.RData"))
prefs2$per <- factor(prefs2$per, levels=c("all", names(phases)))
prefs2$group[prefs2$group == "FALSE"] <- "F"

load(file.path("output", "mig_overall_stage.RData"))
mig.tot$per <- factor(mig.tot$per, levels=c("all", names(phases)))
mig.tot$group[mig.tot$group == "FALSE"] <- "F"
# Plot --------------------------------------------------------------------

# * a. Origination preferences --------------------------------------------

prefs_all <- prefs2 %>%  filter(group!="all" & method=="sqs" & env=="env1")


prefs_ori <- prefs_all %>% 
  dplyr::select(nont, not, ont, ot, per, group) %>% 
  group_by(per, group) %>% 
  summarise_if(is.numeric, mean, na.rm=TRUE) %>% 
  group_by(per, group) %>% 
  do(oddsratio(.$ot, .$not, .$ont, .$nont))

p1 <- ggplot(prefs_ori, aes(x=per, y=OR, ymin=LowerCI, ymax=UpperCI, col=per)) +
  geom_hline(yintercept = 0, col="darkgrey", linetype="dashed")+
  geom_point() +
  geom_errorbar(width=0.05) +
  scale_x_discrete(labels=c("all"="Overall", per_names)) +
  scale_color_manual(values=c("black", clim_clr), guide=FALSE) +
  scale_y_continuous(trans="reverse") +
  labs(y="Origination Preferences") +
  coord_cartesian(xlim=c(1,length(phases)+1), clip="off")+
  facet_grid(~group, labeller=labeller(group=group_names))+
  theme_hce +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1))
# * b. Extinction preferences ---------------------------------------------
prefs_ext <- prefs_all %>% 
  dplyr::select(nent, net, ent, et, per, group) %>% 
  group_by(per, group) %>% 
  summarise_if(is.numeric, mean, na.rm=TRUE) %>% 
  group_by(per, group) %>% 
  do(oddsratio(.$et, .$net, .$ent, .$nent))

p2 <- ggplot(prefs_ext, aes(x=per, y=OR, ymin=LowerCI, ymax=UpperCI, col=per)) +
  geom_hline(yintercept = 0, col="darkgrey", linetype="dashed")+
  geom_point() +
  geom_errorbar(width=0.05) +
  scale_x_discrete(labels=c("all"="Overall", per_names)) +
  scale_color_manual(values=c("black", clim_clr), guide=FALSE) +
  scale_y_continuous(trans="reverse") +
  labs(y="Extinction Preferences") +
  coord_cartesian(xlim=c(1,length(phases)+1), clip="off")+
  facet_grid(~group, labeller=labeller(group=group_names))+
  theme_hce +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1))

# * Migration patterns ----------------------------------------------------
mig.tot2 <- mig.tot %>% 
  filter(method=="sqs" & env=="env1" & group != "all") %>% 
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
  geom_point() +
  geom_errorbarh(height=0)+
  geom_errorbar(width=0) +
  geom_abline(slope=1, intercept = 0, linetype="dashed", col="darkgrey") +
  coord_cartesian(xlim=c(0, 0.5), ylim=c(0,0.5)) +
  scale_color_manual(values=c("black", clim_clr), 
                     labels=c("all"="Overall", per_names), 
                     guide = guide_legend(nrow=1)) +
  facet_grid(~group, labeller=labeller(group=group_names))+
  labs(x="Proportion exported from the tropics", y="Proportion exported \nfrom the extratropics",
       col="Climate Phase") +
  theme_hce

# Save plot ---------------------------------------------------------------

svg(file.path("figs", "Supplement", "Fig_S_preferences_migration_group.svg"), 
    w=10, h=12)
p1+p2 + p3 +
  plot_annotation(tag_prefix = "(", tag_levels = "a", tag_suffix = ")") +
  plot_layout(nrow=3, heights=c(0.3,0.3,0.4), guides = "collect") &
  theme(legend.position='bottom')
dev.off()

