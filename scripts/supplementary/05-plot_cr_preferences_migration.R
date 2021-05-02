# Plots diversfication and dispersal from subsampled data using 
# classical rarefaction.

# load libraries ---------------------------------------------------------
library(tidyverse)
library(divDyn)
library(magrittr)
library(reshape2)

source(file.path("scripts", "utils","functions.R"))
source(file.path("scripts","utils", "plot_parameters.R"))

# Load results ---------------------------------------------------------------
load(file.path("output", "prefs_sp.RData"))
prefs2$per <- factor(prefs2$per, levels=c("all", "G", "T1", "T2", "I"))

load(file.path("output", "mig_overall_stage.RData"))
# Plot --------------------------------------------------------------------

# * a. Origination preferences --------------------------------------------
prefs_all <- prefs2 %>%  filter(group=="all" & method=="cr" & env=="env1")

hist(prefs_all$orig[prefs_all$per=="all"]) #look at normality


prefs_all <- prefs_all %>% group_by(per) %>% 
  do(boot.summary(.$orig))  

p1 <- ggplot(prefs_all, aes(x=per, y=mean, ymin=lwr, ymax=upr, col=per)) +
  geom_hline(yintercept = 0, col="darkgrey", linetype="dashed")+
  geom_point() +
  geom_errorbar(width=0.05) +
  scale_x_discrete(labels=c("O"="Overall", "G"="Greenhouse", "T1"="Climate transition: \ncooling", 
                            "T2"="Climate transition: \nstable", "I"="Icehouse")) +
  scale_color_manual(values=c("black", clim_clr), guide=FALSE) +
  scale_y_continuous(trans="reverse")+
  labs(y="Origination Preferences") +
  theme_hce +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1))


# * Migration patterns ----------------------------------------------------
mig.tot2 <- mig.tot %>% 
  filter(method=="cr" & env=="env1") %>% 
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

p2 <- ggplot(mig.tot2 %>% filter(env=="env1" & group=="all"), 
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
  coord_cartesian(xlim=c(0, 0.3), ylim=c(0,0.3)) +
  scale_color_manual(values=c("black", clim_clr), 
                     labels=c("all"="Overall", "G"="Greenhouse", "T1"="Climate transition: \ncooling", 
                              "T2"="Climate transition: \nstable", "I"="Icehouse"), 
                     guide = guide_legend(nrow=3)) +
  labs(x="Proportion exported from the tropics", y="Proportion exported from the extratropics",
       col="Climate Period") +
  theme_hce


svg(file.path("figs", "Supplement", "Fig_S_preferences_migration_overall_cr.svg"), 
    w=6, h=12)
p1+p2 + 
  plot_annotation(tag_prefix = "(", tag_levels = "a", tag_suffix = ")") +
  plot_layout(nrow=2, guides = "collect") &
  theme(legend.position='bottom')
dev.off()
