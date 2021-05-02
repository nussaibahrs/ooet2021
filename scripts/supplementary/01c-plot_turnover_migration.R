# Plots turnover rates and migration for each plankton group

# load libraries ---------------------------------------------------------
library(tidyverse)
library(divDyn)
library(magrittr)
library(reshape2)
library(forecast)
library(nlme)

source(file.path("scripts", "utils","functions.R"))
source(file.path("scripts","utils", "plot_parameters.R"))

# Load data ---------------------------------------------------------------
load(file.path("data", "wrangled_data.RData")) # sp at species level
load(file.path("output", "turnover_rates.RData"))
load(file.path("output", "mig_stage.RData"))


# Net origination ---------------------------------------------------------

orig.rate <- orig.rate %>% 
  filter(mid <67 & method =="sqs")  %>% 
  left_join(stages[,c("mid", "per")]) %>% 
  filter(!is.na(per))

orig.rate <- orig.rate %>% 
  bind_rows(orig.rate %>%  
              mutate(per = "all")) 

orig.rate <- orig.rate %>% 
  dplyr::select(-method, -mid) 

orig.rate <- orig.rate %>% 
  reshape2::melt(id.vars=c("per", "env", "group")) %>% 
  group_by(per, env, group) %>% 
  do(boot.summary(.$value)) %>% 
  filter(group !="all")

# Net Extinction ----------------------------------------------------------
ext.rate <- ext.rate %>% 
  filter(mid <67 & method =="sqs")  %>% 
  left_join(stages[,c("mid", "per")])

ext.rate <- ext.rate %>% 
  bind_rows(ext.rate %>%  
              mutate(per = "all")) 


ext.rate <- ext.rate %>% 
  dplyr::select(-method, -mid) 


ext.rate <- ext.rate %>% 
  reshape2::melt(id.vars=c("per", "env", "group")) %>% 
  group_by(per, env, group) %>% 
  do(boot.summary(.$value)) %>% 
  filter(group !="all")

orig.rate$per <- factor(orig.rate$per, levels=c("all", names(phases)))
ext.rate$per <- factor(ext.rate$per, levels=c("all", names(phases)))


# Plot --------------------------------------------------------------------
p1 <- ggplot(orig.rate,
             aes(x=per, y=mean, ymin=lwr, ymax=upr, 
                 col=per)) +
  geom_point() +
  geom_errorbar(width=0.05) +
  facet_grid(rows=vars(env), cols=vars(group), 
             labeller=labeller(group=group_names, env=env_names), 
             scales = "free_y") +
  scale_x_discrete(labels=c("all"="Overall", per_names)) +
  scale_color_manual(values=c("black", clim_clr), guide=FALSE) +
  labs(y="Origination rate") +
  coord_cartesian(xlim=c(1,6), clip="off")+
  theme_hce +
  theme(panel.grid.major.y = element_line(colour = "darkgrey", size = 0.2, linetype="dashed"),
         axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1),
        plot.margin = unit(c(1,1,1,3), "lines"),
        axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)))

p2 <- ggplot(ext.rate,
             aes(x=per, y=mean, ymin=lwr, ymax=upr, 
                 col=per)) +
  geom_point() +
  geom_errorbar(width=0.05) +
  facet_grid(rows=vars(env), cols=vars(group), labeller=labeller(group=group_names, env=env_names), scales="free_y") +
  scale_x_discrete(labels=c("all"="Overall", per_names)) +
  scale_color_manual(values=c("black", clim_clr), guide=FALSE) +
  labs(y="Extinction rate") +
  coord_cartesian(xlim=c(1,6), clip="off")+
  theme_hce +
  theme(panel.grid.major.y = element_line(colour = "darkgrey", size = 0.2, linetype="dashed"),
    axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1),
        plot.margin = unit(c(1,1,1,3), "lines"),
        axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)))

# save
svg(file.path("figs", "Supplement", "Fig_S_turnover_rates_sqs.svg"), 
    w=10, h=8)
p1+p2 +
  plot_annotation(tag_prefix = "(", tag_levels = "a", tag_suffix = ")") +
  plot_layout(ncol=1, guides = "collect") &
  theme(legend.position='top')
dev.off()



