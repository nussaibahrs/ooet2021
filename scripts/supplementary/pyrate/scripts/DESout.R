# Plots the results from Pyrate DES

# Load libraries and functions --------------------------------------------

library(tidyverse)
path <- file.path("scripts", "supplementary", "pyrate")
source(file.path(path, "scripts", "functions_DES.R"))
source(file.path("scripts", "utils","functions.R"))
source(file.path("scripts","utils", "plot_parameters.R"))

# Function to get marginal rates ------------------------------------------


plot_file <- file.path(path, "pyrate_output", 
                       "DESin_all_rep1_0_q_0.0_0.0117_2.588_5.333_11.62_15.97_23.03_28.1_33.9_38.0_41.3_47.8_56.0_61.6_66.0_DivdD_DivdE_ML_marginal_rates.log")

file.exists(plot_file)

library(divDyn)
data(stages)
stages$phases <- categorize(stages$stage, phases)
stages <- stages[!is.na(stages$phases),]

mr <- DESout(plot_file, raw=TRUE)

# Migration ---------------------------------------------------------------

nt_to_t <- cbind(reshape2::melt(mr$d12, variable.name="bottom"), env="nt")
t_to_nt <- cbind(reshape2::melt(mr$d21, variable.name="bottom"), env="t")


disp <- rbind(nt_to_t, t_to_nt) 
disp$bottom <- as.numeric(as.character(disp$bottom))

disp <- disp %>% 
  left_join(stages[,c("bottom","phases")])

disp[is.na(disp$phases),]$phases <- "I"
disp <- disp %>% bind_rows(disp %>%  mutate(phases="all"))

mig.tot2 <- disp %>%  group_by(phases, env) %>% 
  summarise(mean=mean(value, na.rm=TRUE), lwr=calcHPD(value)[1], upr=calcHPD(value)[2]) %>% 
  ungroup() 

mig.tot2 <- mig.tot2 %>% 
  pivot_wider(
    names_from = env,
    names_glue = "{env}_{.value}",
    values_from = c(mean, lwr, upr)
  )

mig.tot2$phases <- factor(mig.tot2$phases, levels=c("all", names(phases)))

p1 <- ggplot(mig.tot2 %>%  arrange(phases), 
             aes(
               #Proportion imported into the extratropics, i.e. exported from the tropics
               x=t_mean, 
               xmin=t_lwr,
               xmax=t_upr,
               #Proportion imported into the tropics, i.e. exported from the extratropics
               y=nt_mean,
               ymin=nt_lwr,
               ymax=nt_upr,
               
               col=phases
               
             )
)+ 
  geom_point(size=3) +
  geom_errorbarh(height=0, size=1)+
  geom_errorbar(width=0, size=1) +
  geom_abline(slope=1, intercept = 0, linetype="dashed", col="darkgrey") +
  coord_cartesian(xlim=c(0, 0.3), ylim=c(0,0.3)) +
  scale_color_manual(values=c("all" = "black", clim_clr), 
                     labels=c("all"="All", per_names), 
                     guide = guide_legend(nrow=2)) +
  labs(x="Dispersal rate from the tropics", 
       y="Dispersal rate from the extratropics",
       col="Climate Phase") +
  theme_hce

ggsave(file.path("figs", "Supplement", "Fig_S_pyrate_migration.svg"), p1, w=6, h=6)
