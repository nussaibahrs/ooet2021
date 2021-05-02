# Computes the change in magnitude between rates

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
load(file.path("data", "wrangled_data.RData")) # sp at genus level
load(file.path("output", "turnover_rates.RData"))
load(file.path("output", "mig_stage.RData"))

data(stages)
stages$per <- categorize(stages$stage, phases) 
stages <- stages[!is.na(stages$per),]

# Calculate --------------------------------------------------------------------

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

# significant difference
orig.sig <- orig.rate %>% 
  reshape2::melt(id.vars=c("per", "env", "group")) %>% 
  group_by(per, env, group) 

orig.sig <- orig.sig %>% 
  pivot_wider(id_cols = c(per, group), names_from=env, values_from=value)

orig.sig$p <- NA

for(i in 1:nrow(orig.sig)){
  temp <- orig.sig[i,]
  orig.sig$p[i] <- tryCatch(wilcox.test(temp$t[[1]], temp$nt[[1]])$p.value, 
           error=function(error) return (NA))
}
  
  
orig.rate <- orig.rate %>% 
  reshape2::melt(id.vars=c("per", "env", "group")) %>% 
  group_by(per, env, group) %>% 
  do(boot.summary(.$value)) 

orig.rate <- orig.rate %>%  pivot_wider(id_cols=c(per, group), names_from=c(env), values_from=mean) %>% 
  group_by(per, group) %>% 
  summarise(mag=nt/t)

orig.rate[orig.rate$group == "all",]$mag <- tapply(orig.rate[orig.rate$group != "all",]$mag, orig.rate[orig.rate$group != "all",]$per, mean, na.rm=TRUE)

# Net Extinction ----------------------------------------------------------
ext.rate <- ext.rate %>% 
  filter(mid <67 & method =="sqs")  %>% 
  left_join(stages[,c("mid", "per")])

ext.rate <- ext.rate %>% 
  bind_rows(ext.rate %>%  
              mutate(per = "all")) 


ext.rate <- ext.rate %>% 
  dplyr::select(-method, -mid) 

# significant difference
ext.sig <- ext.rate %>% 
  reshape2::melt(id.vars=c("per", "env", "group")) %>% 
  group_by(per, env, group) 

ext.sig <- ext.sig %>% 
  pivot_wider(id_cols = c(per, group), names_from=env, values_from=value)

ext.sig$p <- NA

for(i in 1:nrow(ext.sig)){
  temp <- ext.sig[i,]
  ext.sig$p[i] <- tryCatch(wilcox.test(temp$t[[1]], temp$nt[[1]])$p.value, 
                            error=function(error) return (NA))
}

ext.rate <- ext.rate %>% 
  reshape2::melt(id.vars=c("per", "env", "group")) %>% 
  group_by(per, env, group) %>% 
  do(boot.summary(.$value)) 

ext.rate <- ext.rate %>%  pivot_wider(id_cols=c(per, group), names_from=c(env), values_from=mean) %>% 
  group_by(per, group) %>% 
  summarise(mag=nt/t)

# ext.rate[ext.rate$group == "all",]$mag <- tapply(ext.rate[ext.rate$group != "all",]$mag, ext.rate[ext.rate$group != "all",]$per, mean, na.rm=TRUE)


# Migration ---------------------------------------------------------------

mig.tot2 <- mig.tot %>% 
  filter(env=="env1" & method=="sqs") %>% 
  reshape2::melt(id.vars=c("group", "env", "method", "bin")) %>% 
  left_join(stages %>% dplyr::select(mid, stg, per), by=c("bin"="stg")) 

#significance test
mig.sig <- mig.tot2 %>% 
  bind_rows(mig.tot2 %>% mutate(per = "all"))  %>% 
  filter(!is.na(per)) %>% 
  pivot_wider(id_cols = c(per, group), names_from=variable, values_from=value)

mig.sig$p <- NA

for(i in 1:nrow(mig.sig)){
  temp <- mig.sig[i,]
  mig.sig$p[i] <- tryCatch(wilcox.test(temp$t[[1]], temp$nt[[1]])$p.value, 
                           error=function(error) return (NA))
}

mig.tot2 <- mig.tot2 %>% 
  bind_rows(mig.tot2 %>% mutate(per = "all"))  %>% 
  group_by(group,per, variable) %>% 
  summarise(mean=mean(value, na.rm=TRUE))

mig.tot2 <-  mig.tot2 %>% 
  ungroup() %>% 
  #proportion imported to proportion exported
  mutate(variable =ifelse(variable=="nt", "t", "nt") ) %>% 
  pivot_wider(names_from = variable,
              names_glue = "{variable}_{.value}",
              values_from = mean) %>% 
  group_by(group, per) %>% 
  summarise(mag=nt_mean/t_mean)


# Compare rates between zones ---------------------------------

comp <- rbind(
  cbind(orig.rate %>% left_join(orig.sig[,c("per", "group", "p")]), rate="origination"),
  cbind(ext.rate %>% left_join(ext.sig[,c("per", "group", "p")]), rate="extinction"),
  cbind(mig.tot2 %>% left_join(mig.sig[,c("per", "group", "p")]), rate="migration")) %>% 
  filter(!is.na(per)) %>% 
  mutate(per=factor(per, levels=c("all", names(phases))),
         rate=factor(rate, levels=c("origination", "extinction", "migration")), 
         mag2 = case_when(p < 0.05 ~ paste0(round(mag,2), "*"),
                         TRUE ~ as.character(round(mag, 2)))
  )



comp <- comp %>%  ungroup %>% 
  dplyr::select(per, group, mag=mag2, rate) %>% 
  mutate(group=recode(group, !!!group_names), 
         per = factor(per, levels=c("all", names(phases)))) 

reshape2::dcast(data=comp, formula=group+rate~per, value.var = "mag") %>% 
  arrange(rate, group) %>% 
  setNames(c("Group", "Rate", "All", per_names)) %T>% 
  write.csv(file.path("output", "magnitude_comparisons.csv"), row.names = FALSE)
