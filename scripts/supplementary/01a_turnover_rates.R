# This script calculate the turnover rates per plankton group

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
data(stages, package="divDyn")

load(file.path("data", "wrangled_data.RData")) # sp at genus level

dat[dat$fossil.group == "FALSE",]$fossil.group <- "F"
dat <- dat %>% filter(FAD != LAD)


# Calculating turnover rates ---------------------------------------------------

# *a. Creating data.frames without migration ----------------------------------

# data.frame to calculate origination

orig.df <- dat
orig.sp <- unique(dat[,c("sp", "bin", "mig1")])
orig.sp <- orig.sp[orig.sp$mig1 == 1,]

#first need to remove migration
for (i in 1:nrow(orig.sp)){
  cat("\r", i, "out of", nrow(orig.sp))
  
  temp <- orig.sp[i,]
  start <- temp$bin
  
  orig.df <- orig.df[!(orig.df$sp == temp$sp & orig.df$bin %in% c(start:95)),]
}

nrow(dat) - nrow(orig.df) #number of observations that counted as migration removed
nrow(orig.df)

# calculate origination rates for all groups
orig.df <- orig.df %>% 
  bind_rows(orig.df %>%  mutate(fossil.group="all"))

# data.frame to calculate extinction
ext.df <- dat
ext.sp <- unique(dat[,c("sp", "bin", "mig1")])
ext.sp <- ext.sp[ext.sp$mig1 == 1,]

#first need to remove migration
for (i in 1:nrow(ext.sp)){
  cat("\r", i, "out of", nrow(ext.sp))
  
  temp <- ext.sp[i,]
  end <- temp$bin
  
  ext.df <- ext.df[!(ext.df$sp == temp$sp & ext.df$bin %in% c(70:end)),]
}

nrow(dat) - nrow(ext.df) #number of observations that counted as migration removed
nrow(ext.df)

# calculate extinction rates for all groups
ext.df <- ext.df %>% 
  bind_rows(ext.df %>%  mutate(fossil.group="all"))


# * b. Origination --------------------------------------------------------
it <- 1000 #number of iterations
q=0.7 #sqs quorum

orig.rate <- list()

ss <- orig.df %>% filter(bin > 78) %>% 
  group_by(fossil.group, env1, bin) %>% 
  tally() %>% 
  filter(n > 100) %>% 
  summarise(n=round(min(n)*0.9))

for(g in names(group_names)){
  for(e in c("t", "nt")){
    temp <- orig.df %>%  filter(fossil.group == g & env1==e)
    
    for(me in c("sqs", "cr")){
    #calculate quota
    dd <- switch(me,
           sqs=divDyn::subsample(temp, q=q, tax="sp", bin="bin", coll="sample_id", type="sqs", iter=it, output="dist"),
           cr=divDyn::subsample(temp, q=ss$n[ss$fossil.group == g & ss$env1==e], 
                        tax="sp", bin="bin", type="cr", iter=it, output="dist"))
      
    orig.rate[[paste0(g, e,me)]] <- dd$ori2f3 %>% 
      as.data.frame() %>% 
      mutate(mid=stages$mid, group=g, env=e, method=me) %>% 
      filter(mid < 67)
  }
}
}

orig.rate <- do.call(rbind.data.frame, orig.rate)


# * c. Extinction --------------------------------------------------------
ext.rate <- list()

ss <- ext.df %>% filter(bin > 78) %>% 
  group_by(fossil.group, env1, bin) %>% 
  tally() %>% 
  filter(n > 100) %>% 
  summarise(n=round(min(n)*0.9))


for(g in names(group_names)){
  for(e in c("t", "nt")){
    temp <- ext.df %>%  filter(fossil.group == g & env1==e)
    
    for(me in c("sqs", "cr")){
      #calculate quota
      
      dd <- switch(me,
                   sqs=divDyn::subsample(temp, q=q, tax="sp", bin="bin", coll="sample_id", type="sqs", iter=it, output="dist"),
                   cr=divDyn::subsample(temp, q=ss$n[ss$fossil.group == g & ss$env1==e], 
                                      tax="sp", bin="bin", type="cr", iter=it, output="dist"))
      
    ext.rate[[paste0(g, e,me)]] <- dd$ext2f3 %>% 
      as.data.frame() %>% 
      mutate(mid=stages$mid, group=g, env=e, method=me) %>% 
      filter(mid < 67)
    }
  }
}

ext.rate <- do.call(rbind.data.frame, ext.rate)

save(orig.rate, ext.rate, file=file.path("output", "turnover_rates.RData"))

# Occurrence table --------------------------------------------------------
dat %>%  filter(bin > 79) %>% 
  group_by(bin, stage, fossil.group, env1) %>% 
  tally() %>% 
  pivot_wider(names_from=c(fossil.group, env1), 
              names_glue="{fossil.group}_{env1}", 
              values_from=n ) %T>%
  write.csv(file.path("output", "occurrence_table.csv"), row.names = FALSE)



