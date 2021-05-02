# This script calculates the proportion of taxa migrating from one zone
# to the other, per plaknton group.

# Load functions ----------------------------------------------------------
library(tidyverse)
library(divDyn)

source(file.path("scripts", "utils","functions.R"))
source(file.path("scripts","utils", "plot_parameters.R"))

# Load data ---------------------------------------------------------------
load(file.path("data", "wrangled_data.RData")) # sp at genus level

#dat$sp <- droplevels(dat$sp) #old levels kept
dat[dat$fossil.group == "FALSE",]$fossil.group <- "F"
dat <- dat %>% filter(FAD != LAD)

dat<- dat %>% mutate(per = case_when(mid > 50 ~ "G",
                                     mid < 50 & mid > 32 ~ "T1",
                                     mid < 32 & mid > 14 ~"T2",
                                     mid < 14 ~ "I")) 


dat <- dat %>% bind_rows(dat %>% mutate(fossil.group = "all"))

dat <- dat %>% bind_rows(dat %>% mutate(per = "all"))

dat$per <- factor(dat$per, levels=c("all","G", "T1", "T2", "I"))

# Bootstrap migration -----------------------------------------------------

it <- 1000

mig.tot <- list()
per <- levels(dat$per)

chunk <- 20
n <- it/chunk

pb <- txtProgressBar(min = 0, max = chunk, style = 3)

#cr quota 
ss <- dat %>% 
  filter(bin > 78) %>% 
  group_by(bin,fossil.group) %>% 
  tally() %>% 
  filter(n > 100) %>% 
  group_by(fossil.group) %>% 
  summarise(n=round(min(n)*0.9))

for(f in names(group_names)){
  df <- dat %>%  filter(fossil.group==f) %>%  dplyr::select(sample_id, sp, bin, env1, env2, mig1, mig2)
  
  for(me in c("sqs", "cr")){
  #do in chunks so as not to overload
  for(i in 1:chunk){
    #create replicates of subsampled data based on number of iterations
    setTxtProgressBar(pb, i)
    sub.df <- replicate(n, bootstr(df,me,ss$n[ss$fossil.group==f]), simplify = FALSE)
    
    
    for(s in c("env1", "env2")){
      m <- gsub("env", "mig", s)
      
      
      #grouping variables for calculation
      grouping_vars <- quos(s, m)
      
      finalCoefs <- lapply(sub.df, function(x) {
        #calculate proportion imported in each environment
        df[x,] %>% dplyr::select(sp, bin, {{s}}, {{m}}) %>%
          distinct() %>%
          group_by(.dots = list(s, m), bin) %>% tally() %>%
          setNames(c("bin", "env", "mig", "n")) %>%
          group_by(env) %>%
          mutate(prop = n/sum(n)) %>%
          ungroup() %>%
          filter(mig > 0) %>% 
          dplyr::select(-n, -mig) %>% 
          spread(env, prop)
      }
      )
      
      #save data
      mig.tot[[paste0(f,s,i,me)]] <- data.table::rbindlist(finalCoefs, fill=TRUE) %>%  
        mutate(group=f, env=s, method=me)
    }
  }
  }
}

mig.tot <-as.data.frame(data.table::rbindlist(mig.tot, fill=TRUE))

save(mig.tot, file=file.path("output", "mig_stage.RData"))

