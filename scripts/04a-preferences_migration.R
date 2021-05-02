# This script calculates origination, extinction and diversitfication 
# preferences, and dispersal proportions across climate zones over the 
# last 66 million years, divided into climate phases.

# Load functions ----------------------------------------------------------
library(tidyverse)
library(divDyn)

source(file.path("scripts","utils", "functions.R"))
source(file.path("scripts","utils", "plot_parameters.R"))
# Load data ---------------------------------------------------------------
load(file.path("data", "wrangled_data.RData")) # sp at genus level

#dat$sp <- droplevels(dat$sp) #old levels kept
dat[dat$fossil.group == "FALSE",]$fossil.group <- "F"
dat <- dat %>% filter(FAD != LAD)

dat$per <- categorize(dat$stage, 
           phases)


dat <- dat %>% bind_rows(dat %>% mutate(fossil.group = "all"))

dat <- dat %>% bind_rows(dat %>% mutate(per = "all"))

dat$per <- factor(dat$per, levels=c("all","W1", "H", "W2", "C", "I"))

# Bootstrap odds ratio -------------------------------------------
it <- 1000
chunk <- 20
n <- it/chunk

pb <- txtProgressBar(min = 0, max = chunk, style = 3)

geog.ev <- levels(dat$per)
prefs <- list()
mig.tot <- list()

#cr quota 
ss <- dat %>% 
  filter(bin > 78) %>% 
  group_by(bin,fossil.group,per) %>% 
  tally() %>% 
  filter(n > 100) %>% 
  group_by(fossil.group,per) %>% 
  summarise(n=round(min(n)*0.9))
ss <- na.omit(ss)

for(p in geog.ev){
  for(f in unique(dat$fossil.group)){
    df <- dat %>%  filter(fossil.group==f, per==p)
    df$ori <- as.factor(ifelse(df$ori == 1, "o", "no"))
    df$ext <- as.factor(ifelse(df$ext == 1, "e", "ne"))
    for(me in c("sqs", "cr")){

      #do in chunks so as not to overload
      for(i in 1:chunk){
        #create replicates of subsampled data based on number of iterations
        setTxtProgressBar(pb, i)
        sub.df <- replicate(n, bootstr(df,me,q=ss$n[ss$fossil.group==f & ss$per == p]
                                       ), simplify = FALSE)
        
        for(s in c("env1")){
          m <- gsub("env", "mig", s)
          
          finalCoefs <- lapply(sub.df, function(x) {
            tryCatch({
              temp <- df[x, ]
              
              temp_ori <- temp[,c("ori", s, "sp")]
              temp_ori <- unique(temp_ori)
              
              temp_ori <- setNames(caret::upSample(temp[,s], temp$ori),
                                   c(s, "ori"))
              
              temp_ori$varo <- paste0(temp_ori$ori, unlist(temp_ori[,s]))
              
              temp_ext <- temp[,c("ext", s, "sp")]
              temp_ext <- unique(temp_ext)
              temp_ext <- setNames(caret::upSample(temp[,s], temp$ext),
                                   c(s, "ext"))
              temp_ext$vare <- paste0(temp_ext$ext, unlist(temp_ext[,s]))
              
              cbind(
                temp_ori%>% group_by(varo) %>% tally() %>% spread(varo,n),
                temp_ext%>% group_by(vare) %>% tally()%>% spread(vare,n)
              )
            }, error=function(e) message(e)
            )
          }
          )
          
          #grouping variables for calculation
          grouping_vars <- quos(s, m)
          
          finalProps <- lapply(sub.df, function(x) {
            #calculate proportion imported in each environment
            df[x,] %>% dplyr::select(sp, {{s}}, {{m}}) %>%
              distinct() %>%
              group_by(.dots = list(s, m)) %>% tally() %>%
              setNames(c("env", "mig", "n")) %>%
              group_by(env) %>%
              mutate(prop = n/sum(n)) %>%
              ungroup() %>%
              filter(mig > 0) %>% 
              dplyr::select(-n, -mig) %>% 
              spread(env, prop)
          }
          )
          
          prefs[[paste(p,f,s,me,i)]] <- data.table::rbindlist(finalCoefs, fill=TRUE) %>% 
            mutate(group=f, per=p, method=me, env=s)
          
          #save data
          mig.tot[[paste(p,f,s,me,i)]] <- data.table::rbindlist(finalProps, fill=TRUE) %>% 
            mutate(group=f, per=p, method=me, env=s)
        }
      }
    }
  }
}

prefs2 <- data.table::rbindlist(prefs, fill=TRUE)
save(prefs2, file=file.path("output", "prefs_sp.RData"))

mig.tot <-as.data.frame(data.table::rbindlist(mig.tot, fill=TRUE))
save(mig.tot, file=file.path("output", "mig_overall_stage.RData"))

# Diversification rate ----------------------------------------------------
df <- dat %>%  filter(fossil.group == "all")

dd <- divDyn::subsample(df, tax="sp", bin="bin", q=0.7, type="sqs", coll = "sample_id", 
                        output="dist", it=it)

#calculate diversification, i.e. origination - extinction rate per iteration
div <- dd$ori2f3 - dd$ext2f3

save(div, file=file.path("output", "div_overall_stage.RData"))
