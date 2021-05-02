# Calculates and plots diversification preferences at genus level

# Load functions ----------------------------------------------------------
library(here)
library(tidyverse)
library(divDyn)

source(here("code", "functions.R"))
source(here("code", "plot_parameters.R"))
# Load data ---------------------------------------------------------------
load(here("data", "wrangled_data_genus.RData")) # sp at genus level

data(stages, package="divDyn")

stages$per <- categorize(stages$stage, phases) 
stages <- stages[!is.na(stages$per),]

dat[dat$fossil.group == "FALSE",]$fossil.group <- "F"
dat <- dat %>% filter(FAD != LAD)

dat<- dat %>% left_join(stages[,c("mid", "per")])

dat <- dat %>% bind_rows(dat %>% mutate(per = "all"))

dat$per <- factor(dat$per, levels=c("all",names(phases)))

# Bootstrap odds ratio -------------------------------------------
it <- 1000
chunk <- 20
n <- it/chunk

pb <- txtProgressBar(min = 0, max = chunk, style = 3)

geog.ev <- levels(dat$per)
prefs <- list()
mig.tot <- list()

for(p in geog.ev){
  cat("\n", p)
  df <- dat %>%  filter(per==p)
  df$ori <- as.factor(ifelse(df$ori == 1, "o", "no"))
  df$ext <- as.factor(ifelse(df$ext == 1, "e", "ne"))
  
  me="sqs"      #do in chunks so as not to overload
  for(i in 1:chunk){
    #create replicates of subsampled data based on number of iterations
    setTxtProgressBar(pb, i)
    sub.df <- replicate(n, bootstr(df,me), simplify = FALSE)
    
    s="env1"
    m <- gsub("env", "mig", s)
    
    finalCoefs <- lapply(sub.df, function(x) {
      tryCatch({
        temp <- df[x, ]
    
      temp_ori <- temp[,c("ori", s, "sp")]
      temp_ori <- unique(temp_ori)
      
      temp_ori <- setNames(caret::upSample(temp[,s], temp$ori),
                           c("env1", "ori"))
      
      temp_ori$varo <- paste0(temp_ori$ori, unlist(temp_ori[,s]))
      
      temp_ext <- temp[,c("ext", s, "sp")]
      temp_ext <- unique(temp_ext)
      temp_ext <- setNames(caret::upSample(temp[,s], temp$ext),
                           c("env1", "ext"))
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
    
    prefs[[paste(p,i)]] <- cbind(
      as.data.frame(data.table::rbindlist(finalCoefs, fill=TRUE)),
      per=p)
    
    #save data
    mig.tot[[paste(p,i)]] <- cbind(as.data.frame(
      data.table::rbindlist(finalProps, fill=TRUE)),
      per=p)
  }
}

prefs2 <- data.table::rbindlist(prefs, fill=TRUE)
save(prefs2, file=here("output", "prefs_gen.RData"))

mig.tot <-as.data.frame(data.table::rbindlist(mig.tot, fill=TRUE))
save(mig.tot, file=here("output", "mig_overall_genus.RData"))

# Plot --------------------------------------------------------------------
load(here("output", "prefs_gen.RData"))
load(here("output", "mig_overall_genus.RData"))

prefs2$per <- factor(prefs2$per, levels=c("all", names(phases)))

prefs_ori <- prefs2 %>% 
  dplyr::select(nont, not, ont, ot, per) %>% 
  group_by(per) %>% 
  summarise_if(is.numeric, mean, na.rm=TRUE) %>% 
  group_by(per) %>% 
  do(oddsratio(.$ot, .$not, .$ont, .$nont))

 
p1 <- ggplot(prefs_ori, aes(x=per, y=OR, ymin=LowerCI, ymax=UpperCI, col=per)) +
  geom_hline(yintercept = 0, col="darkgrey", linetype="dashed")+
  geom_point(size=3) +
  geom_errorbar(width=0.05, size=1) +
  scale_x_discrete(labels=c("all"="All", per_names)) +
  scale_color_manual(values=c("black", clim_clr), guide=FALSE) +
  scale_y_continuous(trans="reverse") +
  labs(y="Origination Preferences") +
  coord_cartesian(xlim=c(1,length(phases)+1), clip="off")+
  theme_hce +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1),
        plot.margin = unit(c(1,1,1,3), "lines"),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

prefs_ext <- prefs2 %>% 
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
  scale_x_discrete(labels=c("all"="All", per_names)) +
  scale_color_manual(values=c("black", clim_clr), guide=FALSE) +
  scale_y_continuous(trans="reverse") +
  labs(y="Extinction Preferences") +
  coord_cartesian(xlim=c(1,length(phases)+1), clip="off")+
  theme_hce +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1),
        plot.margin = unit(c(1,1,1,3), "lines"),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

mig.tot2 <- mig.tot %>% 
  reshape2::melt(id.vars=c("per")) %>%  
  mutate(value=as.numeric(value)) %>% 
  group_by(per, variable) %>% 
  do(boot.summary2(.$value)) %>% 
  ungroup() 

mig.tot2 <- mig.tot2 %>% 
  pivot_wider(
    names_from = variable,
    names_glue = "{variable}_{.value}",
    values_from = c(mean, lwr, upr)
  )

mig.tot2$per <- factor(mig.tot2$per, levels=c("all", names(phases)))

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
  coord_cartesian(xlim=c(0, 0.6), ylim=c(0,0.6)) +
  scale_color_manual(values=c("black", clim_clr), 
                     labels=c("all"="All", per_names)) +
  labs(x="Proportion exported from the tropics", y="Proportion exported from the extratropics",
       col="Climate Period") +
  theme_hce

prefs_div <- prefs2 %>% 
  dplyr::select(nont, not, ont, ot, per) %>% 
  group_by(per) %>% 
  summarise_if(is.numeric, mean, na.rm=TRUE) %>% 
  left_join(
    prefs2 %>% 
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
  geom_errorbar(width=0.05, size=1) +
  geom_point(size=3) +
  scale_x_discrete(labels=c("all"="All", per_names)) +
  scale_color_manual(values=c("black", clim_clr), guide=FALSE) +
  scale_y_continuous(trans="reverse") +
  labs(y="Diversification Preferences") +
  coord_cartesian(xlim=c(1,length(phases)+1), clip="off")+
  theme_hce +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1),
        plot.margin = unit(c(1,1,1,3), "lines"),
        axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)))

svg(here("figs", "Supplement", "Fig_03_preferences_migration_genus.svg"), 
    w=12, h=12)
p1+p2+p4 + p3 +
  plot_annotation(tag_prefix = "(", tag_levels = "a", tag_suffix = ")") +
  plot_layout(nrow=2, guides = "collect") &
  theme(legend.position='bottom')
dev.off()

