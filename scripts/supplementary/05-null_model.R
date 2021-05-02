# This script generates the null dataset for which all analyses for 
# diversification and dispersal are carried out.

#### load libraries and functions ----
library(tidyverse)
library(divDyn)
library(doParallel)

source(file.path("scripts", "utils","functions.R"))
source(file.path("scripts","utils", "plot_parameters.R"))

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

#### load data ####

load(file.path("data", "wrangled_data.RData")) # species level

dat[dat$fossil.group == "FALSE",]$fossil.group <- "F"
dat <- dat %>% filter(FAD != LAD & FAD > 81)
dat$per <- categorize(dat$stage, 
                      phases)
dat <- dat[!is.na(dat$per),]
dat <- dat[,c("sp", "env1", "mig1", "per", "FAD", "LAD", "ori", "ext", "bin", "sample_id")]

# Shuffle data ------------------------------------------------------------
geog.ev <- c("all","W1", "H", "W2", "C", "I")
prefs <- list()
mig.tot <- list()

tr=200

for (j in 1:tr){
  cat("\r", j)
  temp.g <- shuffle(.data = dat, perm_cols = c("env1"))  #shuffle columns 
  temp.g$ori <- as.factor(ifelse(temp.g$ori == 1, "o", "no"))
  temp.g$ext <- as.factor(ifelse(temp.g$ext == 1, "e", "ne"))
  temp.g$mig1 <- 0
  temp.g <- temp.g %>% bind_rows(temp.g %>% mutate(per = "all"))
  
  temp.g$per <- factor(temp.g$per, levels=geog.ev)
  
  tx.list <- unique(temp.g$sp)
  
  temp.g <- foreach(i=1:length(tx.list), .packages = "dplyr", .combine = rbind) %dopar% {
    tx.temp <- tx.list[i]
    
    #get original environment: majority taken if double
    ori_env <- temp.g %>% filter(sp == tx.temp & bin == FAD) %>% 
      group_by(env1) %>% tally() %>% filter(n == max(n)) %>% pull(env1)
    
    temp2 <- temp.g[temp.g$sp == tx.temp,] %>% as.data.frame() 
    temp2[temp2$bin == temp2$FAD,]$env1 <- ori_env #replace 
    
    if (nrow(temp2) > 0){
      x <- table(temp2$bin, temp2$env1) %>% as.data.frame() %>% 
        setNames(c("bin", "env", "n")) %>%
        mutate_if(is.factor, as.character) %>%
        arrange(bin)
      
      if(length(unique(x$env)) > 1) { #if two environments found, then migration occurred
        #get time bins during which taxa lived
        tbin <- x$bin %>% as.numeric() %>% sort() %>% unique()
        
        #compare with previous bin
        for (t in 2:length(tbin)){
          tbin.temp <- x[x$bin == tbin[t],] %>%
            filter(n >0)
          
          if ((nrow(tbin.temp) > 1 & length(ori_env) == 1) || #expansion if previous environment is single and current is double
              (tbin.temp$env != ori_env && length(ori_env) == 1)){ #if taxa left previous environment
            
            temp.g[temp.g$sp == tx.temp & temp.g$bin == tbin[t] & temp.g$env1 != ori_env,]$mig1 <- 1
            ori_env <- tbin.temp$env
            
          }
        }
      }
    }
    temp.g[temp.g$sp == tx.temp,]
  }
  
  s="env1"
  m <- gsub("env", "mig", s)
  
  for(p in geog.ev){
    #subsampling
    temp.p <- subset(temp.g, per==p)
    
    temp <- divDyn::subtrialSQS(temp.p, tax="sp", q=0.7, bin="bin")
    temp <- temp.p[temp,]
    
    
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
    
    finalCoefs <- cbind(
      temp_ori%>% group_by(varo) %>% tally() %>% spread(varo,n),
      temp_ext%>% group_by(vare) %>% tally()%>% spread(vare,n) %>% 
        mutate(it=j)
    )
    
    #grouping variables for calculation
    grouping_vars <- quos(s, m)
    
    finalProps <- temp %>% dplyr::select(sp, {{s}}, {{m}}) %>%
      distinct() %>%
      group_by(.dots = list(s, m)) %>% tally() %>%
      setNames(c("env", "mig", "n")) %>%
      group_by(env) %>%
      mutate(prop = n/sum(n)) %>%
      ungroup() %>%
      filter(mig > 0) %>% 
      dplyr::select(-n, -mig) %>% 
      spread(env, prop)
    
    prefs[[paste(p,j)]] <- finalCoefs %>% 
      mutate(per=p)
    
    #save data
    mig.tot[[paste(p,j)]] <- finalProps %>% 
      mutate(per=p)
  }
}

prefs2 <- data.table::rbindlist(prefs, fill=TRUE)
mig.tot <-as.data.frame(data.table::rbindlist(mig.tot, fill=TRUE))

save(prefs2, mig.tot, file=file.path("output", "null_results.RData"))


# Plot --------------------------------------------------------------------
load(file.path("output", "null_results.RData"))

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
  scale_x_discrete(labels=c("all"="Overall", per_names)) +
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
  scale_x_discrete(labels=c("all"="Overall", per_names)) +
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
                     labels=c("all"="Overall", per_names)) +
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

svg(file.path("figs", "Supplement", "Fig_03_preferences_migration_null.svg"), 
    w=12, h=12)
p1+p2+p4 + p3 +
  plot_annotation(tag_prefix = "(", tag_levels = "a", tag_suffix = ")") +
  plot_layout(nrow=2, guides = "collect") &
  theme(legend.position='bottom')
dev.off()