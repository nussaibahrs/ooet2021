# Calculates and plots net origination & extinction rates, and migration proportions
# for each ocean basin

# load libraries ---------------------------------------------------------
library(tidyverse)
library(divDyn)
library(magrittr)
library(reshape2)
library(foreach)
library(doParallel)
library(raster)
library(sp)

source(file.path("scripts", "utils","functions.R"))
source(file.path("scripts","utils", "plot_parameters.R"))

# Load data ---------------------------------------------------------------
data(stages, package="divDyn")

load(file.path("output", "ocean_basins.RData"))
load(file.path("data", "wrangled_data.RData")) # sp at species level

#dat$sp <- droplevels(dat$sp) #old levels kept
dat[dat$fossil.group == "FALSE",]$fossil.group <- "F"
dat <- dat %>% filter(FAD != LAD)

dat <- dat %>% dplyr::select(sp, fossil.group, sample_id, age, paleolat, paleolong, bin, new_lat, new_long)

#merge east and west pacific
basins$pacific <- raster::union(basins$westpacific, basins$eastpacific)
basins <- basins[-grep("pacific", names(basins)[-length(basins)])]

# create data for basin ---------------------------------------------------
oceans <- names(basins)
ocean.df <- list()

for(o in oceans){
  cat(o)
  #filter occurrences by ocean
  dat.sp <- dat
  coordinates(dat.sp) <- ~paleolong + paleolat
  
  dat.ocean <- sp::over(dat.sp, basins[[o]])
  
  dat.ocean <- dat[names(dat.ocean[!is.na(dat.ocean)]),]
  
  #### generate taxa list with FAD and LAD ----
  tx.list <- dat.ocean %>% group_by(sp, fossil.group) %>%
    summarise(FADage = max(age), FAD = min(bin), 
              LADage = min(age), LAD= max(bin), 
              dur = FADage-LADage) %>%
    filter(FADage <= thres_age) %>%
    filter(FAD != LAD) #remove single interval occurring taxa
  
  dat.ocean <- dat.ocean %>% filter(sp %in% tx.list$sp) #only keep the ones within the tx.list
  dat.ocean <- droplevels(dat.ocean)
  
  #### add origination and extinction data ----
  dat.ocean <- dat.ocean %>% 
    left_join(tx.list, by=c("sp", "fossil.group")) %>%
    mutate(ori = ifelse(FAD == bin, 1, 0), #add origination 
           ext = ifelse(LAD == bin, 1, 0)) #add extinction
  
  head(dat.ocean)
  
  ### Assigning environment
  dat.ocean <- dat.ocean %>%
    mutate(new_lat = abs(new_lat), 
           env1 = ifelse(new_lat < 30, "t", "nt"))
  
  #assigning migration
  dat.ocean$mig1 <- 0
  
  for (i in 1:nrow(tx.list)){
    tx.temp <- tx.list$sp[i]
    
    #get original environment: majority taken if double
    ori_env <- dat.ocean %>% filter(sp == tx.temp & bin == FAD) %>% 
      group_by(env1) %>% tally() %>% filter(n == max(n)) %>% pull(env1)
    
    temp <- dat.ocean[dat.ocean$sp == tx.temp,] %>% as.data.frame() 
    temp[temp$bin == temp$FAD,]$env1 <- ori_env #replace 
    
    x <- table(temp$bin, temp$env1) %>% as.data.frame() %>% 
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
          
          dat.ocean[dat.ocean$sp == tx.temp & dat.ocean$bin == tbin[t] & dat.ocean$env1 != ori_env,]$mig1 <- 1
          ori_env <- tbin.temp$env
          
        } 
        
      }
    }
  }
  
  ocean.df[[o]] <- dat.ocean
}

save(ocean.df, file=file.path("output", "data_per_ocean.RData"))


# Calculate rates ---------------------------------------------------------
load(file.path("output", "data_per_ocean.RData"))

# Calculating turnover rates ---------------------------------------------------
it <- 1000 #number of iterations
q=0.7 #sqs quorum

ocean.rate <- list()

for(o in 1:length(ocean.df)){
  lab <- names(ocean.df)[o]
  dat <- ocean.df[[o]]
  
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
  orig.rate <- list()
  
  for(g in names(group_names)){
    for(e in c("t", "nt")){
      temp <- orig.df %>%  filter(fossil.group == g & env1==e)
      
      
      #calculate quota
      dd <- divDyn::subsample(temp, q=q, tax="sp", bin="bin", coll="sample_id", type="sqs", iter=it, output="dist")
      
      orig.rate[[paste0(g, e)]] <- apply(dd$ori2f3, 1, boot.summary2) %>% 
        bind_rows() %>% 
        mutate(mid=stages$mid[1:nrow(.)], group=g, env=e, ocean=lab)
    }
  }
  
  
  orig.rate <- do.call(rbind.data.frame, orig.rate)
  
  
  # * c. Extinction --------------------------------------------------------
  ext.rate <- list()
  
  for(g in names(group_names)){
    for(e in c("t", "nt")){
      temp <- ext.df %>%  filter(fossil.group == g & env1==e)
      
      for(me in c("sqs", "cr")){
        #calculate quota
        
        dd <- divDyn::subsample(temp, q=q, tax="sp", bin="bin", coll="sample_id", type="sqs", iter=it, output="dist")
        
        ext.rate[[paste0(g, e)]] <- apply(dd$ext2f3, 1, boot.summary2) %>% 
          bind_rows() %>% 
          mutate(mid=stages$mid[1:nrow(.)], group=g, env=e, ocean=lab)
      }
    }
  }
  
  ext.rate <- do.call(rbind.data.frame, ext.rate)
  
  ocean.rate[[lab]] <- rbind(
    cbind(orig.rate, rate="origination"),
    cbind(ext.rate, rate="extinction")
  )
}

mig.tot <- list()
per <- levels(dat$per)

chunk <- 20
n <- it/chunk

pb <- txtProgressBar(min = 0, max = chunk, style = 3)

for(o in 1:length(ocean.df)){
  lab <- names(ocean.df)[o]
  dat <- ocean.df[[o]]
  
  dat <- dat %>%  bind_rows(dat %>%  mutate(fossil.group="all"))
  
  for(f in names(group_names)){
    df <- dat %>%  filter(fossil.group==f) %>%  dplyr::select(sample_id, sp, bin, env1, mig1)
    
    #do in chunks so as not to overload
    for(i in 1:chunk){
      #create replicates of subsampled data based on number of iterations
      
      sub.df <- replicate(n, bootstr(df,"sqs"), simplify = FALSE)
      
      
      s="env1"
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
      mig.tot[[paste0(f,i,o)]] <- data.table::rbindlist(finalCoefs, fill=TRUE) %>%  
        mutate(group=f, ocean=lab)
      
      setTxtProgressBar(pb, i)
    }
  }
}


save(ocean.rate, mig.tot, file=file.path("output", "ocean_turnover_rates.RData"))


# Plot --------------------------------------------------------------------
load(file.path("output", "ocean_turnover_rates.RData"))

data(stages, package="divDyn")

stages$per <- categorize(stages$stage, phases) 
stages <- stages[!is.na(stages$per),]


mig.tot2 <- mig.tot %>% 
  bind_rows() %>% 
  filter(group == "all") %>% 
  reshape2::melt(id.vars=c("group", "bin", "ocean")) %>% 
  left_join(int.tbl, by=c("bin"="stg")) %>% 
  left_join(stages[,c("mid", "per")]) 
  
mig.tot2 <- mig.tot2 %>% 
  bind_rows(mig.tot2 %>%  mutate(per="all")) %>% 
group_by(per,variable, ocean) %>% 
  do(boot.summary2(.$value)) %>% 
  ungroup() %>% 
  #proportion imported to proportion exported
  mutate(variable =ifelse(variable=="nt", "t", "nt") )

mig.tot2 <- mig.tot2 %>% 
  pivot_wider(
    names_from = variable,
    names_glue = "{variable}_{.value}",
    values_from = c(mean, lwr, upr)
  )

mig.tot2 <- mig.tot2 %>%  filter(!is.na(per))

ocean.rate2 <- ocean.rate %>% 
  bind_rows() %>% 
  filter(group == "all" & mid < 67) %>% 
  left_join(stages[,c("mid", "per")]) 


ocean.rate2 <- ocean.rate2 %>%
  bind_rows(ocean.rate2 %>%  mutate(per="all")) 

ocean.rate2 <- ((ocean.rate2 %>% filter(env == "t") %>% 
    dplyr::select(mean, lwr, upr)
) - (ocean.rate2 %>% filter(env == "nt") %>% 
       dplyr::select(mean, lwr, upr)
)) %>% bind_cols(
  ocean.rate2 %>% filter(env == "t") %>% 
    dplyr::select(-mean, -lwr, -upr)) %>% 
  group_by(ocean, rate, per) %>% 
  do(boot.summary2(.$mean))


mig.tot2$per <- factor(mig.tot2$per, levels=c("all", names(phases)))
ocean.rate2$per <- factor(ocean.rate2$per, levels=c("all", names(phases)))

# Atlantic
p1 <- ggplot(ocean.rate2 %>% filter(ocean=="atlantic" & rate=="origination"),
             aes(x=per, 
                 y=mean, ymin=lwr, ymax=upr, 
                 col=per)) +
  geom_point() +
  geom_errorbar(width=0.05) +
  facet_grid(~ocean, labeller = labeller(env=env_names, 
                                            ocean=c(atlantic="Atlantic", indian="Indian", pacific="Pacific")))+
  scale_x_discrete(labels=c("all"="All", per_names)) +
  scale_color_manual(values=c("black", clim_clr), guide=FALSE) +
  scale_y_continuous(trans="reverse", limits=c(2,-2)) +
  labs(y="Net origination rate") +
  coord_cartesian(xlim=c(1,length(phases)+1), clip="off")+
  geom_hline(yintercept = 0, linetype="dashed", col="darkgrey")+
  theme_hce +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1),
        plot.margin = unit(c(1,1,1,3), "lines"),
        axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)))

lims <- c(-2,2)
w <- 1
y1 <- c(lims[1] + w, lims[2] - w)
x1 <- 0.2

p1 <- p1 + annotate("segment", x=-x1, xend=-x1, y=-y1[1], yend=-lims[1],  
                    arrow = arrow(length = unit(0.03, "npc")), 
                    col="darkgrey") +
  annotate("segment", x=-x1, xend=-x1, y=-y1[2], yend=-lims[2],  
           arrow = arrow(length = unit(0.03, "npc")), 
           col="darkgrey") +
  annotate("text", x=-(x1+0.05), y=-y1[1], label="Tropics", angle=90, 
           vjust=-0.2, hjust=1,
           col="darkgrey") +
  annotate("text", x=-(x1+0.05), y=-y1[2], label="Extratropics", angle=90, 
           vjust=-0.2,hjust=0,
           col="darkgrey") 


p2 <- ggplot(ocean.rate2 %>% filter(ocean=="atlantic" & rate=="extinction"),
             aes(x=per, 
                 y=mean, ymin=lwr, ymax=upr, 
                 col=per)) +
  geom_point() +
  geom_errorbar(width=0.05) +
  facet_grid(~ocean, labeller = labeller(env=env_names, 
                                         ocean=c(atlantic="Atlantic", indian="Indian", pacific="Pacific")))+
  scale_x_discrete(labels=c("all"="All", per_names)) +
  scale_color_manual(values=c("black", clim_clr), guide=FALSE) +
  scale_y_continuous(trans="reverse", limits=c(2,-2)) +
  labs(y="Net extinction rate") +
  coord_cartesian(xlim=c(1,length(phases)+1), clip="off")+
  geom_hline(yintercept = 0, linetype="dashed", col="darkgrey")+
  theme_hce +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1),
        plot.margin = unit(c(1,1,1,3), "lines"),
        axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)))

lims <- c(-2,2)
w <- 1
y1 <- c(lims[1] + w, lims[2] - w)
x1 <- 0.2

p2 <- p2 + annotate("segment", x=-x1, xend=-x1, y=-y1[1], yend=-lims[1],  
                    arrow = arrow(length = unit(0.03, "npc")), 
                    col="darkgrey") +
  annotate("segment", x=-x1, xend=-x1, y=-y1[2], yend=-lims[2],  
           arrow = arrow(length = unit(0.03, "npc")), 
           col="darkgrey") +
  annotate("text", x=-(x1+0.05), y=-y1[1], label="Tropics", angle=90, 
           vjust=-0.2, hjust=1,
           col="darkgrey") +
  annotate("text", x=-(x1+0.05), y=-y1[2], label="Extratropics", angle=90, 
           vjust=-0.2,hjust=0,
           col="darkgrey") 



p3 <- ggplot(mig.tot2 %>%  filter(ocean == "atlantic"), 
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
  geom_errorbarh(height=0, alpha=0.3, size=1)+
  geom_errorbar(width=0, alpha=0.3, size=1) +
  geom_abline(slope=1, intercept = 0, linetype="dashed", col="darkgrey") +
  coord_cartesian(xlim=c(0, 0.03), ylim=c(0,0.03)) +
  scale_color_manual(values=c("black", clim_clr), 
                     labels=c("all"="All", per_names), 
                     guide = guide_legend(nrow=1)) +
  labs(x="Proportion exported from the tropics", y="Proportion exported \nfrom the extratropics",
       col="Climate Phase") +
  facet_grid(~ocean, labeller=labeller(ocean=c(atlantic="Atlantic", indian="Indian", pacific="Pacific")))+
  theme_hce +
  theme(axis.title = element_text(size=10, face="bold"),
        legend.title = element_text(size=10, face="bold"))

# Indian
p4 <- ggplot(ocean.rate2 %>% filter(ocean=="indian" & rate=="origination"),
             aes(x=per, 
                 y=mean, ymin=lwr, ymax=upr, 
                 col=per)) +
  geom_point() +
  geom_errorbar(width=0.05) +
  facet_grid(~ocean, labeller = labeller(env=env_names, 
                                         ocean=c(atlantic="Atlantic", indian="Indian", pacific="Pacific")))+
  scale_x_discrete(labels=c("all"="All", per_names)) +
  scale_color_manual(values=c("black", clim_clr), guide=FALSE) +
scale_y_continuous(trans="reverse", limits=c(2,-2)) +
  labs(y="Net origination rate") +
  coord_cartesian(xlim=c(1,length(phases)+1), clip="off")+
  geom_hline(yintercept = 0, linetype="dashed", col="darkgrey")+
  theme_hce +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1))

p5 <- ggplot(ocean.rate2 %>% filter(ocean=="indian" & rate=="extinction"),
             aes(x=per, 
                 y=mean, ymin=lwr, ymax=upr, 
                 col=per)) +
  geom_point() +
  geom_errorbar(width=0.05) +
  facet_grid(~ocean, labeller = labeller(env=env_names, 
                                         ocean=c(atlantic="Atlantic", indian="Indian", pacific="Pacific")))+
  scale_x_discrete(labels=c("all"="All", per_names)) +
  scale_color_manual(values=c("black", clim_clr), guide=FALSE) +
scale_y_continuous(trans="reverse", limits=c(2,-2)) +
  labs(y="Net extinction rate") +
  coord_cartesian(xlim=c(1,length(phases)+1), clip="off")+
  geom_hline(yintercept = 0, linetype="dashed", col="darkgrey")+
  theme_hce +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1))

p6 <- ggplot(mig.tot2 %>%  filter(ocean == "indian"), 
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
  geom_errorbarh(height=0, alpha=0.3, size=1)+
  geom_errorbar(width=0, alpha=0.3, size=1) +
  geom_abline(slope=1, intercept = 0, linetype="dashed", col="darkgrey") +
  coord_cartesian(xlim=c(0, 0.03), ylim=c(0,0.03)) +
  scale_color_manual(values=c("black", clim_clr), 
                     labels=c("all"="All", per_names), 
                     guide = guide_legend(nrow=1)) +
  labs(x="Proportion exported from the tropics", y="Proportion exported \nfrom the extratropics",
       col="Climate Phase") +
  facet_grid(~ocean, labeller=labeller(ocean=c(atlantic="Atlantic", indian="Indian", pacific="Pacific")))+
  theme_hce+
  theme(axis.title = element_text(size=10, face="bold"),
        legend.title = element_text(size=10, face="bold"))

# Pacifics
p7 <- ggplot(ocean.rate2 %>% filter(ocean=="pacific" & rate=="origination"),
             aes(x=per, 
                 y=mean, ymin=lwr, ymax=upr, 
                 col=per)) +
  geom_point() +
  geom_errorbar(width=0.05) +
  facet_grid(~ocean, labeller = labeller(env=env_names, 
                                         ocean=c(atlantic="Atlantic", indian="Indian", pacific="Pacific")))+
  scale_x_discrete(labels=c("all"="All", per_names)) +
  scale_color_manual(values=c("black", clim_clr), guide=FALSE) +
scale_y_continuous(trans="reverse", limits=c(2,-2)) +
  labs(y="Net origination rate") +
  coord_cartesian(xlim=c(1,length(phases)+1), clip="off")+
  geom_hline(yintercept = 0, linetype="dashed", col="darkgrey")+
  theme_hce +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1))

p8 <- ggplot(ocean.rate2 %>% filter(ocean=="pacific" & rate=="extinction"),
             aes(x=per, 
                 y=mean, ymin=lwr, ymax=upr, 
                 col=per)) +
  geom_point() +
  geom_errorbar(width=0.05) +
  facet_grid(~ocean, labeller = labeller(env=env_names, 
                                         ocean=c(atlantic="Atlantic", indian="Indian", pacific="Pacific")))+
  scale_x_discrete(labels=c("all"="All", per_names)) +
  scale_color_manual(values=c("black", clim_clr), guide=FALSE) +
scale_y_continuous(trans="reverse", limits=c(2,-2)) +
  labs(y="Net extinction rate") +
  coord_cartesian(xlim=c(1,length(phases)+1), clip="off")+
  geom_hline(yintercept = 0, linetype="dashed", col="darkgrey")+
  theme_hce +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1))

p9 <- ggplot(mig.tot2 %>%  filter(ocean == "pacific"), 
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
  geom_errorbarh(height=0, alpha=0.3, size=1)+
  geom_errorbar(width=0, alpha=0.3, size=1) +
  geom_abline(slope=1, intercept = 0, linetype="dashed", col="darkgrey") +
  coord_cartesian(xlim=c(0, 0.03), ylim=c(0,0.03)) +
  scale_color_manual(values=c("black", clim_clr), 
                     labels=c("all"="All", per_names), 
                     guide = guide_legend(nrow=1)) +
  labs(x="Proportion exported from the tropics", y="Proportion exported \nfrom the extratropics",
       col="Climate Phase") +
  facet_grid(~ocean, labeller=labeller(ocean=c(atlantic="Atlantic", indian="Indian", pacific="Pacific")))+
  theme_hce+
  theme(axis.title = element_text(size=10, face="bold"),
        legend.title = element_text(size=10, face="bold"))


svg(file.path("figs", "Supplement", "Fig_S_turnover_rates_ocean.svg"),
    w=12, h=8)
p1+p4 + p7 +#origination 
  p2+p5+ p8 + #extinction
  plot_annotation(tag_prefix = "(", tag_levels = "a", tag_suffix = ")") +
  plot_layout(ncol=3, guides = "collect") &
  theme(legend.position='bottom')
dev.off()

svg(file.path("figs", "Supplement", "Fig_S_migration_ocean.svg"),
    w=12, h=4)
p3+p6+p9+
  plot_annotation(tag_prefix = "(", tag_levels = "a", tag_suffix = ")") +
  plot_layout(ncol=3, guides = "collect") &
  theme(legend.position='bottom')
dev.off()

