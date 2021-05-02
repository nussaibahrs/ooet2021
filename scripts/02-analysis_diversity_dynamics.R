# This script calculates the diversity of marine plankton over the last 
# 66 million years per climate zone and hemisphere.

# load libraries ---------------------------------------------------------
library(tidyverse)
library(divDyn)
library(magrittr)
library(reshape2)

source(file.path("scripts", "utils","functions.R"))
source(file.path("scripts","utils", "plot_parameters.R"))

# Load data ---------------------------------------------------------------
data(stages, package="divDyn")

load(file.path("data", "wrangled_data.RData")) # sp at genus level

dat[dat$fossil.group == "FALSE",]$fossil.group <- "F"
dat <- dat %>% filter(FAD != LAD)


# Calculating diversity ---------------------------------------------------
it <- 1000 #number of iterations
q=0.7

# * a. For all zones and hemispfile.path  --------------------------------------
dd_all <- divDyn::subsample(dat, q=q, tax="sp", bin="bin",  
                            type="sqs", coll="sample_id", it=it, output="dist")

dd_all <- apply(dd_all$divCSIB, 1, boot.summary) %>% 
  bind_rows()


# * b. per climate zone ---------------------------------------------------
dd_zone <- list()

for (e in c("t", "nt")){
  temp<- divDyn::subsample(subset(dat, env1==e), q=q, tax="sp", bin="bin",  
                           type="sqs", coll="sample_id", it=it, output="dist")
  
  dd_zone[[e]] <- apply(temp$divCSIB, 1, boot.summary) %>% 
    bind_rows() %>% 
    mutate(env=e)
}


dd_zone <- lapply(dd_zone, function(x) cbind(x, mid=stages$mid)) %>% 
  bind_rows()



# * c. Northern hemispfile.path: per clime zone ---------------------------------
dd_zoneN <- list()

for (e in c("t", "nt")){
  temp<- divDyn::subsample(subset(dat, env1==e & paleolat > 0), q=q, tax="sp", bin="bin",  
                           type="sqs", coll="sample_id", it=it, output="dist")
  
  dd_zoneN[[e]] <- apply(temp$divCSIB, 1, boot.summary) %>% 
    bind_rows() %>% 
    mutate(env=e)
}

dd_zoneN <- lapply(dd_zoneN, function(x) cbind(x, mid=stages$mid)) %>% 
  bind_rows()




# * d. Southern Hemispfile.path: per climate zone --------------------------------

dd_zoneS <- list()

for (e in c("t", "nt")){
  temp<- divDyn::subsample(subset(dat, env1==e & paleolat < 0), q=q, tax="sp", bin="bin",  
                           type="sqs", coll="sample_id", it=it, output="dist")
  
  dd_zoneS[[e]] <- apply(temp$divCSIB, 1, boot.summary) %>% 
    bind_rows() %>% 
    mutate(env=e)
}

dd_zoneS <- lapply(dd_zoneS, function(x) cbind(x, mid=stages$mid)) %>% 
  bind_rows()


# Save data ---------------------------------------------------------------
save(dd_all, dd_zone, dd_zoneN, dd_zoneS, file=file.path("output", "div_time_results.RData"))


# Plots -------------------------------------------------------------------
load(file.path("output", "div_time_results.RData"))
p1 <- ggplot(cbind(stages, dd_all), aes(x=mid, y=mean, ymin=lwr, ymax=upr)) +
  geom_vline(xintercept = ints$bottom, col="darkgrey", linetype="dashed") +
  geom_ribbon(fill="grey80")+
  geom_line(lwd=1) +
  scale_x_continuous(trans="reverse", limits=c(66,0)) +
  labs(x="Age (Ma)", y="Corrected sampled-in-bin diversity") +
  theme_hce

p2 <- ggplot(dd_zone, aes(x=mid, y=mean, ymin=lwr, ymax=upr, group=env)) +
  geom_vline(xintercept = ints$bottom, col="darkgrey", linetype="dashed") +
  geom_ribbon(aes(fill=env), alpha=0.3)+
  geom_line(aes(col=env), lwd=1) +
  scale_color_manual(values=pal, labels=env_names) +
  scale_fill_manual(values=pal, labels=env_names)+
  scale_x_continuous(trans="reverse", limits=c(66,0)) +
  labs(x="Age (Ma)", y="Corrected sampled-in-bin diversity", col="Climate Zone", fill="Climate Zone") +
  theme_hce

p3 <- ggplot(dd_zoneN, aes(x=mid, y=mean, ymin=lwr, ymax=upr, group=env)) +
  geom_vline(xintercept = ints$bottom, col="darkgrey", linetype="dashed") +
  geom_ribbon(aes(fill=env), alpha=0.3)+
  geom_line(aes(col=env), lwd=1) +
  scale_color_manual(values=pal, labels=env_names) +
  scale_fill_manual(values=pal, labels=env_names)+
  scale_x_continuous(trans="reverse", limits=c(66,0)) +
  labs(x="Age (Ma)", y="Corrected sampled-in-bin diversity", col="Climate Zone", fill="Climate Zone") +
  theme_hce

p4 <- ggplot(dd_zoneS, aes(x=mid, y=mean, ymin=lwr, ymax=upr, group=env)) +
  geom_vline(xintercept = ints$bottom, col="darkgrey", linetype="dashed") +
  geom_ribbon(aes(fill=env), alpha=0.3)+
  geom_line(aes(col=env), lwd=1) +
  scale_color_manual(values=pal, labels=env_names) +
  scale_fill_manual(values=pal, labels=env_names)+
  scale_x_continuous(trans="reverse", limits=c(66,0)) +
  labs(x="Age (Ma)", y="Corrected sampled-in-bin diversity", fill="Climate Zone", col="Climate Zone") +
  theme_hce +
  theme(legend.position = "bottom")
# Save plot: Fig_01 -------------------------------------------------------

svg(file.path("figs", "Fig_01_CSIB_diversity.svg"), w=10, h=10)
add_scale(p1) +
  add_scale(p2)+
  add_scale(p3)+
  add_scale(p4) +
  plot_annotation(tag_prefix = "(", tag_levels = "a", tag_suffix = ")") +
  plot_layout(guides = "collect") &
  theme(legend.position='bottom')
dev.off()

