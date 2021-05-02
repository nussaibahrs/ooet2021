# This script calculates the diversity dynamics of each plankton group

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

dd_zone <- list()

for(f in names(group_names)[2:5]){
  for (e in c("t", "nt")){
    temp<- divDyn::subsample(subset(dat, env1==e & fossil.group==f), 
                             q=q, tax="sp", bin="bin",  
                             type="sqs", coll="sample_id", 
                             it=it, output="dist")
    
    dd_zone[[paste(f,e)]] <- apply(temp$divCSIB, 1, boot.summary2) %>% 
      bind_rows() %>% 
      mutate(env=e, group=f)
  }
}

dd_zone <- lapply(dd_zone, function(x) cbind(x, mid=stages$mid)) %>% 
  bind_rows()

p1 <- ggplot(dd_zone, aes(x=mid, y=mean, ymin=lwr, ymax=upr, group=env)) +
  geom_vline(xintercept = ints$bottom, col="lightgrey", linetype="dashed") +
  geom_ribbon(aes(fill=env), alpha=0.3)+
  geom_line(aes(col=env), lwd=1) +
  scale_color_manual(values=pal, labels=env_names) +
  scale_fill_manual(values=pal, labels=env_names)+
  scale_x_continuous(trans="reverse", limits=c(66,0)) +
  labs(x="Age (Ma)", y="Corrected sampled-in-bin diversity", fill="Climate Zone", col="Climate Zone") +
  theme_hce +
  facet_wrap(~group, labeller=labeller(group=group_names))

y=-15
for(i in 1:3){
p1 <- p1 +
  annotate("rect", xmin=syscol$bottom[i], xmax=syscol$top[i], 
           ymin=y, ymax=y*2, col="white", fill="grey95") +
  annotate("text", x=syscol$mid[i], y=y+y/2, label=syscol$system[i], 
           hjust=0.5, fontface="bold", col="grey20") 
}

for(i in 1:5){
  p1 <- p1 +
    annotate("rect", xmin=ints$bottom[i], xmax=ints$top[i], ymin=0, ymax=y, 
             fill=clim_clr[c(3,2,3,4,5)][i], col="white") +
    annotate("text", x=ints$mid[i], y=y/2, 
             label=names(phases)[i],  size=3,
             hjust=0.5, fontface="bold", col="white") 
}

ggsave(file.path("figs", "Supplement", "Fig_S_diversity_group.svg"), p1, w=8,h=8)
