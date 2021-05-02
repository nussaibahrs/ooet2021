# This script calculates the latitudinal diversity gradients for marine 
# plankton over the last 66 million years, divided into climate phases. 


# load libraries ---------------------------------------------------------
library(tidyverse)
library(divDyn)
library(magrittr)
library(reshape2)
library(foreach)
library(doParallel)

source(file.path("scripts","utils", "functions.R"))
source(file.path("scripts","utils", "plot_parameters.R"))

# Load data ---------------------------------------------------------------
data(stages, package="divDyn")

load(file.path("data", "wrangled_data.RData")) # sp at genus level

dat[dat$fossil.group == "FALSE",]$fossil.group <- "F"
dat <- dat %>% filter(FAD != LAD)


##### resolution lat
geog_res <- 10

dat <-   dat %>% select(-new_lat, -new_long) %>% 
  bind_cols(
    
    #latitudinal resolution (cellsize) = 1
    assign_grid_points(dat$paleolong, dat$paleolat, cellsize=c(geog_res,geog_res))
  )  %>%
  rename(new_lat = y, new_long = x) %>%
  filter(!is.na(new_lat))


# create groupings --------------------------------------------------------
dat$fossil.group[dat$fossil.group %in% c("F", "N")] <- "C" #calcareous
dat$fossil.group[dat$fossil.group %in% c("D", "R")] <- "S" #siliceous

dat <- dat %>% 
  bind_rows(
    dat %>%  mutate(fossil.group = "all")
  )

#### 
it = 1000

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)

LDG_div <- list()

for (g in unique(dat$fossil.group)){
  temp.g <- dat %>% filter(fossil.group==g & bin > 80)
  
  tbin <- sort(unique(temp.g$bin))
  div <- list()
  
  temp.t <- list()
  
  for (t in tbin){
    temp <- temp.g %>% filter(bin == t)
    
    temp.res <- foreach(i=1:it, .combine=rbind) %dopar% {
      temp.sqs <- divDyn::subtrialSQS(temp, tax="sp", q=0.5, bin="new_lat")
      cbind(temp[temp.sqs,], tr = i)
    }
    

    
    temp.t[[t]] <- temp.res %>% group_by(new_lat) %>% summarise(n=n_distinct(sp)) %>% mutate(bin = t)
  }
  LDG_div[[g]] <- do.call(rbind.data.frame, temp.t) %>% mutate(group = g)
}

#stop cluster
stopCluster(cl)

LDG_div <- do.call(rbind.data.frame, LDG_div) 
save(LDG_div, file=file.path("output", "LDG_div_overall.RData"))

#### plot data ####
load(file.path("output", "LDG_div_overall.RData"))

int.tbl$per<- categorize(int.tbl$stage, 
           phases)

div_geog_ev <- LDG_div %>% 
  left_join(int.tbl, by=c("bin" ="int")) %>%
  filter(bin != 95)

div_geog_ev <- div_geog_ev[!is.na(div_geog_ev$per),]

div_geog_ev$per <- factor(div_geog_ev$per, levels=c("W1", "H", "W2", "C", "I"))

div_geog_ev$hemispfile.path <- ifelse(div_geog_ev$new_lat <0, "S", "N")

#### getting model summary
library(mgcv)

do_gam <- function(x, data){
  tryCatch({g1 <- gam(n ~ s(x, k=4), data=data)
  cbind.data.frame(rsq=summary.gam(g1)$r.sq, p=summary.gam(g1)$s.pv)
  }, error=function(e) {message(e); cbind.data.frame(rsq=NA, p=NA)}
  )
}

res <- div_geog_ev %>% 
  group_by(group, per, hemispfile.path) %>% 
  do(do_gam(.$new_lat, .))

head(res)

theme_hce <- ggthemes::theme_hc(base_size = 14) + theme(panel.grid.major.y = element_blank(), 
                                                        axis.title = element_text(face="bold"), 
                                                        strip.background = element_rect(fill=alpha("lightgrey", 0.2)), 
                                                        legend.title = element_text(face="bold"),
                                                        legend.position = "top")

group_names <- c(all="Overall", C="Calcareous", S="Siliceous")

p <- ggplot(div_geog_ev, aes(x=new_lat, y=n, group=hemispfile.path)) +
  geom_point(shape=21, aes(fill=group), col="transparent", size=1, alpha=0.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, k=4), se = TRUE, aes(col=group, fill=group) ) +
  facet_grid(group~per, labeller=labeller(group=group_names, per=per_names)) +
  geom_vline(xintercept = 0, linetype="dashed", col="darkgrey")+
  scale_color_manual(values= ldg_clr) +
  scale_fill_manual(values= ldg_clr) +
  labs(x=bquote(bold("Latitude ("*degree*")")), y="Species richness") +
  geom_text(data=res, aes(x=ifelse(hemispfile.path=="N", 10, -80), y=550, 
                          label=ifelse(is.na(rsq), NA, paste0("italic(R) ^ 2 ==", format(rsq, digits=1)))), size=3, inherit.aes = FALSE, parse=TRUE, hjust = 0,
            col="darkgrey")+
  geom_text(data=res, aes(x=ifelse(hemispfile.path=="N", 10, -80), y=475, 
                          label=ifelse(is.na(p), NA, paste0("italic(p)==", format(p, digits  = 1)))), size=3, inherit.aes = FALSE, parse=TRUE, hjust = 0,
            col="darkgrey") +
  theme_hce +
  theme(legend.position = "none")

ggsave(file.path("figs", "Fig_02_LDG_overall.svg"), p, w=8, h=6)

