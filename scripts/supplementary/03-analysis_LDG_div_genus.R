# Computes and plots LDG at genus level for each plankton group

# load libraries ---------------------------------------------------------
library(tidyverse)
library(divDyn)
library(magrittr)
library(reshape2)
library(foreach)
library(doParallel)

source(file.path("code", "functions.R"))
source(file.path("code", "plot_parameters.R"))

# Load data ---------------------------------------------------------------
data(stages, package="divDyn")

load(file.path("data", "wrangled_data_genus.RData")) # sp at genus level

#dat$sp <- droplevels(dat$sp) #old levels kept
dat[dat$fossil.group == "FALSE",]$fossil.group <- "F"
dat <- dat %>% filter(FAD != LAD)


##### resolution lat
geog_res <- 10

dat <-   dat %>% dplyr::select(-new_lat, -new_long) %>% 
  bind_cols(
    
    #latitudinal resolution (cellsize) = 1
    assign_grid_points(dat$paleolong, dat$paleolat, cellsize=c(geog_res,geog_res))
  )  %>%
  rename(new_lat = y, new_long = x) %>%
  filter(!is.na(new_lat))

#### 
it = 1000

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)

LDG_div <- list()

for (g in names(group_names)[2:5]){
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
save(LDG_div, file=file.path("output", "LDG_div_genus.RData"))

#### plot data ####
load(file.path("output", "LDG_div_genus.RData"))

data(stages, package="divDyn")

stages$per <- categorize(stages$stage, phases) 
stages <- stages[!is.na(stages$per),]

LDG_div <- LDG_div %>% 
  left_join(int.tbl, by=c("bin" ="int")) %>%
  filter(bin != 95)

div_geog_ev <- LDG_div %>% left_join(stages[,c("mid", "per")]) %>%
  left_join(LDG_div) %>% 
  filter(!is.na(per))

div_geog_ev$per <- factor(div_geog_ev$per, levels=names(phases))
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

p <- ggplot(div_geog_ev, aes(x=new_lat, y=n, group=hemispfile.path)) +
  geom_point(shape=21, aes(fill=group), col="transparent", size=1, alpha=0.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, k=4), se = TRUE, aes(col=group, fill=group) ) +
  facet_grid(group~per, labeller=labeller(group=group_names, per=per_names)) +
  geom_vline(xintercept = 0, linetype="dashed", col="darkgrey")+
  scale_color_manual(values=ldg_clr[c(3,2,2,3)], labels=group_names) +
  scale_fill_manual(values=ldg_clr[c(3,2,2,3)], labels=group_names) +
  labs(x=bquote(bold("Latitude ("*degree*")")), y="Genus richness") +
  geom_text(data=res, aes(x=ifelse(hemispfile.path=="N", 10, -80), y=95, 
                          label=ifelse(is.na(rsq), NA, paste0("italic(R) ^ 2 ==", format(rsq, digits=2)))), size=3, inherit.aes = FALSE, parse=TRUE, hjust = 0,
            col="darkgrey")+
  geom_text(data=res, aes(x=ifelse(hemispfile.path=="N", 10, -80), y=85, 
                          label=ifelse(is.na(p), NA, paste0("italic(p)==", format(p, digits  = 2)))), size=3, inherit.aes = FALSE, parse=TRUE, hjust = 0,
            col="darkgrey") +
  theme_hce +
  theme(legend.position = "none")

ggsave(file.path("figs", "Supplement", "Fig_02_LDG_genus.svg"), p, w=8, h=10)

