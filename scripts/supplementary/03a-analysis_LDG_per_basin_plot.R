# Plots LDG per ocean basin

# load libraries ---------------------------------------------------------
library(tidyverse)
library(divDyn)
library(magrittr)
library(reshape2)

source(file.path("scripts","utils", "plot_parameters.R"))

# Load data ---------------------------------------------------------------
data(stages, package="divDyn")

stages$per <- categorize(stages$stage, phases) 
stages <- stages[!is.na(stages$per),]

load(file.path("data", "wrangled_data.RData")) # sp at species level

#dat$sp <- droplevels(dat$sp) #old levels kept
dat[dat$fossil.group == "FALSE",]$fossil.group <- "F"
dat <- dat %>% filter(FAD != LAD)

#### plot data ####
load(file.path("output", "LDG_div_ocean.RData"))

LDG_div <- LDG_div %>% 
  left_join(int.tbl, by=c("bin" ="int")) %>%
  filter(bin != 95)

div_geog_ev <- LDG_div %>% left_join(stages[,c("mid", "per")]) %>%
  left_join(LDG_div) %>% 
  filter(!is.na(per))

div_geog_ev$per <- factor(div_geog_ev$per, levels=names(phases))
div_geog_ev$hemispfile.path <- ifelse(div_geog_ev$new_lat <0, "S", "N")

head(div_geog_ev)

#### getting model summary
library(mgcv)

oceans <- unique(div_geog_ev$ocean)

do_gam <- function(x, data){
  tryCatch({g1 <- gam(n ~ s(x, k=3), data=data)
  cbind.data.frame(rsq=summary.gam(g1)$r.sq, p=summary.gam(g1)$s.pv)
  }, error=function(e) {message(e); cbind.data.frame(rsq=NA, p=NA)}
  )
}


res <- div_geog_ev %>% 
  group_by(per, hemispfile.path, ocean) %>% 
  do(do_gam(.$new_lat, .))

head(res)

p <- ggplot(div_geog_ev, aes(x=new_lat, y=n, group=hemispfile.path)) +
  geom_point(shape=21, fill="black", col="transparent", size=1, alpha=0.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, k=3), se = TRUE, 
              col="black", fill="black") +
  facet_grid(ocean~per, 
             labeller=labeller(ocean=c(atlantic="Atlantic",
                                       indian = "Indian", 
                                       pacific = "Pacific"),
                               per=per_names)) +
  geom_vline(xintercept = 0, linetype="dashed", col="darkgrey")+
  coord_cartesian(xlim=c(-90,90))+
  # scale_color_manual(values=viridis(6)[c(2,4,5)]) +
  # scale_fill_manual(values=viridis(6)[c(2,4,5)]) +
  labs(x=bquote(bold("Latitude ("*degree*")")), y="Species richness") +
  geom_text(data=res, aes(x=ifelse(hemispfile.path=="N", 10, -80), y=400, 
                          label=ifelse(is.na(rsq), NA, paste0("italic(R) ^ 2 ==", format(rsq, digits=2)))), size=3, inherit.aes = FALSE, parse=TRUE, hjust = 0,
            col="darkgrey")+
  geom_text(data=res, aes(x=ifelse(hemispfile.path=="N", 10, -80), y=350, 
                          label=ifelse(is.na(p), NA, paste0("italic(p)==", format(p, digits  = 2)))), size=3, inherit.aes = FALSE, parse=TRUE, hjust = 0,
            col="darkgrey") +
  theme_hce +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5, face=2))

ggsave(file.path("figs", "Supplement", "Fig_02_LDG_ocean.svg"), p, w=8, h=7.5)
