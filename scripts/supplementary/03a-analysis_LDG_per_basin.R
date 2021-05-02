# Computes LDG at species level for each plankton group 
# per ocean basin.

# load libraries ---------------------------------------------------------
library(tidyverse)
library(divDyn)
library(magrittr)
library(reshape2)
library(foreach)
library(doParallel)
library(raster)

source(file.path("scripts", "utils","functions.R"))

# Load data ---------------------------------------------------------------
data(stages, package="divDyn")

load(file.path("output", "ocean_basins.RData"))
load(file.path("data", "wrangled_data.RData")) # sp at species level

#dat$sp <- droplevels(dat$sp) #old levels kept
dat[dat$fossil.group == "FALSE",]$fossil.group <- "F"
dat <- dat %>% filter(FAD != LAD)

#merge pacific
basins[["pacific"]] <- raster::union(basins$eastpacific, 
                                    basins$westpacific)

basins <- basins[-grep("pacific", names(basins))[1:2]]

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

oceans <- names(basins)

LDG_div <- list()

for(o in oceans){
  cat(o, "\n")
  #filter occurrences by ocean
  dat.sp <- dat
  coordinates(dat.sp) <- ~paleolong + paleolat
  
  dat.ocean <- sp::over(dat.sp, basins[[o]])
  
  dat.ocean <- dat[names(dat.ocean[!is.na(dat.ocean)]),]
  
  #for (g in names(group_names)[2:5]){
  #  cat(g)
    temp.g <- dat.ocean %>% filter(bin > 80)
    
    tbin <- sort(unique(temp.g$bin))
    div <- list()
    
    temp.t <- list()
    
    for (t in tbin){
      temp <- temp.g %>% filter(bin == t)
      
      temp.res <- foreach(i=1:it, .combine=rbind) %dopar% {
        temp.sqs <- divDyn::subtrialSQS(temp, tax="sp", q=0.7, bin="new_lat")
        
        if(nrow(temp[temp.sqs,]) > 0){
          cbind(temp[temp.sqs,], tr = i)
        } else{
          NA
        }
      }
      
      
      if(length(na.omit(temp.res)) > 0){
        temp.t[[t]] <- temp.res %>% 
           group_by(new_lat) %>% summarise(n=n_distinct(sp)) %>% mutate(bin = t)
      }
    }
    
    LDG_div[[paste(o)]] <- temp.t %>% 
      bind_rows() %>% 
      mutate(ocean=o)
  #}
}
#stop cluster
stopCluster(cl)

LDG_div <- do.call(rbind.data.frame, LDG_div) 
save(LDG_div, file=file.path("output", "LDG_div_ocean.RData"))
