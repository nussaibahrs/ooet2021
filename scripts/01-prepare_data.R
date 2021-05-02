# This script prepared the data files to be used in the analyses from 
# the original download from the Nepture Database. It calculates the 
# first appearance datum and last appearance datum of each species. 
# It also categorises each occurrence into the tropics and extratropics, 
# and calculate when dispersal from one zone to the other may have 
# happened. 

#### load libraries & functions ----
library(tidyverse)
library(magrittr)
library(divDyn)

# load functions
source(file.path("scripts", "utils", "functions.R"))

# Create folder if doesn't exist ----
dir.create(file.path("output"))

#### load data ----
path = file.path("data", "original")

dat <- multMerge(path, sep="\t") %>%
  janitor::clean_names()

colnames(dat)[1:8] <- c("taxon.ID", "genus", "gq", "species", "sq", 
                        "subspecies", "status", "fossil.group")
colnames(dat)[17:19] <- c("lat", "long", "age")
colnames(dat)[23:24] <- c("paleolat", "paleolong")

nrow(dat) #number of species occurrences

#working at species level
dat$sp <- paste(dat$genus, dat$species)
dat$sp <- droplevels(dat$sp)

#### Set latitude and time bins ----

#### set resolution and recalculate coordinates
geog_res <- 1
temp_res <- 5

dat <-   dat %>% 
  bind_cols(
    
    #latitudinal resolution (cellsize) = 1
    assign_grid_points(dat$paleolong, dat$paleolat, cellsize=c(geog_res,geog_res))
  )  %>%
  rename(new_lat = y, new_long = x) %>%
  filter(!is.na(new_lat)) #%>%


#### binning at stage level ####
data(stages)

#create array to store bin numbers
stg <- array(dim=nrow(dat))

for (i in 1:nrow(stages)){
  temp <- which(dat$age <= stages$bottom[i] & dat$age >= stages$top[i])
  stg[temp] <- stages$stg[i]
}

dat$bin <- stg

dat <- dat %>% inner_join(stages %>% 
                            select(stg, bottom, mid, top, stage, short), by=c("bin"="stg")) %>% #add stages age for plotting later 
  filter(!is.na(bin))

#### generate taxa list with FAD and LAD ----
tx.list <- dat %>% group_by(sp, fossil.group) %>%
  summarise(FADage = max(age), FAD = min(bin), 
            LADage = min(age), LAD= max(bin), 
            dur = FADage-LADage) %>%
  filter(FADage <= 83.6) %>%
  filter(FAD != LAD) #remove single interval occurring taxa

dat <- dat %>% filter(sp %in% tx.list$sp) #only keep the ones within the tx.list
dat <- droplevels(dat)

#### add origination and extinction data ----
dat <- dat %>% 
  left_join(tx.list, by=c("sp", "fossil.group")) %>%
  mutate(ori = ifelse(FAD == bin, 1, 0), #add origination 
         ext = ifelse(LAD == bin, 1, 0)) #add extinction

#### set climate zone ----
dat <- dat %>%
  mutate(new_lat = abs(new_lat), 
         env1 = ifelse(new_lat < 30, "t", "nt")
  )  


##### migration

dat$mig1 <- 0

for (i in 1:nrow(tx.list)){
  tx.temp <- tx.list$sp[i]
  
  #get original environment: majority taken if double
  ori_env <- dat %>% filter(sp == tx.temp & bin == FAD) %>% 
    group_by(env1) %>% tally() %>% filter(n == max(n)) %>% pull(env1)
  
  temp <- dat[dat$sp == tx.temp,] %>% as.data.frame() 
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
        
        dat[dat$sp == tx.temp & dat$bin == tbin[t] & dat$env1 != ori_env,]$mig1 <- 1
        ori_env <- tbin.temp$env
        
      } 
      
    }
  }
}

#### save file ####
save(dat, tx.list, int.tbl, geog_res, temp_res, thres_age, file=file.path("data", "wrangled_data.RData"))
