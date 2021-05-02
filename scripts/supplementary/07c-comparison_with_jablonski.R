# Compares rates in Jablonski (2006) to rates in this study

#### load libraries & functions ----
library(tidyverse)
library(magrittr)
library(divDyn)

#load functions
source(file.path("scripts", "utils","functions.R"))

# Jablonski ----
total=c(38,73,52)
nt_ori <- c(25,53,39)
(39+53+25)/(13+21+13)
sum(nt_ori)/(sum(total-nt_ori))


#### load data ----
path = file.path("data", "original")

dat <- multMerge(path, sep="\t") %>%
  janitor::clean_names()

colnames(dat)[1:8] <- c("taxon.ID", "genus", "gq", "species", "sq", "subspecies", "status", "fossil.group")
colnames(dat)[17:19] <- c("lat", "long", "age")
colnames(dat)[23:24] <- c("paleolat", "paleolong")

nrow(dat) #number of species occurrences

#working at genus level
dat$sp <- paste(dat$genus)
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
  filter(!is.na(new_lat))

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

nrow(dat)

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

##### Assigning environment ----

#if using stage binning
int.tbl <- stages
int.tbl$int <- int.tbl$stg

dat <- dat %>%
  mutate(new_lat = abs(new_lat), 
         env1 = ifelse(new_lat < 30, "t", "nt")) 


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

#calculating migration
it <- 1000
chunk <- 20
n <- it/chunk

mig.tot <- list()
for(i in 1:chunk){
sub.df <- replicate(10, bootstr(dat,"sqs"), simplify = FALSE)

#grouping variables for calculation
grouping_vars <- quos("env1", "mig1")

finalProps <- lapply(sub.df, function(x) {
  #calculate proportion imported in each environment
  df[x,] %>% dplyr::select(sp, env1, mig1) %>%
    distinct() %>%
    group_by(env1, mig1) %>% tally() %>%
    setNames(c("env", "mig", "n")) %>%
    group_by(env) %>%
    mutate(prop = n/sum(n)) %>%
    ungroup() %>%
    filter(mig > 0) %>% 
    dplyr::select(-n, -mig) %>% 
    spread(env, prop)
}
)

mig.tot[[i]] <- data.table::rbindlist(finalProps, fill=TRUE) 
}

mig <- data.table::rbindlist(mig.tot, fill=TRUE)  %>% 
  reshape2::melt() %>%  
  mutate(value=as.numeric(value)) %>% 
  group_by(variable) %>% 
  do(boot.summary2(.$value)) %>% 
  ungroup()

mig$mean[2]/mig$mean[1]
