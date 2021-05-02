# Runs the Pyrate DES model

# Load libraries ----------------------------------------------------------

library(dplyr)
path <- file.path("scripts", "supplementary", "pyrate")

# Fossil data -------------------------------------------------------------
load(file.path("data", "wrangled_data.RData"))
dat[dat$fossil.group == "FALSE",]$fossil.group <- "F"

dat <- dat %>% filter(FAD != LAD)

df <- dat[dat$FAD > 80,] # including Maastrictian

df <- df %>% 
  # group_by(sp, new_long, new_lat) %>% 
  # summarise(earliestAge = max(age), latestAge=min(age), age=mean(age)) %>% 
  rename(decimalLongitude = new_long, decimalLatitude=new_lat, scientificName=sp) %>% 
  select(scientificName, decimalLatitude, decimalLongitude, higherGeography=env1, age)

write.table(df, file=file.path(path, "DES_input_data", "fossil_data_all.txt"), sep="\t", row.name=F)

# Create DES file ---------------------------------------------------------

source(file.path(path, "scripts", "functions_DES.R"))

data(stages, package="divDyn")
timeint <- c(stages$bottom[81:95],0)

#create DESin files
foss <- read.table(file.path(path, "DES_input_data", "fossil_data_all.txt"), sep = "\t", header = T)
foss <- foss[complete.cases(foss),]

DES_input_data <- DESin(foss, timeint = timeint, age="age")
writeDESin(DES_input_data, "ooet_pyrate/pyrate_output/DESin_all")

# Run DES -----------------------------------------------------------------
phases <- list(W1=c("Danian", "Selandian-Thanetian"),
               H = c("Ypresian"),
               W2 = c("Lutetian", "Bartonian", "Priabonian"),
               C = c("Rupelian", "Chattian", "Lower Miocene", 
                     "Middle Miocene", "Upper Miocene", "Pliocene"),
               I = c("Pleistocene", "Holocene")
)

timeper <- stages[82:95,]
timeper$phases <- divDyn::categorize(timeper$stage, phases)

# run des
runDES("pyrate_output/DESin_foram_rep1.txt", 
       qtimes=c(sort(tapply(timeper$bottom, timeper$phases, max), decreasing = T), 0), type="diversity", run=TRUE)


