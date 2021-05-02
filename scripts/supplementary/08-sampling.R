# Determine gappiness (adapted from Kiessling et al. 2010)

# Load functions and data -------------------------------------------------

# functions
source(file.path("code", "functions.R"))

# data
load(file.path("data", "wrangled_data.RData")) # sp at genus level

# Gapiness ----------------------------------------------------------------
# * Raw -------------------------------------------------------------------

raw.prop <- gappiness(dat, "env1", c("t", "nt"))

# * Subsampled ------------------------------------------------------------

it <- 1000
chunk <- 20
n <- it/chunk

pb <- txtProgressBar(min = 0, max = chunk, style = 3)
finalSamp <- list()

#do in chunks so as not to overload
for(i in 1:chunk){
   #create replicates of subsampled data based on number of iterations
   setTxtProgressBar(pb, i)
   
#subsampling
sub.df <- replicate(n, bootstr(dat,"sqs"), simplify = FALSE)

finalSamp[[i]] <- lapply(sub.df, function(x){
   df <- dat[x,]
   
   gappiness(df, "env1", c("t", "nt"))
})

}

finalSamp <- unlist(finalSamp)

pps <- cbind(t=finalSamp[grep("prob\\.t", names(finalSamp))],
             nt=finalSamp[grep("prob\\.nt", names(finalSamp))])
wilcox.test(pps[,1],pps[,2])
pps <- apply(pps, 2, mean)

error <- cbind(t=finalSamp[grep("error\\.t", names(finalSamp))],
               nt=finalSamp[grep("error\\.nt", names(finalSamp))])
error <- apply(error, 2, mean)


# Combine and save results ------------------------------------------------
res <- cbind(matrix(paste(round(raw.prop$prob, 3), "±", round(raw.prop$error, 3)), ncol=1),
matrix(paste(round(pps, 3), "±", round(error, 3)), ncol=1))

colnames(res) <- c("Raw Data", "Subsampled Data")
rownames(res) <- c("Tropics", "Extratropics")

write.csv(res, file.path("output", "sampling_probs.csv"))
   