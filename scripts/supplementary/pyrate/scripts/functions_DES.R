#' Create DES input, adapted from the package speciesgeocodeR
#'
#' @param x data file of fossils occurrences containing coordinates and ages.
#' @param recent data file for recent occurrences containing coordinates (optional). 
#' @param taxon column name for taxon, "scientificName" by default.
#' @param area column name for region, "higherGeography" by default.
#' @param age column name for age, if actual ages of fossil used.
#' @param age1 column name for minimum age, "earliestAge" by default
#' @param age2 column name for maximum age, "latestAge" by default
#' @param bin.size bin size for binning data in equal time intervals, in Myr.
#' @param timeint custom time bins if unequal time intervals required.
#' @param reps number of replicates - only works when age1 and age2 are provided
#' @param verbose print output?
#'
#' @return
#' @export
#'
#' @examples
DESin <- function(x, 
                  recent=NULL, 
                  taxon = "scientificName",
                  area = "higherGeography",
                  age=NULL,
                  age1 = "earliestAge",
                  age2 = "latestAge",
                  bin.size = 5, 
                  timeint = NULL,
                  reps = 3, 
                  verbose = FALSE) {
  
  # load fossil data
  if (is.data.frame(x)) {
    dat <- x
  } else {
    dat <- read.table(x, sep = "\t", header = TRUE, row.names = NULL)
  }
  
  # CHECK IF this is still necessary, and why the summary method still uses
  # midpoints
  if(is.null(age)){ # only if age is not provided
    if (!"midpointage" %in% names(dat)) {
      dat$midpointage <- (dat[[age1]] + dat[[age2]])/2
      warning("column midpointage not found, calculating from earliestage and latestage")
    }
  } else {
    dat$midpointage <- dat$age
    reps <- 1
    cat("Actual ages provided, no randomisation will be carried out")
  }
  
  if(!is.null(timeint)) bin.size <- NULL
  
  if(!is.null(recent)){
    # load and prepare recent data
    if (is.data.frame(recent)) {
      rece <- recent
    } else {
      rece <- read.table(recent, header = TRUE, sep = "\t", stringsAsFactors = FALSE, 
                         row.names = NULL)
    }
    
    # Calculate no. of taxon per region
    rece[[area]] <- as.character(rece[[area]])
    regs <- unique(rece[[area]])
    rece[rece[[area]] == regs[1], area] <- 1
    rece[rece[[area]] == regs[2], area] <- 2
    
    rece <- unique(rece)
    rece[[area]] <- as.numeric(rece[[area]])
    rece <- aggregate(rece[[area]] ~ rece[[taxon]], FUN = sum)
    names(rece) <- c(taxon, area)
  } else rece <- NA
  
  # code fossil data
  outp <- list()
  
  for (i in 1:reps) {
    if (verbose) {
      print(sprintf("producing replicate %s of %s", i, reps))
    }
    
    if(is.null(age)){
      # simulate random age between min and max
      dat$age <- sapply(seq(1, nrow(dat)), function(x) stats::runif(1, max = dat[[age1]][x], 
                                                                    min = dat[[age2]][x]))
    } else {
      #assign actual ages
      dat$age <- dat[,age]
    }
    
    #bin data according to bin size or custom time interval
    if (!is.null(bin.size)){
      # define age class cutter and cut ages into timebins
      cutter <- timeint <-  seq(0, ceiling(max(dat$age)), by = bin.size)
    } else {
      cutter <- sort(timeint)
    }
    
    dat$timeint <- as.numeric(as.character(cut(dat$age, breaks = cutter, 
                                               digits = 5, labels = cutter[-length(cutter)]),
                                           right=FALSE
    ))
    dat <- dat[complete.cases(dat),]
    
    # following only works if data is sorted. Make changes for the future.
    dat <- dat[order(dat$age),]
    
    # code the presence in each regions per species
    dat.list <- split(dat, dat[[taxon]])
    
    binned <- lapply(dat.list, function(k) {
      dat.out <- data.frame(timebin = cutter, 
                            area1 = rep(0, length(cutter)), 
                            area2 = rep(0, length(cutter)))
      
      if (length(k[[area]] == sort(unique(dat[[area]])[1])) > 0) {
        dat.out[dat.out$timebin %in% unlist(k[k[[area]] == unique(dat[[area]])[1], "timeint"]), "area1"] <- 1
      }
      if (length(k[[area]] == sort(unique(dat[[area]])[2])) > 0) {
        dat.out[dat.out$timebin %in% unlist(k[k[[area]] == unique(dat[[area]])[2], "timeint"]), "area2"] <- 2
      }
      
      presence <- rowSums(dat.out[, 2:3])
      return(presence)
    })
    
    # set timebins before first appearance to NaN
    out <- lapply(binned, function(k) {
      if (!any(k > 0)) {
        k <- rep("nan", length(k))
        return(as.numeric(k))
      } else {
        if (max(which(k > 0)) < length(k)) {
          k[(max(which(k > 0)) + 1):length(k)] <- "nan"
          return(as.numeric(k))
        } else {
          return(k)
        }
      }
    })
    
    # output format
    out <- do.call("rbind.data.frame", out)
    
    if(is.null(bin.size)){
      bin.size=0
      
    } 
    
    names(out) <- (cutter + bin.size / 2)
    
    out <- rev(out)
    out[[taxon]] <- names(dat.list)
    outp[[i]] <- out[,c(taxon,(cutter + bin.size / 2) )]
  }
  
  if(!is.null(recent)){
    
    # combine recent and fossil data
    outp2 <- lapply(outp, function(k) {
      outo <- merge(k, rece, by = taxon, all.x = TRUE)
      outo[[area]][is.na(outo[[area]])] <- 0
      #rownames(outo) <- outo[, 1]
      #outo <- outo[, -1]
      names(outo)[ncol(outo)] <- 0
      return(outo)
    })
    
  } else outp2 <- outp
  
  # make sure all replicates cover the same time spann, i.e. add additional
  # columns before the first time column
  meas <- sapply(outp2, "ncol")
  
  if (max(meas) != min(meas)) {
    numb <- which(meas < max(meas))
    for (i in numb) {
      dat.int <- outp2[[i]]
      repl <- nrow(dat.int) * (max(meas) - meas[i])  # how many NaNs are needed
      dat.comb <- c(rep(NaN, times = repl), unlist(dat.int[,-1]))
      dat.int <- data.frame(matrix(dat.comb, 
                                   nrow = nrow(dat.int), 
                                   ncol = max(meas)-1, 
                                   byrow = FALSE))
      dat.int <- data.frame(outp[[i]][taxon],
                            dat.int)
      names(dat.int) <- names(outp2[[which(meas == max(meas))[1]]])
      #rownames(dat.int) <- rownames(outp2[[i]])
      outp2[[i]] <- dat.int
    }
  }
  
  # create output object
  outp <- list(input_fossils = dat, 
               input_recent = rece, 
               DES_replicates = outp2,
               bin_size = ifelse(bin.size==0,NA, bin.size),
               timeint = timeint,
               area = area,
               regions = c(`1`=unique(dat$higherGeography)[1], `2`=unique(dat$higherGeography)[2]),
               taxon = taxon)
  
  names(outp) <- c("input_fossils", "input_recent", "DES_replicates", "bin_size", "timeint", 
                   "area", "regions", "taxon")
  class(outp) <- c("DESin", "list")
  return(outp)
}

#' Exports DES input file to be used for PyDES. Adapted from https://github.com/dsilvestro/PyRate
#'
#' @param x DESin object
#' @param file filename (without extension)
#' @param save_metadata (Logical) Should the DESin be saved as an RData file?
#'
#' @return
#' @export
#'
#' @examples
writeDESin <- function(x, file, save_metadata=T) {
  for (i in 1:length(x[["DES_replicates"]])) {
    write.table(x[["DES_replicates"]][[i]], paste(file, "_rep", i, ".txt", 
                                                  sep = ""), na = "NaN", sep = "\t", row.names = F, quote = FALSE)
  }
  
  if(save_metadata){
    save(x, file=sprintf("%s_metadata.RData", file))
  }
}

#' Summary of DESin object. Adapted from https://github.com/dsilvestro/PyRate
#'
#' @param object 
#'
#' @return
#' @export
#'
#' @examples
summary.DESin <- function(object, ...) {
  ares <- split(object[["input_fossils"]], f = object[["input_fossils"]][[object$area]])
  
  outp.nams <- c("Minimum_age", 
                 "Maximum_age", 
                 "Number of records",
                 "Mean record age", 
                 "Number of taxa",
                 "Mean taxon age")
  
  area.1 <- c(round(min(ares[[1]]$midpointage), 1),
              round(max(ares[[1]]$midpointage),1),
              nrow(ares[[1]]), 
              round(mean(ares[[1]]$midpointage), 1), 
              length(unique(ares[[1]][[object$taxon]])),
              round(mean(aggregate(ares[[1]]$midpointage, by = list(ares[[1]][[object$taxon]]),min)$x), 1))
  
  area.2 <- c(round(min(ares[[2]]$midpointage), 1),
              round(max(ares[[2]]$midpointage),  1), 
              nrow(ares[[2]]),
              round(mean(ares[[2]]$midpointage), 1),
              length(unique(ares[[2]][[object$taxon]])),
              round(mean(aggregate(ares[[2]]$midpointage, by = list(ares[[2]][[object$taxon]]), min)$x), 1))
  
  list(Number_of_areas = length(ares), 
       Input_Data = setNames(data.frame(row.names=outp.nams, 
                                        area.1, area.2), object$regions), 
       Number_of_Replicates = length(object[["DES_replicates"]]))
}

# Copied from https://github.com/dsilvestro/PyRate

is.DESin <- function(x) {
  inherits(x, "DESin")
}

#' Plot DESin data. Adapted from https://github.com/dsilvestro/PyRate
#'
#' @param x 
#' @param ribbon 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
plot.DESin <- function(x, ribbon = TRUE, ...) {
  require(ggplot2)
  # species in all areas
  area1 <- lapply(x[["DES_replicates"]], function(k) {
    k[is.na(k)] <- 0
    k[k != 1] <- 0
    colSums(k[,-1])
  })
  area1 <- do.call("rbind.data.frame", area1)
  
  area2 <- lapply(x[["DES_replicates"]], function(k) {
    k[is.na(k)] <- 0
    k[k != 2] <- 0
    k[k == 2] <- 1
    colSums(k[,-1])
  })
  area2 <- do.call("rbind.data.frame", area2)
  
  areaB <- lapply(x[["DES_replicates"]], function(k) {
    k[is.na(k)] <- 0
    k[k != 3] <- 0
    k[k == 3] <- 1
    colSums(k[,-1])
  })
  areaB <- do.call("rbind.data.frame", areaB)
  
  times <- as.numeric(as.character(names(x[["DES_replicates"]][[1]][-1])))
  
  dat.plo <- data.frame(time = rep(times, 3), 
                        mean = c(round(colMeans(area1), 1),
                                 round(colMeans(area2), 1), 
                                 round(colMeans(areaB), 1)), 
                        lwr = c(do.call(pmin, data.frame(t(area1))), 
                                do.call(pmin, data.frame(t(area2))),
                                do.call(pmin, data.frame(t(areaB)))), 
                        upr = c(do.call(pmax, data.frame(t(area1))), 
                                do.call(pmax, data.frame(t(area2))),
                                do.call(pmax, data.frame(t(areaB)))), 
                        area = c(rep(x$regions[1], length(times)), 
                                 rep(x$regions[2], length(times)), 
                                 rep("both", length(times))))
  
  plo <- ggplot() + 
    geom_line(data = dat.plo, 
              aes_string(x = "time", y = "mean",
                         group = "area", col = " area")) +
    scale_x_reverse() + 
    xlab("Time") + 
    ylab("Species") + 
    theme_bw() + 
    theme(legend.title = element_blank())
  
  if (ribbon) {
    plo <- plo + 
      geom_ribbon(data = dat.plo, 
                  aes_string(x = "time", 
                             ymax = "upr",
                             ymin = "lwr", 
                             group = "area",
                             fill = "area"), 
                  alpha = 1/5)
  }
  plo
}


#' Basic function to run PyRateDES2
#' Need to update to add looping through file, maybe parse DESin file
#' @param d file name of DESin file
#' @param qtimes shift time (Q)
#' @param type type of analysis, currently just Time-Dependent for both
#' @param run 
#'
#' @return
#' @export
#'
#' @examples
runDES <- function(d,
                   qtimes = NULL, 
                   type = "time",
                   algo=2,
                   run=FALSE){
  
  if(!is.null(qtimes)) qtimes <- paste("-q", paste(qtimes, collapse=" "))
  
 
    type <- switch(type, 
                   time="-TdD -TdE",
                   diversity="-DivdD -DivdE")
  
  args <- c("python", "PyRateDES2.py -d", d, "-A", algo, qtimes, type)
  
  cmdargs <-paste(args, collapse = " ")
  
  if(run){
    system(cmdargs)
  } else {
    return(cmdargs)
  }
}

#' Calculate mean and confidence interval of marginal rates from PyDES output (currently only works for equal bin sizes)
#'
#' @param file filename of PyDES output (usually ends with marginal_rates.log)
#' @param burnin Burning (no of rows to exclude)
#' @param raw (Logical) returns raw values
#' @param time_intervals custom time interval for which to calculate rates (not sure this works yet)
#'
#' @return
#' @export
#'
DESout <- function(file, burnin=0, raw=FALSE, time_intervals=NULL){
  # message if burnin = 0
  if(burnin == 0) cat("\rBurnin was set to 0. Specyfying a higher burnin, e.g. 100 will exclude the first 100 samples")
  
  # read file and exclude header and burnin rows
  rtt <- read.table(file, skip=1+burnin)
  rtt <- rtt[,-1] # remove iterations column
  hds <- strsplit(readLines(file)[1], "\t")[[1]][-1] #headers
  rtt <- rtt[,1:length(hds)]
  colnames(rtt) <- hds
  
  #extract time intervals from ages
  time <- as.numeric(gsub(".+_(.+)$", "\\1", hds))
  
  if(!raw){
    hpd <- matrix(0, nrow=2, ncol=length(time))
    
    colnames(hpd) <- colnames(rtt)
    
    #calculate posterior rates
    for(i in 1:ncol(rtt)){
      par <- rtt[,i]
      hpd[,i] <- calcHPD(par, 0.95) #python function from utilities
    }
    
    rate_mean = apply(rtt, 2, mean) #mean of rates
    
    # saving the results, imported from python, need to find a better way
    d12_index = grep("d12", colnames(hpd))
    d21_index = grep("d21", colnames(hpd))
    e1_index = grep("e1", colnames(hpd))
    e2_index = grep("e2", colnames(hpd))
    d12_mean = rate_mean[d12_index]
    d21_mean = rate_mean[d21_index]
    e1_mean = rate_mean[e1_index]
    e2_mean = rate_mean[e2_index]
    
    d12_hpd = hpd[,d12_index]
    d21_hpd = hpd[,d21_index]
    e1_hpd = hpd[,e1_index]
    e2_hpd = hpd[,e2_index]
    
    # # caluclate time intervals again for results output
    # ints <- data.frame(bottom=time_intervals[-length(time_intervals)],
    #                    top=time_intervals[2:length(time_intervals)])
    # ints$mid <- (ints$top + ints$bottom)/2
    # 
    mr <- cbind.data.frame(time=unique(time),
                           d12_mean, d21_mean, e1_mean, e2_mean,
                           d12_lwr=d12_hpd[1,], d21_lwr=d21_hpd[1,], 
                           e1_lwr=e1_hpd[1,], e2_lwr=e2_hpd[1,],
                           d12_upr=d12_hpd[2,], d21_upr=d21_hpd[2,], 
                           e1_upr=e1_hpd[2,], e2_upr=e2_hpd[2,])
    rownames(mr) <- NULL
  } else {
    d12_index = grep("d12", colnames(rtt))
    d21_index = grep("d21", colnames(rtt))
    e1_index = grep("e1", colnames(rtt))
    e2_index = grep("e2", colnames(rtt))
    
    mr <- list(d12 = rtt[,d12_index],
               d21 = rtt[,d21_index],
               e1 = rtt[,e1_index],
               e2 = rtt[,e2_index])
    
    mr <- lapply(mr, function(x) setNames(x, unique(time)))
  }
  return(mr)
}

#' Calculate Posterior rates- Translated from python code from: https://github.com/dsilvestro/PyRate
#'
#' @param data  vector of values
#' @param level significance level
#'
#' @return
#' @export
#'
#' @examples
calcHPD <- function(data, level=0.95){
  
  if(level <0 | level > 1){
    stop("Level should be between 0 and 1")
  }
  
  d <- sort(data)
  nData = length(data)
  
  nIn <- floor(level * nData)
  
  if (nIn < 2) stop ("not enough data")	
  
  i <- 1
  r <- d[i+nIn-1] - d[i]
  
  for(k in sort(1:abs(length(data)- nIn-1))){
    rk <- d[k+nIn-1] - d[k]
  }
  
  if (rk < r) r <- rk
  i = k
  
  # if( 0 <= i & i <= (i+nIn-1) & (i+nIn-1) < length(d)	) stop()
  
  return (c(d[i], d[i+nIn-1]))
}

