#merge csv files
multMerge = function(mypath, sep=" "){
  filenames = list.files(path = mypath, full.names = TRUE)
  datalist = lapply(filenames, 
                    function(x){read.csv(file = x,
                                         header = TRUE,
                                         stringsAsFactors = FALSE, 
                                         sep= sep)})
  Reduce(function(x,y) {merge(x, y, all = TRUE)}, datalist)
}


### assign grid points based on origin and cell size -----
assign_grid_points <- function(x,y, origin=c(0,0), cellsize=c(5,5)) {
  xy <- cbind(x,y)
  temp <- t(apply(xy, 1, function(z) cellsize/2+origin+cellsize*(floor((z - origin)/cellsize))))
  return(as.data.frame(temp))
}

#function for boostrapping
bootstr <- function(df,method,q=NULL){
  
  ss.t<- switch(method,
                sqs=  divDyn::subtrialSQS(df, tax="sp", q = 0.7, bin = "bin", coll="sample_id"),
                cr = divDyn::subtrialCR(df, q, bin = "bin")
  )
}

# Function for random shuffling -------------------------------------------
shuffle <- function(.data, n=1, perm_cols){
  cols_ids <- match(perm_cols, colnames(.data))
  ids <- seq_len(nrow(.data))
  n_ids <- rerun(n, sample(ids))
  
  purrr::map_dfr(n_ids, function(x){
    .data[ids, cols_ids] <- .data[x, cols_ids]
    .data 
  })
  
}

#calculate mean and normal ci of bootstrapped data
boot.summary <- function(x, na.rm=TRUE, conf.level = 0.95) {
  if(na.rm) x <- na.omit(x)
  
  z = (1 - conf.level)/2  
  
  xbar = mean(x)
  error <- qnorm(1-z)*sd(x)
  
  return(
    setNames(data.frame(matrix(c(xbar, xbar-error, xbar+error), nrow=1)), 
             c("mean", "lwr","upr"))
  )
}

boot.summary2 <- function(x, na.rm=TRUE, conf.level = 0.95) {
  if(na.rm) x <- na.omit(x)
  
  z = (1 - conf.level)/2  

  xbar = mean(x)
  quanx <- quantile(x, probs=c(z, 1-z))
  
  return(
    setNames(data.frame(matrix(c(xbar, quanx), nrow=1)), 
                        c("mean", "lwr","upr"))
    )
}

oddsratio <- function(n00, n01, n10, n11, alpha = 0.05){
  #
  #  Compute the odds ratio between two binary variables, x and y,
  #  as defined by the four numbers nij:
  #
  #    n00 = number of cases where x = 0 and y = 0
  #    n01 = number of cases where x = 0 and y = 1
  #    n10 = number of cases where x = 1 and y = 0
  #    n11 = number of cases where x = 1 and y = 1
  #
  OR <- (n00 * n11)/(n01 * n10)
  #OR <- (n00/(n00+n01))/(n10/(n11+n10))
  #
  #  Compute the confidence intervals:
  #
  siglog <- sqrt((1/n00) + (1/n01) + (1/n10) + (1/n11))
  zalph <- qnorm(1 - alpha/2)
  logOR <- log(OR)
  loglo <- logOR - zalph * siglog
  loghi <- logOR + zalph * siglog
  
  pval <- fisher.test(matrix(c(n00, n01, n10, n11), nrow=2))$p.value
  
  oframe <- data.frame(LowerCI = loglo, OR = logOR, UpperCI = loghi, alpha = alpha, pval=pval)
  oframe
}

# Gappiness from Wolfgang Kiessling: Kiessling (2010) Table S2

gappiness <- function(data, categ, x, tax="sp", bin="bin"){
  # simple version just determine pps
  pps <- numeric(2)
  sd.pps <- numeric(2)
  
  for (m in 1:length(x)) {
    temp <- subset(data, eval(parse(text = paste0(categ, '=="', x[m], '"'))))
    
    # new ranges by ordering this dataframe
    temp <- temp[order(temp[,tax], temp[,bin]),]
    
    temp <- temp[,c(tax, bin)]
    FAD <- subset(temp, duplicated(temp[,tax])==F)
    
    temp <- temp[order(temp[,tax], -temp[,bin]),]
    LAD <- subset(temp, duplicated(temp[,tax])==F)
    
    range <- LAD[,bin]-FAD[,bin]
    ranges <- data.frame(FAD[,tax], range)
    colnames(ranges)[1] <- tax
    
    
    temp <- merge(temp, ranges, by.x="sp", by.y="sp")
    temp <- subset(temp, range>1)
    
    ranges <- subset(temp, duplicated(temp[,tax])==F)
    
    # opportunities for sampling (excluding range endpoints)
    s.opp <- sum(ranges$range-1)
    
    # sampled within range
    temp2 <- unique(temp)
    gens <- table(factor(temp2[,tax]))
    s <- nrow(temp2)-2*length(gens)
    
    # sampling probability
    pps[m] <- s/s.opp
    # standard error 
    sd.pps[m] <- sqrt(pps[m]*(1-pps[m])/s.opp)
    
  }
  return(list(prob=setNames(pps, x), error=setNames(sd.pps, x)))
}