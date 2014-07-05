# Project: SmartTab Anonymization
# Author: Tran Quoc Hoan
# Start date: 2013-Nov-14

# File: ComplexityNonlinear.R
# Do: get nonlinear statistics of data
# Last edited on: 2014-January-11

#######################################################################################
# function to get time lag
#   (1) the first time lag of zero-crossing auto-correlation or first minimum correlation
#   if value in (1) > lagMax then return lagMax
#######################################################################################
getTimeLag <- function(x,lagMax) {
  my.ts = ts(x)
  if (min(x) == max(x)) return (0)
  acfCo <- acf(my.ts,lag.max=lagMax, plot=F)$acf
  result <- which.min(acfCo)-1
  if (acfCo[result+1] < 0) {
    negative <- which(acfCo < 0)
    prev <- negative[1]-1
    after <- negative[1]
    if (abs(acfCo[prev]) < abs(acfCo[after]) ) lag2 <- prev
    else lag2 <- after
    result <- min(result, lag2)
  }  
  return (result)
}

getTimeLagPart <- function(x) {
  len <- length(x)
  epLen <- len%/%5
  x1 <- x[1:epLen]
  x2 <- x[epLen:(2*epLen)]
  x3 <- x[(2*epLen):(3*epLen)]
  x4 <- x[(3*epLen):(4*epLen)]
  x5 <- x[(4*epLen):len]
  
  lag1 <- getTimeLag(x1,500)
  lag2 <- getTimeLag(x2,500)
  lag3 <- getTimeLag(x3,500)
  lag4 <- getTimeLag(x4,500)
  lag5 <- getTimeLag(x5,500)
  
  return (c(lag1,lag2,lag3,lag4,lag5))
}

####################################################################
# distance function using in appEn and samEn
####################################################################
distance <- function(x, len, lag, index1, index2) {
  u <- x[index1+(0:(len-1))*lag] - x[index2+(0:(len-1))*lag]
  Max <- max(abs(u))
  return (Max)
}

###########################################################################
# Calculate Sample Entropy
###########################################################################
sampleEn <- function(x, Num, len, lag, tolerance) {
  if (lag == 0) return (0)
  a <- NULL
  b <- NULL
  
  for (i in 1:(Num-len*lag)) {
    b[i] <- 0
    a[i] <- 0
    for (j in 1:(Num-len*lag)) {
      if ( (j != i) && (distance(x, len, lag, i, j) <= tolerance) ) b[i] <- b[i]+1
      if ( (j != i) && (distance(x, len+1, lag, i, j) <= tolerance) ) a[i] <- a[i]+1
      # Notice b[i] >= a[i]
    }
    b[i] <- b[i]/(Num-(len-1)*lag)
    a[i] <- a[i]/(Num-(len-1)*lag)
  }
  
  A <- sum(a[1:(Num-len*lag)])/(Num-len*lag)
  B <- sum(b[1:(Num-len*lag)])/(Num-len*lag)
  
  if (B == 0) sampEnResult <- 1000
  else sampEnResult <- -log(A/B)/log(2)
  return (sampEnResult)
}

##################################################################################################
# Calculate Approximate Entropy (including self-matched)
##################################################################################################
appEn <- function(x, Num, len, lag, tolerance) {
  if (lag == 0) return (0)
  a <- NULL
  b <- NULL
  c <- NULL
  for (i in 1:(Num-len*lag)) {
    b[i] <- 0
    a[i] <- 0
    for (j in 1:(Num-len*lag)) {
      if ( distance(x, len, lag, i, j) <= tolerance ) b[i] <- b[i]+1
      if ( distance(x, len+1, lag, i, j) <= tolerance ) a[i] <- a[i]+1
    }
    c[i] <- -log(a[i]/b[i])/log(2)
  }
  
  appEnResult <- sum(c[1:(Num-len*lag)])/(Num-len*lag)
  return (appEnResult)
}


##################################################################################################
# Calculate Approximate Entropy (including self-matched) & Sample Entropy (exclusing self-matched)
##################################################################################################
AppSamEn <- function(x, Num, len, lag, tolerance) {
  if (lag == 0) return (0)
  a <- NULL
  b <- NULL
  c <- NULL
  for (i in 1:(Num-len*lag)) {
    b[i] <- 0
    a[i] <- 0
    for (j in 1:(Num-len*lag)) {
      if ( (j != i) && (distance(x, len, lag, i, j) <= tolerance) ) b[i] <- b[i]+1
      if ( (j != i) && (distance(x, len+1, lag, i, j) <= tolerance) ) a[i] <- a[i]+1
    }
    c[i] <- -log((a[i]+1)/(b[i]+1))/log(2)
    b[i] <- b[i]/(Num-(len-1)*lag)
    a[i] <- a[i]/(Num-(len-1)*lag)
  }
  
  A <- sum(a[1:(Num-len*lag)])/(Num-len*lag)
  B <- sum(b[1:(Num-len*lag)])/(Num-len*lag)
  
  if (B==0) sampEnResult <- 1000 # NA value
  else sampEnResult <- -log(A/B)/log(2)
  
  appEnResult <- sum(c[1:(Num-len*lag)])/(Num-len*lag)
  
  return (c(appEnResult,sampEnResult))
}
###########################################################################
# Euclid distance function
###########################################################################
sqrtDist <- function(x, len, lag, index1, index2) {
  u <- x[index1+(0:(len-1))*lag] - x[index2+(0:(len-1))*lag]
  result <- sum(u*u)
  return (result)
}

###########################################################################
# Find the index of nearest neighbor
###########################################################################
indexNN <- function(x, len, lag, index) {
  n <- length(x)
  MIN <- 10000
  minIndex <- index
  for (k in 1:(n-lag*len)) {
    if (k != index) {
      tmp <- sqrtDist(x, len, lag, index, k)
      if (tmp < MIN) {
        MIN <- tmp
        minIndex <- k
      }
    }
  }
  return (minIndex)
}

###########################################################################
#
# Select parameter m (for calculating entropy)
#        determines the length of the sequences
# Method: False Nearest Neighbor
#
###########################################################################
getFNNrate <- function(x, len, lag, tolerance) {
  N <- length(x)
  numFNN <- 0
  std <- sd(x)
  
  for (i in 1:(N-lag*len)) {
    j <- indexNN(x, len, lag, i)
    rd <- sqrtDist(x, len, lag, i, j)
    rd1 <- sqrtDist(x, len+1, lag, i, j)
    diff <- abs(x[i+len*lag]-x[j+len*lag])
    if (diff > rd*tolerance) numFNN <- numFNN+1
  }
  result <- numFNN/(N-lag*len)
  return (result)
}

# next time : from 108
for (i in 1:target_length) {
  
  target <- target_names[i]
  
  # Find id name
  splitID <- strsplit(target,"/",fixed=TRUE)
  splitLength <- length(splitID[[1]])
  idName <- splitID[[1]][splitLength-1]
  
  
  fileName <- paste("../Data/IREF/November/",idName, "-",fromDay,"-",toDay,".csv",sep="")
  # open file and output it plot
  ?try
  devData <- try(read.csv(fileName))
  if (inherits(devData, 'try-error')) next
  
  numberDataPoints <- length(devData$Value)
  
  samplePoints <- 1000
  if (numberDataPoints > samplePoints) {
    sData <- devData$Value[1:samplePoints]
    numberDataPoints <- length(sData)
    minData <- min(sData)
    maxData <- max(sData)
    meanData <- mean(sData)
    std <- sd(sData)
    
    sData2 <- devData$Value[2:samplePoints]-devData$Value[1:(samplePoints-1)]
    std2 <- sd(sData2)
    rmax2 <- 0
    rmax3 <- 0
    rmax4 <- 0
    rmax5 <- 0
    rmax6 <- 0
    rmax7 <- 0
    
    if (std2 > 0) {
      rmax2 <- (-0.036+0.26*sqrt(std/std2))/((numberDataPoints/1000)^(1/4))
      rmax3 <- (-0.08+0.46*sqrt(std/std2))/((numberDataPoints/1000)^(1/4))
      rmax4 <- (-0.12+0.62*sqrt(std/std2))/((numberDataPoints/1000)^(1/4))
      rmax5 <- (-0.16+0.78*sqrt(std/std2))/((numberDataPoints/1000)^(1/4))
      rmax6 <- (-0.19+0.91*sqrt(std/std2))/((numberDataPoints/1000)^(1/4))
      rmax7 <- (-0.2+1.0*sqrt(std/std2))/((numberDataPoints/1000)^(1/4))
    }
    
    rpre1 <- 0.2
    rpre2 <- 0.01
    rpre3 <- 0.25
    rpre4 <- 0.1
    
    m <- 2
    # print linear metric
    # print(c(i, idName, minData, maxData, meanData, std, std2, numberDataPoints))
    
    # prin non-linear metric
        if (minData < maxData) {
          timeLag1 <- getTimeLag(sData,numberDataPoints%/%10) 
          
    #      noiseData <- minData + (maxData-minData)*sample(0:(numberDataPoints-1))/(numberDataPoints-1)
    #      nstd <- sd(noiseData)
          
    #       nsData2 <- noiseData[2:samplePoints]-noiseData[1:(samplePoints-1)]
    #       nstd2 <- sd(nsData2)
    #       nrmax2 <- (-0.036+0.26*sqrt(nstd/nstd2))/((numberDataPoints/1000)^(1/4))
          
    #      timeLag2 <- getTimeLag(noiseData,100)
                                      
          ApSaMax <- AppSamEn(sData, numberDataPoints, m, timeLag1, rmax2*std)
          ApSa1 <- AppSamEn(sData, numberDataPoints, m, timeLag1, rpre1*std)
          ApSa2 <- AppSamEn(sData, numberDataPoints, m, timeLag1, rpre2*std)
          ApSa3 <- AppSamEn(sData, numberDataPoints, m, timeLag1, rpre3*std)
          ApSa4 <- AppSamEn(sData, numberDataPoints, m, timeLag1, rpre4*std)
          
          print(c(i, idName, timeLag1, rmax2, ApSaMax, rpre1, ApSa1, rpre2, ApSa2, rpre3, ApSa3, rpre4, ApSa4))
          
        }
    
    #     # print predicted time lag for all devices
    #     timeLag <- getTimeLag(sData,numberDataPoints%/%2-1)
    #     print(c(i,idName,numberDataPoints, timeLag))
    
    #   if (minData < maxData) {
    
    #      tsData <- ts(sData)    
    #       acfFileName <- paste("../Data/IREF/November_acf2/",idName, "-",fromDay,"-",toDay,".png",sep="")
    #       dev.copy(png,acfFileName,width=1200,height=600)
    #       acf(tsData,lag.max=numberDataPoints/2)
    #       dev.off()
    
    #       timeLag <- getTimeLag(sData,100)
    #       saEn <- sampleEn(sData, numberDataPoints, 2, timeLag, r)
    #       print(c(timeLag,saEn))
    
    #    }
    
    #test for time lag
    #     r <- std*0.2
    #     for (lag in 1:499) {
    #          saEn <- AppSamEn(sData, numberDataPoints, 2, lag, r)
    #          print(c(lag,std,saEn))
    #     }
    
#     timeLag <- getTimeLag(sData,numberDataPoints%/%10) 
#     
#     print(c(i, idName, numberDataPoints, minData, maxData, meanData))
#     print (c(std, std2, timeLag))
#     print(c("rmax:", rmax2, rmax3, rmax4, rmax5, rmax6, rmax7))
#     print(c("emDim", "Tolerance", "AppEn", "SampleEn"))
#     # test for tolerance
#     for (m in 3:7) {
#       for (k in 1:100) {
#         r <- k*std/100
#         ApSa <- AppSamEn(sData, numberDataPoints, m, timeLag, r)
#         print(c(m, k/100, ApSa))
#       }
#     }
    
    
    #         # embedding dimension variability
    #         for (m in 2:20) {
    #           r <- 0.2*std
    #           ApSa <- AppSamEn(sData, numberDataPoints, m, timeLag, r)
    #           print(c(m, ApSa))
    #         }
    
    # test for embedding dimension    
    #     for (tol in 1:10) {
    #       for (dim in 1:50) {
    #         rate <- getFNNrate(sData, dim, timeLag, 5*tol)
    #         if (rate < 0.005) break
    #         print(c(5*tol,dim,timeLag,rate))
    #       }
    #     }
    
  }
  
}

#dev.off()