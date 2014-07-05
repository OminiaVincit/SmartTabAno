# Project: SmartTab Anonymization
# Author: Tran Quoc Hoan
# Start date: 2013-Nov-14

# File: ComplexityEntropy.R
# Do: get complexity : ApEntropy and SaEntropy 
# Last edited on: 2014-January-01

library(ggplot2)
library(tseriesChaos)
library(scatterplot3d)
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
  
  if (A*B == 0) sampEnResult <- 1000
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
  if (lag == 0) return (c(0,0))
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
  
  if (A*B==0) sampEnResult <- 1000 # NA value
  else sampEnResult <- -log(A/B)/log(2)
  
  appEnResult <- sum(c[1:(Num-len*lag)])/(Num-len*lag)
  
  return (c(appEnResult,sampEnResult))
}

#################################################################################################
# Fast algorithm of ApEn & SaEn
#################################################################################################


repmat <- function(X,m,n) {
  # R equivalent of repmat (matlab)
  dim(X) <- c(length(X),1)
  mx = nrow(X)
  nx = ncol(X)
  matrix(X,mx*m,nx*n)
}

ApEnSaFast <- function(data, N, dim, lag, tolerance) {
  if (lag == 0) return (c(0,0))
  apResult <- c(0,0)
  saResult <- c(0,0)
  for (j in 1:2) {
    m <- dim+j-1
    phi <- c(1:(N-(m-1)*lag))
    teta <- c(1:(N-dim*lag))
    dataMat <- mat.or.vec(m,(N-(m-1)*lag))
    
    # setting up data matrix
    for (i in 1:m) {
      dataMat[i,] <- data[(i-1)*lag + (1:(N-(m-1)*lag))]
    }
    
    # counting similar patterns using distance calculation
    for (i in 1:(N-(m-1)*lag)) {
      tempMat <- abs(dataMat - repmat(dataMat[,i],1,N-(m-1)*lag))
      boolMat <- colSums(tempMat > tolerance)
      phi[i] <- sum(boolMat<1)/(N-(m-1)*lag)
      if (i <= N-dim*lag) teta[i] <- sum(boolMat<1)-1
    }
    
    # summing over the counts
    apResult[j] <- sum(log(phi))/(N-(m-1)*lag)
    saResult[j] <- sum(teta)
    
  }
  apen <- (apResult[1]-apResult[2])/log(2)
  saen <- log(saResult[1]/saResult[2])/log(2)
  
  return (c(apen,saen))
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

laptop <- c(6,57,72,78,111,116,121,66,112,125)
refri <- c(21, 87, 103)
#knownDv <- c(3, 6, 7, 9, 21, 24, 25, 28, 47, 57, 72, 74, 87, 89, 92, 93, 95, 103, 109, 111, 115, 122, 123, 131)
set.seed(983423)
N <- 25000
ranData <- rnorm(N)
#raEn1 <- ApEnSaFast(ranData, N, 2, 1, 0.1)
raEn2 <- ApEnSaFast(ranData, N, 2, 1, 0.2)
print(raEn2)
#raEn3 <- ApEnSaFast(ranData, N, 2, 1, 0.3)
#raEn4 <- ApEnSaFast(ranData, N, 2, 1, 0.4)
#print(c(raEn1, raEn2, raEn3, raEn4))

traceLap <- c("dev_5AE2CA.csv","dev_B1B603.csv","dev_B7BEA1.csv","dev_B8148C.csv", "dev_D320C8.csv")
traceRef <- c("dev_76C07F2.csv","dev_98C08A.csv","dev_98C08A2.csv","dev_599393.csv","dev_B7BE29.csv","dev_B7E6F4.csv","dev_B83B9E.csv","dev_D325D9.csv","dev_D325D9B.csv","dev_D331DA.csv", "dev_D331DA2.csv", "dev_D32131A.csv", "dev_D32131A.csv", "dev_D32131B.csv", "dev_D32131C.csv", "dev_undef1.csv", "dev_undef2.csv")
traceDesk <- c("dev_11F01E.csv","dev_59AC8E.csv","dev_59AC89.csv","dev_59ADA7.csv","dev_7296D7.csv","dev_B7E6FA.csv","dev_D35C05A.csv","dev_D35C05B.csv","dev_D337C9.csv","dev_D337C9B.csv","dev_D337C9C.csv","dev_schalli.csv", "dev_schalli2.csv", "dev_denis.csv")
traceTV <- c("dev_330A3.csv", "dev_B80E51.csv", "dev_B80E51B.csv", "dev_B8121D.csv", "dev_B81116.csv", "dev_B83416.csv", "dev_C3E6D1.csv", "dev_D35F73.csv", "dev_D369E0.csv", "dev_D33097.csv")
laptopad <- c("dev_D320C8.csv", "dev_B8198B1.csv", "dev_B8198B2.csv", "testData.csv")
for (i in 1:target_length) {
  
  target <- target_names[i]
  
  # Find id name
  splitID <- strsplit(target,"/",fixed=TRUE)
  splitLength <- length(splitID[[1]])
  idName <- splitID[[1]][splitLength-1]
  
  fileName <- paste("../Data/IREF/November_splines/",idName, "-",fromDay,"-",toDay,".csv",sep="")
  #fileName <- paste("../Data/Tracebase/Laptop/", i,sep="")
  # open file and output it plot
  
  ?try
  devData <- try(read.csv(fileName))
  if (inherits(devData, 'try-error')) next
  
  numberDataPoints <- length(devData[,2])

  if (numberDataPoints >= N) {
    # sData <- devData$Value[1:numberDataPoints]
    k <- numberDataPoints%/%N
    sData <- devData[,2][1:N]
    std <- sd(sData)
#    print("out")
#     maxData <- max(sData)
#     minData <- min(sData)
#    meanData <- mean(sData)
     timeLag <- getTimeLag(sData, N%/%10)
#    print(c(i, std, timeLag, numberDataPoints))
    

    for (m in 2:2) {
      for (r in 2:2) {
        result <- ApEnSaFast(sData, N, m, timeLag, r*std/10)
#         if (r == 1) rate <- result/raEn1
#         if (r == 2) rate <- result/raEn2
#         if (r == 3) rate <- result/raEn3
#         if (r == 4) rate <- result/raEn4
        rate <- result/raEn2
        print(c(i, idName, std, timeLag, r/10, result, rate))
      }
    }
    
    # plot two dimensions data
#     timeLag <- 20
#     x <- sData[1:25000]
#     y <- sData[timeLag+(1:25000)]
#     z <- sData[2*timeLag+(1:25000)]
#     #plot (1:25200, sData[1:25200], type="l")
#     plot3d(x,y, z, type="l")
    
#     # filter data
#     minData <- min(sData)
#     sData <- sData[sData > minData]
#     numberNonzeroDataPoints <- length(sData)
#     std <- sd(sData)
 
    #   print predicted time lag for all devices
    # timeLag <- getTimeLag(sData,1000)
    # print(c(i, idName, timeLag, std))
    
    #mtout <- mutual(sData, 50, lag.max=2000, plot=TRUE)
    
    # x <- window(rossler.ts, start=90)
    # xyz <- embedd(sData, m=3, d=20)
    # plot(xyz,type="l",main="")
    # scatterplot3d(xyz,type="l")
    # recurr(sData, m=3, d=2, start.time=1, end.time=500)
    #stplot(sData, m=3, d=8, idt=1, mdt=250)
    #fn.out <- false.nearest(sData,m=10,d=timeLag,t=180,rt=10)
    #plot(fn.out)
    
#     m <- 2
#     for (k in 1:4) {
#       saEn <- ApEnSaFast(sData, N, m, timeLag, k*std/10)
#       print(c(i, idName, timeLag, std, k/10, saEn))
#     }
    
    
#     r <- std*0.2
#     for (lag in 1:(N%/%10)) {
#               saEn <- ApEnSaFast(sData, N, 2, lag, r)
#               print(c(lag,saEn,timeLag,N))
#     }
    
    # minData <- min(sData)
    # filter data
    # sData <- sData[sData > minData]
    
    # numberDataPointsNonzero <- length(sData)
    # print(c(i,idName,numberDataPoints, numberDataPointsNonzero))
#     if (numberDataPointsNonzero >= N) {
#       sData <- sData[1:N]
#       timeLag <- getTimeLag(sData,N%/%10)
#       std <- sd(sData)
#       ApSa1 <- AppSamEn(sData, N, 2, timeLag, std*0.01)
#       ApSa2 <- AppSamEn(sData, N, 2, timeLag, std*0.1)
#       ApSa3 <- AppSamEn(sData, N, 2, timeLag, std*0.2)
#       ApSa4 <- AppSamEn(sData, N, 2, timeLag, std*0.3)
#       ApSa5 <- AppSamEn(sData, N, 2, timeLag, std*0.4)
#       print(c(i,idName, N, std, timeLag, ApSa1, ApSa2, ApSa3, ApSa4, ApSa5))
#       #print(c(i,idName,timeLag, numberDataPoints, numberDataPointsNonzero))
#     }
    
    #minData <- min(sData)
    #maxData <- max(sData)
    #meanData <- mean(sData)

    # print(c(i,idName,minData, numberDataPoints, numberDataPointsNonzero))
#     sData2 <- devData$Value[2:N]-devData$Value[1:(N-1)]
#     std2 <- sd(sData2)
#     rmax2 <- 0
#     rmax3 <- 0
#     rmax4 <- 0
#     rmax5 <- 0
#     rmax6 <- 0
#     rmax7 <- 0
#     
#     if (std2 > 0) {
#       rmax2 <- (-0.036+0.26*sqrt(std/std2))/((numberDataPoints/1000)^(1/4))
#       rmax3 <- (-0.08+0.46*sqrt(std/std2))/((numberDataPoints/1000)^(1/4))
#       rmax4 <- (-0.12+0.62*sqrt(std/std2))/((numberDataPoints/1000)^(1/4))
#       rmax5 <- (-0.16+0.78*sqrt(std/std2))/((numberDataPoints/1000)^(1/4))
#       rmax6 <- (-0.19+0.91*sqrt(std/std2))/((numberDataPoints/1000)^(1/4))
#       rmax7 <- (-0.2+1.0*sqrt(std/std2))/((numberDataPoints/1000)^(1/4))
#     }
    
    
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
#    r <- std*0.2
#     lag <- getTimeLag(sData,numberDataPoints%/%5)
#     saEn <- AppSamEn(sData, numberDataPoints, 2, lag, r)
#     print(c(i, idName, lag,saEn))
#     
#     for (lag in 1:499) {
#          saEn <- AppSamEn(sData, numberDataPoints, 2, lag, r)
#          print(c(lag,saEn))
#     }
    
    # print to csv file
#     resultFile <- paste("../Data/IREF/November_result/",idName, "-",fromDay,"-",toDay,"-",N,".csv",sep="")
#     timeLag <- getTimeLag(sData,N%/%10) 
# 
#     sink(file=resultFile, type="output")
#     print(c(i,idName))
#         # test for tolerance
#         for (m in 2:9) {
#           for (k in 1:100) {
#             r <- k*std/100
#             ApSa <- ApEnSaFast(sData, N, m, timeLag, r)
#             print(c(m, timeLag, k/100, ApSa))
#           }
#         }
#     sink()

    
#     # embedding dimension variability
#     for (m in 2:50) {
#       r <- 0.2*std
#       ApSa <- AppSamEn(sData, numberDataPoints, m, timeLag, r)
#       print(c(m, ApSa))
#     }
    
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