# Project: SmartTab Anonymization
# Author: Tran Quoc Hoan
# Start date: 2013-Nov-13

# File: GetPattern.R
# Do: get pattern from time series data
# Last edited on: 2013-Nov-13

# Read data first
# rawData <- read.csv("../Data/selab/000D6F0001A44A42_20131012.csv")
# rawData <- read.csv("../Data/selab_compl/000D6F00039D9CDF_20131018_compl.csv")
# rawData <- read.csv("../Data/selab/dk_1105_1112.txt")
sData <- rawData$value
sampleTimeSeriesData <- ts(sData)
plot.ts(sampleTimeSeriesData, type="l",lwd=1,col="blue", col.main="red", main="Original Time Series Data")

# Find local-minimum points: X[t0] <= X[t0+-1]
# Find local-maximum points: X[t0] >= X[t0+-1]

# Store in vector (t,v,p) with
# t = time, v = value, p = min(-1) max(+1) both(0)
numberDataPoints <- length(sData)
localPoints <- list()
numberLocalPoints <- 0
epsilon <- 10
-epsilon
  
# get stable points
stablePoints <- list()
stablePointsTwo <- list()
numberStablePoints <- 0

# Get local points
for (i in 1:numberDataPoints) {
  # mode for local-min, local-max 
  localMode <- 10
  
  # get backward and forward difference
  backwardDiff1 <- 0
  forwardDiff1 <- 0
  
  if (i==1) backwardDiff1 <- 0
  else backwardDiff1 <- (sData[i] - sData[i-1])
  
  if (i == numberDataPoints) forwardDiff1 <- 0
  else forwardDiff1 <- (sData[i+1] - sData[i])
  
  if ( (backwardDiff1 <= 0) && (forwardDiff1 >= 0) ) localMode <- -1
  if ( (backwardDiff1 >= 0) && (forwardDiff1 <= 0) ) localMode <- 1
  if ( (backwardDiff1 == 0)&&(forwardDiff1 == 0) ) localMode <- 0
  
  if (localMode < 10) {
    # write to local array
    numberLocalPoints <- numberLocalPoints+1
    localPoints[[numberLocalPoints]] <- c(i, sData[i], localMode)
  }
  
  
}

# put localPoints (time,value) into new data set
localData <- NULL
localData$time <- 1:numberLocalPoints
localData$value <- 1:numberLocalPoints
for (i in 1:numberLocalPoints) {
  localData$time[i] <- localPoints[[c(i,1)]]
  localData$value[i] <- localPoints[[c(i,2)]]
}

# create a new set stable points
numberStablePoints <- 0
# Get local points
for (i in 1:numberLocalPoints) {
  # mode for local-min, local-max 
  localMode <- 10
  
  # get backward and forward difference
  if (i == 1) backwardDiff1 <- 0 
  else backwardDiff1 <- (localPoints[[c(i,2)]] - localPoints[[c(i-1,2)]])
  
  if (i == numberLocalPoints) forwardDiff1 <- 0
  else forwardDiff1 <- (localPoints[[c(i+1,2)]] - localPoints[[c(i,2)]])
  
  if (i <= 2) backwardDiff2 <- backwardDiff1 
  else backwardDiff2 <- (localPoints[[c(i,2)]] - localPoints[[c(i-2,2)]])
  
  if (i >= numberLocalPoints-1) forwardDiff2 <- forwardDiff1
  else forwardDiff2 <- (localPoints[[c(i+2,2)]] - localPoints[[c(i,2)]])
  
  sumBackwardDiff <- backwardDiff1 + backwardDiff2
  sumForwardDiff <- forwardDiff1 + forwardDiff2
  
#   if ( (backwardDiff1 <= -epsilon) && (forwardDiff1 >= 0) ) localMode <- -1
#   if ( (backwardDiff1 <= 0) && (forwardDiff1 >= epsilon) ) localMode <- -1
#   if ( (backwardDiff1 >= epsilon) && (forwardDiff1 <= 0) ) localMode <- 1
#   if ( (backwardDiff1 >= 0) && (forwardDiff1 <= -epsilon) ) localMode <- 1
  # if ( (backwardDiff1 == 0)&&(forwardDiff1 == 0) ) localMode <- 0
  
  if ( (backwardDiff1 <= -epsilon) || (backwardDiff1 >= epsilon) ) localMode <- -1
  if ( (forwardDiff1 >= epsilon) || (forwardDiff1 <= -epsilon) ) localMode <- 1
  
#     # vector backwardDiff mode
#     backwardDiffMode <- 0
#     if ( (backwardDiff1 >= 0) && (backwardDiff2 >= 0) && (sumBackwardDiff >= 0) ) {
#       backwardDiffMode <- 1
#     }
#     if ( (backwardDiff1 <= 0) && (backwardDiff2 <= 0) && (sumBackwardDiff <= 0) ) {
#       backwardDiffMode <- -1
#     }
#     
#     # vector forwardDiff mode
#     forwardDiffMode <- 0
#     if ( (forwardDiff1 >= 0) && (forwardDiff2 >=0) && (sumForwardDiff >= 0) ) 
#       forwardDiffMode <- 1
#     if ( (forwardDiff1 <= 0) && (forwardDiff2 <= 0) && (sumForwardDiff <= 0) ) 
#       forwardDiffMode <- -1
#     
#     stableMode <- 10
#     
#     # stableMin
#     # if ((backwardDiff2 <= 0) && (forwardDiff2 >= 0)) stableMode <- -1
#     if ((backwardDiffMode == -1) && (forwardDiffMode == 1)) stableMode <- -1
#     
#     # stableMax
#     # if ((backwardDiff2 >= 0) && (forwardDiff2 <= 0)) stableMode <- 1
#     if ((backwardDiffMode == 1) && (forwardDiffMode == -1)) stableMode <- 1
  
    stableMode <- 1
    if ((localMode < 10) && (stableMode < 10)) {
      # write to stablePoints list
      numberStablePoints <- numberStablePoints+1
      stablePoints[[numberStablePoints]] <- c(localPoints[[c(i,1)]], localPoints[[c(i,2)]], stableMode)
    }
  
}

# divide again
epsilon <- 20

# create a new set stable points
numberStablePointsTwo <- 0
# Get local points
for (i in 1:numberStablePoints) {
  # mode for local-min, local-max 
  localMode <- 10
  
  # get backward and forward difference
  if (i == 1) backwardDiff1 <- 0
  else backwardDiff1 <- (stablePoints[[c(i,2)]] - stablePoints[[c(i-1,2)]])
  
  if (i == numberStablePoints) forwardDiff1 <- 0
  else forwardDiff1 <- (stablePoints[[c(i+1,2)]] - stablePoints[[c(i,2)]])
  
  if (i <= 2) backwardDiff2 <- backwardDiff1
  else backwardDiff2 <- (stablePoints[[c(i,2)]] - stablePoints[[c(i-2,2)]])
  
  if (i >= numberStablePoints-1) forwardDiff2 <- forwardDiff1
  else forwardDiff2 <- (stablePoints[[c(i+2,2)]] - stablePoints[[c(i,2)]])
  
  sumBackwardDiff <- backwardDiff1 + backwardDiff2
  sumForwardDiff <- forwardDiff1 + forwardDiff2
  
  if ( (backwardDiff1 <= -epsilon) || (backwardDiff1 >= epsilon) ) localMode <- -1
  if ( (forwardDiff1 >= epsilon) || (forwardDiff1 <= -epsilon) ) localMode <- 1
  # if ( (backwardDiff1 == 0)&&(forwardDiff1 == 0) ) localMode <- 0
  
  
  # vector backwardDiff mode
  backwardDiffMode <- 0
  if ( (backwardDiff1 >= 0) && (backwardDiff2 >= 0) && (sumBackwardDiff > 0) ) {
    backwardDiffMode <- 1
  }
  if ( (backwardDiff1 <= 0) && (backwardDiff2 <= 0) && (sumBackwardDiff < 0) ) {
    backwardDiffMode <- -1
  }
  
  # vector forwardDiff mode
  forwardDiffMode <- 0
  if ( (forwardDiff1 >= 0) && (forwardDiff2 >=0) && (sumForwardDiff > 0) ) 
    forwardDiffMode <- 1
  if ( (forwardDiff1 <= 0) && (forwardDiff2 <= 0) && (sumForwardDiff < 0) ) 
    forwardDiffMode <- -1
  
  stableMode <- 10
  
  # stableMin
  if ((backwardDiff2 <= 0) && (forwardDiff2 >= 0)) stableMode <- -1
  if ((backwardDiffMode == -1) && (forwardDiffMode == 1)) stableMode <- -1
  
  # stableMax
  if ((backwardDiff2 >= 0) && (forwardDiff2 <= 0)) stableMode <- 1
  if ((backwardDiffMode == 1) && (forwardDiffMode == -1)) stableMode <- 1
  
  
  stableMode <- 1
  if ((localMode < 10) && (stableMode < 10)) {
    # write to stablePoints list
    numberStablePointsTwo <- numberStablePointsTwo+1
    stablePointsTwo[[numberStablePointsTwo]] <- c(stablePoints[[c(i,1)]], stablePoints[[c(i,2)]], stableMode)
  }
  
}


# draw new data set
# plot(localData$time, localData$value, type="o", col="blue", col.main="red",main="Local Points Graph")

# put stablePointsTwo (time,value) into new data set
nData <- NULL
nData$time <- 1:numberStablePoints
nData$value <- 1:numberStablePoints
for (i in 1:numberStablePoints) {
  nData$time[i] <- stablePoints[[c(i,1)]]
  nData$value[i] <- stablePoints[[c(i,2)]]
}

# put stablePointsTwo (time,value) into new data set
nDataTwo <- NULL
nDataTwo$time <- 1:numberStablePointsTwo
nDataTwo$value <- 1:numberStablePointsTwo

for (i in 1:numberStablePointsTwo) {
  nDataTwo$time[i] <- stablePointsTwo[[c(i,1)]]
  nDataTwo$value[i] <- stablePointsTwo[[c(i,2)]]
}

# 
# extract convex & concave


# draw new data set
#plot (1, type="n", xlim=c(1,numberDataPoints), ylim=c(0,150), xlab="time", ylab="Value", main="")

lines(localData$time, localData$value, lwd=1, type="p",col="red")
lines(nData$time, nData$value, lwd=3, type="p", col="green")
lines(nDataTwo$time, nDataTwo$value, lwd=3, type="p", col="purple")
#lines(nDataTwo$time, nDataTwo$value, lwd=3, type="l", col="red")
#lines(nDataTwo$time, nDataTwo$value, type="l", col="red")