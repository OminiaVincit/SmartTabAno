# Project: SmartTab Anonymization
# Author: Tran Quoc Hoan
# Start date: 2013-Nov-14

# File: SplinesComplement.R
# Do: complements data using splines 
# Last edited on: 2014-January-24

library(akima)

for (i in 1:target_length) {

  target <- target_names[i]
  
  # Find id name
  splitID <- strsplit(target,"/",fixed=TRUE)
  splitLength <- length(splitID[[1]])
  idName <- splitID[[1]][splitLength-1]
  
  # read from input
  inputFile <- paste("../Data/IREF/November/",idName, "-",fromDay,"-",toDay,".csv",sep="")
  ?try
  devData <- try(read.csv(inputFile))
  
  if (inherits(devData, 'try-error')) next
  
  dt <- devData[,1]
  y <- devData[,2]
  Min <- min(y)
  Max <- max(y)
  
  len <- length(dt)
  
  if ((Min < Max) && (len > 1000)) {
    # convert x into numeric data
    dt <- strptime(dt, format="%Y-%m-%dT%H:%M") 
  
    begin <- strptime("2013-10-28 00:00:00", format="%Y-%m-%d %H:%M")
    end <- strptime("2013-12-01 23:58:00", format="%Y-%m-%d %H:%M")
  
    nBegin <- as.numeric(begin)
    nEnd <- as.numeric(end)
    indexEnd <- (nEnd-nBegin)%/%120+1
    x <- (as.numeric(dt)-nBegin)%/%120+1
    xnew <- seq(1,indexEnd,1)
    ynew <- aspline(x,y,xnew)
  
    outDt <- NULL
    outDt$Date <- ynew$x
    outDt$Value <- ynew$y
  
    # write complements file to ouput
    print(c(i,idName))
    outputFile <- paste("../Data/IREF/November_splines/",idName, "-",fromDay,"-",toDay,".csv",sep="")
    write.csv(outDt, outputFile,row.names=FALSE)
  }
}

#dev.off()