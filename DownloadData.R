# Project: SmartTab Anonymization
# Author: Tran Quoc Hoan
# Start date: 2013-Nov-14

# File: DownloadData.R
# Do: download smart-tab data for i-ref-storage.gutp.ic.i.u-tokyo.ac.jp.
#     load auto-correlation graph
# Last edited on: 2013-December-20

###################################################################################
# You should also type the command below
#
# install.packages("ggplot2", dep=T)
# install.packages("xtsExtra", dep=T) or install.packages("xts", dep=T) instead
# 
###################################################################################

#library(xts)
#library(ggplot2)

# to store number of data in each channel (device)
# linearMetric <- NULL
# linearMetric$Id <- NULL
# linearMetric$Num <- NULL
# linearMetric$Mean <- NULL
# linearMetric$Std <- NULL
# linearMetric$Max <- NULL
# linearMetric$Min <- NULL



for (i in 2:target_length) {
  print(i)
  target <- target_names[i]
  
  # Find id name
  splitID <- strsplit(target,"/",fixed=TRUE)
  splitLength <- length(splitID[[1]])
  idName <- splitID[[1]][splitLength-1]
  
#   linearMetric$Id[i] <- idName
#   linearMetric$Num[i] <- 0
#   linearMetric$Mean[i] <- 0
#   linearMetric$Std[i] <- 0
#   linearMetric$Max[i] <- 0
#   linearMetric$Min[i] <- 0
  
  # download file from URL
  targets_link <- paste("http://iref-storage.gutp.ic.i.u-tokyo.ac.jp/fiapapi/get?",
                   "from=", fromstr, "&to=", tostr,
                   "&pointid=", target, "&head=true",
                   sep="")
  
  # before download it, set options(timeout=300)
  fileName <- paste("../Data/IREF/ThreeMonths/",idName, "-",fromDay,"-",toDay,".csv",sep="")
#   
#   # fileName <- paste("testData",idName,".csv",sep="")
  download.file(targets_link,fileName)
#   
  # open file and output it plot
#   ?try
#   rawData <- try(read.csv(fileName))
#   if (inherits(rawData, 'try-error')) next
# 
# #   data.df <- data.frame(index(rawData),coredata(rawData),stringsAsFactors=FALSE)
# #   colnames(data.df) <- c("Index","Date","Value")
#   
#   numberDataPoints <- length(rawData$Value)
# # 
#   if (numberDataPoints > 0) {
# 
#     linearMetric$Num[i] <- numberDataPoints
#     linearMetric$Mean[i] <- mean(rawData$Value)
#     linearMetric$Std[i] <- sd(rawData$Value)
#     linearMetric$Max[i] <- max(rawData$Value)
#     linearMetric$Min[i] <- min(rawData$Value)
#     MAX <- max(rawData$Value)
    
#     if (linearMetric$Min[i] != linearMetric$Max[i]) {
#       # autocorrelation image file
#       imageFileName <- paste("../Data/IREF/November_acf/",idName, "-",from_day,"-",to_day,"acf1000.png",sep="")
#       dev.copy(png,imageFileName,width=1200,height=600)
#       my.ts <- ts(rawData$Value[1:1000])
#       acf(my.ts,lag.max=500, main="Auto correlation function (Data set 2)")
#       dev.off()
#     }
#     # write to image file
#     imageFileName <- paste("../Data/IREF/November_plot/",idName, "-",from_day,"-",to_day,".png",sep="")
#     dev.copy(png,imageFileName,width=1200,height=600)
#     plot (1, type="n", xlim=c(1,numberDataPoints), ylim=c(0,MAX), xlab="time", ylab="Value", main="")
#     lines(rawData$Date, rawData$Value, lwd=1, type="l",col="blue")
#     #plot(x=data.df$Date, y=data.df$Value, type="l",xlab="time",ylab="Value",main="")
#     dev.off()
#   }
}
# 
#outFile <- paste("../Data/IREF/linearMetric","-",from_day,"-",to_day,".csv",sep="")
#write.csv(linearMetric, outFile)