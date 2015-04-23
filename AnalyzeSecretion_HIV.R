rm(list=ls())

# Analyze Secretion Template File
library(stringr)
library(foreign)
library(data.table)
reorganizeTable <- function(data, baseName=NA, convertToNumeric=TRUE, nameCol='Measurement', valueCol='Value')
{
     library(plyr)
     idCols <- names(data)
     idCols <- idCols[-which(idCols %in% c(nameCol,valueCol))]
     newData <- data.frame(stringsAsFactors=FALSE)
     measurements <- unique(data[,nameCol])
     #     m <- measurements[1]
     #     browser()
     for(m in measurements)
     {
          if(is.na(baseName))
          {
               newColName <- m
               newColName <- gsub(' ','.',newColName, fixed=TRUE) # Get rid of extraneous spaces
          }else
          {
               newColName <- paste(baseName,'.',m, sep='')
               newColName <- gsub(' ','.',newColName, fixed=TRUE) # Get rid of extraneous spaces
          }

          temp <- data[data[,nameCol]==m,]
          temp2 <- temp[,c(idCols,valueCol)]
          names(temp2)[names(temp2)==valueCol] <- newColName
          if(nrow(newData) == 0)
          {
               newData <- temp2
          }else
          {
               newData <- merge(newData, temp2, by=idCols)
          }
     }

     if(convertToNumeric)
     {
          for(n in idCols)
          {
               newData[,n] <- as.numeric(as.character(newData[,n]))
          }
     }

     return(newData)
}

path <- '/Users/jaywarrick/Desktop/A Sandbox/Test/20140805/File - Results/x0_y0.txt'
path <- '/Users/jaywarrick/Desktop/A Sandbox/Test/Dataset Name/File - Results/x0_y0.txt'

data <- reorganizeTable(read.arff(path), convertToNumeric=FALSE);

# Filter weird data
data <- subset(data, ratio > 0)

data$HIV <- as.character(data$HIV)
HIVLevels <- unique(data$HIV)
data$HIVNumber <- as.numeric(str_extract(as.character(data$HIV), "\\d+\\.*\\d*")) # Just grabs the number portion of the string, leaving out the units
data$logRatio <- 0
data$logRatio <- log(data$ratio)

# # Grab whatever portion of the data you are interested in...
# calData <- subset(data, Cells=='None')
# LNCaPData <- subset(data, Cells=='LNCaP')
# PC3Data <- subset(data, Cells=='PC3')

# Summarize the data
summary <- data.frame()
temp <- data.table(data)
data$signal <- log(data$secretionSig)-log(data$secretionMed)
bg <- subset(data, HIV=='0')
#summary <- temp[,lapply(.SD,mean),by=c('HIV')][,c('HIV','ratio','logRatio'),with=FALSE]
summary <- temp[,list(mean=mean(signal), mad=mad(signal), sd=sd(signal), p=t.test(x=signal, y=bg$signal)$p.value, standard.score=(mean(signal)-mean(bg$signal))/mad(bg$signal), n.microwells=length(signal)), by=c('HIV')]
data$signal <- log(data$secretionSig)-log(data$secretionMed)

# Plot the histograms
plot(c(), c(), xlim=c(0,10), ylim=range(0,15), xlab='HIV Concentration [virus/well]', ylab='Signal per Bead [au]')
for(hiv in unique(data$HIVNumber))
{
     temp <- subset(data, HIVNumber==hiv)
     if(hiv=='0')
     {
          points(temp$HIVNumber + 0.001, temp$ratio, pch=20, col=rgb(0,0,0,0.2))
          points(temp$HIVNumber[1] + 0.001, mean(temp$ratio), pch=20, col='red')
     }
     else
     {
          points(temp$HIVNumber, temp$ratio, pch=20, col=rgb(0,0,0,0.2))
          points(temp$HIVNumber[1], mean(temp$ratio), pch=20, col='red')
     }
}

# Plot the histograms of signal ratios
plot(density(data$ratio), col=rgb(1,1,1,0.5), xlim=c(0,25), ylim=range(0,0.3), xlab='Signal per Bead [au]', ylab='Density', main='')
cols <- seq(0,1,length.out=length(HIVLevels)+1)
i <- 1
for(hiv in HIVLevels)
{
     temp <- subset(data, HIV==hiv)
     lines(density(temp$ratio, adjust=0.75), col=gray(cols[i]), lwd=2)
     print(hiv)
     print(cols[i])
     i = i + 1;
}
legend('topright', legend=HIVLevels, col=gray(cols[1:length(HIVLevels)]), lwd=2)

# Plot the histograms of signals instead of ratios
# plot(density(data$signal), col=rgb(1,1,1,0.5), xlim=c(0,3.5), ylim=range(0,1.5), xlab='Signal per Bead [au]', ylab='Density', main='')
cols <- seq(0,1,length.out=length(HIVLevels)+1)
i <- 1
for(hiv in HIVLevels)
{
     temp <- subset(data, HIV==hiv)
     if(i == 1)
     {
          # If first time through loop, call plot
          plot(density(data$signal), col=rgb(1,1,1,0.5), xlim=c(0,1.1), ylim=range(0,5), xlab='Microwell Signal [au]', ylab='Density', main='')
     }
     lines(density(temp$signal, adjust=0.75), col=gray(cols[i]), lwd=2)
     print(hiv)
     print(cols[i])
     i = i + 1;
}
legend('topright', legend=HIVLevels, col=gray(cols[1:length(HIVLevels)]), lwd=2)

# temp <- subset(data, HIV=='3')
# hist(temp$Ratio, col=gray(cols[3]), lwd=2, breaks=40)


