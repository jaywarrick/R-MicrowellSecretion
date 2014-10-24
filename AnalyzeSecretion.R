# Analyze Secretion Template File
library(stringr)
library(foreign)
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

data <- reorganizeTable(read.arff(path), convertToNumeric=FALSE);

data$PSA <- as.character(data$PSA)
PSALevels <- unique(data$PSA)
data$PSANumber <- as.numeric(str_extract(as.character(data$PSA), "\\d+\\.*\\d*")) # Just grabs the number portion of the string, leaving out the units
data$Ratio[data$Ratio > 0] <- log(data$Ratio[data$Ratio > 0])

# Grab whatever portion of the data you are interested in...
calData <- subset(data, Cells=='None')
LNCaPData <- subset(data, Cells=='LNCaP')
PC3Data <- subset(data, Cells=='PC3')

# Summarize the data
summary <- data.frame()
for(psa in PSALevels)
{
     psa <- "0"
     temp <- subset(calData, PSA==psa)
     summary <- rbind(summary, data.frame(PSA=psa, RatioMean=median(temp$Ratio), RatioMad=mad(temp$Ratio)))
}

# Plot the histograms
plot(c(), c(), xlim=c(0.001,100), ylim=range(0,8), xlab='PSA Concentration [ng/mL]', ylab='Secretion per Bead [au]', log='x')
for(psa in PSALevels)
{
     temp <- subset(calData, PSA==psa)
     PSANumber <- as.numeric(str_extract(as.character(psa), "\\d+\\.*\\d*"))
     if(psa=='0')
     {
          points(temp$PSANumber + 0.001, temp$Ratio, pch=20, col=rgb(0,0,0,0.02))
          points(PSANumber + 0.001, mean(temp$Ratio), pch=20, col='red')
     }
     else
     {
          points(temp$PSANumber, temp$Ratio, pch=20, col=rgb(0,0,0,0.02))
          points(PSANumber, mean(temp$Ratio), pch=20, col='red')
     }
}

# Plot the histograms
plot(density(calData$Ratio), col=rgb(1,1,1,0.5), xlim=c(0,5), ylim=range(0,2), xlab='Secretion per Bead [au]', ylab='Density', main='')
cols <- seq(0,1,length.out=length(PSALevels)+1)
i <- 1
for(psa in PSALevels)
{
     temp <- subset(calData, PSA==psa)
     lines(density(temp$Ratio), col=gray(cols[i]), lwd=2)
     print(psa)
     print(cols[i])
     i = i + 1;
}
legend('topright', legend=PSALevels, col=gray(cols[1:4]), lwd=2)

# Plot the cell data histograms
temp <- LNCaPData
cellData <- subset(LNCaPData, cellCount > 0)
emptyData <- subset(LNCaPData, cellCount == 0)
lines(density(cellData$Ratio), col=rgb(1,0,0,1), lwd=2)
lines(density(emptyData$Ratio), col=rgb(0,0,1,1), lwd=2)


weird <- subset(data, Ratio < 0)
rs <- sort(unique(temp$ImRow))
cs <- sort(unique(temp$ImCol))
for(r in rs)
{
     x <- subset(temp, ImRow==r & ImCol==c)
     lines(density(x$Ratio), col=gray(cols[i]), lwd=2)
     print(psa)
     print(cols[i])
     i = i + 1;
}

