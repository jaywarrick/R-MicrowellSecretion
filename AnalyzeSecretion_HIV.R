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
path <- '/Users/jaywarrick/Desktop/A Sandbox/Test/20140923 - HIV Test/File - Results/x0_y0.txt'

data <- reorganizeTable(read.arff(path), convertToNumeric=FALSE);

data$HIV <- as.character(data$HIV)
HIVLevels <- unique(data$HIV)
data$HIVNumber <- as.numeric(str_extract(as.character(data$HIV), "\\d+\\.*\\d*")) # Just grabs the number portion of the string, leaving out the units
data$Ratio[data$Ratio > 0] <- log(data$Ratio[data$Ratio > 0])
data <- subset(data, Ratio > 0)

# # Grab whatever portion of the data you are interested in...
# calData <- subset(data, Cells=='None')
# LNCaPData <- subset(data, Cells=='LNCaP')
# PC3Data <- subset(data, Cells=='PC3')

# Summarize the data
summary <- data.frame()
for(hiv in HIVLevels)
{
     temp <- subset(data, HIV==hiv)
     summary <- rbind(summary, data.frame(HIV=hiv, RatioMean=median(temp$Ratio), RatioMad=mad(temp$Ratio)))
}

# Plot the histograms
plot(c(), c(), xlim=c(0.001,100), ylim=range(0,8), xlab='HIV Concentration [virus/well]', ylab='Signal per Bead [au]', log='x')
for(hiv in HIVLevels)
{
     temp <- subset(data, HIV==hiv)
     HIVNumber <- as.numeric(str_extract(as.character(hiv), "\\d+\\.*\\d*"))
     if(hiv=='0')
     {
          points(temp$HIVNumber + 0.001, temp$Ratio, pch=20, col=rgb(0,0,0,0.02))
          points(HIVNumber + 0.001, mean(temp$Ratio), pch=20, col='red')
     }
     else
     {
          points(temp$HIVNumber, temp$Ratio, pch=20, col=rgb(0,0,0,0.02))
          points(HIVNumber, mean(temp$Ratio), pch=20, col='red')
     }
}

# Plot the histograms
plot(density(data$Ratio), col=rgb(1,1,1,0.5), xlim=c(0,5), ylim=range(0,2), xlab='Signal per Bead [au]', ylab='Density', main='')
cols <- seq(0,1,length.out=length(HIVLevels)+1)
i <- 1
for(hiv in HIVLevels)
{
     temp <- subset(data, HIV==hiv)
     lines(density(temp$Ratio), col=gray(cols[i]), lwd=2)
     print(hiv)
     print(cols[i])
     i = i + 1;
}
legend('topright', legend=HIVLevels, col=gray(cols[1:length(HIVLevels)]), lwd=2)



weird <- subset(data, Ratio < 0)
print(weird)

