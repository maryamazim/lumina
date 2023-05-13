# Script to structure trial 1 raw data coming from the machine into R DFs
# Jul 20, 2022
# Oct 30, 2022

# Read in Illumination Data 
projectFolder=" "  # include path to the project folder

illumData <- read.csv(paste0(projectFolder,"illumDataPC.csv"), header=F)
odData <- read.csv(paste0(projectFolder,"odDataPC.csv"), header=F)
delMatrix <- as.matrix(read.csv(paste0(projectFolder,"design_matrix.csv"), header=F,
                                 check.names=F,  col.names=2:6, row.names=c('B','C','D','E','F', 'G', 'H', 'I','J','K')))

# Function to return a DF with illum data of all replicates and types at a time point!  
get_1TT <- function(d, delMatrix) {
  illum <- d[-(1:2)]
  illum[illum=="OVER"] <- max(as.numeric(illum[illum != "OVER"]))
  df <- data.frame(replicate=rep(row.names(delMatrix), each=10), 
                   type=c(t(delMatrix)),
                   time = as.numeric(d[1]),
                   temperature =as.numeric(d[2]),
                   illumination = as.numeric(illum))
  df
}
# extract and merge illum and od data, one data frame at a time!
illumDF <- NULL
for(i in 2:ncol(illumData)) {
  df <- get_1TT(illumData[,i], delMatrix=delMatrix)
  illumDF <- rbind(illumDF, df)
}
odDF <- NULL
for(i in 2:ncol(odData)) {
  df <- get_1TT(odData[,i], delMatrix=delMatrix)
  odDF <- rbind(odDF, df)
}
names(odDF)[5] <- 'growth'
illumDF$growth <- odDF$growth

for (i in 1:nrow(illumDF)) {
  x <- unlist(strsplit(as.character(illumDF$type[i]), split= ' '))
  illumDF$plmr[i] = x[1]
  illumDF$conc[i] = as.numeric(strtrim(x[2], width=nchar(x[2])-1))
}

# Look at data using poly()
# illum data
illumDF$type <- factor(illumDF$type)[,drop=T]
illumDF$plmr <- factor(illumDF$plmr)[,drop=T]
illumDF$hrs <- illumDF$time/(60^2)

## The model below shows that Mucus at .5% depresses QS the most
plotType <- function(iData, whatPlmr, pred=0) {
  x <- tapply(iData$illumination, iData$type, mean)
  x <- x[grep(names(x), pattern=whatPlmr)]
  y <- as.numeric(gsub("%", "", unlist(strsplit(names(x), split= ' '))[seq(2, 2*length(x), 2)]))
  z = data.frame(y=sort(y), x = x[order(y)])
  mnConc = z[z$x == min(z$x), 'y']; mnIllum = round(z[z$x == min(z$x), 'x'])
  if(whatPlmr != 'Mucus') {
    plot(iData[iData$plmr==whatPlmr, 'conc'], iData[iData$plmr==whatPlmr, 'illumination'], 
         pch='.', xlab='Concentration', ylab='Fluorescence', 
         main=paste0(whatPlmr, ': min at (', mnConc, ', ', mnIllum, ')'))
  } else {
    plot(iData[iData$plmr==whatPlmr, 'conc'], iData[iData$plmr==whatPlmr, 'illumination'], 
         pch='.', xlab='Concentration', ylab='Fluorescence', 
         main=paste0('Mucin', ': min at (', mnConc, ', ', mnIllum, ')'))
  }
  points(z$y, z$x, col='orange', pch=16, type='b')
  points(z[z$x == min(z$x), 'y'], z[z$x == min(z$x), 'x'], pch=16, col='white')
  points(z[z$x == min(z$x), 'y'], z[z$x == min(z$x), 'x'], pch=25)
}
plotType(illumDF, 'Mucus')
plotType(illumDF, 'Agarose')
plotType(illumDF, 'Alginate')
plotType(illumDF, 'CMC')
plotType(illumDF, 'eDNA')

# Plot minimums of all polymers in the study --
illum_min <- data.frame(plmr = c('Mucin', 'Agarose', 'Alginate', 'CMC', 'eDNA'), 
                        conc = c(.5, .25, .25, 0, .002), 
                        illum = c(9071, 16042, 13530, 17113, 16383))

# Model for each polymer-concentration combination as a type 
illumDF.poly <- lm(illumination ~ type + poly(growth,3) + temperature + poly(hrs, 5), 
                   data = illumDF[illumDF$conc != 0,])

# Model of choice  
illumDF.poly <- lm(illumination ~ plmr+ poly(conc, 2) %in% plmr + temperature + poly(growth,3) + poly(hrs, 5), 
                  data = illumDF)

summary(illumDF.poly)
anova(illumDF.poly)

newD <- illumDF
newD$temperature <- mean(newD$temperature)
newD$growth <- mean(newD$growth)
newD[,c('pred', 'lower', 'upper')] <- predict(illumDF.poly, newdata=newD, interval = 'confidence')

plot(newD$time, newD$illumination)  # not adjusted for anything
plot(newD$time, newD$pred)  # adjusted for OD

# Minimums differed after full-model adjustment
newD2 = newD; newD2$illumination = newD2$pred
plotType(newD2, 'Mucus')
plotType(newD2, 'Agarose')
plotType(newD2, 'Alginate')
plotType(newD2, 'CMC')
plotType(newD2, 'eDNA')

