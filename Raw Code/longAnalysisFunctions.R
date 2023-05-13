# Script with functions from QS4 longAnalysis scripts
# Aug 13, 2022
#

# Function to return a DF with illum data of all replicates and types at a time point!  
get_1TT <- function(d, delMatrix, num) {
  illum <- d[-(1:2)]
  df <- data.frame(replicate=rep(row.names(delMatrix), each=num), 
                   type=c(t(delMatrix)),
                   time = d[1],
                   temperature = d[2],
                   illumination = illum)
  df
}

# plot data with specific adjustments
plot_type <- function(illumDF, type, y='illumination', y.lab="", cex=1, pch=16, average=FALSE, points=FALSE, col='black', main="", cex.main=1, ulMargin=0) {
  typeDF <- illumDF[illumDF$type == type,]
  if(y.lab == "") { 
    y.lab = y
    if(average == TRUE) y.lab = paste0(y.lab, ' (Average)')
  }
  if(main == "") main=paste0('type = ', type)
  y.idx = grep(paste0("^", y, "$"), names(typeDF)) # index of variable to plot
  if(average == TRUE) { 
    x <- tapply(typeDF[,y.idx], typeDF$time, mean)
    typeDFAvg <- data.frame(time = as.numeric(names(x)), illumAvg = x)
    plotRange = c(0, max(typeDFAvg$illumAvg)+ulMargin)
    if(points == TRUE) {
      points(typeDFAvg$time/(60^2), typeDFAvg$illumAvg, col=col, pch=pch, cex=cex, ylim=plotRange)
    } else {
      plot(typeDFAvg$time/(60^2), typeDFAvg$illumAvg, xlab='Time (Hours)', 
           ylab=y.lab, col=col, pch=pch, cex=cex, ylim=plotRange,
           main=main, cex.main=cex.main)
    }
  } else {
    plotRange = c(0, max(typeDF[,y.idx])+ulMargin)
    if(points == TRUE) {
      points(typeDF$time/(60^2), typeDF[,y.idx], col=col, pch=pch, cex=cex, ylim=plotRange)
    } else {
      plot(typeDF$time/(60^2), typeDF[,y.idx], xlab='Time (Hours)', 
           ylab=y.lab, col=col, pch=pch, cex=cex, ylim=plotRange,
           main=main, cex.main=cex.main)
    }
  }
}

# This function plot data as is and connects averages at each time point
# send a data frame with columns 'type' and 'illumination' and a type (whatPlmr) to plot
plotType <- function(iData, whatPlmr) {
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
# To use: plotType(illumDF, 'Mucus')


