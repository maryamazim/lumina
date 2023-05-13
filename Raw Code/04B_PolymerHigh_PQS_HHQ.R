# Script to structure raw data of higher concentrations of mucus with pqsBC and PQS with WT strain 
# coming from the machine into R DFs
# Dec 20, 2022
#
library(reshape2)

# Read in Illumination Data 
projectFolder=" "  # include path to the project folder 

illumData <- read.csv(paste0(projectFolder,"illumData.csv"), header=F)
odData <- read.csv(paste0(projectFolder,"odData.csv"), header=F)
delMatrix <- as.matrix(read.csv(paste0(projectFolder,"design_matrix.csv"), header=F,
                                check.names=F,  col.names=1:12, row.names=c('A','B','C','D', 'E','F','G','H')))

# functions to create data frame and plot data
source(paste0(projectFolder, "longAnalysisFunctions.R"))

# extract and merge illum and od data, one data frame at a time!
illumDF <- NULL
for(i in 2:ncol(illumData)) {
  df <- get_1TT(illumData[,i], delMatrix=delMatrix, num=12)
  illumDF <- rbind(illumDF, df)
}
odDF <- NULL
for(i in 2:ncol(odData)) {
  df <- get_1TT(odData[,i], delMatrix=delMatrix, num=12)
  odDF <- rbind(odDF, df)
}
names(odDF)[5] <- 'growth'
illumDF$growth <- odDF$growth

# create mucin and HHQ concentrations
illumDF$PA = 1; 
illumDF[illumDF$type=="MUC(.005%)"|illumDF$type=="MUC(.05%)"|illumDF$type=="MUC(.5%)"|illumDF$type =="LB", 'PA'] = 0;
illumDF$mcnConc = 0; 
illumDF[grep(illumDF$type, pattern="\\(.005%\\)"), 'mcnConc'] = .005
illumDF[grep(illumDF$type, pattern="\\(.05%\\)"), 'mcnConc'] = .05; 
illumDF[grep(illumDF$type, pattern="\\(.5%\\)"), 'mcnConc'] = .5

illumDF$hhqConc = 0; 
illumDF[grep(illumDF$type, pattern="+HHQ .01%"), 'hhqConc'] = .01
illumDF[grep(illumDF$type, pattern="+HHQ .02%"), 'hhqConc'] = .02; 
illumDF[grep(illumDF$type, pattern="+HHQ .05%"), 'hhqConc'] = .05
illumDF[grep(illumDF$type, pattern="+HHQ .1%"), 'hhqConc'] = .1; 
illumDF[grep(illumDF$type, pattern="+HHQ .5%"), 'hhqConc'] = .5
illumDF[grep(illumDF$type, pattern="+HHQ 1%"), 'hhqConc'] = 1
illumDF$hrs <- illumDF$time/(60^2)

mucIllum = illumDF[illumDF$PA == 0, ]  # culture not included
illumDF <- illumDF[illumDF$PA == 1, ]  # culture included

# Set illumination/absorbance covariates - 
# Adjustment 1 (adj1) is to subtract illumination of Muc w/t culture across time
illumCov <- tapply(mucIllum$illumination, list(mucIllum$time, mucIllum$mcnConc), mean)
absCov <- tapply(mucIllum$growth, list(mucIllum$time, mucIllum$mcnConc), mean)
illumDF$illumCov <- 0
illumDF$absCov <- 0
for(i in 1:nrow(illumDF)) {
  illumDF[i, 'illumCov'] <- illumCov[as.character(illumDF[i, 'time']), as.character(illumDF[i, 'mcnConc'])]
  illumDF[i, 'absCov'] <- absCov[as.character(illumDF[i, 'time']), as.character(illumDF[i, 'mcnConc'])]
}
illumDF$illumination.adj1 = illumDF$illumination - illumDF$illumCov
illumDF$growth.adj1 = illumDF$growth - illumDF$absCov

# Another way to adjust for substance illumination is to subtract base illumination of the same media 
# Base illumination per mcnConc and hhqConc combination
basetimerange = c(1, 2) # min and max base time in hours
`%inrange%` <- function(value, range) (value <= range[2] & value >= range[1])
xtrct = illumDF$hrs %inrange% basetimerange
# adj.2 illumination
base_illum <- tapply(
  illumDF[xtrct,'illumination'], 
  list(illumDF[xtrct,'mcnConc'], illumDF[xtrct,'hhqConc']), mean)

for(i in 1:nrow(illumDF)) 
  illumDF[i, 'base_illum'] <- base_illum[as.character(illumDF$mcnConc[i]), as.character(illumDF$hhqConc[i])]
illumDF$illumination.adj2 <- illumDF$illumination - illumDF$base_illum

# adj.2 growth
base_growth <- tapply(
  illumDF[xtrct,'growth'], 
  list(illumDF[xtrct,'mcnConc'], illumDF[xtrct,'hhqConc']), mean)

for(i in 1:nrow(illumDF)) 
  illumDF[i, 'base_growth'] <- base_growth[as.character(illumDF$mcnConc[i]), as.character(illumDF$hhqConc[i])]
illumDF$growth.adj2 <- illumDF$growth - illumDF$base_growth


mcn=0.5
iul = max(mucIllum[mucIllum$mcnConc==mcn, 'illumination'])+500
x <- tapply(mucIllum[mucIllum$mcnConc==mcn, 'illumination'], mucIllum[mucIllum$mcnConc==mcn, 'hrs'], mean)
png(file="/Users/mazim/Documents/RahmeLab/Summer2022/Analysis/QS4/manuscript/high_res_plots/Figure_2.png",
    res = 1200, units = 'in', width = 6.2, height = 5.56)
plot(names(x), x, pch=16, ylab='Observed Fluorescence', xlab='Time (hours)', 
     ylim=c(0,iul), main="Fluorescence of Mucin - no culture")
grid(NULL, NULL)
mcn=.05
x <- tapply(mucIllum[mucIllum$mcnConc==mcn, 'illumination'], mucIllum[mucIllum$mcnConc==mcn, 'hrs'], mean)
points(names(x), x, pch=17, col='dark red') 
mcn=0
x <- tapply(mucIllum[mucIllum$mcnConc==mcn, 'illumination'], mucIllum[mucIllum$mcnConc==mcn, 'hrs'], mean)
points(names(x), x, pch=18, col='orange') 
# (.5) stable at around 5000
# (.05) stable at around 3000
# (.005) upward trend with possible contamination for one of the replicates.
legend(0, 1500, legend=c('MUC 0.5%', 'MUC 0.05%', 'MUC 0%'), pch=16:18, col=c('black', 'dark red', 'orange'))
dev.off()

# Plot absorbance 
gul = 1
mcn=.5
plot(mucIllum[mucIllum$mcnConc==mcn, 'hrs'], 
     mucIllum[mucIllum$mcnConc==mcn, 'growth'], 
     pch=16, ylab='Absorbance', xlab='Time (hours)', ylim=c(0,gul),
     main="Absorbance of Mucin - no culture")
mcn=.05
points(mucIllum[mucIllum$mcnConc==mcn, 'hrs'], 
       mucIllum[mucIllum$mcnConc==mcn, 'growth'],
       pch = 17, col = "dark red")
mcn=0
points(mucIllum[mucIllum$mcnConc==mcn, 'hrs'], 
       mucIllum[mucIllum$mcnConc==mcn, 'growth'],
       pch = 18, col = "orange")
legend(0, 1, legend=c('MUC 0.5%', 'MUC 0.05%', 'MUC 0%'), pch=16:18, col=c('black', 'dark red', 'orange'))

# (.5) slight trend down (at .4 - .3) - Mucin shouldn't degrade
# (.05) stable at around 0.1
# (.005) noticeable growth - obvious contamination 

plot_type(illumDF, 'MUC(0%)+HHQ 0%', y="illumination.adj2", col='bisque1', ulMargin=2000)
plot_type(illumDF, 'MUC(0%)+HHQ .01%', y="illumination.adj2", points=T, col='burlywood3')
plot_type(illumDF, 'MUC(0%)+HHQ .02%', y="illumination.adj2", points=T, col='coral3')
plot_type(illumDF, 'MUC(0%)+HHQ .05%', y="illumination.adj2", points=T, col='darkgray')
plot_type(illumDF, 'MUC(0%)+HHQ .1%', y="illumination.adj2", points=T, col='darkred')
plot_type(illumDF, 'MUC(0%)+HHQ .5%', y="illumination.adj2", points=T, col='darkorchid4')
 
plot_type(illumDF, 'MUC(0%)+HHQ 1%', y="illumination.adj2", points=T, col='darkblue')
# plots show that HHQ >= .01 were enough for causing QS activity
# Probably greater HHQ concentrations clouds the solution and partially blocks light

# Check 0 mucin concentrations 
plot_type(illumDF, 'MUC(0%)+HHQ 0%', y="illumination.adj2", y.lab="Fluorescence", cex=.7, pch=17, average=T, col='darkgray', 
          main="C. Mucin-downregulated QS at HHQ of 0.01 (observed data)", cex.main=1, ulMargin=2000)
plot_type(illumDF, 'MUC(.005%)+HHQ 0%', y="illumination.adj2", y.lab="Fluorescence", cex=.7, pch=16, average=T, points=T, col='blue')
plot_type(illumDF, 'MUC(.05%)+HHQ 0%', y="illumination.adj2", y.lab="Fluorescence", cex=.7, pch=18, average=T, points=T, col='black')
plot_type(illumDF, 'MUC(.5%)+HHQ 0%', y="illumination.adj2", y.lab="Fluorescence", cex=.7, pch=15, average=T, points=T, col='coral3')
legend(0, 4100, legend=c('Muc 0.005', 'Muc 0.05', 'Muc 0.5', 'Baseline'), col=c('blue','black','coral3', 'darkgray'), pch=c(16,18,15,17), cex=.8)

# Look at impact of Muc 0.5 on varying levels of hhq
d <- illumDF[illumDF$mcnConc == 0.5 & illumDF$hhqConc != 1,]
d$hhqConc <- factor(d$hhqConc)[,drop=T]
m <- lm(illumination.adj2 ~ hhqConc + poly(hrs, 3) + temperature, data = d)

summary(m)
newD <- d
newD$temperature <- 0 #mean(newD$temperature)
newD$pred <- predict(m, newdata=newD)


# Look at activity of hhq concentrations at 0.05 and 0.5 levels of Mucin
# hhq of { 0, .01, .02 }
d <- illumDF[illumDF$hhqConc <= .02 & illumDF$mcnConc >= 0,]
d$mcnConc <- factor(d$mcnConc)[,drop=T]; d$hhqConc <- factor(d$hhqConc)[,drop=T]
m <- lm(illumination.adj2 ~ poly(hrs, 3) * mcnConc + hhqConc + temperature, data = d)
summary(m)
newD <- d
newD$temperature <- 0 
newD$pred <- predict(m, newdata=newD)
ypl = range(newD$pred)
xtrct = newD$hhqConc==.01 & newD$mcnConc==.05
plot(newD[xtrct, 'hrs'], newD[xtrct,'pred'],
     pch=18, cex=.7, col='black',ylab='Fluorescence', ylim=ypl,
     xlab='Time (Hours)', main="A. Mucin-downregulated QS at HHQ of 0.01 (model-adjusted)", cex.main=1.3)
xtrct = newD$hhqConc==.01 & newD$mcnConc==.5
points(newD[xtrct, 'hrs'], newD[xtrct,'pred'], pch=15, cex=.7, col='coral3')
xtrct = newD$hhqConc==.01 & newD$mcnConc==.005
points(newD[xtrct, 'hrs'], newD[xtrct,'pred'], pch=16, cex=.7, col='blue')
xtrct = newD$hhqConc==0 & newD$mcnConc==0
points(newD[xtrct, 'hrs'], newD[xtrct,'pred'], pch=17, cex=.7, col='darkgray')
legend(0, max(newD$pred), legend=c('Muc 0.005', 'Muc 0.05', 'Muc 0.5', 'Baseline'), col=c('blue','black','coral3', 'darkgray'), pch=c(16,18,15,17), cex=.8)

xtrct = newD$hhqConc==.02 & newD$mcnConc==.05
plot(newD[xtrct, 'hrs'], newD[xtrct,'pred'],
     pch=18, cex=.7, col='black',ylab='Fluorescence', ylim=ypl,
     xlab='Time (Hours)', main="B. Mucin-downregulated QS at HHQ of 0.02 (model-adjusted)", cex.main=1)
xtrct = newD$hhqConc==.02 & newD$mcnConc==.5
points(newD[xtrct, 'hrs'], newD[xtrct,'pred'], pch=15, cex=.7, col='coral3')
xtrct = newD$hhqConc==.02 & newD$mcnConc==.005
points(newD[xtrct, 'hrs'], newD[xtrct,'pred'], pch=16, cex=.7, col='blue')
xtrct = newD$hhqConc==0 & newD$mcnConc==0
points(newD[xtrct, 'hrs'], newD[xtrct,'pred'], pch=17, cex=.7, col='darkgray')
legend(0, max(newD$pred), legend=c('Muc 0.005', 'Muc 0.05', 'Muc 0.5', 'Baseline'), col=c('blue','black','coral3', 'darkgray'), pch=c(16,18,15,17), cex=.8)

# Check 0 HHQ concentration
xtrct = newD$hhqConc==0 & newD$mcnConc==.05
plot(newD[xtrct, 'hrs'], newD[xtrct,'pred'],
     pch=18, cex=.7, col='black',ylab='Fluorescence', ylim=ypl,
     xlab='Time (Hours)', main="B. Mucin-downregulated QS at HHQ of 0.02 (observed data)", cex.main=1)
xtrct = newD$hhqConc==0 & newD$mcnConc==.5
points(newD[xtrct, 'hrs'], newD[xtrct,'pred'], pch=15, cex=.7, col='coral3')
xtrct = newD$hhqConc==0 & newD$mcnConc==.005
points(newD[xtrct, 'hrs'], newD[xtrct,'pred'], pch=16, cex=.7, col='blue')
xtrct = newD$hhqConc==0 & newD$mcnConc==0
points(newD[xtrct, 'hrs'], newD[xtrct,'pred'], pch=17, cex=.7, col='darkgray')
legend(0, max(newD$pred), legend=c('Muc 0.005', 'Muc 0.05', 'Muc 0.5', 'Baseline'), col=c('blue','black','coral3', 'darkgray'), pch=c(16,18,15,17), cex=.8)
