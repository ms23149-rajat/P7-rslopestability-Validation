#!/usr/bin/R

##############################################################################
#
# MODULE:       r.slope.stability.evolution.R
#
# AUTHOR:       Martin Mergili
#
# PURPOSE:      The slope stability model - R script for visualizing
#               the evolution of areas with factor of safety < 1 or of the
#               average slope failure probability
#
# COPYRIGHT:    (c) 2013 - 2021 by the author
#               (c) 2020 - 2021 by the University of Graz
#               (c) 2013 - 2020 by the BOKU University, Vienna
#               (c) 2015 - 2020 by the University of Vienna
#               (C) 1993 - 2021 by the R Development Core Team
#
# VERSION:      20210525 (25 May 2021)
#
#               This program is free software under the GNU General Public
#               License (>=v2).
#
##############################################################################


# Loading library
library(stats)

# Importing arguments defined in r.slope.stability.py
args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]
numtiles <- as.integer(args[2])
pflag <- args[3]
status_max <- as.numeric(args[4])

# Defining data file name
intable <- paste(prefix, "_results/", prefix, "_files/", prefix, "_evolution.txt", 
    sep = "")

# Creating vectors from data
elldens <- scan(intable, nlines = 1)  #ellipsoid density
numsim <- scan(intable, skip = 1, nlines = 1)  #number of ellipsoids
elldens_max <- length(elldens)

# Creating plot file
evolplot = paste(prefix, "_results/", prefix, "_plots/", prefix, "_plot_evolution.png", 
    sep = "")
png(filename = evolplot, width = 12, height = 10, units = "cm", res = 300)

# Defining margins
par(mar = c(3.2, 3.2, 0.5, 1))

# Preparing plot
plot(NULL, axes = F, xlab = NA, ylab = NA, xlim = c(1, elldens_max), ylim = c(-0.025 * 
    status_max, 1.025 * status_max), type = "l", lty = 1, lwd = 3, col = "blue", 
    log = "x")

# Plotting lines
sumpix <- scan(intable, skip = 2, nlines = 1)  #number of pixels
sumstatus <- scan(intable, skip = 3, nlines = 1)  #indicator for instability, sum
if (numtiles > 1) {
    # for multi-core processing:
    sumpix <- sumpix * 0  #resetting number of pixels
    sumstatus <- sumstatus * 0  #resetting indicator for instability, sum
    j <- 1  #initializing counter for loop over all tiles
    repeat {
        # starting loop over all tiles
        if (j > numtiles) 
            break  #break if last tile was reached
 else {
            numpix <- scan(intable, skip = 2 + 4 * (j - 1), nlines = 1)  #number of pixels
            statusabs <- scan(intable, skip = 3 + 4 * (j - 1), nlines = 1)  #indicator for instability, sum
            status <- scan(intable, skip = 4 + 4 * (j - 1), nlines = 1)  #indicator for instability, average
            sumpix <- sumpix + numpix  #updating sum number of pixels
            sumstatus <- sumstatus + statusabs  #updating sum of indicator for instability
            evol <- scan(intable, skip = 5 + 4 * (j - 1), nlines = 1)  #evolution of indicator
            lines(elldens, status, lty = 3, lwd = 1, col = "blue")  #plotting data points as line
            j <- j + 1
        }
    }
}
lines(elldens, sumstatus/sumpix, lty = 1, lwd = 3, col = "red")  #plotting data points of average over all tiles or for single result as line

# Defining requirement for zero point (number of simulations divided by ellipsoid
# density without change) req_zeropoint=100

# Plotting data points j<-1 #initializing counter for loop over all tiles repeat
# { #starting loop over all tiles if ( j > numtiles ) break #break if last tile
# was reached else {

# status<-scan(intable, skip=4+4*(j-1), nlines=1) #indicator for instability,
# average evol<-scan(intable, skip=5+4*(j-1), nlines=1) #evolution of indicator

# Vector for zero point i2<-1 #initial counter for new vector ctrl_evol<-1
# #initial control for new vector elldenssel<-vector() #selected ellipsoid
# density numsimsel<-vector() #selected number of ellipsoids statussel<-vector()
# #selected number of critical cells evolsel<-vector() #selected change of number
# of critical cells i<-2 #initial control for loop status_prev<--1 #initial
# number of critical cells of previous loop while ( i <= elldens_max ) { #loop
# over all ellipsoid densities: if ( status[i] == status_prev ) { #if change is
# zero: if ( ctrl_evol == 1 ) { #if control is still positive:
# elldenssel[i2]<-elldens[i] #writing values to new vector
# numsimsel[i2]<-numsim[i] statussel[i2]<-status[i] evolsel[i2]<-evol[i] i2<-i2+1
# #updating counter for new vector ctrl_evol<-0 #updating control for new vector
# } } status_prev<-status[i] #updating number of critical cells of previous loop
# i<-i+req_zeropoint #increasing control for loop by 10 }
# points(elldenssel,statussel, pch=16, col='green') #plotting zero data point

# Labelling zero data point if reached if ( ctrl_evol == 0 ) {

# if ( numtiles == 1 ) { #if non-parallelized mode applies:

# dell=vector('expression',1) dell[1]=substitute(expression(italic(d)[ell] ==
# dellips), list(dellips = elldenssel[1]))[2] #label text for ellipsoid density

# nell=vector('expression',1) nell[1]=substitute(expression(italic(n)[ell] ==
# nellips), list(nellips = numsimsel[1]))[2] #label text for number of ellipsoids

# text(x=elldenssel, y=statussel+0.06, labels=dell, col='green') #Labelling zero
# point - ellipsoid density text(x=elldenssel, y=statussel+0.12, labels=nell,
# col='green') #Labelling zero point - number of ellipsoids } } j <- j + 1 } }

# Plotting bounding box and axes
box()
axis(side = 1, tck = -0.02, labels = NA)
axis(side = 2, tck = -0.02, labels = NA)
axis(side = 1, lwd = 0, line = -0.4)
axis(side = 2, lwd = 0, line = -0.4)

# Legend legend('topleft', pch=c(16,NA), pt.cex=c(1,NA), lty=c(NA,1),
# lwd=c(NA,3), col=c('green','blue'), legend=c('Zero point','Evolution'),
# text.col='black', bty='n', inset=c(0.05,0))
if (numtiles == 1) {
    legend("topleft", pch = NA, pt.cex = NA, lty = 1, lwd = 3, col = "red", legend = c("Evolution"), 
        text.col = "black", bty = "n", inset = c(0.05, 0))
} else {
    legend("topleft", pch = NA, pt.cex = NA, lty = c(1, 3), lwd = c(3, 1), col = c("red", 
        "blue"), legend = c("Average", "Tiles"), text.col = "black", bty = "n", inset = c(0.05, 
        0))
}

# Plotting axis labels
mtext(side = 1, "Ellipsoid density", line = 2)
if (pflag == "0") {
    mtext(side = 2, "Fraction of area with factor of safety < 1", line = 2)
} else {
    mtext(side = 2, "Average of slope failure probability", line = 2)
}

# Closing plot file
dev.off()
