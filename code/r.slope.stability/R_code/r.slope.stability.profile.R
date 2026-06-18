#!/usr/bin/R

##############################################################################
#
# MODULE:       r.slope.stability.map.R
#
# AUTHOR:       Martin Mergili
#
# PURPOSE:      The slope stability model - R script for the preparation
#               of a longitudinal profile through an ellipsoid
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
prefix <- args[1]  #prefix for file names
numlayers <- as.integer(args[2])  #number of layers
layercrit <- as.integer(args[3])  #number of critical layer
foscrit <- as.numeric(args[4])  #factor of safety of critical layer

# Defining data file name
intable <- paste(prefix, "_results/", prefix, "_files/", prefix, "_profile.txt", 
    sep = "")

# Creating vectors from data
horpos <- scan(intable, nlines = numlayers)  #horizontal position
elev <- scan(intable, skip = numlayers, nlines = numlayers)  #surface elevation
layers <- scan(intable, skip = 2 * numlayers, nlines = numlayers)  #elevation of layers
slips <- scan(intable, skip = 3 * numlayers, nlines = numlayers)  #elevation of slip surfaces
horposcrit <- scan(intable, skip = layercrit - 1, nlines = 1)  #horizontal position of critical slip surface
slipcrit <- scan(intable, skip = 3 * numlayers + layercrit - 1, nlines = 1)  #elevation of critical slip surface
seepres <- scan(intable, skip = 4 * numlayers, nlines = 1)  #contribution of seepage force to shear resistance
seepforce <- scan(intable, skip = 4 * numlayers + 1, nlines = 1)  #contribution of deepage force to shear force
shearres <- scan(intable, skip = 4 * numlayers + 2, nlines = 1)  #shear resistance
shearforce <- scan(intable, skip = 4 * numlayers + 3, nlines = 1)  #shear force
groundwatertable <- scan(intable, skip = 4 * numlayers + 4, nlines = 1)  #elevation of groundwater table

# Computing difference between shear resistance and shear force
difference <- shearforce + shearres

# Setting parameters needed for plot layout
maxhorpos <- max(horpos, na.rm = TRUE)  #maximum horizontal position
maxelev <- max(elev, na.rm = TRUE)  #maximum terrain elevation
minelev <- min(layers, na.rm = TRUE)  #minimum elevation of deepest layer
minforce <- min(shearforce, na.rm = TRUE)  #minimum of all forces
maxforce <- max(shearres, na.rm = TRUE)  #maximum of all forces
diffel <- maxelev - minelev  #difference between maximum and minimum elevation
diffforce <- maxforce - minforce  #difference between maximum and minimum force
minelev <- minelev - 0.04 * diffel  #correcting minimum elevation (for neat display)
maxelev <- minelev + diffel * 1.24  #correcting maximum elevation (for neat display)
minforce <- maxforce - diffforce * 1.2  #correcting minimum force (for neat display)
intel <- maxelev - minelev  #updating difference between maximum and minimum elevation
intforce <- maxforce - minforce  #updating difference between maximum and minimum force
ratio <- intforce/intel  #computing ratio between force and elevation differeences (for neat display)

# Formatting legend labels
fos = vector("expression", 1)  #initializing vectors
rleg <- vector("expression", 1)
tleg <- vector("expression", 1)
fsleg <- vector("expression", 1)
diffleg <- vector("expression", 1)
fos[1] <- substitute(expression(FoS[crit] == fosx), list(fosx = format(round(foscrit, 
    digits = 2), nsmall = 2)))[2]  #fator of safety
rleg <- expression(italic(R))  #shear resistance
tleg <- expression(italic(T))  #shear force
fsleg <- expression(italic(T)[s])  #seepage force
diffleg <- expression(italic(R) - italic(T))  #difference between shear resistance and shear force

# Creating plot file
profileplot = paste(prefix, "_results/", prefix, "_plots/", prefix, "_plot_profile.tif", 
    sep = "")
tiff(filename = profileplot, width = 12, height = 14, units = "cm", res = 300)

# Defining margins
par(mar = c(3.2, 3.2, 0.5, 3.2))

# Plotting bars
barplot(shearres, axes = F, xlab = NA, ylab = NA, ylim = c(minforce, maxforce + ratio * 
    intel), border = NA, col = rgb(0.5, 0.85, 0.5))  #shear resistance

par(new = TRUE)
barplot(shearforce, axes = F, xlab = NA, ylab = NA, ylim = c(minforce, maxforce + 
    ratio * intel), border = NA, col = rgb(0.85, 0.5, 0.5))  #shear force

par(new = TRUE)
barplot(seepres, axes = F, xlab = NA, ylab = NA, ylim = c(minforce, maxforce + ratio * 
    intel), border = NA, col = rgb(0.7, 0.7, 0.9))  #contribution of seepage force to shear resistance

par(new = TRUE)
barplot(seepforce, axes = F, xlab = NA, ylab = NA, ylim = c(minforce, maxforce + 
    ratio * intel), border = NA, col = rgb(0.7, 0.7, 0.9))  #contribution of seepage force to shear force

par(new = TRUE)
barplot(difference, axes = F, xlab = NA, ylab = NA, ylim = c(minforce, maxforce + 
    ratio * intel), border = rgb(0, 0, 0), col = NA)  #difference between shear resistance and shear force

axis(side = 4, tck = -0.02, labels = NA)  #y axis for bars
axis(side = 4, lwd = 0, line = -0.4)

# Plotting lines
par(new = TRUE)
plot(x = horpos, y = slips, axes = F, xlab = NA, ylab = NA, xlim = c(0, maxhorpos), 
    ylim = c(minelev - intforce/ratio, maxelev), type = "l", lty = 1, lwd = 1, col = "brown")  #elevation of slip surfaces
lines(x = horposcrit, y = groundwatertable, lty = 3, lwd = 1.5, col = "blue")  #elevation of groundwater table
lines(x = horposcrit, y = slipcrit, lty = 1, lwd = 2, col = "red")  #elevation of critical slip surface
lines(x = horpos, y = layers, lty = 3, lwd = 1, col = "green")  #elevation of layers
lines(x = horpos, y = elev, lty = 1, lwd = 2, col = "black")  #terrain elevation

# Legends legend('topleft', legend='Terrain', lty=1, lwd=2, col='black',
# text.col='black', bty='n', inset=c(0.0,-0.01), horiz=TRUE, x.intersp=0.5)
# #terrain elevation
legend("topleft", legend = "GW", lty = 3, lwd = 1.5, col = "blue", text.col = "black", 
    bty = "n", inset = c(0, -0.01), horiz = TRUE, x.intersp = 0.5)  #terrain elevation
legend("topleft", legend = "Layer", lty = 3, lwd = 1, col = "green", text.col = "black", 
    bty = "n", inset = c(0.29, 0 - 0.01), horiz = TRUE, x.intersp = 0.5)  #elevation of layers
legend("topleft", legend = "Ellipsoid", lty = 1, lwd = 1, col = "brown", text.col = "black", 
    bty = "n", inset = c(0.56, -0.01), horiz = TRUE, x.intersp = 0.5)  #elevation of slip surfaces
legend("top", legend = fos[1], lty = 1, lwd = 2, col = "red", text.col = "red", bty = "n", 
    inset = c(0, 0.04))  #elevation of critical slip surface, factor of safety
legend("bottomleft", legend = rleg, border = NA, fill = rgb(0.5, 0.85, 0.5), text.col = "black", 
    bty = "n", inset = c(0.05, 0), horiz = TRUE)  #shear resistance
legend("bottomleft", legend = tleg, border = NA, fill = rgb(0.85, 0.5, 0.5), text.col = "black", 
    bty = "n", inset = c(0.25, 0), horiz = TRUE)  #shear force
legend("bottomleft", legend = fsleg, border = NA, fill = rgb(0.7, 0.7, 0.9), text.col = "black", 
    bty = "n", inset = c(0.45, 0), horiz = TRUE)  #seepage force
legend("bottomleft", legend = diffleg, border = "black", fill = NA, text.col = "black", 
    bty = "n", inset = c(0.65, 0), horiz = TRUE)  #difference between shear resistance and shear force

box()  #bounding box
axis(side = 1, tck = -0.02, labels = NA)  #x axis
axis(side = 1, lwd = 0, line = -0.4)
axis(side = 2, tck = -0.02, labels = NA)  #y axis for lines
axis(side = 2, lwd = 0, line = -0.4)

# Plotting axis labels
mtext(side = 1, "Horizontal distance (m)", line = 2)  #x axis
mtext(side = 2, "Elevation (m)", line = 2)  #y axis for lines
mtext(side = 4, "Forces (kN)", line = 2)  #y axis for bars

# Closing plot file
dev.off()
