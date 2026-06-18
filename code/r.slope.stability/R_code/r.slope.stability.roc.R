#!/usr/bin/R

##############################################################################
#
# MODULE:       r.slope.stability.roc.R
#
# AUTHOR:       Martin Mergili
#
# PURPOSE:      The slope stability model - R script for the preparation
#               of an ROC plot relating the computed slope failure proba-
#               bility or susceptibility to observed landslide areas
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


# Loading libraries
# library(rgdal)
library(ROCR)

# Importing arguments (defined in r.slope.stability.py)
args <- commandArgs(trailingOnly = TRUE)
temppath <- args[1]
prefix <- args[2]

# Importing raster maps
observ <- paste(temppath, "/observed.tif", sep = "")
sind <- paste(temppath, "/sindex.tif", sep = "")

# Generating data frames from raster maps
observed <- readGDAL(fname = observ)  #observation
sindex <- readGDAL(fname = sind)  #susceptibility index

# Preprocessing data frames
sindex$band1[sindex$band1 < 0] <- NaN
observed$band1[observed$band1 > 0] <- 1
observed$band1[observed$band1 < 0] <- NaN
sindex$band1[is.na(observed$band1) == TRUE] <- NaN
observed$band1[is.na(sindex$band1) == TRUE] <- NaN

sindexf <- sindex$band1[!is.na(sindex$band1)]
observedf <- observed$band1[!is.na(observed$band1)]

# ROC curve
roc <- prediction(sindexf, observedf)  #ROC prediction
proc <- performance(roc, "tpr", "fpr")  #performance
arocv <- performance(roc, "auc")  #area under curve
thresholds <- data.frame(cut = proc@alpha.values[[1]], fpr = proc@x.values[[1]], 
    tpr = proc@y.values[[1]])  #threshold levels

# Formatting area under curve value for display
aroc <- vector("expression", 1)
aroc[1] <- substitute(expression(italic(A)[ROC] == aroca), list(aroca = format(round(as.numeric(arocv@y.values), 
    digits = 3), nsmall = 3)))[2]

# Selection of threshold levels for display according to true positives (y-axis)
# in order to avoid overlap of labels
len <- length(thresholds[, 1])  #number of levels
ctrlval <- 0.01  #initial control value
i2 <- 1  #initial counter for new vector
thrvect <- vector()  #threshold levels
fprvect <- vector()  #false positive rate
tprvect <- vector()  #true positive rate
for (i in 1:len) {
    # loop over all levels:
    tprvectc <- thresholds[i, 3]  #true positive rate
    if (tprvectc >= ctrlval) {
        # if true positive rate is higher than control:
        thrvect[i2] <- thresholds[i, 1]  #writing values to new vector
        fprvect[i2] <- thresholds[i, 2]
        tprvect[i2] <- thresholds[i, 3]
        ctrlval <- ctrlval + 0.07  #new control value
        i2 <- i2 + 1  #updating counter for new vector
    }
}

# Creating plot file
rocplot = paste(prefix, "_results/", prefix, "_plots/", prefix, "_plot_roc.png", 
    sep = "")
png(filename = rocplot, width = 12, height = 10, units = "cm", res = 300)

# Defining margins
par(mar = c(3.2, 3.2, 0.5, 0.5))

# Plotting ROC and data points
plot(x = fprvect, y = tprvect, axes = F, xlab = NA, ylab = NA, pch = 16, col = "green", 
    xlim = c(-0.05, 1), ylim = c(0, 1.02))  #threshold levels as points
plot(proc, add = TRUE, axes = F, xlab = NA, ylab = NA, lwd = 3, col = "blue")  #ROC curve
clip(0, 1, 0, 1)  #constraining drawing area
abline(a = 0, b = 1, col = "brown", lty = 2, lwd = 1)  #random prediction line
clip(-1, 2, -1, 2)  #removing constraint of drawing area
points(x = fprvect, y = tprvect, pch = 16, col = "green")  #again, threshold levels as points (to show them in the foreground)
text(x = fprvect - 0.035, y = tprvect + 0.03, labels = format(round(thrvect, 2), 
    nsmall = 1), cex = 0.75, col = "green")  #labels of threshold levels

# Plotting bounding box and axes
box()
axis(side = 1, tck = -0.02, labels = NA)
axis(side = 2, tck = -0.02, labels = NA)
axis(side = 1, lwd = 0, line = -0.4)
axis(side = 2, lwd = 0, line = -0.4)

# Plotting axis labels
mtext(side = 1, "rFP/rON", line = 1.9)
mtext(side = 2, "rTP/rOP", line = 2)

# Legend
legend("bottomright", pch = c(NA, NA, 16), lty = c(1, 3, NA), lwd = c(3, 1, NA), 
    col = c("blue", "brown", "green"), legend = c(aroc[1], "Random prediction", "Selected threshold levels"), 
    text.col = c("blue", "black", "black"), bty = "n")

# Closing plot file
dev.off()

write(round(as.numeric(arocv@y.values), 3), file = "sumtable_aucroc.txt", append = TRUE)
