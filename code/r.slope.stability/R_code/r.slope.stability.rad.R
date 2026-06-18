#!/usr/bin/R

##############################################################################
#
# MODULE:       r.slope.stability.rad.R
#
# AUTHOR:       Martin Mergili
#
# PURPOSE:      The slope stability model - R script for the preparation
#               of a radar chart illustrating true positve, true negative,
#               false positive and false negative model predictions
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
library(fmsb)

# Importing arguments defined in r.slope.stability.py
args <- commandArgs(trailingOnly = TRUE)
temppath <- args[1]
prefix <- args[2]

# Defining data file name and creating data frame
intable <- paste(temppath, "/radtable.txt", sep = "")
radtable <- read.table(intable, row.names = c("Max", "Min", "Model", "Perfect", "Average"), 
    col.names = c("TP", "TN", "FN", "FP"))

# Extracting confirmation parameters
tp <- radtable[3, 1]  #true positive rate
tn <- radtable[3, 2]  #true negative rate
fn <- radtable[3, 3]  #false negative rate
fp <- radtable[3, 4]  #false positive rate
op <- radtable[1, 1]  #observed positive rate
on <- radtable[1, 2]  #observed negative rate

# Creating plot file
radplot <- paste(prefix, "_results/", prefix, "_plots/", prefix, "_plot_rad.png", 
    sep = "")
png(filename = radplot, width = 15, height = 12, units = "cm", res = 300)

# Defining margins
par(mar = c(0, 0, 0, 0))

# Radar chart
radarchart(radtable, axistype = 0, pcol = c("blue", "black", "brown"), plty = c(1, 
    1, 5), plwd = c(3, 1, 1), pdensity = c(0, 0, 0), pfcol = c(NA, NA, NA), maxmin = TRUE, 
    centerzero = TRUE, cglcol = 200, vlabels = c(" ", " ", " ", " "))

# Plotting tp, tn, fn and pf points in different colour
points(x = c(0, -tn/on, 0, fp/on), y = c(-fn/op, 0, tp/op, 0), pch = 16, col = "green")

# Plotting legend
legend("bottomleft", legend = c("Model prediction", "Random prediction", "Perfect prediction"), 
    lty = c(1, 3, 1), lwd = c(3, 1, 1), col = c("blue", "brown", "black"), bty = "n", 
    pch = c(NA, 16, 16))  #main legend
points(x = -1.513, y = -0.965, pch = 16, col = "green")  #additional point (different colour)

# Plotting labels
text(x = c(0, -1.3, 0, 1.3), y = c(1.2, 0.05, -1.1, 0.05), label = c(paste("rTP: ", 
    format(round(tp, 1), nsmall = 1), "%", sep = ""), paste("rTN: ", format(round(tn, 
    1), nsmall = 1), "%", sep = ""), paste("rFN: ", format(round(fn, 1), nsmall = 1), 
    "%", sep = ""), paste("rFP: ", format(round(fp, 1), nsmall = 1), "%", sep = "")), 
    col = "green", font = 1)  #model
text(x = c(0, -1.3, 0, 1.3), y = c(1.1, -0.05, -1.2, -0.05), label = c(paste("rOP: ", 
    format(round(op, 1), nsmall = 1), "%", sep = ""), paste("rON: ", format(round(on, 
    1), nsmall = 1), "%", sep = ""), paste("rOP: ", format(round(op, 1), nsmall = 1), 
    "%", sep = ""), paste("rON: ", format(round(on, 1), nsmall = 1), "%", sep = "")))  #observation

# Plotting additional labels characterizing the prediction quality
text(x = -0.6, y = 0.6, label = paste("successful", "prediction", sep = "\n"), srt = 45)
text(x = 0.6, y = 0.6, label = paste("conservative", "prediction", sep = "\n"), srt = 315)
text(x = 0.6, y = -0.6, label = paste("failed", "prediction", sep = "\n"), srt = 45)
text(x = -0.6, y = -0.6, label = paste("non-conservative", "prediction", sep = "\n"), 
    srt = 315)

# Closing plot file
dev.off()
