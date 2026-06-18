#!/usr/bin/R

##############################################################################
#
# MODULE:       r.slope.stability.reg.R
#
# AUTHOR:       Martin Mergili
#
# PURPOSE:      The slope stability model - R script for the preparation
#               of a scatterplot with linear regression relating model
#               results and observations on the basis of slope units
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
temppath <- args[1]
prefix <- args[2]
pflag <- args[3]

# Defining data file name
intable <- paste(temppath, "/regtable.txt", sep = "")

# Creating vectors from data
observation <- scan(intable, nlines = 1)  #observation
model <- scan(intable, skip = 1, nlines = 1)  #model
weight <- scan(intable, skip = 2, nlines = 1)  #weight for regression
whsize <- scan(intable, skip = 3, nlines = 1)  #weight for point size

# Linear regression
lmodel <- lm(formula = model ~ observation, weights = weight)

# R-squared and p values of the regression
lsummary <- summary(lmodel)
r2 <- lsummary$adj.r.squared
pval <- lsummary$coefficients[2, 4]
r2d <- vector("expression", 2)

# Formatting R-squared and p values for line legend
r2d[1] <- substitute(expression(italic(R)^2 == val1), list(val1 = format(r2, digits = 2)))[2]
r2d[2] <- substitute(expression(italic(p) == val2), list(val2 = format(pval, digits = 2)))[2]

# Formatting labels for point legend
val3a <- 0.1
val3b <- 0.05
val3c <- 0.01
val3 <- vector("expression", 3)
val3[1] <- substitute(expression(italic(A)[u] == val3d ~ km^2), list(val3d = format(val3a, 
    digits = 2)))[2]
val3[2] <- substitute(expression(italic(A)[u] == val3e ~ km^2), list(val3e = format(val3b, 
    digits = 2)))[2]
val3[3] <- substitute(expression(italic(A)[u] == val3f ~ km^2), list(val3f = format(val3c, 
    digits = 2)))[2]

# Creating plot file
regplot <- paste(prefix, "_results/", prefix, "_plots/", prefix, "_plot_reg.png", 
    sep = "")
png(filename = regplot, width = 12, height = 10, units = "cm", res = 300)

# Defining margins
par(mar = c(3.2, 3.2, 0.5, 0.5))

# Plotting data points
plot(observation, model, cex = whsize, axes = F, xlab = NA, ylab = NA, xlim = c(0, 
    1), ylim = c(-0.14, 1.14), col = "green")

# Plotting lines
clip(0, 1, 0, 1)  #constraining start and end of line
abline(a = 0, b = 1, lty = 3, lwd = 1, col = "brown")  #perfect fit
abline(lmodel, col = "blue", lty = 1, lwd = 3)  #regression
clip(-1, 2, -1, 2)  #removing constraint of drawing area

# Plotting legend
legend("top", legend = val3, pt.cex = c(1.581, 1.118, 0.5), pch = 1, bty = "n", col = "green", 
    horiz = TRUE, inset = c(0, -0.015))  #data points
legend("bottomleft", legend = r2d[2], lty = 1, lwd = 3, col = "blue", text.col = "blue", 
    bty = "n", horiz = TRUE, inset = c(0, -0.025))  #regression, p value
legend("bottomleft", legend = r2d[1], text.col = "blue", bty = "n", horiz = TRUE, 
    inset = c(0.31, -0.015))  #regression, R-squared value
legend("bottomleft", legend = "Perfect fit", lty = 3, lwd = 1, col = "brown", bty = "n", 
    horiz = TRUE, inset = c(0.62, -0.02))  #perfect fit

# Plotting bounding box and axes
box()
axis(side = 1, tck = -0.02, labels = NA)
axis(side = 2, tck = -0.02, labels = NA)
axis(side = 1, lwd = 0, line = -0.4)
axis(side = 2, lwd = 0, line = -0.4)

# Plotting axis labels
mtext(side = 1, "Fraction of observed landslide area", line = 2)
if (pflag == "0") {
    mtext(side = 2, "Fraction of area with factor of safety < 1", line = 2)
} else {
    mtext(side = 2, "Average of slope failure probability", line = 2)
}

# Closing plot file
dev.off()
