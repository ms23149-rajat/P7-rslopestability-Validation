#!/usr/bin/env python3

##############################################################################
#
# MODULE:       r.slope.stability.py
#
# AUTHORS:      Martin Mergili, Massimiliano Alvioli, Ivan Marchesini
# CONTRIBUTORS: Wolfgang Fellin, Fausto Guzzetti, Markus Metz,
#               Markus Neteler, and Mauro Rossi
#
# PURPOSE:      The slope stability model
#
# COPYRIGHT:    (c) 2013 - 2021 by the author
#               (c) 2020 - 2021 by the University of Graz
#               (c) 2013 - 2021 by the BOKU University, Vienna
#               (c) 2015 - 2020 by the University of Vienna
#               (c) 2013 - 2021 by the CNR-IRPI Perugia
#               (c) 2000 - 2021 by the GRASS Development Team
#
#               This program is free software under the GNU General Public
#               License (>=v2). Read the file COPYING that comes with GRASS
#               for details.
#
##############################################################################

#%module
#% description: The slope stability model
#% keywords: raster
#% keywords: slope stability
#% keywords: factor of safety
#% keywords: limit equilibrium model
#% keywords: siding surface
#% keywords: slope failure probability
#% keywords: slope stability
#%end

#%flag
#% key: a
#% description: Additional output rasters
#% guisection: flags
#%end

#%flag
#% key: d
#% description: Documentation
#% guisection: flags
#%end

#%flag
#% key: m
#% description: Multi-core processing
#% guisection: flags
#%end

#%flag
#% key: p
#% description: Slope failure probability
#% guisection: flags
#%end

#%flag
#% key: s
#% description: Excessive use of segment library
#% guisection: flags
#%end

#%flag
#% key: t
#% description: Truncation of ellipsoids
#% guisection: flags
#%end

#%flag
#% key: v
#% description: empirical confirmation and visualization
#% guisection: flags
#%end

#%flag
#% key: z
#% description: Test mode for multi-core processing (results from tiles are not put together)
#% guisection: flags
#%end

#%option
#% key: prefix
#% type: string
#% description: Prefix for output files and folders
#% required: yes
#% multiple: no
#%end

#%option
#% key: cellsize
#% type: string
#% description: Pixel size (-c, -i, -l)
#% required: no
#% multiple: no
#%end

#%option
#% key: adm
#% type: string
#% description: Area of detail map (N, S, W, E; -v)
#% required: no
#% multiple: yes
#%end

#%option
#% key: model
#% type: string
#% description: Model to be used
#% required: no
#% multiple: no
#%end

#%option
#% key: elldens
#% type: string
#% description: Ellipsoid density (-c, -l)
#% required: no
#% multiple: no
#%end

#%option
#% key: observed
#% type: string
#% gisprompt: old,raster,cell
#% description: Name of observed landslides raster map (-v)
#% required: no
#% multiple: no
#%end

#%option
#% key: elevation
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of elevation raster map (-c, -i, -l)
#% required: no
#% multiple: no
#%end

#%option
#% key: soilclass
#% type: string
#% gisprompt: old,raster,cell
#% description: Name of soil class raster map (-c, -i, -l)
#% required: no
#% multiple: no
#%end

#%option
#% key: gwdepth
#% type: string
#% gisprompt: old,raster,cell
#% description: Depth of groundwater table raster map (-c, -i, -l)
#% required: no
#% multiple: no
#%end

#%option
#% key: pseudostatic
#% type: string
#% description: Name of peak ground acceleration raster map, seisimic coefficient (-c, -i, -l)
#% required: no
#% multiple: yes
#%end

#%option
#% key: newmark
#% type: string
#% description: Name of distance from epicentre raster map, moment magnitude (-c, -i, -l)
#% required: no
#% multiple: yes
#%end

#%option
#% key: newmarkcoef
#% type: string
#% description: Coefficients for Newmark displacement (-c, -i, -l)
#% required: no
#% multiple: yes
#%end

#%option
#% key: newmarkref
#% type: string
#% gisprompt: old,raster,cell
#% description: Name of critical displacement raster map (-c, -i, -l)
#% required: no
#% multiple: no
#%end

#%option
#% key: arias
#% type: string
#% description: parameters for computing Arias intensity (-c, -i, -l)
#% required: no
#% multiple: yes
#%end

#%option
#% key: slopeunits
#% type: string
#% gisprompt: old,raster,cell
#% description: Name of slope units raster map (-v)
#% required: no
#% multiple: no
#%end

#%option
#% key: dxl
#% type: string
#% gisprompt: old,raster,cell
#% description: Name of average dx slope of layers raster map (-l)
#% required: no
#% multiple: no
#%end

#%option
#% key: dyl
#% type: string
#% gisprompt: old,raster,cell
#% description: Name of average dy slope of layers raster map (-l)
#% required: no
#% multiple: no
#%end

#%option
#% key: numlayers
#% type: string
#% description: Number of layers for each soil class, comma separated (-c, -i)
#% required: no
#% multiple: yes
#%end

#%option
#% key: depthmaps
#% type: string
#% gisprompt: old,raster,cell
#% description: Names of depth of layer bottom raster maps, comma separated (-c, -i)
#% required: no
#% multiple: yes
#%end

#%option
#% key: depthvals
#% type: string
#% description: Depth of layer bottom values, comma separated (-c, -i)
#% required: no
#% multiple: yes
#%end

#%option
#% key: depthstats
#% type: string
#% description: Statistical properties of top layer depth (-p)
#% required: no
#% multiple: yes
#%end

#%option
#% key: prelayers
#% type: string
#% description: Prefix for depth of layer bottom raster maps (-l)
#% required: no
#% multiple: no
#%end

#%option
#% key: classlayers
#% type: string
#% description: Soil class to be assigned to each layer, comma separated (-l)
#% required: no
#% multiple: yes
#%end

#%option
#% key: maxlayers
#% type: string
#% description: Maximum number of layers relevant for one ellipsoid (-l)
#% required: no
#% multiple: no
#%end

#%option
#% key: geotech
#% type: string
#% description: Geotechnical parameters for each soil class and/or layer, comma separated (-c, -i, -l)
#% required: no
#% multiple: yes
#%end

#%option
#% key: seepage
#% type: string
#% description: Seepage direction (-c, -i, -l)
#% required: no
#% multiple: no
#%end

#%option
#% key: ellips
#% type: string
#% description: Ellipsoid parameters, comma separated (-c, -l)
#% required: no
#% multiple: yes
#%end

#%option
#% key: sampling
#% type: string
#% description: Sampling parameters, comma separated (-p)
#% required: no
#% multiple: yes
#%end

#%option
#% key: regpar
#% type: string
#% description: Regression parameters, comma separated (-p)
#% required: no
#% multiple: yes
#%end

#%option
#% key: segsize
#% type: string
#% description: Size of one segment in rows or columns (-c, -i, -l)
#% required: no
#% multiple: no
#%end

#%option
#% key: nsegs
#% type: string
#% description: Number of segments to be held in memory (-c, -i, -l)
#% required: no
#% multiple: no
#%end

#%option
#% key: tilesx
#% type: string
#% description: Number of tiles in x direction (-m)
#% required: no
#% multiple: no
#%end

#%option
#% key: tilesy
#% type: string
#% description: Number of tiles in y direction (-m)
#% required: no
#% multiple: no
#%end

#%option
#% key: overlap
#% type: string
#% description: Overlap of tiles (-m, -v)
#% required: no
#% multiple: no
#%end

#%option
#% key: ncores
#% type: string
#% description: Number of cores to use (-m)
#% required: no
#% multiple: no
#%end


# Importing libraries

import grass.script as grass
from grass.script import core as grasscore
import os
import queue
import re
import shutil
import subprocess
import sys
import threading
import time


# Defining fundamental variables

exitFlag = 0

queueLock = threading.Lock()
workQueue = queue.Queue()

ambvars = grass.gisenv()  # path to grass data base
locpath = ambvars.GISDBASE + "/" + ambvars.LOCATION_NAME  # path to grass location
mapsetpath = (
    ambvars.GISDBASE + "/" + ambvars.LOCATION_NAME + "/" + ambvars.MAPSET
)  # path to grass mapset directory
temppath = "TEMP/rtemp"  # path to temporary directory
scriptpath = "$HOME/.grass7/addons/etc/r.slope.stability.rcode"  # path to R scripts
scriptpath2 = "$HOME/.grass7/addons/scripts"  # path to r.slope.stability.mult


# Function for paths


def which(file):
    for path in os.environ["PATH"].split(":"):
        if os.path.exists(path + "/" + file):
            return path + "/" + file
    return None


# GRASS binary

for i in range(70, 80):
    if not which("grass" + str(i)) is None:
        bingrass = which("grass" + str(i))


# Function for error message


def ErrorMessage(specify):  # function for error message:
    grass.message(" ")
    grass.error("Please revise the " + specify + ".")
    grass.message(" ")
    sys.exit()


# Class for threading


class myThread(threading.Thread):
    def __init__(self, name, q):
        threading.Thread.__init__(self)
        self.name = name
        self.q = q

    def run(self):
        process_data(self.name, self.q)


# Function for start of multi-core processing


def StartBatch(jid):
    print("Processing tile no. %s" % jid)
    esegui = "bash " + "TEMP/rtemp" + str(jid) + "/tmp" + str(jid) + "/batch" + str(jid)
    os.system(esegui + " < /dev/null > " + "TEMP/rtemp" + str(jid) + "/out" + str(jid))
    print("Tile no. %s completed." % jid)
    return


# Function supporting start of multi-core processing


def process_data(threadName, q):
    while not exitFlag:
        queueLock.acquire()
        if not workQueue.empty():
            time.sleep(0.25)
            data = q.get()
            queueLock.release()
            StartBatch(data)
        else:
            queueLock.release()


# Calling main function


def main():

    # Setting flags and options

    aflag = flags["a"]
    dflag = flags["d"]
    mflag = flags["m"]
    pflag = flags["p"]
    sflag = flags["s"]
    tflag = flags["t"]
    vflag = flags["v"]
    zflag = flags["z"]
    pf = options["prefix"]
    cellsize = options["cellsize"]
    adm = options["adm"]
    model = options["model"]
    elldens = options["elldens"]
    observed = options["observed"]
    elevation = options["elevation"]
    pseudostatic = options["pseudostatic"]
    newmark = options["newmark"]
    newmarkcoef = options["newmarkcoef"]
    newmarkref = options["newmarkref"]
    arias = options["arias"]
    soilclass = options["soilclass"]
    gwdepth = options["gwdepth"]
    slopeunits = options["slopeunits"]
    dxl = options["dxl"]
    dyl = options["dyl"]
    numlayers = options["numlayers"]
    depthmaps = options["depthmaps"]
    depthvals = options["depthvals"]
    depthstats = options["depthstats"]
    prelayers = options["prelayers"]
    classlayers = options["classlayers"]
    maxlayers = options["maxlayers"]
    geotech = options["geotech"]
    seepage = options["seepage"]
    ellips = options["ellips"]
    sampling = options["sampling"]
    regpar = options["regpar"]
    segsize = options["segsize"]
    nsegs = options["nsegs"]
    tilesx0 = options["tilesx"]
    tilesy0 = options["tilesy"]
    overlap0 = options["overlap"]
    ncores0 = options["ncores"]

    # Making sure that all necessary parameters are available

    iflag = None
    cflag = None
    lflag = None
    rflag = None

    if model:

        if model == "i":  # infinite slope stability model
            iflag = "1"
        elif model == "c":  # soil class mode with Hovland model
            cflag = "1"
        elif model == "l":  # soil class mode with revised Hovland model
            lflag = "1"
        elif model == "cr":  # layer mode with Hovland model
            cflag = "1"
            rflag = "1"
        elif model == "lr":  # layer mode with revised Hovland model
            lflag = "1"
            rflag = "1"
        else:
            ErrorMessage("valid model (option model)")

        if (
            elldens == "0"
        ):  # no multi-core processing or slope failure probability, but documentation if only one ellipsoid is analyzed
            dflag = 1
            mflag = None
            pflag = None
        if (
            iflag
        ):  # no additional output, documentation or truncation for infinite slope stability mode
            aflag = None
            dflag = None
            tflag = None
        if pflag:  # no documentation with slope failure probability
            dflag = None
        if zflag:  # no documentation or confirmation and visualization for test mode
            dflag = None
            vflag = None

        if not pf:
            pf = "slst"
        if not seepage:
            seepage = "1"
        if not elldens:
            elldens = "2000"
        if not sampling:
            sampling = "2,4,4,0,1,1,1"
        if not regpar:
            regpar = "-583,21711"
        if not maxlayers:
            maxlayers = "-9999"
        if not segsize:
            segsize = "16"
        if not nsegs:
            nsegs = "16"
        if not ncores0:
            ncores0 = "8"

        if not elevation:
            ErrorMessage("elevation raster map (option elevation)")
        if not soilclass:
            ErrorMessage("soil class raster map (option soilclass)")
        if not geotech:
            ErrorMessage("geotechnical parameters (option geotech)")

        if (
            model == "c" or model == "cr" or model == "l" or model == "lr"
        ):  # ellipsoid modes
            if not ellips:
                ErrorMessage("ellipsoid parameters (option ellips)")
            if vflag:  # empirical confirmation:
                if observed and not overlap0:
                    ErrorMessage("width of edge (option overlap)")

        if model == "c" or model == "cr":  # soil class mode
            if not numlayers:
                ErrorMessage("number of layers (option numlayers)")
            if not depthmaps and not depthvals:
                ErrorMessage(
                    "maps or values of depth of layer bottoms (option depthmaps or depthvals)"
                )

        if model == "l" or model == "lr":  # layer mode
            if not prelayers:
                ErrorMessage("prefix for layer bottom maps (option prelayers)")
            if not classlayers:
                ErrorMessage(
                    "geotechnical classes to be assigned to each layer (option classlayers)"
                )

        if mflag:  # multi-core processing
            if not tilesx0:
                ErrorMessage("number of tiles in x direction (option tilesx)")
            if not tilesy0:
                ErrorMessage("number of tiles in y direction (option tilesy)")
            if not overlap0:
                ErrorMessage("overlap between tiles (option overlap)")

        if not cellsize:
            cellsize = grass.raster_info(elevation)["nsres"]  # cell size

    print()
    print()
    print(
        "##############################################################################"
    )
    print("#")
    print("# MODULE:       r.slope.stability 2.0 Pre-release")
    print("#")
    print("# AUTHORS:      Martin Mergili, Massimiliano Alvioli, Ivan Marchesini")
    print("# CONTRIBUTORS: Wolfgang Fellin, Fausto Guzzetti, Markus Metz,")
    print("#               Markus Neteler, and Mauro Rossi")
    print("#")
    print("# PURPOSE:      The slope stability model")
    print("#")
    print("# COPYRIGHT:    (c) 2009 - 2021 by the author")
    print("#               (c) 2020 - 2021 by the University of Graz")
    print("#               (c) 2009 - 2021 by the BOKU University, Vienna")
    print("#               (c) 2015 - 2020 by the University of Vienna")
    print("#               (c) 2013 - 2021 by the CNR-IRPI Perugia")
    print("#               (c) 2000 - 2021 by the GRASS Development Team")
    print("#               (C) 1993 - 2021 by the R Development Core Team")
    print("#")
    print("# VERSION:      20210715 (15 July 2021)")
    print("#")
    print("# PREFIX:       %s" % pf)
    print("#")
    print("#               This program is free software under the GNU General Public")
    print("#               License (>=v2). Read the file COPYING that comes with GRASS")
    print("#               for details.")
    print("#")
    print(
        "##############################################################################"
    )
    print()
    print()
    print("1. PREPARING ENVIRONMENT.")

    # Preparing environment

    jidstr = "0"

    os.system("rm -rf " + "TEMP")  # if necessary, deleting temporary directory
    os.system("mkdir " + "TEMP")  # creating temporary directory
    os.system("mkdir " + "TEMP/rtemp" + jidstr)  # creating temporary directory

    vonly = (
        0
    )  # setting control for confirmation mode only to 0 (negative, initial value)
    d_numlayers = 2  # setting initial value for number of layers

    if pseudostatic:  # if pseudostatic seismic slope stability analysis is applied:
        pseudostatic = list(map(str, pseudostatic.split(",")))  # splitting input

    if newmark:  # if Newmark seismic slope stability analysis is applied:

        newmark = list(map(str, newmark.split(",")))

        if not newmarkcoef:
            newmarkcoef = "0.847,-10.62,6.587,1.84,0.295"
        newmarkcoef = list(map(str, newmarkcoef.split(",")))

        if not arias:
            arias = "4"
        parias = ["0", "0", "0", "0", "0"]
        arias = list(map(str, arias.split(",")))
        for i in range(0, 5):
            for j in arias:
                if str(i + 1) == j:
                    parias[i] = "1"

    if pflag:
        length_all = (
            12
        )  # length of list with geotechnical parameters per layer and soil class
    else:
        length_all = 6

    if mflag:  # converting parameters for multi-core processing to numbers:
        tilesx = int(tilesx0)
        tilesy = int(tilesy0)
        ncores = int(ncores0)

    if overlap0:
        overlap = float(overlap0)  # converting overlap string to number
    elif iflag:
        overlap = 2 * float(
            cellsize
        )  # default for overlap with infinite slope stability model

    if cflag or iflag:  # for soil class mode or infinite slope stability mode:

        if depthmaps:  # if layer depths are defined by maps:
            depthmaps = list(map(str, depthmaps.split(",")))  # splitting input
            j = 0  # initializing counter

            for i in depthmaps:  # loop over all parts of input:
                j = j + 1  # updating counter
                grass.run_command(
                    "g.copy", rast=(i, "xdepth" + str(j)), overwrite=True, quiet=True
                )  # temporary layer depth maps

    if cflag or iflag or lflag:  # for slope stability model:
        p1file = open(
            "TEMP/rtemp" + jidstr + "/param1.txt", "w"
        )  # opening parameter files for writing
        p2file = open("TEMP/rtemp" + jidstr + "/param2.txt", "w")
        p3file = open("TEMP/rtemp" + jidstr + "/param3.txt", "w")
        p4file = open("TEMP/rtemp" + jidstr + "/cellsize.txt", "w")

        # Writing param2.txt file

        print("param2", file=p2file)  # writing header to file

        if cflag or iflag:  # for soil class or infinite slope stability mode:
            numlayers = list(map(str, numlayers.split(",")))  # splitting input
            imaxxx = 0
            for i in numlayers:
                print(
                    i, file=p2file
                )  # writing number of layers in each soil class to file
                if int(i) > imaxxx:
                    imaxxx = int(i)

        else:  # for layer mode:
            geotech_gamma = ["99999"]  # initializing lists of geotechnical parameters
            geotech_c = ["99999"]
            geotech_phi = ["89"]
            geotech_theta = ["0"]
            geotech_cstdev = ["0"]
            geotech_cmin = ["99999"]
            geotech_cmax = ["99999"]
            geotech_phistdev = ["0"]
            geotech_phimin = ["89"]
            geotech_phimax = ["89"]
            geotech = list(map(str, geotech.split(",")))  # splitting input
            numcat = (len(geotech) + 1) / (length_all - 1)  # number of soil classes
            j = 0  # initializing counter
            for i in geotech:  # loop over all parts of input:
                j = j + 1  # updating counter
                if j == length_all:
                    j = 1  # if end of class is reached, go to next class
                if j == 2:
                    geotech_gamma.append(i)  # soil specific weight
                if j == 3:
                    geotech_c.append(i)  # cohesion
                if j == 4:
                    geotech_phi.append(i)  # friction angle
                if j == 5:
                    geotech_theta.append(i)  # saturated water content
                if pflag:
                    if j == 6:
                        geotech_cstdev.append(i)  # cohesion, standard deviation
                    if j == 7:
                        geotech_cmin.append(i)  # cohesion, minimum
                    if j == 8:
                        geotech_cmax.append(i)  # cohesion, maximum
                    if j == 9:
                        geotech_phistdev.append(
                            i
                        )  # angle of internal friction, standard deviation
                    if j == 10:
                        geotech_phimin.append(i)  # angle of internal friction, minimum
                    if j == 11:
                        geotech_phimax.append(i)  # angle of internal friction, maximum
                if j > 1:
                    print(i, file=p2file)  # writing parameters to file

        # Writing param3.txt file

        print("param3", file=p3file)  # writing header to file

        if cflag or iflag:  # for soil class or infinite slope stability mode:
            if depthvals:
                depthvals = list(map(str, depthvals.split(",")))  # splitting input
            geotech = list(map(str, geotech.split(",")))  # splitting input
            length = len(geotech) + 1  # reading length of dataset
            k = 0  # initializing counter for geotechnical parameters
            n = 0  # initializing counter for layer bottom depth values
            for (
                i
            ) in (
                geotech
            ):  # loop over geotechnical parameter and layer bottom depth values:
                if k not in list(range(0, length, length_all)) and k not in list(
                    range(1, length, length_all)
                ):  # for geotechnical parameter values:
                    print(i, file=p3file)  # writing parameter values to file
                k = k + 1  # increasing counter for geotechnical parameters
                if k in range(
                    length_all, length, length_all
                ):  # for layer bottom depth values:
                    if depthmaps:
                        print(
                            -9999, file=p3file
                        )  # if depth map is provided, writing no data value to file
                    else:
                        print(
                            depthvals[n], file=p3file
                        )  # else, writing depth of layer bottom to file
                        n = n + 1  # increasing counter for layers

        else:  # for layer mode:
            classlayers = list(map(str, classlayers.split(",")))  # splitting input
            j = 0  # initializing counter
            for i in classlayers:  # loop over all parts of input:
                j = j + 1  # updating counter
                print(
                    geotech_gamma[int(i)], file=p3file
                )  # writing soil specific weight for layer to file
                print(
                    geotech_c[int(i)], file=p3file
                )  # writing cohesion for layer to file
                print(
                    geotech_phi[int(i)], file=p3file
                )  # writing friction angle for layer to file
                print(
                    geotech_theta[int(i)], file=p3file
                )  # writing saturated water content for layer to file

                if pflag:
                    print(
                        geotech_cstdev[int(i)], file=p3file
                    )  # writing cohesion, standard deviation for layer to file
                    print(
                        geotech_cmin[int(i)], file=p3file
                    )  # writing cohesion, minimum for layer to file
                    print(
                        geotech_cmax[int(i)], file=p3file
                    )  # writing cohesion, maximum for layer to file
                    print(
                        geotech_phistdev[int(i)], file=p3file
                    )  # writing friction angle, standard deviation for layer to file
                    print(
                        geotech_phimin[int(i)], file=p3file
                    )  # writing friction angle, minimum for layer to file
                    print(
                        geotech_phimax[int(i)], file=p3file
                    )  # writing friction angle, maximum for layer to file

                grass.run_command(
                    "g.copy",
                    rast=(prelayers + str(j), "xdepth" + str(j)),
                    overwrite=True,
                    quiet=True,
                )  # temporary layer depth maps
            d_numlayers = j  # number of layers

        # Writing param1.txt file

        print("param1", file=p1file)  # writing header to file
        if cflag or iflag:
            maxcat = int(
                grass.raster_info(soilclass)["max"]
            )  # reading number of soil classes
        print("TEMP/rtemp" + jidstr, file=p1file)  # writing temporary directory to file
        print(ambvars.MAPSET, file=p1file)  # writing mapset to file
        if cflag:
            print("1", file=p1file)  # writing control for soil class mode to file
        else:
            print("0", file=p1file)
        if iflag:
            print(
                "1", file=p1file
            )  # writing control for infinite slope stability mode to file
        else:
            print("0", file=p1file)
        if lflag:
            print("1", file=p1file)  # writing control for layer mode to file
        else:
            print("0", file=p1file)
        if rflag:
            print("1", file=p1file)  # writing control for revised Hovland model
        else:
            print("0", file=p1file)
        if tflag:
            print(
                "1", file=p1file
            )  # writing control for mode with truncated ellipsoids to file
        else:
            print("0", file=p1file)
        if aflag:
            print("1", file=p1file)  # writing control for writing additional output
        else:
            print("0", file=p1file)
        print(seepage, file=p1file)  # writing control for seepage direction to file

        if pseudostatic:
            print(
                "1", file=p1file
            )  # writing control for pseudostatic seismic slope stability analysis
        else:
            print("0", file=p1file)
        if newmark:
            print(
                "1", file=p1file
            )  # writing control for Newmark seismic slope stability analysis
        else:
            print("0", file=p1file)
        if newmarkref:
            print("1", file=p1file)  # writing control for critical displacement map
        else:
            print("0", file=p1file)

        if pseudostatic:
            print(pseudostatic[1], file=p1file)  # writing seismic coefficient
        if newmark:
            print(newmark[1], file=p1file)  # writing earthquake depth
            print(newmark[2], file=p1file)  # writing moment magnitude
            print(newmark[3], file=p1file)  # writing parameter S1
            print(newmark[4], file=p1file)  # writing parameter S2
            for i in range(0, 5):
                print(
                    newmarkcoef[i], file=p1file
                )  # writing coefficients for Newmark displacement
            for i in range(0, 5):
                print(parias[i], file=p1file)  # writing arias equations

        if dflag:
            print("1", file=p1file)  # writing control for documentation to file
        else:
            print("0", file=p1file)

        if cflag or iflag:
            print(maxcat, file=p1file)  # writing number of soil classes to file
            print(
                int(imaxxx), file=p1file
            )  # writing temporary value for number of layers to file
            # !!!BE CAREFUL THIS IS EXPERIMENTAL AND MAY CAUSE AN ERROR

        else:
            print(numcat, file=p1file)  # writing number of soil classes to file
            print(j, file=p1file)  # writing number of layers to file
            print(
                maxlayers, file=p1file
            )  # writing maximum number of layers for one ellipsoid to file

        if sflag:
            print(
                "1", file=p1file
            )  # writing control for excessive use of segment library to file
        else:
            print("0", file=p1file)
        if pflag:
            print(
                "1", file=p1file
            )  # writing control for slope failure probability to file
        else:
            print("0", file=p1file)

        if pflag:
            sampling = list(map(str, sampling.split(",")))  # splitting input
            for isteps in sampling:
                print(isteps, file=p1file)  # writing sampling parameters to file

            if sampling[4] == "4":
                regpar = list(map(str, regpar.split(",")))  # splitting input
                for iregpar in regpar:
                    print(iregpar, file=p1file)  # writing regression parameters to file

        if (cflag or iflag) and pflag and (tflag or iflag):
            if int(sampling[3]) > 1:
                depthstats = list(map(str, depthstats.split(",")))  # splitting input
                for istats in depthstats:
                    print(
                        istats, file=p1file
                    )  # writing statistical properties of top layer to file

        print(segsize, file=p1file)  # writing segment size to file
        print(nsegs, file=p1file)  # writing number of segments to file

        # Completing param3.txt file

        if cflag or lflag:  # for modes with ellipsoids:
            ellpar = list(map(str, ellips.split(",")))  # ellipsoid parameters
            for i in ellpar:
                print(i, file=p3file)  # writing ellipsoid parameters to file
            print(elldens, file=p3file)  # writing ellipsoid density to file

        print(cellsize, file=p4file)  # writing pixel size to file

        p1file.close()  # closing parameter files
        p2file.close()
        p3file.close()
        p4file.close()

    # Setting and storing region

    if cflag or iflag or lflag:
        grass.run_command(
            "g.region", flags="a", res=cellsize
        )  # updating raster cell size

    if mflag or vflag:  # if mode involves extent change, storing values

        region = grasscore.region()  # reading region of study area
        e0 = region["e"]  # reading E boundary of study area
        w0 = region["w"]  # reading W boundary of study area
        n0 = region["n"]  # reading N boundary of study area
        s0 = region["s"]  # reading S boundary of study area

    # Preparing temporary raster maps for slope stability model

    if cflag or iflag or lflag:  # for modes with slope stability model:

        grass.run_command(
            "g.copy", rast=(elevation, "xelev"), overwrite=True, quiet=True
        )  # temporary elevation map
        grass.run_command(
            "g.copy", rast=(soilclass, "xsoilclass"), overwrite=True, quiet=True
        )  # temporary soil class map
        if gwdepth:
            grass.run_command(
                "g.copy", rast=(gwdepth, "xgwdepth"), overwrite=True, quiet=True
            )  # temporary groundwater table map
        else:
            grass.mapcalc('"xgwdepth" = 0', overwrite=True, quiet=True)

        if pseudostatic:
            grass.run_command(
                "g.copy", rast=(pseudostatic[0], "xpga"), overwrite=True, quiet=True
            )  # temporary peak ground acceleration map
        if newmark:
            grass.run_command(
                "g.copy", rast=(newmark[0], "xdistance"), overwrite=True, quiet=True
            )  # temporary distance from epicentre map
        if newmarkref:
            grass.run_command(
                "g.copy", rast=(newmarkref, "xnewmarkref"), overwrite=True, quiet=True
            )  # temporary critical displacement map

    if lflag and seepage == "2":  # for layer mode with layer-parallel seepage:

        print("%i soil classes and %i layers." % (numcat, j))

        if dxl and dyl:  # if average slope of layers is given, reading input rasters:

            grass.run_command(
                "g.copy", rast=(dxl, "xdzdyl"), overwrite=True, quiet=True
            )  # average slope in x direction
            grass.run_command(
                "g.copy", rast=(dyl, "xdzdxl"), overwrite=True, quiet=True
            )  # average slope in y direction

        else:  # if not given, computing average slope of layers:
            grass.run_command(
                "r.slope.aspect",
                elevation="xelev",
                dx="xdx0",
                dy="xdy0",
                quiet=True,
                overwrite=True,
            )
            grass.mapcalc('"sumdxl"=0', overwrite=True, quiet=True)
            grass.mapcalc('"sumdyl"=0', overwrite=True, quiet=True)
            grass.mapcalc('"numdxl"=0', overwrite=True, quiet=True)
            grass.mapcalc('"numdyl"=0', overwrite=True, quiet=True)

            for m in range(1, j + 1):
                grass.mapcalc(
                    '"xelevl"=if("xdepth%s">0,"xelev"-"xdepth%s",1/0)' % (m, m),
                    overwrite=True,
                    quiet=True,
                )
                grass.run_command(
                    "r.slope.aspect",
                    elevation="xelevl",
                    dx="dxl",
                    dy="dyl",
                    quiet=True,
                    overwrite=True,
                )
                grass.mapcalc('"dxtestl"=isnull("dxl")', overwrite=True, quiet=True)
                grass.mapcalc('"dytestl"=isnull("dyl")', overwrite=True, quiet=True)
                grass.mapcalc(
                    '"sumdxl"=if("dxtestl"==0,"sumdxl"+"dxl","sumdxl")',
                    overwrite=True,
                    quiet=True,
                )
                grass.mapcalc(
                    '"sumdyl"=if("dytestl"==0,"sumdyl"+"dyl","sumdyl")',
                    overwrite=True,
                    quiet=True,
                )
                grass.mapcalc(
                    '"numdxl"=if("dxtestl"==0,"numdxl"+1,"numdxl")',
                    overwrite=True,
                    quiet=True,
                )
                grass.mapcalc(
                    '"numdyl"=if("dytestl"==0,"numdyl"+1,"numdyl")',
                    overwrite=True,
                    quiet=True,
                )
                grass.run_command(
                    "g.remove", flags="f", type="rast", name="dxl,dyl", quiet=True
                )

                sys.stdout.write("Computing slope for layer %s of %s ...\r" % (m, j))
                sys.stdout.flush()

            grass.mapcalc(
                '"xdxl0"=if("numdxl"==0,"xdx0","sumdxl"/"numdxl")',
                overwrite=True,
                quiet=True,
            )
            grass.mapcalc(
                '"xdyl0"=if("numdyl"==0,"xdy0","sumdyl"/"numdyl")',
                overwrite=True,
                quiet=True,
            )
            grass.mapcalc('xdzdx="xdy0"', overwrite=True, quiet=True)
            grass.mapcalc('xdzdy="xdx0"', overwrite=True, quiet=True)
            grass.mapcalc('xdzdxl="xdyl0"', overwrite=True, quiet=True)
            grass.mapcalc('xdzdyl="xdxl0"', overwrite=True, quiet=True)
            grass.run_command(
                "g.remove",
                flags="f",
                type="rast",
                name="xdx0,xdy0,xdxl0,xdyl0,xelevl,dxtestl,dytestl,sumdxl,sumdyl,numdxl,numdyl",
                quiet=True,
            )

            grass.run_command(
                "g.copy", rast=("xdzdxl", "%s_dyl" % pf), overwrite=True, quiet=True
            )  # average slope in x direction
            grass.run_command(
                "g.copy", rast=("xdzdyl", "%s_dxl" % pf), overwrite=True, quiet=True
            )  # average slope in y direction

            print("Computing slopes ... completed.                  ")

        print()

    else:  # if average slope of layers is not required, setting to no data value

        grass.mapcalc('"xdzdxl"=-9999', overwrite=True, quiet=True)
        grass.mapcalc('"xdzdyl"=-9999', overwrite=True, quiet=True)

    if cflag or iflag or lflag:  # mode with slope stability model:

        print("2. EVALUATING SLOPE STABILITY.")

        start = time.time()  # storing time (start of main computation)

        # Single-core processing

        if not mflag:  # if the study area is computed as one entity:

            os.environ["XINT"] = "0"  # exporting id
            grass.run_command(
                "r.slope.stability.main"
            )  # executing slope stability model

            if aflag:
                grass.mapcalc(
                    '"r_sindex"=float(round(10000*"r_cumfos_crit"/"r_cumsurf"))/10000',
                    overwrite=True,
                    quiet=True,
                )
            # computing susceptibility index

            # Making result raster maps permanent

            if iflag or not (
                elldens == "0"
            ):  # for ellipsoid modes with more than one ellipsoid or infinite slope stability mode:
                grass.run_command(
                    "g.rename",
                    rast=("r_fos", "%s_fos" % pf),
                    overwrite=True,
                    quiet=True,
                )  # factor of safety (fos)
                # if pseudostatic: grass.run_command('g.rename', rast=('r_fos_pseudo','%s_fos_pseudo' %pf), overwrite=True, quiet=True) #pseudostatic factor of safety (fos)
                if newmark:
                    grass.run_command(
                        "g.rename",
                        rast=("r_dnewmark", "%s_dnewmark" % pf),
                        overwrite=True,
                        quiet=True,
                    )  # Newmark displacement
                if pflag:
                    grass.run_command(
                        "g.rename",
                        rast=("r_pfail", "%s_pfail" % pf),
                        overwrite=True,
                        quiet=True,
                    )  # slope failure probability
                    if newmark and newmarkref:
                        grass.run_command(
                            "g.rename",
                            rast=("r_fos_stdev", "%s_dnewmark_stdev" % pf),
                            overwrite=True,
                            quiet=True,
                        )  # stdev of Newmark displacement
                    else:
                        grass.run_command(
                            "g.rename",
                            rast=("r_fos_stdev", "%s_fos_stdev" % pf),
                            overwrite=True,
                            quiet=True,
                        )  # stdev of fos

            if (cflag or lflag) and not (elldens == "0") and aflag:
                # for ellipsoid modes with more than one ellipsoid and additional output:
                grass.run_command(
                    "g.rename",
                    rast=("r_sindex", "%s_sindex" % pf),
                    overwrite=True,
                    quiet=True,
                )  # susceptibility index
                grass.run_command(
                    "g.rename",
                    rast=("r_lcrit", "%s_lcrit" % pf),
                    overwrite=True,
                    quiet=True,
                )  # length of critical ellipsoid
                grass.run_command(
                    "g.rename",
                    rast=("r_wcrit", "%s_wcrit" % pf),
                    overwrite=True,
                    quiet=True,
                )  # width of critical ellipsoid
                grass.run_command(
                    "g.rename",
                    rast=("r_zbcrit", "%s_zbcrit" % pf),
                    overwrite=True,
                    quiet=True,
                )  # zb of critical ellipsoid

            if not elldens == "0":  # for multiple ellipsoids:
                grass.run_command(
                    "g.rename",
                    rast=("r_ddepth", "%s_depth" % pf),
                    overwrite=True,
                    quiet=True,
                )  # depth of critical slip surface

            if (
                lflag and not elldens == "0" and aflag
            ):  # for layer mode with more than one ellipsoid and additional output:
                grass.run_command(
                    "g.rename",
                    rast=("r_layer", "%s_layer" % pf),
                    overwrite=True,
                    quiet=True,
                )  # id of critical layer

            if elldens == "0":  # for mode with one ellipsoid only:
                grass.run_command(
                    "g.rename",
                    rast=("r_depth", "%s_depth" % pf),
                    overwrite=True,
                    quiet=True,
                )  # depth of critical slip surface

        # Multi-core processing

        else:  # if the study area is split into a defined number of tiles:

            # Creating tiles

            ntiles = tilesy * tilesx  # number of tiles
            dx = region["e"] - region["w"]  # E-W-extent of study area
            dy = region["n"] - region["s"]  # N-S-extent of study area
            tdx = (dx + overlap * (tilesx - 1)) / tilesx  # E-W-extent of 1 tile
            tdy = (dy + overlap * (tilesy - 1)) / tilesy  # N-S-extent of 1 tile

            print()
            print(
                "Splitting study area into %s tiles of %s x %s m." % (ntiles, tdx, tdy)
            )

            jid = 0  # resetting counter for tiles

            for nx in range(
                1, tilesx + 1
            ):  # starting loop over all tiles in E-W-direction:
                for ny in range(
                    1, tilesy + 1
                ):  # starting loop over all tiles in N-S-direction:

                    jid = jid + 1  # updating counter for tiles
                    jidstr = str(jid)

                    tx1 = w0 + (nx - 1) * (tdx - overlap)  # W boundary of tile
                    ty1 = s0 + (ny - 1) * (tdy - overlap)  # S boundary of tile
                    tx2 = tx1 + tdx  # E boundary of tile
                    ty2 = ty1 + tdy  # N boundary of tile

                    tx1 = float(cellsize) * round(
                        tx1 / float(cellsize), 0
                    )  # adjusting W boundary of tile to pixel size
                    ty1 = float(cellsize) * round(
                        ty1 / float(cellsize), 0
                    )  # adjusting S boundary of tile to pixel size
                    tx2 = float(cellsize) * round(
                        tx2 / float(cellsize), 0
                    )  # adjusting E boundary of tile to pixel size
                    ty2 = float(cellsize) * round(
                        ty2 / float(cellsize), 0
                    )  # adjusting N boundary of tile to pixel size

                    grass.run_command(
                        "g.region", n=ty2, s=ty1, w=tx1, e=tx2
                    )  # setting region to extent of tile

                    print()
                    print("Writing batch file for tile %s." % jid)
                    print()
                    print("Extent adjusted to pixel size:")
                    print("West:  %s m" % tx1)
                    print("East:  %s m" % tx2)
                    print("South: %s m" % ty1)
                    print("North: %s m" % ty2)
                    print()

                    # Creating batch file

                    os.mkdir("TEMP/rtemp" + jidstr)
                    os.mkdir(
                        "TEMP/rtemp" + jidstr + "/tmp%s" % jid
                    )  # creating directory for batch file
                    strtmp = (
                        "TEMP/rtemp" + jidstr + "/tmp%s/batch" % jid
                    )  # file name for batch file
                    out = open(strtmp + str(jid), "w")  # creating batch file
                    os.system(
                        "rm -rf " + locpath + "/map%s" % jid
                    )  # removing old mapset for tile
                    grass.run_command(
                        "g.mapset", flags="c", mapset="map%s" % jid
                    )  # creating new mapset for tile
                    grass.run_command(
                        "g.mapsets", mapset=ambvars.MAPSET, operation="add"
                    )  # making original mapset active

                    grass.run_command(
                        "g.mapset", mapset=ambvars.MAPSET
                    )  # switching back to original mapset
                    os.environ["PATH"] += (
                        os.pathsep
                        + os.path.join("TEMP/rtemp" + jidstr + "/tmp%s") % jid
                    )  # adding path to batch file

                    if aflag:
                        aflag_exp = 1  # exporting controls to r.slope.stability.mult
                    else:
                        aflag_exp = 0
                    if pflag:
                        pflag_exp = 1
                    else:
                        pflag_exp = 0
                    if newmark:
                        newmark_exp = 1
                    else:
                        newmark_exp = 0

                    print(
                        """#!/bin/bash
export jid=%s
export ntiles=%s
export SHELL=\"/bin/bash\"
export west=%s
export sout=%s
export east=%s
export nort=%s
export cellsize=%s
export aflag=%s
export pflag=%s
export newmark=%s
export rtemp=%s
export GRASS_BATCH_JOB=%s/r.slope.stability.mult
rm -f ~/.grassrc6
cp %s0/*.txt TEMP/rtemp%s/
%s -text %s/map%s
unset GRASS_BATCH_JOB"""
                        % (
                            jidstr,
                            ntiles,
                            tx1,
                            ty1,
                            tx2,
                            ty2,
                            cellsize,
                            aflag_exp,
                            pflag_exp,
                            newmark_exp,
                            temppath + jidstr,
                            scriptpath2,
                            temppath,
                            jidstr,
                            bingrass,
                            locpath,
                            jidstr,
                        ),
                        file=out,
                    )  # creating batch file

                    out.close()  # closing batch file
                    fd = os.open(strtmp + str(jid), os.O_RDONLY)  # opening batch file
                    os.fchmod(fd, 0o755)  # making batch file executable
                    os.close(fd)  # closing batch file

            # Executing batch processing

            start_batch = time.time()  # storing time (start of multi-core processing)

            neff = min(ncores, ntiles)
            threadList = list(range(1, neff + 1))
            nameList = list(range(1, ntiles + 1))
            threads = []
            for tName in threadList:
                thread = myThread(tName, workQueue)
                thread.start()
                threads.append(thread)
            queueLock.acquire()
            for word in nameList:
                workQueue.put(word)
            queueLock.release()
            while not workQueue.empty():
                pass
            global exitFlag
            exitFlag = 1
            for t in threads:
                t.join()
            print()
            print("Batch processing completed.")
            print()

            stop_batch = time.time()  # storing time (stop of multi-core processing)
            comptime_batch = (
                stop_batch - start_batch
            )  # storing computing time for batch processing in seconds

            # Collecting results from all tiles

            if cflag or lflag:
                lbtmmaxmax = (
                    0
                )  # initializing overall maximum number of layers for one ellipsoid
            for jid in nameList:  # loop over all tiles:
                if not zflag:  # if mode is not test mode:

                    jidstr = str(jid)

                    if iflag or not (elldens == "0"):
                        # for ellipsoid modes with more than one ellipsoid or infinite slope stability mode:
                        grass.run_command(
                            "r.in.gdal",
                            input="TEMP/rtemp" + jidstr + "/r_fos" + str(jid),
                            output="r_fos" + str(jid),
                            overwrite=True,
                            quiet=True,
                        )  # factor of safety (fos)
                        # if pseudostatic: grass.run_command('r.in.gdal', input='TEMP/rtemp'+jidstr+'/r_fos_pseudo'+str(jid),
                        # output='r_fos_pseudo'+str(jid), overwrite=True, quiet=True) #factor of safety (fos)
                        if newmark:
                            grass.run_command(
                                "r.in.gdal",
                                input="TEMP/rtemp" + jidstr + "/r_dnewmark" + str(jid),
                                output="r_dnewmark" + str(jid),
                                overwrite=True,
                                quiet=True,
                            )  # Newmark displacement
                        if pflag:
                            grass.run_command(
                                "r.in.gdal",
                                input="TEMP/rtemp" + jidstr + "/r_pfail" + str(jid),
                                output="r_pfail" + str(jid),
                                overwrite=True,
                                quiet=True,
                            )  # slope failure probability
                            grass.run_command(
                                "r.in.gdal",
                                input="TEMP/rtemp" + jidstr + "/r_fos_stdev" + str(jid),
                                output="r_fos_stdev" + str(jid),
                                overwrite=True,
                                quiet=True,
                            )  # standard deviation of fos

                    if (cflag or lflag) and not (elldens == "0") and aflag:
                        # for ellipsoid modes with more than one ellipsoid and additional output:
                        grass.run_command(
                            "r.in.gdal",
                            input="TEMP/rtemp" + jidstr + "/r_sindex" + str(jid),
                            output="r_sindex" + str(jid),
                            overwrite=True,
                            quiet=True,
                        )  # susceptibility index
                        grass.run_command(
                            "r.in.gdal",
                            input="TEMP/rtemp" + jidstr + "/r_lcrit" + str(jid),
                            output="r_lcrit" + str(jid),
                            overwrite=True,
                            quiet=True,
                        )  # length of critical ellipsoid
                        grass.run_command(
                            "r.in.gdal",
                            input="TEMP/rtemp" + jidstr + "/r_wcrit" + str(jid),
                            output="r_wcrit" + str(jid),
                            overwrite=True,
                            quiet=True,
                        )  # width of critical ellipsoid
                        grass.run_command(
                            "r.in.gdal",
                            input="TEMP/rtemp" + jidstr + "/r_zbcrit" + str(jid),
                            output="r_zbcrit" + str(jid),
                            overwrite=True,
                            quiet=True,
                        )  # zb of critical ellipsoid

                    if not elldens == "0":
                        grass.run_command(
                            "r.in.gdal",
                            input="TEMP/rtemp" + jidstr + "/r_ddepth" + str(jid),
                            output="r_ddepth" + str(jid),
                            overwrite=True,
                            quiet=True,
                        )  # depth of critical slip surface

                    if (
                        lflag and not (elldens == "0") and aflag
                    ):  # for layer mode with more than one ellipsoid and additional output:
                        grass.run_command(
                            "r.in.gdal",
                            input="TEMP/rtemp" + jidstr + "/r_layer" + str(jid),
                            output="r_layer" + str(jid),
                            overwrite=True,
                            quiet=True,
                        )  # id of critical layer

                if cflag or lflag:  # for ellipsoid modes:

                    sumfile = open(
                        "TEMP/rtemp" + jidstr + "/summary3" + str(jid) + ".txt", "r"
                    )
                    # opening summary file for reading
                    lbtmmax = (
                        sumfile.readline()
                    )  # reading maximum number of layers for one ellipsoid
                    lbtmmax = int(lbtmmax.replace("\n", ""))
                    lbtmalw = (
                        sumfile.readline()
                    )  # reading maximum number of layers allowed
                    lbtmalw = int(lbtmalw.replace("\n", ""))
                    sumfile.close()  # closing file
                    if lbtmmax > lbtmmaxmax:
                        lbtmmaxmax = (
                            lbtmmax
                        )  # updating maximum number of layers for one ellipsoid, if necessary

            grass.run_command(
                "g.region", n=n0, s=s0, w=w0, e=e0, quiet=True
            )  # restoring original extent

            # Initializing result maps

            if not zflag:  # if mode is not test mode:

                if iflag or not (
                    elldens == "0"
                ):  # for ellipsoid modes with more than one ellipsoid or infinite slope stability mode:
                    grass.mapcalc('"%s_fos"=9999' % pf, overwrite=True, quiet=True)
                    # if pseudostatic: grass.mapcalc('"%s_fos_pseudo"=9999' %pf, overwrite=True, quiet=True)
                    if newmark:
                        grass.mapcalc(
                            '"%s_dnewmark"=-9999' % pf, overwrite=True, quiet=True
                        )
                    if pflag:
                        grass.mapcalc(
                            '"%s_pfail"=-9999' % pf, overwrite=True, quiet=True
                        )
                        if newmark and newmarkref:
                            grass.mapcalc(
                                '"%s_dnewmark_stdev"=-9999' % pf,
                                overwrite=True,
                                quiet=True,
                            )
                        else:
                            grass.mapcalc(
                                '"%s_fos_stdev"=-9999' % pf, overwrite=True, quiet=True
                            )

                if (cflag or lflag) and not (elldens == "0") and aflag:
                    # for ellipsoid modes with more than one ellipsoid and additional output:
                    grass.mapcalc('"%s_sindex"=-9999' % pf, overwrite=True, quiet=True)
                    grass.mapcalc('"%s_lcrit"=-9999' % pf, overwrite=True, quiet=True)
                    grass.mapcalc('"%s_wcrit"=-9999' % pf, overwrite=True, quiet=True)
                    grass.mapcalc('"%s_zbcrit"=-9999' % pf, overwrite=True, quiet=True)

                grass.mapcalc('"%s_depth"=-9999' % pf, overwrite=True, quiet=True)
                if (
                    lflag and not (elldens == "0") and aflag
                ):  # for layer mode with more than one ellipsoid:
                    grass.mapcalc('"%s_layer"=-9999' % pf, overwrite=True, quiet=True)

                # Superimposing results from all tiles

                for jid in range(1, ntiles + 1):  # loop over all tiles:

                    if iflag or not (elldens == "0"):
                        # for ellipsoid modes with more than one ellipsoid or infinite slope stability mode:

                        grass.mapcalc(
                            '"%s_fos"=if(isnull("r_fos%s")==1,"%s_fos",if("r_fos%s">"%s_fos","%s_fos","r_fos%s"))'
                            % (pf, str(jid), pf, str(jid), pf, pf, str(jid)),
                            overwrite=True,
                            quiet=True,
                        )
                        if newmark:
                            grass.mapcalc(
                                '"%s_dnewmark"=if(isnull("r_dnewmark%s")==1,"%s_dnewmark",if("r_dnewmark%s"<"%s_dnewmark","%s_dnewmark","r_dnewmark%s"))'
                                % (pf, str(jid), pf, str(jid), pf, pf, str(jid)),
                                overwrite=True,
                                quiet=True,
                            )
                        if pflag:
                            grass.mapcalc(
                                '"%s_pfail"=if(isnull("r_pfail%s")==1,"%s_pfail",if("r_pfail%s"<"%s_pfail","%s_pfail","r_pfail%s"))'
                                % (pf, str(jid), pf, str(jid), pf, pf, str(jid)),
                                overwrite=True,
                                quiet=True,
                            )
                            if newmark and newmarkref:
                                grass.mapcalc(
                                    '"%s_dnewmark_stdev"=if(isnull("r_fos_stdev%s")==1,"%s_dnewmark_stdev",if("r_dnewmark%s"<"%s_dnewmark","%s_dnewmark_stdev","r_fos_stdev%s"))'
                                    % (pf, str(jid), pf, str(jid), pf, pf, str(jid)),
                                    overwrite=True,
                                    quiet=True,
                                )
                            else:
                                grass.mapcalc(
                                    '"%s_fos_stdev"=if(isnull("r_fos_stdev%s")==1,"%s_fos_stdev",if("r_fos%s">"%s_fos","%s_fos_stdev","r_fos_stdev%s"))'
                                    % (pf, str(jid), pf, str(jid), pf, pf, str(jid)),
                                    overwrite=True,
                                    quiet=True,
                                )

                    if (cflag or lflag) and not (elldens == "0") and aflag:
                        # for ellipsoid modes with more than one ellipsoid and additional output:
                        grass.mapcalc(
                            '"%s_sindex"=if(isnull("r_sindex%s")==1,"%s_sindex",if("r_sindex%s"<"%s_sindex","%s_sindex","r_sindex%s"))'
                            % (pf, str(jid), pf, str(jid), pf, pf, str(jid)),
                            overwrite=True,
                            quiet=True,
                        )
                        grass.mapcalc(
                            '"%s_lcrit"=if(isnull("r_lcrit%s")==1,"%s_lcrit",if("r_fos%s">"%s_fos","%s_lcrit","r_lcrit%s"))'
                            % (pf, str(jid), pf, str(jid), pf, pf, str(jid)),
                            overwrite=True,
                            quiet=True,
                        )
                        grass.mapcalc(
                            '"%s_wcrit"=if(isnull("r_wcrit%s")==1,"%s_wcrit",if("r_fos%s">"%s_fos","%s_wcrit","r_wcrit%s"))'
                            % (pf, str(jid), pf, str(jid), pf, pf, str(jid)),
                            overwrite=True,
                            quiet=True,
                        )
                        grass.mapcalc(
                            '"%s_zbcrit"=if(isnull("r_zbcrit%s")==1,"%s_zbcrit",if("r_fos%s">"%s_fos","%s_zbcrit","r_zbcrit%s"))'
                            % (pf, str(jid), pf, str(jid), pf, pf, str(jid)),
                            overwrite=True,
                            quiet=True,
                        )

                    grass.mapcalc(
                        '"%s_depth"=if(isnull("r_ddepth%s")==1,"%s_depth",if("r_fos%s">"%s_fos","%s_depth","r_ddepth%s"))'
                        % (pf, str(jid), pf, str(jid), pf, pf, str(jid)),
                        overwrite=True,
                        quiet=True,
                    )

                    if (
                        lflag and not (elldens == "0") and aflag
                    ):  # for layer mode with more than one ellipsoid and additional output:
                        grass.mapcalc(
                            '"%s_layer"=if(isnull("r_layer%s")==1,"%s_layer",if("r_fos%s">"%s_fos","%s_layer","r_layer%s"))'
                            % (pf, str(jid), pf, str(jid), pf, pf, str(jid)),
                            overwrite=True,
                            quiet=True,
                        )

                for jid in range(1, ntiles + 1):
                    os.system(
                        "rm -rf " + locpath + "/map" + str(jid)
                    )  # removing mapsets for all tiles

        stop = time.time()  # storing time (end of main computation)
        comptime = stop - start  # storing computing time in seconds

        # Creating structure for storing results

        if os.path.exists(pf + "_results"):
            os.system(
                "rm -rf " + pf + "_results"
            )  # if result directory already exists, deleting it
        os.mkdir(pf + "_results")  # creating directory for results
        os.mkdir(pf + "_results/" + pf + "_files")  # creating directory for files
        os.mkdir(
            pf + "_results/" + pf + "_hmaps"
        )  # creating directory for maps with hillshade
        os.mkdir(pf + "_results/" + pf + "_maps")  # creating directory for maps
        os.mkdir(pf + "_results/" + pf + "_plots")  # creating directory for plots
        os.mkdir(pf + "_results/" + pf + "_tiffs")  # creating directory for tiffs

        # Exporting result raster maps to tiff files

        if not zflag:  # if mode is not test mode:

            if iflag or not (
                elldens == "0"
            ):  # for ellipsoid modes with more than one ellipsoid or infinite slope stability mode:

                grass.run_command(
                    "r.out.gdal",
                    input=pf + "_fos",
                    output=pf + "_results/" + pf + "_tiffs/" + pf + "_fos.tif",
                    format="GTiff",
                    flags="c",
                    quiet=True,
                    overwrite=True,
                )

                if newmark:
                    grass.run_command(
                        "r.out.gdal",
                        input=pf + "_dnewmark",
                        output=pf + "_results/" + pf + "_tiffs/" + pf + "_dnewmark.tif",
                        format="GTiff",
                        flags="c",
                        quiet=True,
                        overwrite=True,
                    )
                if pflag:
                    grass.run_command(
                        "r.out.gdal",
                        input=pf + "_pfail",
                        output=pf + "_results/" + pf + "_tiffs/" + pf + "_pfail.tif",
                        format="GTiff",
                        flags="c",
                        quiet=True,
                        overwrite=True,
                    )

                    if newmark and newmarkref:
                        grass.run_command(
                            "r.out.gdal",
                            input=pf + "_dnewmark_stdev",
                            output=pf
                            + "_results/"
                            + pf
                            + "_tiffs/"
                            + pf
                            + "_dnewmark_stdev.tif",
                            format="GTiff",
                            flags="c",
                            quiet=True,
                            overwrite=True,
                        )

                    else:
                        grass.run_command(
                            "r.out.gdal",
                            input=pf + "_fos_stdev",
                            output=pf
                            + "_results/"
                            + pf
                            + "_tiffs/"
                            + pf
                            + "_fos_stdev.tif",
                            format="GTiff",
                            flags="c",
                            quiet=True,
                            overwrite=True,
                        )

            if (cflag or lflag) and not (elldens == "0") and aflag:
                # for ellipsoid modes with more than one ellipsoid and additional output:
                grass.run_command(
                    "r.out.gdal",
                    input=pf + "_sindex",
                    output=pf + "_results/" + pf + "_tiffs/" + pf + "_sindex.tif",
                    format="GTiff",
                    flags="c",
                    quiet=True,
                    overwrite=True,
                )
                grass.run_command(
                    "r.out.gdal",
                    input=pf + "_lcrit",
                    output=pf + "_results/" + pf + "_tiffs/" + pf + "_lcrit.tif",
                    format="GTiff",
                    flags="c",
                    quiet=True,
                    overwrite=True,
                )
                grass.run_command(
                    "r.out.gdal",
                    input=pf + "_wcrit",
                    output=pf + "_results/" + pf + "_tiffs/" + pf + "_wcrit.tif",
                    format="GTiff",
                    flags="c",
                    quiet=True,
                    overwrite=True,
                )
                grass.run_command(
                    "r.out.gdal",
                    input=pf + "_zbcrit",
                    output=pf + "_results/" + pf + "_tiffs/" + pf + "_zbcrit.tif",
                    format="GTiff",
                    flags="c",
                    quiet=True,
                    overwrite=True,
                )

            grass.run_command(
                "r.out.gdal",
                input=pf + "_depth",
                output=pf + "_results/" + pf + "_tiffs/" + pf + "_depth.tif",
                format="GTiff",
                flags="c",
                quiet=True,
                overwrite=True,
            )

            if (
                lflag and not (elldens == "0") and aflag
            ):  # for layer mode with more than one ellipsoid and additional output:
                grass.run_command(
                    "r.out.gdal",
                    input=pf + "_layer",
                    output=pf + "_results/" + pf + "_tiffs/" + pf + "_layer.tif",
                    format="GTiff",
                    flags="c",
                    quiet=True,
                    overwrite=True,
                )

            print("Exporting result maps ... completed.")

        # Writing table for evolution plot

        if (cflag or lflag) and not elldens == "0" and not zflag:
            # for ellipsoid-based modes with more than one ellipsoid:

            if mflag:
                jidmin = 1
                nmaxloop = ntiles + 1  # number of loops + 1 for multi-core processing
            else:
                jidmin = 0
                nmaxloop = 1  # number of loops + 1 for single-core processing
            evolfile_linesmin = 99999999
            # initial value of number of lines of evolution file
            evoltable = open(
                pf + "_results/" + pf + "_files/" + pf + "_evolution.txt", "w"
            )  # opening evolution table file for writing

            for jid in range(jidmin, nmaxloop):  # loop over all tiles
                evolfile = open(
                    "TEMP/rtemp" + str(jid) + "/evolution" + str(jid) + ".txt", "r"
                )  # opening temporary evolution file
                evolfile_lines = sum(1 for line in evolfile)  # reading number of lines
                evolfile.close()  # closing file
                evolfile = open(
                    "TEMP/rtemp" + str(jid) + "/evolution" + str(jid) + ".txt", "r"
                )  # again opening temporary evolution file
                cellallraw = evolfile.readline()  # total number of cells
                cellall = float(cellallraw.replace("\n", ""))  # removing newline
                evolfile.close()  # closing file
                if (
                    evolfile_lines > 2
                    and cellall > 0
                    and evolfile_lines < evolfile_linesmin
                ):
                    evolfile_linesmin = evolfile_lines
                # updating minimum number of lines

            status_max = 0.0
            # initial value of maximum status related to all pixels
            jidval = 0  # initial value for number of valid tiles
            for jid in range(jidmin, nmaxloop):  # loop over all tiles

                evolfile = open(
                    "TEMP/rtemp" + str(jid) + "/evolution" + str(jid) + ".txt", "r"
                )  # opening temporary evolution file
                evolfile_lines = sum(1 for line in evolfile)  # reading number of lines
                evolfile.close()  # closing file
                evolfile = open(
                    "TEMP/rtemp" + str(jid) + "/evolution" + str(jid) + ".txt", "r"
                )  # opening temporary evolution file again

                cellallraw = evolfile.readline()  # total number of cells
                cellall = float(cellallraw.replace("\n", ""))  # removing newline

                if (
                    evolfile_lines > 2 and cellall > 0
                ):  # if tile is valid (contains data)
                    jidval = jidval + 1  # updating number of valid tiles

                    if jidval == 1:  # if first valid tile:
                        elldensch = [0]  # array for ellipsoid density
                        numsim = [0]  # array for number of ellipsoids
                    cellayer_abs = [
                        0
                    ]  # array for status (number of pixels with fos smaller than 1 or sum of pfail)
                    cellayer = [
                        0
                    ]  # array for status (ratio of pixels with fos smaller than 1 or sum of pfail)
                    cellchange = [0]  # array for change of status
                    cellayer_pre = 0  # initializing status of previous step

                    for i in range(
                        1, evolfile_linesmin
                    ):  # loop over all lines (steps):
                        changeall = evolfile.readline()  # reading first line
                        changeall = changeall.replace("\n", "")  # removing newline
                        change = list(
                            map(str, changeall.split("\t"))
                        )  # splitting input
                        if jidval == 1:
                            elldensch.append(
                                change[0]
                            )  # appending ellipsoid density to array
                            numsim.append(
                                change[1]
                            )  # appending number of tested ellipsoids to array
                        cellayer_abs.append(change[2])  # appending status to array
                        cellayer_ratio = (
                            float(change[2]) / cellall
                        )  # status related to all pixels
                        cellayer.append(
                            cellayer_ratio
                        )  # appending status related to all pixels to array
                        cellchange_new = float(cellayer_ratio) - float(
                            cellayer_pre
                        )  # change of status related to all pixels
                        cellchange.append(
                            cellchange_new
                        )  # append change of status related to all pixels to array
                        cellayer_pre = (
                            cellayer_ratio
                        )  # store status related to all pixels for use in next step
                        if cellayer_ratio > status_max:
                            status_max = (
                                cellayer_ratio
                            )  # updating maximum status related to all pixels, if necessary

                    evolfile.close()  # closing temporary evolution file

                    if jidval == 1:  # if first valid tile:
                        for i in range(1, evolfile_linesmin):  # loop over all lines:
                            print(
                                "%s\t" % elldensch[i], end=" ", file=evoltable
                            )  # writing ellipsoid density to file
                        print("\n", end=" ", file=evoltable)
                        for i in range(1, evolfile_linesmin):  # loop over all lines:
                            print(
                                "%s\t" % numsim[i], end=" ", file=evoltable
                            )  # writing number of ellipsoids to file
                        print("\n", end=" ", file=evoltable)
                    for i in range(1, evolfile_linesmin):  # loop over all lines:
                        print(
                            "%i\t" % cellall, end=" ", file=evoltable
                        )  # writing status to file
                    print("\n", end=" ", file=evoltable)
                    for i in range(1, evolfile_linesmin):  # loop over all lines:
                        print(
                            "%s\t" % cellayer_abs[i], end=" ", file=evoltable
                        )  # writing status to file
                    print("\n", end=" ", file=evoltable)
                    for i in range(1, evolfile_linesmin):  # loop over all lines:
                        print(
                            "%s\t" % cellayer[i], end=" ", file=evoltable
                        )  # writing status related to all pixels to file
                    print("\n", end=" ", file=evoltable)
                    for i in range(1, evolfile_linesmin):  # loop over all lines:
                        print(
                            "%s\t" % cellchange[i], end=" ", file=evoltable
                        )  # writing change of status to file
                    print("\n", end=" ", file=evoltable)

            evoltable.close()  # closing evolution table file

        # Writing table for profile plot

        if (cflag or lflag) and dflag:  # for ellipsoid modes with documentation:

            if mflag:  # making summary file(s) permanent:
                for jid in range(1, ntiles + 1):
                    shutil.copyfile(
                        "TEMP/rtemp0" + "/summary1" + str(jid) + ".txt",
                        pf
                        + "_results/"
                        + pf
                        + "_files/"
                        + pf
                        + "_summary"
                        + str(jid)
                        + ".txt",
                    )
            else:
                shutil.copyfile(
                    "TEMP/rtemp0" + "/summary10.txt",
                    pf + "_results/" + pf + "_files/" + pf + "_summary.txt",
                )

            if elldens == "0":  # for ellipsoid modes with only one ellipsoid:
                shutil.copyfile(
                    "TEMP/rtemp0" + "/summary20.txt",
                    pf + "_results/" + pf + "_files/" + pf + "_profile.txt",
                )
                # making profile file permanent
                sumfile = open(
                    "TEMP/rtemp0" + "/summary30.txt", "r"
                )  # opening summary file for reading
                sumfile.readline()  # reading first line (not needed here)
                sumfile.readline()  # reading second line (not needed here)
                prof_numlayers = sumfile.readline()
                prof_numlayers = int(
                    prof_numlayers.replace("\n", "")
                )  # extracting number of layers
                prof_slipcrit = sumfile.readline()
                prof_slipcrit = int(
                    prof_slipcrit.replace("\n", "")
                )  # extracting id of critical layer
                prof_foscrit = sumfile.readline()
                prof_foscrit = float(
                    prof_foscrit.replace("\n", "")
                )  # extracting critical factor of safety or slope failure probability
                sumfile.close()  # closing file

        # Reading maximum number of layers for one ellipsoid

        if (
            cflag or lflag
        ) and not mflag:  # for single-core processing ellipsoid modes:
            sumfile = open(
                "TEMP/rtemp0" + "/summary30.txt", "r"
            )  # opening summary file for reading
            lbtmmaxmax = sumfile.readline()  # reading maximum number of layers
            lbtmmaxmax = int(lbtmmaxmax.replace("\n", ""))
            lbtmalw = sumfile.readline()  # reading maximum number of layers allowed
            lbtmalw = int(lbtmalw.replace("\n", ""))
            sumfile.close()  # closing file

        # Writing documentation file

        if (
            not vflag or (model)
        ) and not zflag:  # if mode does not include confirmation only and is not test mode:
            documentation = open(
                pf + "_results/" + pf + "_files/" + pf + "_documentation.txt", "w"
            )  # opening documentation file for writing

            if aflag:
                print("a\t1", file=documentation)
            else:
                print("a\t0", file=documentation)
            if dflag:
                print("d\t1", file=documentation)
            else:
                print("d\t0", file=documentation)
            if mflag:
                print("m\t1", file=documentation)
            else:
                print("m\t0", file=documentation)
            if pflag:
                print("p\t1", file=documentation)
            else:
                print("p\t0", file=documentation)
            if sflag:
                print("s\t1", file=documentation)
            else:
                print("s\t0", file=documentation)
            if tflag:
                print("t\t1", file=documentation)
            else:
                print("t\t0", file=documentation)
            print("cellsize\t%s" % cellsize, file=documentation)

            if adm:  # extracting area of detail map
                adm = list(map(str, adm.split(",")))
                admN = adm[0]
                admS = adm[1]
                admW = adm[2]
                admE = adm[3]
                print("admN\t%s" % admN, file=documentation)
                print("admS\t%s" % admS, file=documentation)
                print("admW\t%s" % admW, file=documentation)
                print("admE\t%s" % admE, file=documentation)
            else:
                print("admN\tNA", file=documentation)
                print("admS\tNA", file=documentation)
                print("admW\tNA", file=documentation)
                print("admE\tNA", file=documentation)

            if model:
                print("model\t%s" % model, file=documentation)
            else:
                print("model\tNA", file=documentation)
            if elldens:
                print("elldens\t%s" % elldens, file=documentation)
            else:
                print("elldens\tNA", file=documentation)
            if observed:
                print("observed\t%s" % observed, file=documentation)
            else:
                print("observed\tNA", file=documentation)
            print("elevation\t%s" % elevation, file=documentation)
            if pseudostatic:
                print("pseudostatic\t1", file=documentation)
            else:
                print("pseudostatic\t0", file=documentation)
            if newmark:
                print("newmark\t1", file=documentation)
            else:
                print("newmark\t0", file=documentation)
            if newmarkcoef:
                print("newmarkcoef\t%s" % newmarkcoef, file=documentation)
            else:
                print("newmarkcoef\tNA", file=documentation)
            if newmarkref:
                print("newmarkref\t%s" % newmarkref, file=documentation)
            else:
                print("newmarkref\tNA", file=documentation)
            if arias:
                print("arias\t%s" % arias, file=documentation)
            else:
                print("arias\tNA", file=documentation)
            print("soilclass\t%s" % soilclass, file=documentation)
            if gwdepth:
                print("gwdepth\t%s" % gwdepth, file=documentation)
            else:
                print("gwdepth\tNA", file=documentation)
            if slopeunits:
                print("slopeunits\t%s" % slopeunits, file=documentation)
            else:
                print("slopeunits\tNA", file=documentation)
            if dxl:
                print("dxl\t%s" % dxl, file=documentation)
            else:
                print("dxl\tNA", file=documentation)
            if dyl:
                print("dyl\t%s" % dyl, file=documentation)
            else:
                print("dyl\tNA", file=documentation)
            if numlayers:
                print("numlayers\t%s" % numlayers, file=documentation)
            else:
                print("numlayers\tNA", file=documentation)
            if depthmaps:
                print("depthmaps\t%s" % depthmaps, file=documentation)
            else:
                print("depthmaps\tNA", file=documentation)
            if depthvals:
                print("depthvals\t%s" % depthvals, file=documentation)
            else:
                print("depthvals\tNA", file=documentation)
            if prelayers:
                print("prelayers\t%s" % prelayers, file=documentation)
            else:
                print("prelayers\tNA", file=documentation)
            if classlayers:
                print("classlayers\t%s" % classlayers, file=documentation)
            else:
                print("classlayers\tNA", file=documentation)
            if maxlayers:
                print("maxlayers\t%s" % maxlayers, file=documentation)
            else:
                print("maxlayers\tNA", file=documentation)
            print("geotech\t%s" % geotech, file=documentation)
            print("seepage\t%s" % seepage, file=documentation)
            if ellips:
                print("ellips\t%s" % ellips, file=documentation)
            else:
                print("ellips\tNA", file=documentation)
            if sampling:
                print("sampling\t%s" % sampling, file=documentation)
            else:
                print("sampling\tNA", file=documentation)
            if regpar:
                print("regpar\t%s" % regpar, file=documentation)
            else:
                print("regpar\tNA", file=documentation)
            print("segsize\t%s" % segsize, file=documentation)
            print("nsegs\t%s" % nsegs, file=documentation)
            if tilesx0:
                print("tilesx\t%s" % tilesx0, file=documentation)
            else:
                print("tilesx\tNA", file=documentation)
            if tilesy0:
                print("tilesy\t%s" % tilesy0, file=documentation)
            else:
                print("tilesy\tNA", file=documentation)
            if overlap0:
                print("overlap\t%s" % overlap0, file=documentation)
            elif iflag:
                print("overlap\t%s" % str(2 * float(cellsize)), file=documentation)
            else:           
                print("overlap\tNA", file=documentation)
            if ncores0:
                print("ncores\t%s" % ncores0, file=documentation)
            else:
                print("ncores\tNA", file=documentation)

            print("number of layers\t%i" % d_numlayers, file=documentation)
            if (not model == "i") and mflag:
                print("number of valid tiles\t%i" % jidval, file=documentation)
            else:
                print("number of valid tiles\t0", file=documentation)

            if (not model == "i") and elldens == "0":
                print("number of layers\t%i" % prof_numlayers, file=documentation)
                print("id of critical layer\t%i" % prof_slipcrit, file=documentation)
                print(
                    "critical factor of safety\t%.4f" % prof_foscrit, file=documentation
                )

            if (not model == "i") and not elldens == "0":
                print(
                    "maximum y value of evolution plot\t%.5f" % status_max,
                    file=documentation,
                )

            documentation.close()  # closing documentation file

    # confirmation and visualization

    if vflag:

        print()
        print("3. CONFIRMING AND VISUALIZING MODEL RESULTS.")
        print()

        if not model:  # if mode only includes confirmation and visualization:

            # Reading documentation file for flags and options

            documentation = open(
                pf + "_results/" + pf + "_files/" + pf + "_documentation.txt", "r"
            )  # opening documentation file

            aflag = documentation.readline()
            aflag = aflag.replace("a\t", "").replace("\n", "")
            documentation.readline()
            mflag = documentation.readline()
            mflag = mflag.replace("m\t", "").replace("\n", "")
            pflag = documentation.readline()
            pflag = pflag.replace("p\t", "").replace("\n", "")
            documentation.readline()
            documentation.readline()
            cellsize = documentation.readline()
            cellsize = cellsize.replace("cellsize\t", "").replace("\n", "")
            admN = documentation.readline()
            admN = admN.replace("admN\t", "").replace("\n", "")
            admS = documentation.readline()
            admS = admS.replace("admS\t", "").replace("\n", "")
            admW = documentation.readline()
            admW = admW.replace("admW\t", "").replace("\n", "")
            admE = documentation.readline()
            admE = admE.replace("admE\t", "").replace("\n", "")
            model = documentation.readline()
            model = model.replace("model\t", "").replace("\n", "")
            elldens = documentation.readline()
            elldens = elldens.replace("elldens\t", "").replace("\n", "")
            observed = documentation.readline()
            observed = observed.replace("observed\t", "").replace("\n", "")
            elevation = documentation.readline()
            elevation = elevation.replace("elevation\t", "").replace("\n", "")
            pseudostatic = documentation.readline()
            pseudostatic = pseudostatic.replace("pseudostatic\t", "").replace("\n", "")
            newmark = documentation.readline()
            newmark = newmark.replace("newmark\t", "").replace("\n", "")
            documentation.readline()
            documentation.readline()
            documentation.readline()
            soilclass = documentation.readline()
            soilclass = soilclass.replace("soilclass\t", "").replace("\n", "")
            documentation.readline()
            slopeunits = documentation.readline()
            slopeunits = slopeunits.replace("slopeunits\t", "").replace("\n", "")
            documentation.readline()
            documentation.readline()
            documentation.readline()
            documentation.readline()
            documentation.readline()
            documentation.readline()
            documentation.readline()
            documentation.readline()
            geotech = documentation.readline()
            geotech = geotech.replace("geotech\t", "").replace("\n", "")
            geotech = list(map(str, geotech.split(",")))
            documentation.readline()
            ellips = documentation.readline()
            ellips = ellips.replace("ellips\t", "").replace("\n", "")
            documentation.readline()
            documentation.readline()
            documentation.readline()
            documentation.readline()
            documentation.readline()
            documentation.readline()
            overlap = documentation.readline()
            overlap = overlap.replace("overlap\t", "").replace("\n", "")
            documentation.readline()

            iflag = None
            cflag = None
            lflag = None

            if model == "i":  # infinite slope stability model
                iflag = "1"
            elif model == "c":  # soil class mode with Hovland model
                cflag = "1"
            elif model == "l":  # soil class mode with revised Hovland model
                lflag = "1"
            elif model == "cr":  # layer mode with Hovland model
                cflag = "1"
            elif model == "lr":  # layer mode with revised Hovland model
                lflag = "1"

            if aflag == "0":
                aflag = None
            if mflag == "0":
                mflag = None
            if pflag == "0":
                pflag = None
            if admN == "NA":
                adm = None
            else:
                adm = "1"
            if elldens == "NA":
                elldens = None
            if observed == "NA":
                observed = None
            if pseudostatic == "0":
                pseudostatic = None
            if newmark == "0":
                newmark = None
            if slopeunits == "NA":
                slopeunits = None
            if overlap == "NA":
                overlap = None

            d_numlayers = documentation.readline()
            d_numlayers = d_numlayers.replace("number of layers\t", "").replace(
                "\n", ""
            )
            d_numlayers = int(d_numlayers)
            jidval = documentation.readline()
            jidval = jidval.replace("number of valid tiles\t", "").replace("\n", "")
            if (cflag or lflag) and mflag:
                jidval = int(jidval)

            if (cflag or lflag) and elldens == "0":
                prof_numlayers = documentation.readline()
                prof_numlayers = int(
                    prof_numlayers.replace("number of layers\t", "").replace("\n", "")
                )
                prof_slipcrit = documentation.readline()
                prof_slipcrit = int(
                    prof_slipcrit.replace("id of critical layer\t", "").replace(
                        "\n", ""
                    )
                )
                prof_foscrit = documentation.readline()
                prof_foscrit = float(
                    prof_foscrit.replace("critical factor of safety\t", "").replace(
                        "\n", ""
                    )
                )

            if (cflag or lflag) and not elldens == "0":
                status_max = documentation.readline()
                status_max = float(
                    status_max.replace(
                        "maximum y value of evolution plot\t", ""
                    ).replace("\n", "")
                )

            documentation.close()  # closing documentation file

            grass.run_command("g.region", res=cellsize)  # updating pixel size
            grass.run_command(
                "g.copy", rast=(elevation, "xelev"), overwrite=True, quiet=True
            )  # temporary elevation map
            grass.run_command(
                "g.copy", rast=(soilclass, "xsoilclass"), overwrite=True, quiet=True
            )  # temporary soil class map

            vonly = 1  # setting control for confirmation mode only to 1 (positive)

        # Preparing spatial data sets for confirmation and visualization

        if observed and overlap and not (elldens == "0"):
            # if observation raster map and overlap were specified and more than one ellipsoid was tested:

            xvalmin = w0 + float(overlap)  # boundaries for confirmation
            yvalmin = s0 + float(overlap)
            xvalmax = e0 - float(overlap)
            yvalmax = n0 - float(overlap)

            grass.mapcalc(
                '"xsoilclass"=if("xsoilclass">-1,"xsoilclass", 1/0)',
                overwrite=True,
                quiet=True,
            )  # cleaned soil class map
            grass.mapcalc(
                '"xsoilclassinv"=isnull("xsoilclass")', overwrite=True, quiet=True
            )  # inverted soil class map
            grass.run_command(
                "r.buffer",
                flags="z",
                input="xsoilclassinv",
                output="xconfirmation000",
                distances=overlap,
                quiet=True,
                overwrite=True,
            )
            grass.mapcalc(
                '"xconfirmation00"=isnull("xconfirmation000")',
                overwrite=True,
                quiet=True,
            )
            grass.mapcalc(
                '"xconfirmation0"=if("xconfirmation00">0,1,1/0)',
                overwrite=True,
                quiet=True,
            )
            grass.mapcalc(
                '"xconfirmation"=if(x()<%s||x()>%s||y()<%s||y()>%s,1/0,"xconfirmation0")'
                % (str(xvalmin), str(xvalmax), str(yvalmin), str(yvalmax)),
                overwrite=True,
                quiet=True,
            )  # confirmation area map

            grass.run_command(
                "r.to.vect",
                input="xconfirmation",
                output="d_confirmation",
                type="area",
                flags="v" + "s",
                quiet=True,
                overwrite=True,
            )  # temporary confirmation area map as vector
            grass.run_command(
                "v.out.ogr",
                flags="c",
                input="d_confirmation",
                output="TEMP/rtemp0/d_confirmation.shp",
                format="ESRI_Shapefile",
                type="boundary",
                quiet=True,
                overwrite=True,
            )  # temporary confirmation area map as shapefile

            grass.run_command(
                "g.copy", rast=(observed, "xobserved0"), overwrite=True, quiet=True
            )
            grass.mapcalc(
                '"xobserved"="xconfirmation"*"xobserved0"', overwrite=True, quiet=True
            )  # temporary observed landslides map
            grass.run_command(
                "r.out.gdal",
                input="xobserved",
                output="TEMP/rtemp0" + "/observed.tif",
                format="GTiff",
                flags="c",
                quiet=True,
                overwrite=True,
            )  # exporting observed landslides raster map to tiff

            grass.mapcalc(
                '"xlobserved"=if("xobserved0">0,1,null())', overwrite=True, quiet=True
            )
            grass.run_command(
                "r.to.vect",
                input="xlobserved",
                output="d_lobserved",
                type="area",
                flags="v" + "s",
                quiet=True,
                overwrite=True,
            )  # temporary observed landslides map as vector
            grass.run_command(
                "v.out.ogr",
                flags="c",
                input="d_lobserved",
                output="TEMP/rtemp0/d_lobserved.shp",
                format="ESRI_Shapefile",
                type="boundary",
                quiet=True,
                overwrite=True,
            )  # temporary observed landslides map as shapefile

        grass.run_command(
            "r.relief", input="xelev", output="xhillshade", quiet=True, overwrite=True
        )  # temporary hillshade map
        grass.run_command(
            "r.out.gdal",
            input="xhillshade",
            output="TEMP/rtemp0" + "/hillshade.tif",
            format="GTiff",
            flags="c",
            quiet=True,
            overwrite=True,
        )  # hillshade as geotiff

        if pflag:
            pflag2 = "1"
        else:
            pflag2 = "0"

        # Map plots

        if aflag and (cflag or lflag):
            ellparam = list(
                map(str, ellips.split(","))
            )  # extracting ellipsoid parameters
        if not pflag:
            maxdepth = float(
                grass.raster_info("%s_depth" % pf)["max"]
            )  # reading maximum critical slip surface depth
        if newmark:
            maxdisp = float(
                grass.raster_info("%s_dnewmark" % pf)["max"]
            )  # reading maximum value of Newmark displacement
        if newmark and pflag:
            maxstdev = (
                maxdisp
            )  # defining maximum value of standard deviation of Newmark displacement

        d_maxval = []  # initializing list of maximum values for map plots
        d_minval = []  # initializing list of minimum values for map plots
        d_maxlab = 0.0  # initializing label of maximum value
        nummaps = 1  # initializing counter for number of maps

        plotmaps = []

        if (
            iflag or not elldens == "0"
        ):  # for infinite slope mode or modes with multiple ellipsoids:
            # factor of safety
            d_maxval.append(2.5)
            d_minval.append(0.5)
            nummaps = nummaps + 1
            plotmaps.append("fos")

        if not pflag:
            d_maxlab = float(
                maxdepth
            )  # maximum depth of critical slip surface for display
            d_maxval.append(float(maxdepth))
            d_minval.append(0.0)
            nummaps = nummaps + 1
            plotmaps.append("depth")

        if pflag:  # slope failure probability and standard deviation of fos
            d_maxval.append(1.0)
            d_minval.append(0.0)
            nummaps = nummaps + 1
            plotmaps.append("pfail")
            if not (newmark and newmarkref):
                plotmaps.append("fos_stdev")
                d_maxval.append(2.5)
                d_minval.append(0.0)
                nummaps = nummaps + 1

        if aflag and (cflag or lflag):  # additional maps with ellipsoid modes:
            d_maxval.append(float(ellparam[0]))  # critical slip surface length
            d_minval.append(float(ellparam[1]))
            d_maxval.append(float(ellparam[2]))  # critical slip surface width
            d_minval.append(float(ellparam[3]))
            d_maxval.append(float(ellparam[6]))  # critical offset of ellipsoid centre
            d_minval.append(float(ellparam[7]))
            nummaps = nummaps + 3
            plotmaps.append("lcrit")
            plotmaps.append("wcrit")
            plotmaps.append("zbcrit")

            if not pflag:  # for factor of safety:
                d_maxval.append(1.0)  # susceptibility index
                d_minval.append(0.0)
                nummaps = nummaps + 1
                plotmaps.append("sindex")

        if lflag and aflag:  # additional map with layer mode:
            d_maxval.append(float(d_numlayers))  # critical layer
            d_minval.append(1.0)
            nummaps = nummaps + 1
            plotmaps.append("layer")

        if newmark and (
            iflag or not elldens == "0"
        ):  # for Newmark seismic slope stability analysis:
            # Newmark displacement
            d_maxval.append(float(maxdisp))
            d_minval.append(0.0)
            nummaps = nummaps + 1
            plotmaps.append("dnewmark")

        if (
            newmark and pflag and (iflag or not elldens == "0")
        ):  # for probabilistic Newmark seismic slope stability analysis:
            # Newmark displacement - standard deviation
            d_maxval.append(float(maxstdev))
            d_minval.append(0.0)
            nummaps = nummaps + 1
            plotmaps.append("dnewmark_stdev")

        i = 0
        for d_map in plotmaps:  # loop over all map plots:

            d_interval = (
                d_maxval[i] - d_minval[i]
            ) / 6  # splitting range of values into six classes
            d_lev0 = str(d_minval[i])
            d_lev1 = str(d_minval[i] + 1 * d_interval)
            d_lev2 = str(d_minval[i] + 2 * d_interval)
            d_lev3 = str(d_minval[i] + 3 * d_interval)
            d_lev4 = str(d_minval[i] + 4 * d_interval)
            d_lev5 = str(d_minval[i] + 5 * d_interval)
            d_lev6 = str(d_maxval[i])
            i += 1

            if observed:
                ctrl_observed = (
                    "1"
                )  # control for availability of observation raster map
            else:
                ctrl_observed = "0"

            if pflag:
                pflagi = 1
            else:
                pflagi = 0

            subprocess.call(
                "Rscript %s/r.slope.stability.map.R %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s"
                % (
                    scriptpath,
                    "TEMP/rtemp0",
                    pf,
                    str(n0),
                    str(s0),
                    str(w0),
                    str(e0),
                    cellsize,
                    d_map,
                    d_lev0,
                    d_lev1,
                    d_lev2,
                    d_lev3,
                    d_lev4,
                    d_lev5,
                    d_lev6,
                    str(i),
                    "0",
                    d_maxlab,
                    "1",
                    ctrl_observed,
                    pflagi,
                ),
                shell=True,
            )  # overview maps with hillshade

            if adm:
                subprocess.call(
                    "Rscript %s/r.slope.stability.map.R %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s"
                    % (
                        scriptpath,
                        "TEMP/rtemp0",
                        pf,
                        admN,
                        admS,
                        admW,
                        admE,
                        cellsize,
                        d_map,
                        d_lev0,
                        d_lev1,
                        d_lev2,
                        d_lev3,
                        d_lev4,
                        d_lev5,
                        d_lev6,
                        str(i),
                        "1",
                        d_maxlab,
                        "1",
                        ctrl_observed,
                        pflagi,
                    ),
                    shell=True,
                )  # detail maps with hillshade

            subprocess.call(
                "Rscript %s/r.slope.stability.map.R %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s"
                % (
                    scriptpath,
                    "TEMP/rtemp0",
                    pf,
                    str(n0),
                    str(s0),
                    str(w0),
                    str(e0),
                    cellsize,
                    d_map,
                    d_lev0,
                    d_lev1,
                    d_lev2,
                    d_lev3,
                    d_lev4,
                    d_lev5,
                    d_lev6,
                    str(i),
                    "0",
                    d_maxlab,
                    "0",
                    ctrl_observed,
                    pflagi,
                ),
                shell=True,
            )  # overview maps without hillshade

            if adm:
                subprocess.call(
                    "Rscript %s/r.slope.stability.map.R %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s"
                    % (
                        scriptpath,
                        "TEMP/rtemp0",
                        pf,
                        admN,
                        admS,
                        admW,
                        admE,
                        cellsize,
                        d_map,
                        d_lev0,
                        d_lev1,
                        d_lev2,
                        d_lev3,
                        d_lev4,
                        d_lev5,
                        d_lev6,
                        str(i),
                        "1",
                        d_maxlab,
                        "0",
                        ctrl_observed,
                        pflagi,
                    ),
                    shell=True,
                )  # detail maps without hillshade

        # ROC plot

        if (
            (cflag or lflag or (iflag and pflag))
            and observed
            and overlap
            and not (elldens == "0")
            and (aflag or pflag)
        ):
            # for ellipsoid modes with observed landslides and overlap specified,
            # with slope failure probability or additional output:

            if not pflag:  # susceptibility index
                grass.mapcalc(
                    '"xsindex111"=isnull("%s_sindex")' % pf, overwrite=True, quiet=True
                )
                grass.mapcalc(
                    '"xsindex11"=if("xsindex111"==1,0,"%s_sindex")' % pf,
                    overwrite=True,
                    quiet=True,
                )
            else:  # slope failure probability
                grass.mapcalc(
                    '"xsindex111"=isnull("%s_pfail")' % pf, overwrite=True, quiet=True
                )
                grass.mapcalc(
                    '"xsindex11"=if("xsindex111"==1,0,"%s_pfail")' % pf,
                    overwrite=True,
                    quiet=True,
                )

            grass.mapcalc(
                '"xsindex1"=if("xsindex11"==-9999,1/0,"xsindex11")',
                overwrite=True,
                quiet=True,
            )
            grass.mapcalc(
                '"xsindex"="xconfirmation"*"xsindex1"', overwrite=True, quiet=True
            )
            grass.mapcalc(
                '"xsindex"=if(isnull("xobserved")==1,1/0,"xsindex")',
                overwrite=True,
                quiet=True,
            )  # temporary pfail/sindex map
            grass.mapcalc(
                '"xobserved1"=if(isnull("xsindex")==1,1/0,"xobserved")',
                overwrite=True,
                quiet=True,
            )  # temporary observation map

            grass.run_command(
                "r.out.gdal",
                input="xsindex",
                output="TEMP/rtemp0" + "/sindex.tif",
                format="GTiff",
                flags="c",
                quiet=True,
                overwrite=True,
            )  # exporting temporary sindex/pfail map to tiff
            grass.run_command(
                "r.out.gdal",
                input="xobserved1",
                output="TEMP/rtemp0" + "/observed.tif",
                format="GTiff",
                flags="c",
                quiet=True,
                overwrite=True,
            )  # exporting temporary observation map to tiff

            subprocess.call(
                "Rscript %s/r.slope.stability.roc.R %s %s --slave --quiet"
                % (scriptpath, "TEMP/rtemp0", pf),
                shell=True,
            )  # ROC plot relating sindex/pfail to observation

        # Radar chart

        if observed and overlap and not (elldens == "0") and not pflag:
            # if observation and overlap were specified and more than one ellipsoid was tested for the factor of safety:

            grass.mapcalc(
                '"xfosbin"=if("%s_fos"==9999,null(),if("%s_fos"<1,1,0))' % (pf, pf),
                overwrite=True,
                quiet=True,
            )  # binary fos map
            grass.mapcalc(
                '"xobsbin"=if("xobserved">0,1,0)', overwrite=True, quiet=True
            )  # binary observation map
            grass.mapcalc(
                '"xtp"=if("xfosbin"==1&&"xobsbin"==1,4,if("xfosbin"==0&&"xobsbin"==0,3, \
                if("xfosbin"==1&&"xobsbin"==0,2,if("xfosbin"==0&&"xobsbin"==1,1,null()))))',
                overwrite=True,
                quiet=True,
            )  # combination of fos and observation map
            grass.run_command(
                "r.stats",
                flags="c" + "n",
                separator="x",
                input="xtp",
                output="%s/confirmation.txt" % "TEMP/rtemp0",
                quiet=True,
                overwrite=True,
            )  # scores tp=4, tn=3, fp=2, fn=1

            valfile = open(
                "TEMP/rtemp0" + "/confirmation.txt", "r"
            )  # opening raw file with confirmation scores
            valfile_lines = sum(1 for line in valfile)  # reading number of lines
            valfile.close()  # closing file

            if (
                valfile_lines < 4
            ):  # if number of lines is smaller than 4 (at least one of the scores is 0):
                valfile = open(
                    "TEMP/rtemp0" + "/confirmation.txt", "r"
                )  # opening file again
                lr = valfile.readline()  # reading first line
                l1 = lr.find(
                    "1x"
                )  # checking for fn score (if not found, -1 is returned)
                if l1 > -1:
                    lr = valfile.readline()  # if found, reading next line
                l2 = lr.find("2x")  # checking for fp score
                if l2 > -1:
                    lr = valfile.readline()  # if found, reading next line
                l3 = lr.find("3x")  # checking for tn score
                if l3 > -1:
                    lr = valfile.readline()  # if found, reading next line
                l4 = lr.find("4x")  # checking for tp score
                valfile.close()  # closing file again
            else:
                l1 = (
                    1
                )  # if number of lines is 4, setting all controls to positive (any value larger than -1 would be fine)
                l2 = 1
                l3 = 1
                l4 = 1

            valfile = open(
                "TEMP/rtemp0" + "/confirmation.txt", "r"
            )  # once more opening temporary file with confirmation scores
            if l1 == -1:
                fn = float(0)  # if no fn score was identified before, setting to 0
            else:
                fn = valfile.readline()  # else, reading line
                fn = float(
                    fn.replace("1x", "").replace("\n", "")
                )  # extracting fn score
            if l2 == -1:
                fp = float(0)  # if no fp score was identified before, setting to 0
            else:
                fp = valfile.readline()  # else, reading line
                fp = float(
                    fp.replace("2x", "").replace("\n", "")
                )  # extracting fp score
            if l3 == -1:
                tn = float(0)  # if no tn score was identified before, setting to 0
            else:
                tn = valfile.readline()  # else, reading line
                tn = float(
                    tn.replace("3x", "").replace("\n", "")
                )  # extracting tn score
            if l4 == -1:
                tp = float(0)  # if no tp score was identified before, setting to 0
            else:
                tp = valfile.readline()  # else, reading line
                tp = float(
                    tp.replace("4x", "").replace("\n", "")
                )  # extracting tp score
            valfile.close()  # closing file for the last time

            total = tp + tn + fp + fn  # sum of scores
            tpr = 100 * tp / total  # tp rate in percent
            tnr = 100 * tn / total  # tn rate in percent
            fpr = 100 * fp / total  # fp ratte in percent
            fnr = 100 * fn / total  # fn rate in percent
            tlr = (
                100 * tp / (tp + fn)
            )  # rate of true positive predictions out of all positive observations
            tnlr = (
                100 * tn / (tn + fp)
            )  # rate of true negative predictions out of all negative observations
            mpr = tpr + fpr  # total rate of positive predictions
            mnr = tnr + fnr  # total rate of negative predictions
            tr = tnr + tpr  # total rate of true predictions
            fr = fnr + fpr  # total rate of false predictions
            opr = tpr + fnr  # rate of positive observations
            onr = tnr + fpr  # rate of negative observations

            print(
                "True positive rate\trTP\t%.1f%%" % tpr
            )  # printing tp, tn, fp and fn rates
            print("True negative rate\trTN\t%.1f%%" % tnr)
            print("False positive rate\trFP\t%.1f%%" % fpr)
            print("False negative rate\trFN\t%.1f%%" % fnr)
            print()

            resultfile = open(
                pf + "_results/" + pf + "_files/" + pf + "_confirmation.txt", "w"
            )  # opening file for confirmation scores for writing
            print("Criterion\tShortcut\tValue", file=resultfile)  # writing file header
            print("Rate of positive observations\trOP\t%.1f%%" % opr, file=resultfile)
            print("Rate of negative observations\trON\t%.1f%%" % onr, file=resultfile)
            print(
                "True positive rate\trTP\t%.1f%%" % tpr, file=resultfile
            )  # writing rates to file
            print("True negative rate\trTN\t%.1f%%" % tnr, file=resultfile)
            print("False positive rate\trFP\t%.1f%%" % fpr, file=resultfile)
            print("False negative rate\trFN\t%.1f%%" % fnr, file=resultfile)
            print(
                "True positive rate out of all positive observations\trTP/rOP\t%.1f%%"
                % tlr,
                file=resultfile,
            )
            print(
                "True negative rate out of all negative observations\trTN/rON\t%.1f%%"
                % tnlr,
                file=resultfile,
            )
            print(
                "Rate of total positive predictions\trMP\t%.1f%%" % mpr, file=resultfile
            )
            print(
                "Rate of total negative predictions\trMN\t%.1f%%" % mnr, file=resultfile
            )
            print("Rate of true predictions\trT\t%.1f%%" % tr, file=resultfile)
            print("Rate of false predictions\trF\t%.1f%%" % fr, file=resultfile)
            resultfile.close()  # closing file

            radiotable = open(
                "TEMP/rtemp0" + "/radtable.txt", "w"
            )  # opening temporary input file for radioplot of confirmation scores for writing
            print(
                "%.1f\t%.1f\t%.1f\t%.1f" % (opr, onr, opr, onr), file=radiotable
            )  # writing theoretical maxima to file
            print("0\t0\t0\t0", file=radiotable)  # writing theoretical minima to file
            print(
                "%.1f\t%.1f\t%.1f\t%.1f" % (tpr, tnr, fnr, fpr), file=radiotable
            )  # writing scores obtained by model to file
            print(
                "%.1f\t%.1f\t0\t0" % (opr, onr), file=radiotable
            )  # writing theoretical scores for perfect prediction to file as reference
            print(
                "%.1f\t%.1f\t%.1f\t%.1f"
                % (mpr * opr / 100, mnr * onr / 100, mnr * opr / 100, mpr * onr / 100),
                file=radiotable,
            )  # writing theoretical random scores to file as reference
            radiotable.close()  # closing file

            subprocess.call(
                "Rscript %s/r.slope.stability.rad.R %s %s --slave --quiet"
                % (scriptpath, "TEMP/rtemp0", pf),
                shell=True,
            )
            # radar chart with confirmation scores

            sumtable = open("sumtable_scores.txt", "a")
            sumtable.write(
                geotech[3]
                + "\t"
                + geotech[4]
                + "\t"
                + str(tpr)
                + "\t"
                + str(tnr)
                + "\t"
                + str(fpr)
                + "\t"
                + str(fnr)
                + "\n"
            )
            sumtable.close()

        # Evolution plot

        if (
            cflag or lflag
        ) and not elldens == "0":  # for ellipsoid modes with more than one ellipsoid:

            if mflag:
                numtiles = jidval  # number of tiles (for multi-core mode)
            else:
                numtiles = 1  # number of tiles (for single-core mode)
            subprocess.call(
                "Rscript %s/r.slope.stability.evolution.R %s %s %s %s --slave --quiet"
                % (scriptpath, pf, numtiles, pflag2, status_max),
                shell=True,
            )

        # Profile plot

        if (
            (cflag or lflag) and not pflag and elldens == "0"
        ):  # for ellipsoid-based factor of safety with only one ellipsoid:
            subprocess.call(
                "Rscript %s/r.slope.stability.profile.R %s %i %i %.3f --slave --quiet"
                % (scriptpath, pf, prof_numlayers, prof_slipcrit, prof_foscrit),
                shell=True,
            )

        # Regression and scatterplot based on slope units

        if slopeunits and observed and overlap and not elldens == "0":
            # if slope units, observation and overlap are provided and more than one ellipsoid was tested:

            grass.run_command(
                "g.copy", rast=(slopeunits, "xslopeunits"), overwrite=True, quiet=True
            )  # temporary slope units map

            if not pflag:
                grass.mapcalc(
                    '"xmodint"=if(isnull("xfosbin")==1||isnull("xobserved")==1,1/0,int("xfosbin"*1000))',
                    overwrite=True,
                    quiet=True,
                )  # binary fos integer raster map
            else:
                grass.mapcalc(
                    '"xpfail"=if("%s_pfail"==-9999,null(),"%s_pfail")' % (pf, pf),
                    overwrite=True,
                    quiet=True,
                )  # cleaned pfail map
                grass.mapcalc(
                    '"xmodint"=if(isnull("xpfail")==1||isnull("xobserved")==1,1/0,int("xpfail"*1000))',
                    overwrite=True,
                    quiet=True,
                )  # integer raster map of slope failure probability
            grass.mapcalc(
                '"xobsint"=if(isnull("xmodint")==1||isnull("xobserved")==1,1/0,int("xobserved"*1000))',
                overwrite=True,
                quiet=True,
            )  # integer raster map of observed landslides
            grass.mapcalc(
                '"xall"=if(isnull("xmodint")==1||isnull("xobserved")==1,1/0,int(1))',
                overwrite=True,
                quiet=True,
            )
            # integer raster map for cell count
            grass.mapcalc(
                '"xunint"=int("xslopeunits")', overwrite=True, quiet=True
            )  # integer raster map of slope units

            grass.run_command(
                "r.statistics",
                base="xunint",
                cover="xobsint",
                method="average",
                output="xobsstat",
                overwrite=True,
                quiet=True,
            )  # average of observed landslides map per slope unit
            grass.run_command(
                "r.statistics",
                base="xunint",
                cover="xmodint",
                method="average",
                output="xmodstat",
                overwrite=True,
                quiet=True,
            )  # average of prediction map per slope unit
            grass.run_command(
                "r.statistics",
                base="xunint",
                cover="xall",
                method="sum",
                output="xallstat",
                overwrite=True,
                quiet=True,
            )  # cell count per slope unit

            grass.run_command(
                "r.stats",
                flags="l" + "n",
                separator="x",
                input="xobsstat",
                output="TEMP/rtemp0" + "/obsstat.txt",
                quiet=True,
                overwrite=True,
            )  # file with average of observed landslides map
            grass.run_command(
                "r.stats",
                flags="l" + "n",
                separator="x",
                input="xmodstat",
                output="TEMP/rtemp0" + "/modstat.txt",
                quiet=True,
                overwrite=True,
            )  # file with average of prediction map
            grass.run_command(
                "r.stats",
                flags="l" + "n",
                separator="x",
                input="xallstat",
                output="TEMP/rtemp0" + "/allstat.txt",
                quiet=True,
                overwrite=True,
            )  # file with cell count

            valfile = open(
                "TEMP/rtemp0" + "/obsstat.txt", "r"
            )  # opening file with average of observed landslides per slope unit
            valfile_lines = sum(1 for line in valfile)  # counting number of lines
            valfile.close()  # closing file

            regtable = open(
                "TEMP/rtemp0" + "/regtable.txt", "w"
            )  # opening input table file for regression for writing

            obsfile = open(
                "TEMP/rtemp0" + "/obsstat.txt", "r"
            )  # again, opening file with average of observed landslides
            for val in range(1, valfile_lines):  # loop over all lines:
                obs = obsfile.readline()  # reading line
                obs = re.sub(
                    "[0-9]*x", "", obs
                )  # extracting desired value (second column)
                obs = (
                    float(obs.replace("\n", ""))
                ) / 1000  # removing newline and converting to float
                print(
                    "%.3f\t" % obs, end=" ", file=regtable
                )  # writing to input file for regression (first line)
            obsfile.close()  # closing file

            print("\n", end=" ", file=regtable)

            modfile = open(
                "TEMP/rtemp0" + "/modstat.txt", "r"
            )  # opening file with average of prediction map per slope unit
            for val in range(1, valfile_lines):  # loop over all lines:
                mod = modfile.readline()  # reading line
                mod = re.sub(
                    "[0-9]*x", "", mod
                )  # extracting desired value (second column)
                mod = (
                    float(mod.replace("\n", ""))
                ) / 1000  # removing newline and converting to float
                print(
                    "%.3f\t" % mod, end=" ", file=regtable
                )  # writing to input file for regression (second line)
            modfile.close()  # closing file

            print("\n", end=" ", file=regtable)

            alcmax = 0  # initializing maximum slope unit area
            allfile = open(
                "TEMP/rtemp0" + "/allstat.txt", "r"
            )  # opening file with cell count per slope unit
            for val in range(1, valfile_lines):  # loop over all lines:
                alc = allfile.readline()  # reading line
                alc = re.sub(
                    "[0-9]*x", "", alc
                )  # extracting desired value (second column)
                alc = (
                    float(alc.replace("\n", ""))
                ) / 1000  # removing newline and converting to float
                if alc > alcmax:
                    alcmax = alc  # updating maximum slope unit area, if necessary
                print(
                    "%.3f\t" % alc, end=" ", file=regtable
                )  # writing to input file for regression (third line)
            allfile.close()  # closing file

            print("\n", end=" ", file=regtable)

            whtfile = open(
                "TEMP/rtemp0" + "/allstat.txt", "r"
            )  # again opening file with cell count per slope unit
            for val in range(1, valfile_lines):  # loop over all lines:
                wht = whtfile.readline()  # reading line
                wht = re.sub(
                    "[0-9]*x", "", wht
                )  # extracting desired value (second column)
                wht = (
                    float(wht.replace("\n", ""))
                ) / 1000  # removing newline and converting to float
                # wht1 = 10 * pow(wht / alcmax, 0.5)  # relative slope unit area
                wht2 = 5 * pow(
                    wht * float(cellsize) * float(cellsize) / 1000, 0.5
                )  # absolute slope unit area
                print(
                    "%.3f\t" % wht2, end=" ", file=regtable
                )  # writing to input file for regression (fourth line)
            whtfile.close()  # closing file

            print("", file=regtable)
            regtable.close()  # closing input table file for regression

            subprocess.call(
                "Rscript %s/r.slope.stability.reg.R %s %s %s --slave --quiet"
                % (scriptpath, "TEMP/rtemp0", pf, pflag2),
                shell=True,
            )  # regression and scatterplot

    # Writing computing time to file

    if not vonly:  # if mode does not only include confirmation and visualization:
        if not mflag:
            comptime_batch = 0
        timefile = open(
            pf + "_results/" + pf + "_files/" + pf + "_time.txt", "w"
        )  # opening time file for writing
        print(
            "%s\t%.5f\t%.5f" % (pf, comptime, comptime_batch), file=timefile
        )  # writing computing time to file
        timefile.close()  # closing time file

    # Cleaning file system, printing message of completion and exiting

    if not zflag:
        os.system("rm -rf " + "TEMP")  # removing temporary directory
    grass.run_command(
        "g.remove", flags="f", type="rast", pattern="r_*", quiet=True
    )  # removing raw result maps
    if not zflag:
        grass.run_command(
            "g.remove", flags="f", type="rast", pattern="x*", quiet=True
        )  # removing temporary raster maps
    grass.run_command(
        "g.remove", flags="f", type="vect", pattern="d_*", quiet=True
    )  # removing temporary vector maps

    grass.run_command("g.region", flags="d")  # resetting default region

    print()
    if not vonly:
        print("Completed in %.1f seconds (net computing time)." % comptime)
    else:
        print("Completed.")

    if (
        cflag or lflag
    ) and not vonly:  # for ellipsoid modes with multi-core processing:

        if lflag:
            lbtmalw = lbtmalw - 1
            # reducing number of allowed layers by one for layer mode
        if lbtmalw >= lbtmmaxmax:
            print("Maximum number of layers for one ellipsoid: %i" % lbtmmaxmax)

        else:  # if too many layers for at least one ellipsoid, print warning message:

            print()
            print(
                "NUMBER OF LAYERS FOR AT LEAST ONE ELLIPSOID EXCEEDS POSSIBLE MAXIMUM:"
            )
            print("MODEL RESULTS ARE INVALID!")
            print()
            print("Maximum number of layers: %i" % lbtmmaxmax)
            print("Maximum number allowed: %i" % lbtmalw)

        print()

    print("Please find the collected results in the directory %s_results." % pf)
    print()

    sys.exit()  # exit


if __name__ == "__main__":
    options, flags = grass.parser()
    main()
