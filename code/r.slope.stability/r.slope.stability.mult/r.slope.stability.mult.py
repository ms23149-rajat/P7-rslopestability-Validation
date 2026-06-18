#!/usr/bin/env python3

##############################################################################
#
# MODULE:       r.slope.stability.mult.py
#
# AUTHOR:       Martin Mergili
# CONTRIBUTORS: Massimiliano Alvioli, Ivan Marchesini
#
# PURPOSE:      The slope stability model - Script for multi-core processing
#
# COPYRIGHT:    (c) 2013 - 2021 by the author
#               (c) 2020 - 2021 by the University of Graz
#               (c) 2013 - 2021 by the BOKU University, Vienna
#               (c) 2015 - 2020 by the University of Vienna
#               (c) 2013 - 2021 by the CNR-IRPI Perugia
#               (c) 2000 - 2021 by the GRASS Development Team
#
# VERSION:      20210625 (25 June 2021)
#
#               This program is free software under the GNU General Public
#               License (>=v2). Read the file COPYING that comes with GRASS
#               for details.
#
##############################################################################

#%module
#% description: the slope stability model: multi-core processing
#% keywords: raster
#% keywords: factor of safety
#% keywords: limit equilibrium model
#% keywords: sliding surface
#% keywords: slope failure probability
#% keywords: slope stability
#%end

import grass.script as grass  # importing libraries
import os


def main():

    n = os.environ["nort"]  # importing N boundary of tile
    s = os.environ["sout"]  # importing S boundary of tile
    e = os.environ["east"]  # importing E boundary of tile
    w = os.environ["west"]  # importing W boundary of tile
    cellsize = os.environ["cellsize"]  # importing raster cell size
    aflag = os.environ["aflag"]  # importing flag for writing additional output
    pflag = os.environ["pflag"]  # importing flag for slope failure probability
    newmark = os.environ["newmark"]  # importing control for Newmark displacement
    jid = os.environ["jid"]  # importing tile number
    os.environ["XINT"] = jid  # exporting tile number
    rtemp = os.environ["rtemp"]  # importing path to temporary directory

    grass.run_command("g.region", n=n, s=s, e=e, w=w, res=cellsize)  # setting region

    print()
    print("Executing slope stability model.")
    print()

    grass.run_command(
        "r.slope.stability.main"
    )  # executing r.slope.stability

    print()
    print("Completed.")
    print()
    print("Exporting result.")
    print()

    grass.run_command(
        "r.out.gdal",
        input="r_fos",
        output=rtemp + "/r_fos" + str(jid),
        format="GTiff",
        flags="c",
        quiet=True,
        overwrite=True,
    )  # exporting factor of safety (fos) to tiff
    grass.run_command(
        "r.out.gdal",
        input="r_ddepth",
        output=rtemp + "/r_ddepth" + str(jid),
        format="GTiff",
        flags="c",
        quiet=True,
        overwrite=True,
    ) # exporting depth of deepest slip surface with fos lower than 1 to tiff
    if pflag == "1":
        grass.run_command(
            "r.out.gdal",
            input="r_pfail",
            output=rtemp + "/r_pfail" + str(jid),
            format="GTiff",
            flags="c",
            quiet=True,
            overwrite=True,
        )  # exporting slope failure probability to tiff
        grass.run_command(
            "r.out.gdal",
            input="r_fos_stdev",
            output=rtemp + "/r_fos_stdev" + str(jid),
            format="GTiff",
            flags="c",
            quiet=True,
            overwrite=True,
        )  # exporting standard deviation of fos to tiff
    if aflag == "1":
        grass.run_command(
            "r.out.gdal",
            input="r_dfos",
            output=rtemp + "/r_dfos" + str(jid),
            format="GTiff",
            flags="c",
            quiet=True,
            overwrite=True,
        )
        # exporting fos for the deepest slip surface with fos lower than 1 to tiff
        grass.run_command(
            "r.out.gdal",
            input="r_depth",
            output=rtemp + "/r_depth" + str(jid),
            format="GTiff",
            flags="c",
            quiet=True,
            overwrite=True,
        )  # exporting depth of slip surface with lowest fos to tiff
        grass.mapcalc(
            '"r_sindex"=float(round(10000*"r_cumfos_crit"/"r_cumsurf"))/10000',
            overwrite=True,
            quiet=True,
        )
        # computing susceptibility index
        grass.run_command(
            "r.out.gdal",
            input="r_sindex",
            output=rtemp + "/r_sindex" + str(jid),
            format="GTiff",
            flags="c",
            quiet=True,
            overwrite=True,
        )  # exporting susceptibility index to tiff
        grass.run_command(
            "r.out.gdal",
            input="r_layer",
            output=rtemp + "/r_layer" + str(jid),
            format="GTiff",
            flags="c",
            quiet=True,
            overwrite=True,
        )  # exporting id of critical layer to tiff
        grass.run_command(
            "r.out.gdal",
            input="r_lcrit",
            output=rtemp + "/r_lcrit" + str(jid),
            format="GTiff",
            flags="c",
            quiet=True,
            overwrite=True,
        )  # exporting length of critical ellipsoid to tiff
        grass.run_command(
            "r.out.gdal",
            input="r_wcrit",
            output=rtemp + "/r_wcrit" + str(jid),
            format="GTiff",
            flags="c",
            quiet=True,
            overwrite=True,
        )  # exporting width of critical ellipsoid to tiff
        grass.run_command(
            "r.out.gdal",
            input="r_zbcrit",
            output=rtemp + "/r_zbcrit" + str(jid),
            format="GTiff",
            flags="c",
            quiet=True,
            overwrite=True,
        )  # exporting zb of critical ellipsoid to tiff

    if newmark == "1":
        grass.run_command(
            "r.out.gdal",
            input="r_dnewmark",
            output=rtemp + "/r_dnewmark" + str(jid),
            format="GTiff",
            flags="c",
            quiet=True,
            overwrite=True,
        )  # exporting Newmark displacement to tiff

    print()
    print("Completed.")
    print()


if __name__ == "__main__":
    options, flags = grass.parser()
    main()
