#!/bin/bash
#
# reproduce_idukki.sh -- Reproduce the Idukki, Kerala real-terrain
# application of r.slope.stability documented in ../README.md
# ("Real-terrain application: Idukki, Kerala") and ../notes.md
# (issues 6-8).
#
# PREREQUISITES (not automated by this script):
#
#   1. Run ../install.sh first. This script assumes r.slope.stability is
#      already installed AND patched, including the issue-8 ellipsoid-
#      sampling fix -- without it, step 8 below will crash intermittently.
#
#   2. A DEM covering your area of interest, already reprojected to a
#      real metric CRS (this project used UTM zone 43N / EPSG:32643 for
#      Kerala). This project's DEM was a ~30m Copernicus GLO-30 tile
#      clipped to a ~22km x 22km area of Idukki district; reproducing
#      that exact clip is outside this script's scope. Set DEM_PATH and
#      DEM_EPSG below to match your own DEM, or substitute your own AOI
#      entirely -- nothing past step 1 is Idukki-specific by necessity,
#      only by the choices made for this project.
#
#   3. A landslide inventory (point or polygon) for empirical confirmation,
#      with a usable type/class attribute field if you want to filter by
#      failure type the way this project did (SS = shallow slide). This
#      project used Hao et al. (2020), the Kerala 2018 inventory, manually
#      downloaded from https://doi.org/10.17026/dans-x6c-y7x2 (DANS
#      repository -- not freely scriptable, fetch and unzip it yourself).
#
# Adjust the variables below to your own paths/CRS/filters before running.
#
set -e

WORKDIR="$HOME/landslide-toolkit/P7"
DEM_PATH="${DEM_PATH:-$WORKDIR/idukki_dem.tif}"
DEM_EPSG="${DEM_EPSG:-EPSG:32643}"
INVENTORY_SHP="${INVENTORY_SHP:-$WORKDIR/kerala_inventory/Kerela landslide.shp}"
TYPE_FIELD="${TYPE_FIELD:-Type_of_sl}"
TYPE_VALUE="${TYPE_VALUE:-SS}"
LOCATION="$WORKDIR/idukki"
WGS84_LOC="$WORKDIR/kerala_wgs84"

if [ ! -f "$DEM_PATH" ]; then
    echo "ERROR: DEM not found at $DEM_PATH"
    echo "See prerequisites above -- set DEM_PATH, or place your DEM there."
    exit 1
fi
if [ ! -f "$INVENTORY_SHP" ]; then
    echo "ERROR: landslide inventory not found at $INVENTORY_SHP"
    echo "See prerequisites above -- download and unzip it, then set INVENTORY_SHP."
    exit 1
fi

echo "== 1. Creating GRASS location ($DEM_EPSG) and importing DEM =="
rm -rf "$LOCATION"
grass -c "$DEM_EPSG" "$LOCATION" --exec g.region -p
grass "$LOCATION/PERMANENT" --exec r.in.gdal input="$DEM_PATH" output=dem --overwrite
grass "$LOCATION/PERMANENT" --exec g.region raster=dem -p

echo "== 2. Importing landslide inventory (native CRS) and reprojecting =="
# v.in.ogr refuses cross-CRS import by design (it would silently mislabel
# coordinates, not transform them) -- import into a location matching the
# inventory's own CRS first, then use v.proj to do the real reprojection.
rm -rf "$WGS84_LOC"
grass -c "$INVENTORY_SHP" "$WGS84_LOC" --exec v.in.ogr \
    input="$INVENTORY_SHP" output=landslides_all --overwrite
grass "$LOCATION/PERMANENT" --exec v.proj location="$(basename "$WGS84_LOC")" \
    mapset=PERMANENT input=landslides_all output=landslides_all --overwrite

echo "== 3. Spatially filtering to events within the DEM extent =="
grass "$LOCATION/PERMANENT" --exec v.in.region output=dem_extent --overwrite
grass "$LOCATION/PERMANENT" --exec v.select ainput=landslides_all atype=point \
    binput=dem_extent operator=overlap output=landslides_in_tile --overwrite

echo "== 4. Filtering by failure type ($TYPE_FIELD = $TYPE_VALUE) =="
grass "$LOCATION/PERMANENT" --exec v.extract input=landslides_in_tile \
    where="${TYPE_FIELD}='${TYPE_VALUE}'" output=landslides_filtered --overwrite

echo "== 5. Rasterizing observed landslides (explicit 0/1, no NULLs) =="
# r.slope.stability's confirmation script reads >0 as landslide, <0 as
# excluded, and (critically) leaves 0 as "confirmed stable" -- the whole
# region must be 0 or 1, never NULL, or confirmation stats silently break.
grass "$LOCATION/PERMANENT" --exec v.to.rast input=landslides_filtered \
    output=landslides_pts use=val value=1 type=point --overwrite
grass "$LOCATION/PERMANENT" --exec r.mapcalc \
    "observed = if(isnull(dem), null(), if(isnull(landslides_pts), 0, 1))" --overwrite

echo "== 6. Downloading and reclassifying ESA WorldCover land cover =="
# Adjust the tile name below if your AOI falls in a different 3x3-degree
# WorldCover tile -- see https://esa-worldcover.s3.eu-central-1.amazonaws.com/readme.html
LANDCOVER_TILE="${LANDCOVER_TILE:-N09E075}"
LANDCOVER_TIF="$WORKDIR/landcover_${LANDCOVER_TILE}.tif"
if [ ! -f "$LANDCOVER_TIF" ]; then
    wget "https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/ESA_WorldCover_10m_2021_v200_${LANDCOVER_TILE}_Map.tif" \
        -O "$LANDCOVER_TIF"
fi
grass "$WGS84_LOC/PERMANENT" --exec r.in.gdal input="$LANDCOVER_TIF" output=landcover --overwrite
grass "$LOCATION/PERMANENT" --exec g.region raster=dem
# method=nearest is mandatory for categorical data -- bilinear/cubic would
# blend adjacent class codes into meaningless intermediate values.
grass "$LOCATION/PERMANENT" --exec r.proj location="$(basename "$WGS84_LOC")" \
    mapset=PERMANENT input=landcover output=landcover method=nearest --overwrite

echo "== 7. Building 4-class soilclass (forest / grassland / cleared / water) =="
# WorldCover codes: 10=tree cover, 30=grassland, 80=water (-> dummy class 4,
# "never fails", so the engine always has a defined value); everything
# else (cropland/built-up/bare) -> class 3, cleared/disturbed.
grass "$LOCATION/PERMANENT" --exec r.mapcalc \
    "soilclass = if(landcover == 80, 4, if(landcover == 10, 1, if(landcover == 30, 2, 3)))" --overwrite

echo "== 8. Building per-class soil-depth and bedrock-floor rasters =="
# depth1 = depth to bottom of the real (variable) soil layer, by class.
# depth2 = depth to a uniform, non-failing bedrock floor (bounds how deep
# ellipsoids can search) -- shared across all classes, see geotech= below.
# Adjust these to match the geotech= parameters for your own classes.
grass "$LOCATION/PERMANENT" --exec r.mapcalc \
    "depth1 = if(soilclass==1,2.5,if(soilclass==2,1.5,1.0))" --overwrite
grass "$LOCATION/PERMANENT" --exec r.mapcalc "depth2 = 5.0" --overwrite

echo "== 9. Running r.slope.stability (probabilistic, full confirmation) =="
# NOTE: do not pass gwdepth= alongside seepage=1 -- the authors' own
# working examples never combine these, and doing so crashes the engine
# (notes.md, issue 6). cellsize must match your DEM's actual resolution.
grass "$LOCATION/PERMANENT" --exec g.region raster=dem
grass "$LOCATION/PERMANENT" --exec r.slope.stability -a -p -t -v \
    prefix=idukki_ss_pf cellsize=30 model=c elevation=dem soilclass=soilclass \
    seepage=1 observed=observed numlayers=2,2,2,2 depthmaps=depth1,depth2 \
    geotech=1,1,18000,15000,32,42,4000,5000,25000,4,24,40,1,2,27000,97500,88,0,1,95000,100000,1,87,89,2,1,18000,10000,30,38,3000,3000,17000,4,22,38,2,2,27000,97500,88,0,1,95000,100000,1,87,89,3,1,18000,8000,28,35,2500,2000,14000,4,20,36,3,2,27000,97500,88,0,1,95000,100000,1,87,89,4,1,18000,99000,89,50,1,98000,100000,1,88,89,4,2,27000,97500,88,0,1,95000,100000,1,87,89 \
    sampling=2,3,3,0,1,1,1 elldens=100 \
    ellips=100,20,30,8,3,1,0.5,0.0,3,1,-9999,-9999,-9999,-9999 overlap=300

echo "== 10. Computing confirmation statistics (bypasses rgdal-blocked R script) =="
# r.slope.stability's own R confirmation step needs readOGR/readGDAL
# (rgdal), retired from CRAN in 2023 -- see notes.md issue 3. Computing
# the same rTP/rTN statistics directly in GRASS sidesteps it entirely.
grass "$LOCATION/PERMANENT" --exec g.region raster=dem
grass "$LOCATION/PERMANENT" --exec r.mapcalc \
    "predicted_unstable = if(idukki_ss_pf_pfail >= 0.5, 1, 0)" --overwrite
grass "$LOCATION/PERMANENT" --exec r.stats -c input=predicted_unstable,observed

echo "== Done. =="
echo "Result rasters (FoS, Pf, susceptibility index, etc.) are at:"
echo "  $(find "$HOME" -maxdepth 1 -iname 'idukki_ss_pf_results' 2>/dev/null)"
