#!/bin/bash
#
# install.sh -- P-7 (r.slope.stability) setup on the SURGE 2026 lab workstation
#
# Assumes: GRASS GIS 7.8+ already installed system-wide (lab has 8.3.2),
# R 4.x already installed system-wide, no sudo / apt-get access.
# Run from inside ~/landslide-toolkit/P7 (or adjust WORKDIR below).
#
set -e

WORKDIR="$HOME/landslide-toolkit/P7"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

echo "== 1. Downloading authors' source code =="
SRC_ZIP="r.slope.stability20pre_20210715.zip"
if [ ! -f "$SRC_ZIP" ]; then
    wget "https://www.landslidemodels.org/r.slope.stability/software/$SRC_ZIP"
fi
rm -rf r_slope_stability
unzip -q "$SRC_ZIP" -d r_slope_stability

echo "== 2. Compiling and installing via g.extension =="
# --exec runs non-interactively and sidesteps the GRASS_BATCH_JOB issue (see notes.md, issue 5)
grass --tmp-location EPSG:4326 --exec g.extension \
    extension=r.slope.stability \
    url="$WORKDIR/r_slope_stability/r.slope.stability"

ADDON_SCRIPT="$HOME/.grass8/addons/scripts/r.slope.stability"
if [ ! -f "$ADDON_SCRIPT" ]; then
    echo "ERROR: expected addon script not found at $ADDON_SCRIPT"
    echo "If your GRASS major version isn't 8, addons install under ~/.grassN/ instead -- adjust ADDON_SCRIPT above."
    exit 1
fi

echo "== 3. Patching known GRASS 7 -> GRASS 8 incompatibilities (see notes.md) =="

# Issue 2: -text -> --text
sed -i 's/-text/--text/g' "$ADDON_SCRIPT"

# Issue 4: scriptpath2 hardcoded to .grass7/addons/scripts
sed -i 's#\.grass7/addons/scripts#\.grass8/addons/scripts#g' "$ADDON_SCRIPT"

# Issue 1: bingrass fallback
python3 - "$ADDON_SCRIPT" <<'PYEOF'
import sys, pathlib
p = pathlib.Path(sys.argv[1])
text = p.read_text()
old = 'for i in range(70, 80):\n    if not which("grass" + str(i)) is None:\n        bingrass = which("grass" + str(i))'
new = ('bingrass = which("grass")\n\n'
       'if bingrass is None:\n'
       '    for i in range(70, 80):\n'
       '        if not which("grass" + str(i)) is None:\n'
       '            bingrass = which("grass" + str(i))')
if old in text:
    p.write_text(text.replace(old, new))
    print("Patched bingrass fallback.")
elif 'bingrass = which("grass")' in text:
    print("bingrass fallback already patched, skipping.")
else:
    print("WARNING: bingrass pattern not found -- addon version may differ; check manually.")
PYEOF

# Issue 3: library(rgdal) -- retired from CRAN, comment out (stopgap; see notes.md caveat)
RCODE_DIR="$HOME/.grass7/addons/etc/r.slope.stability.rcode"
if [ -d "$RCODE_DIR" ]; then
    sed -i 's/^library(rgdal)/# library(rgdal)/' "$RCODE_DIR/r.slope.stability.map.R"
    sed -i 's/^library(rgdal)/# library(rgdal)/' "$RCODE_DIR/r.slope.stability.roc.R"
else
    echo "WARNING: $RCODE_DIR not found. Find the real R-script directory with:"
    echo "  grep -n 'scriptpath = ' $ADDON_SCRIPT"
    echo "and comment out 'library(rgdal)' in r.slope.stability.map.R and r.slope.stability.roc.R there."
fi

echo "== 4. Personal R library + packages =="
mkdir -p "$HOME/R/library"
grep -qxF 'R_LIBS_USER=~/R/library' "$HOME/.Renviron" 2>/dev/null || \
    echo 'R_LIBS_USER=~/R/library' >> "$HOME/.Renviron"

Rscript -e 'install.packages(c("sp","Rcpp","fmsb","ROCR","raster"), repos="https://cloud.r-project.org")'
# NOTE: rgdal, rgeos, maptools (listed in the official manual) are retired from CRAN
# and cannot be installed any more -- step 3 above patches around them instead.

echo "== 5. Downloading official test datasets =="
for f in slideslope.zip slideland.zip; do
    if [ ! -f "$f" ]; then
        echo "Downloading $f ..."
        wget "https://www.landslidemodels.org/r.slope.stability/data/$f"
    fi
    name="${f%.zip}"
    [ -d "$name" ] || unzip -q "$f"
done
wget -nc https://www.landslidemodels.org/r.slope.stability/data/so.slope.stability.start.sh
wget -nc https://www.landslidemodels.org/r.slope.stability/data/sl.slope.stability.start.sh

echo "== Done. =="
echo "Smoke test:"
echo "  grass $WORKDIR/slideslope/PERMANENT"
echo "  r.slope.stability --help"
echo "See README.md for the validated example commands."
