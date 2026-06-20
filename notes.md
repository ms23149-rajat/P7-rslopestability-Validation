# Notes — P-7 (r.slope.stability)

## Environment

- Workstation: labws2, IP 172.26.248.92 (hostname `labws2.iitk.ac.in` doesn't resolve reliably — connect via IP).
- GRASS GIS 8.3.2, R 4.3.3, Ubuntu 24.04.4 LTS.
- Conda env `p7_rslopestability` (python=3.11) created for general housekeeping. The model itself is C + the GRASS Python API + R, not a Python package, so this env mostly just keeps a clean shell — it isn't actually a dependency of the model.
- R packages installed into a personal library (`~/R/library`, set via `R_LIBS_USER` in `~/.Renviron`), since the lab account has no write access to the system R library and no sudo.

## Methodology followed

Per the project's required loop: read the source papers first (GMD 2014a and Geomorphology 2014b), then locate and download the authors' own code and training data rather than reimplementing anything from the papers' equations, get it running on the authors' own test case, and document failures honestly rather than smoothing them over.

## The five GRASS 7 → GRASS 8 compatibility issues

The addon (`r.slope.stability` v2.0 Pre-release, 20210715) was written and tested against GRASS 7.8. The lab machine runs GRASS 8.3.2 — a two-major-version jump — and every issue below traces back to that gap, not to anything wrong with the model's science.

### 1. `bingrass` NameError

`r.slope.stability` looks for the GRASS binary by scanning for version-suffixed names (`grass70` … `grass79`), a GRASS-7-era naming convention. On this machine the binary is just `/usr/bin/grass`, no suffix, so the scan never matches and `bingrass` is never assigned — but the multi-core tiling code uses it anyway, raising a `NameError` the first time `-m` is used.

Fix, in `~/.grass8/addons/scripts/r.slope.stability`, before the version-suffix loop:

```python
bingrass = which("grass")

if bingrass is None:
    for i in range(70, 80):
        if not which("grass" + str(i)) is None:
            bingrass = which("grass" + str(i))
```

### 2. `-text` vs `--text`

GRASS 8's launcher is argparse-based and only accepts double-dash long flags; the addon's tile-launch template hardcodes the old single-dash `-text`. Every tile failed instantly with `grass: error: unrecognized arguments: -text`, though the wrapper's threading code doesn't check subprocess exit codes, so it printed `Tile no. N completed` regardless of whether the tile actually did anything.

Fix:
```bash
sed -i 's/-text/--text/g' ~/.grass8/addons/scripts/r.slope.stability
```

### 3. `library(rgdal)`

`rgdal` was retired from CRAN in 2023 and can no longer be installed at all. Two of the addon's R helper scripts (`r.slope.stability.map.R`, `r.slope.stability.roc.R`) load it. These actually live at `~/.grass7/addons/etc/r.slope.stability.rcode/` — note `.grass7`, even though the rest of the addon installed under `.grass8` — because that path is hardcoded in the addon's own Makefile, independent of which GRASS version's `g.extension` performs the install.

Fix (stopgap, not complete — see caveat):
```bash
sed -i 's/^library(rgdal)/# library(rgdal)/' ~/.grass7/addons/etc/r.slope.stability.rcode/r.slope.stability.map.R
sed -i 's/^library(rgdal)/# library(rgdal)/' ~/.grass7/addons/etc/r.slope.stability.rcode/r.slope.stability.roc.R
```

Caveat: this silences the `library()` import error but not later calls to `readOGR()` / `readGDAL()` (functions `rgdal` provided), which still fail with `could not find function`. In practice this only affects the landslide-boundary overlay on the PNG maps and the ROC-curve plot — the actual FoS/Pf computation, the exported TIFFs, and the printed confirmation statistics (`rTP`/`rTN`/`rFP`/`rFN`) are all written out successfully before that point.

### 4. `scriptpath2` hardcoded to `.grass7`

Separately from issue 3, the path used to find `r.slope.stability.mult` (the script multi-core processing runs per tile) is hardcoded as `$HOME/.grass7/addons/scripts`, but `g.extension` under GRASS 8.3.2 actually installs it at `~/.grass8/addons/scripts/`. Unlike issue 3, no file exists at the old path here, so this one needed a real fix rather than a workaround.

Fix:
```bash
sed -i 's#\.grass7/addons/scripts#\.grass8/addons/scripts#g' ~/.grass8/addons/scripts/r.slope.stability
```

### 5. `GRASS_BATCH_JOB` not auto-executing

With issues 1–4 fixed, each tile's `grass --text <location>/mapN` call succeeded and even auto-created the `mapN` mapset on first use (confirmed: `map1` through `map25` all exist on disk with proper GRASS mapset structure) — but the batch job script never actually ran inside them. Manually entering one tile's session and checking confirmed `GRASS_BATCH_JOB` was correctly set and the tile's environment variables (`jid`, `nort`, `west`, etc.) were correctly exported, yet `g.region -p` still showed the full untiled extent rather than the tile's sub-extent that the batch job is supposed to set via `g.region`.

This points to GRASS 8.3.2's `--text` launcher no longer auto-executing `GRASS_BATCH_JOB` on startup the way GRASS 7 did. Its replacement appears to be the explicit `--exec` flag visible in the launcher's own usage banner. This wasn't patched in the addon; instead, the underlying engine was confirmed directly by invoking `r.slope.stability.main` by hand inside a manually-launched session (see README, "What was confirmed to run"). A clean fix would likely involve rewriting the wrapper to call `grass --exec` instead of relying on the environment variable — out of scope for this reproduction, but the natural next step if multi-core processing turns out to be genuinely necessary for a larger real-terrain study area later.

### 6. Redundant `gwdepth` parameter conflicting with `seepage=1`

The authors' own working examples (e.g. `sl_ell_pf`) never pass a `gwdepth=` raster when using `seepage=1` — that flag already tells the model to treat the water table as following the ground surface (fully saturated), using the per-layer saturated-water-content values supplied in `geotech=`. Supplying an additional, redundant `gwdepth=0` raster on top of this caused an immediate crash in pre-loop setup (signal SIGFPE), almost certainly because the engine fed a literal zero into a formula expecting a water-table-depth ratio rather than an absolute value.

Fix: when `seepage=1` is used, omit `gwdepth=` entirely.

### 7. NULL soil-class cells are not handled gracefully

When applying the model to real terrain, areas excluded from analysis (e.g. open water, masked out of `soilclass` as NULL) caused unpredictable crashes once randomly-generated ellipsoids happened to land on or near them. The engine appears to assume every cell in the computational region has a valid, defined soil class.

Fix: rather than leaving excluded areas as NULL, assign them their own dummy soil class with deliberately extreme "never fails" geotechnical parameters (very high cohesion and friction angle). Their FoS/Pf results are physically meaningless and can be masked out during interpretation, but the engine itself needs every cell defined.

### 8. `rand() % 0` in ellipsoid geometry sampling (genuine bug, patched)

Confirmed via direct GDB debugging (see below) rather than source-reading alone. In the ellipsoid randomization code (`main.c`, ellipsoid length/width/depth/zb sampling, originally around lines 1644–1793), each dimension is sampled via:

\`\`\`c
randmax = ( int ) ( 1000 * ( XMAX - XMIN ) + 0.5 );
val = ( ( float ) ( rand()%randmax ) ) / 1000 + XMIN;
\`\`\`

When the `RAT_GEOM_MIN`/`RAT_GEOM_MAX` (length:width ratio) constraints clamp `ELLBMAX` and `ELLBMIN` to nearly identical values (which happens legitimately for some random ellipsoid orientations, confirmed via GDB to differ by as little as 0.0003), `1000 * (ELLBMAX - ELLBMIN)` rounds down to 0, making `randmax = 0`. `rand() % 0` is undefined behavior in C and reliably raised `SIGFPE`. This crash was data-dependent and occurred at inconsistent points across runs (immediately, or after testing hundreds of thousands of ellipsoids), since it depends on the specific random values drawn.

Interestingly, the authors had already patched this exact issue for the depth dimension only (`if (ELLCMAX > ELLCMIN + 0.001)`, a small safety margin), but never extended the same protection to the zb, length, and width dimensions.

Fix: added a `randmax > 0` guard to all three remaining vulnerable spots (zb, length, width), falling back to the boundary value when the valid sampling range rounds to zero width:

\`\`\`c
if ( randmax > 0 ) ell_b = ( ( float ) ( ( rand()%randmax ) ) ) / 1000 + ELLBMIN;
else ell_b = ELLBMIN;
\`\`\`

Diagnosed using GDB directly against the compiled binary (`gdb --args ~/.grass8/addons/bin/r.slope.stability.main`, reproducing the crash live and inspecting variable values with `print`), after parameter-file capture from a crashed run (`~/TEMP/rtemp0/param*.txt`) confirmed real inputs could be reproduced outside the Python wrapper.

## What this means in practice

The scientific engine (`r.slope.stability.main`) and every single-process invocation mode are fully functional and have been exercised across all of the authors' original example commands except the tiled one, including full empirical confirmation against the bundled landslide inventory. Only the automatic multi-core orchestration layer is broken, and its root cause is understood precisely rather than just observed. For study areas in the range tested here (under ~1 million cells), single-process runtime (30 s–130 s) is fast enough that tiling isn't actually needed in practice — worth keeping in mind before spending more time patching issue 5 for the real Kerala-scale work.

Following this debugging, the model was successfully applied to a real ~22km × 22km tile of Idukki district (30m DEM, real Kerala 2018 landslide inventory, land-cover-differentiated soil parameters) — see README for results.

## AI assistance

Claude was used throughout as an interactive, step-by-step mentor: explaining both source papers page by page before any coding began, walking through GRASS concepts and commands one at a time, and helping diagnose each of the five compatibility issues above by reasoning through the GRASS 7→8 changes together with the actual error messages and source lines. All commands were run by hand on the lab workstation; Claude has no direct access to the lab machine and never executed anything on it.
