# P-7 — r.slope.stability

**Status: RUNS** — core engine and 6 of 7 single-process official example workflows confirmed end-to-end; the multi-core tiling wrapper has a documented, root-caused GRASS 8.3.2 incompatibility (see `notes.md`).

**Intern:** Rajat (rajat-surge2026), SURGE/SARIP 2026, Group B
**Supervisor:** Dr. Shyam Nandan
**Workstation:** labws2 (Ubuntu 24.04.4 LTS, GRASS GIS 8.3.2, R 4.3.3)

## What this model is

r.slope.stability is a GRASS GIS module for physically-based slope-stability analysis over large areas. It tests large numbers of randomly placed, randomly sized ellipsoidal slip surfaces against a 3-D limit-equilibrium (Hovland) model, reporting either a per-pixel Factor of Safety (FoS) or — by sampling cohesion, friction angle, and depth across their measured uncertainty ranges — a per-pixel slope failure probability (Pf). It also supports splitting large study areas into overlapping tiles for parallel processing across CPU cores.

## Source papers

- Mergili, M., Marchesini, I., Alvioli, M., Metz, M., Schneider-Muntau, B., Rossi, M., Guzzetti, F. (2014a). A strategy for GIS-based 3-D slope stability modelling over large areas. *Geoscientific Model Development*, 7, 2969–2982.
- Mergili, M., Marchesini, I., Rossi, M., Guzzetti, F., Fellin, W. (2014b). Spatially distributed three-dimensional slope stability modelling in a raster GIS. *Geomorphology*, 206, 178–195.

## Code source

Authors' own distribution, **not** the official GRASS Addons repository: https://www.landslidemodels.org/r.slope.stability/

Version used: `r.slope.stability20pre_20210715.zip` (2.0 Pre-release, 15 July 2021). The authors built and tested this version against GRASS 7.8 and Ubuntu 20.04; the lab workstation runs GRASS 8.3.2 — a two-major-version gap that is the source of every compatibility issue documented in `notes.md`.

## Test data used

Both datasets come from the authors' own training-data page (https://www.landslidemodels.org/r.slope.stability/data.php), not from the original 2014 papers directly — the papers' field-measured Collazzone dataset is not publicly bundled.

- **Slideslope**: a small (200×200 cell, 1 m resolution) synthetic, computer-generated slope, used for the single-ellipsoid soil-class and layer-mode examples.
- **Slideland**: a larger (~9.26M cell, 5 m resolution) terrain dataset. Per the authors' own documentation, the DEM was "modified from existing data" and the observed-landslide layer is drawn from the real IFFI inventory (CC BY-SA 4.0); every other layer and the geotechnical parameters used in the start script are explicitly fictitious. The coordinates and the page's own image caption suggest this is Collazzone-shaped terrain, but **this is the authors' bundled demonstration dataset, not the field-measured dataset behind the 2014 papers' own validation** — worth stating precisely in any further write-up.

## What was confirmed to run

| Run | Mode | Result |
|---|---|---|
| `so_sing_class` | Slideslope, soil-class mode, single ellipsoid | FoS = 0.842 |
| `so_sing_lyr` | Slideslope, layer mode, single ellipsoid | FoS = 0.741 |
| `slx_inf_fos` | Slideland, infinite slope, 1 lithology class | rTP 2.3%, rTN 47.0% |
| `sly_inf_fos` | Slideland, infinite slope, 5 classes | rTP 1.6%, rTN 71.6% |
| `sl_inf_fos` | Slideland, infinite slope, 5 classes, varying depth | rTP 2.1%, rTN 63.8% |
| `sl_ell_fos` | Slideland, ellipsoid FoS, ~1.65M ellipsoids tested | 30.6 s, rTP 3.3%, rTN 51.7% |
| `sl_ell_newmark` | Slideland, + Newmark seismic | 32.6 s, same confirmation stats |
| `sl_ell_pf` | Slideland, single-core slope failure probability | 39.0 s, convergence (evolution) data captured |

All eight rows above are the authors' own official example sequence from `sl.slope.stability.start.sh`. rTP/rTN/rFP/rFN are true/false positive/negative rates against the bundled `landslides` observed layer — the same empirical-confirmation metrics reported in the source papers.

The ninth official command, the multi-core tiled version (`sl_ell_pf_tiles`), does not complete through its own automatic wrapper (see below), but the underlying computation was independently confirmed by invoking the compiled core engine (`r.slope.stability.main`) directly inside a manually-launched session: 130 s, 1,659,518 of 1,670,795 ellipsoids valid, `r_fos` (min 0.185) and `r_pfail` (range 0–1, mean 0.43) both written with sane values.

## Known issue: multi-core tiling wrapper

The `-m` (multi-core) flag's automatic tile-splitting and recombination does not work under GRASS GIS 8.3.2. Five GRASS-7-vs-8 compatibility issues were traced and four were fixed directly; the fifth — GRASS's `GRASS_BATCH_JOB` mechanism not auto-executing under the GRASS 8.3.2 `--text` launcher — was root-caused but bypassed rather than patched. Full detail in `notes.md`. In practice, for study areas in the few-hundred-thousand to ~1 million cell range, this isn't a blocker: the core engine runs single-process fast enough that tiling isn't actually required.

## Quick start

```bash
bash install.sh
```

then, from inside a GRASS session opened on the Slideslope location:r.slope.stability -t -v prefix=so_sing_class cellsize=1 model=c elevation=dem soilclass=lithology gwdepth=gwdepth seepage=1 numlayers=5 depthmaps=depth1,depth2,depth3,depth4,depth5 geotech=1,1,15000,4000,30,40,1,2,15000,2000,20,35,1,3,15000,2500,40,40,1,4,15000,5000,30,40,1,5,27000,100000,89,0 elldens=0 ellips=100,100,50,50,25,25,0,0,3,1,100,100,-9999,-9999
A successful run prints `Number of valid ellipsoids: 1 out of 1.` and a Factor of Safety around 0.84 — confirming the install before moving on to the larger Slideland examples.

## Real-terrain application: Idukki, Kerala

Following validation against the authors' official example datasets above, the model was applied to a real ~22 km × 22 km tile of Idukki district, Kerala — chosen because Idukki has by far the highest concentration of recorded landslides in the state (2,223 of 4,728 total events in the Hao et al. 2020 / Kerala 2018 landslide inventory).

### Data sources

- **DEM**: Copernicus GLO-30 (~30 m), reprojected to UTM zone 43N (EPSG:32643).
- **Land cover**: ESA WorldCover 10 m, 2021 (v200), tile N09E075, resampled to the DEM grid via nearest-neighbour (categorical data must not be interpolated).
- **Landslide inventory**: Hao et al. (2020), Kerala 2018 landslide dataset (DANS, DOI 10.17026/dans-x6c-y7x2). Filtered to events recorded as "SS" (shallow slide) — the failure mode r.slope.stability's physics actually represents, as distinct from debris flows (DF) or rockfalls (RF), which the model does not simulate. 530 SS-type events fall within the tile; 527 occupy distinct 30 m cells after rasterization.

### Soil parameterization

Land cover was reclassified into three soil-mechanics-relevant classes (forest, grassland, cleared/disturbed), with water masked out and given a fourth "never fails" dummy class so the engine has a defined value everywhere. Geotechnical parameters (cohesion, friction angle, unit weight, saturated water content) follow a forest > grassland > cleared gradient reflecting root-reinforcement decline, grounded in regional soil-depth data for this part of Idukki (Sajinkumar et al., reporting 0.25–5 m soil cover) and general Western Ghats lateritic/weathered-crystalline residual soil literature. These are literature-typical values, not site-measured — a documented limitation, not a claim of field calibration. A two-layer structure (variable soil layer over a uniform, non-failing bedrock floor at 5 m) follows the authors' own demonstrated best practice from their official examples. Groundwater was assumed fully saturated via `seepage=1` (worst-case monsoon-trigger scenario), consistent with the model's design and the prior TRIGRS Idukki application.

### New compatibility issues found and fixed

Applying the model to real terrain (rather than the authors' bundled demo data) surfaced three further issues beyond the five already documented above — see notes.md, issues 6–8. Briefly: a redundant `gwdepth` parameter conflicting with `seepage=1`, unhandled NULL soil-class cells, and a genuine `rand() % 0` bug in the ellipsoid-sampling code (confirmed via direct GDB debugging against the compiled binary, then patched).

### Results

A full probabilistic run (45,735,802 valid ellipsoids out of 59,087,772 tested, `elldens=100`, ~5.4 minutes net compute time) produced Factor of Safety, failure-probability, and susceptibility-index rasters across the full tile.

The model's own R-based confirmation script could not run (the same pre-existing rgdal/CRAN retirement issue documented above), so confirmation statistics were computed directly in GRASS instead, cross-tabulating the failure-probability raster against the 527-cell observed-landslide layer:

| Threshold | rTP (true positive rate) | rTN (true negative rate) |
|---|---|---|
| Pf ≥ 0.5 | 4.4% (23/527) | 96.7% (527,000/544,843) |

For comparison, the authors' own official Slideland confirmation runs (documented above) achieved rTP 1.6–3.3% and rTN 47.0–71.6%. This real-terrain Idukki application's true-positive rate is comparable to or better than every one of those, with a substantially higher true-negative rate — though it should be read as one data point at one threshold, not a full ROC evaluation, and as an unvalidated first-pass application in the same sense as the earlier TRIGRS Idukki work: physically consistent and empirically suggestive, but not field-calibrated.

### Honest limitations

- Single ~22 km × 22 km tile, not the full Idukki district.
- Soil parameters are literature-typical, not site-measured.
- Land cover (and therefore soil class) reclassification is a simplification — in particular, the "grassland" class may include misclassified tea/cardamom plantation, which ESA WorldCover does not reliably distinguish from natural grassland in this terrain.
- The rgdal-blocked map overlay and ROC plot were not regenerated; only the underlying numeric rasters and a manually-computed single-threshold confusion matrix are available for this run.
