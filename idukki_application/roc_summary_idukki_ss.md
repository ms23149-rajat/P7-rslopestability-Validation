# P-7 r.slope.stability — Idukki ROC/AUC Results

**Generated:** 2026-06-26  
**Script:** Python (rasterio + sklearn), replacing broken rgdal overlay path

## Setup

| Item | Value |
|---|---|
| Pf raster | `idukki_application/tiffs/idukki_ss_pf_pfail.tif` |
| Raster CRS | EPSG:32643 (UTM 43N) |
| Raster extent | 736×743 pixels, ~22×22 km Idukki tile |
| Pf range | 0.000 – 0.750 |
| Inventory | Kerala Landslide Inventory, SS (shallow slide) events only |
| SS points in tile | 530 of 1,760 total SS records |
| Negative samples | 2,650 random stable pixels (5× positives, excluding SS locations) |

## Results

| Metric | Value |
|---|---|
| **AUC-ROC** | **0.685** |
| Positive samples | 530 SS crown points |
| Negative samples | 2,650 random stable pixels |

## Performance at key Pf thresholds

| Threshold | TPR (rTP) | FPR | TNR (rTN) |
|---|---|---|---|
| 0.1 | 0.613 | 0.341 | 0.659 |
| 0.2 | 0.275 | 0.158 | 0.842 |
| 0.3 | 0.185 | 0.113 | 0.887 |
| 0.4 | 0.036 | 0.031 | 0.969 |
| 0.5 | 0.036 | 0.031 | 0.969 |
| 0.6 | 0.036 | 0.031 | 0.969 |
| 0.7 | 0.000 | 0.006 | 0.994 |

## Interpretation

AUC = 0.685 indicates weak-to-moderate discriminatory skill, better than random
(0.5) but below the 0.7 threshold commonly cited as acceptable in landslide
susceptibility literature.

**Why the AUC is modest — honest assessment:**

1. **Literature parameters, not field-measured** — cohesion, friction angle, and
   unit weight were taken from published ranges for Kerala laterite/regolith, not
   site-measured. Incorrect soil parameters shift Pf uniformly across the tile,
   which reduces discrimination without necessarily reducing spatial pattern skill.

2. **Single rainfall scenario** — Pf was computed at one fixed rainfall input
   (Kerala 2018-style). ROC skill depends on the modelled rainfall matching the
   actual triggering intensity for each event; a single scenario cannot match all
   530 events simultaneously.

3. **Point inventory vs pixel raster** — SS crown points mark initiation locations;
   the Pf raster represents per-pixel stability. A crown point falling on a
   Pf=0.10 pixel does not mean the model is wrong — the failure may have initiated
   at a pixel boundary or just upslope.

4. **Random negatives may include unrecorded failures** — "stable" pixels drawn
   randomly from the tile include terrain that may have failed before the inventory
   period or failed in unreported events. This inflates apparent false-negative rate
   and deflates AUC.

5. **Conservative parameter set** — Pf_max = 0.75 (model never reaches certainty).
   Calibrating C or φ to site-measured values would likely shift the Pf
   distribution and improve AUC.

**Context from threshold table:** The previously reported rTP = 4.4% at Pf ≥ 0.5
was an artefact of using too conservative a threshold. At Pf ≥ 0.1 the model
captures 61.3% of SS events (TPR = 0.613), at the cost of a 34.1% false-positive
rate. The optimal operating threshold (balancing TPR and TNR) lies near Pf = 0.2–0.3,
where TPR = 27.5% / FPR = 15.8%.

**Comparison with literature:** Physically-based slope stability models applied
to regional-scale inventories typically report AUC 0.60–0.75 without site
calibration (e.g., Mergili et al. 2014, Tran et al. 2018). AUC = 0.685 is
within this expected range for an uncalibrated first-pass application.
