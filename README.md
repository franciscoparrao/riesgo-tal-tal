# Reach, not initiation — code

Code for the manuscript *"Reach, not initiation: an event-based test of the hazard term in composite
risk indices for a fan-built city in the Atacama Desert"* (Parra, Gil-Costa, Bonacic and Marín;
submitted to *Natural Hazards*).

The pipeline builds a conventional composite risk index (initiation susceptibility, vulnerability,
exposure and lack of coping capacity) for the commune of Taltal (Chile), replaces the hazard term
by a terrain-based reach index (Height Above Nearest Drainage, HAND) selected on five other communes
affected by the same event, and validates both against the damage of the aluvión (sediment-laden
flow) of 25–26 March 2015.

## Contents

| Path | Content |
|---|---|
| `_targets.R`, `run.R`, `config.yml` | [`targets`](https://docs.ropensci.org/targets/) pipeline and its configuration |
| `R/` | Pipeline functions (one file per component: exposure, vulnerability, coping capacity, reach hazard, risk, validation, robustness, selection, official zoning benchmark, data deposit) |
| `bin/` | Helper scripts: r5r isochrones (`isocrona.R`), receiving-channel area for the GFI (`area_cauce_destino.py`), Flow-Py runs (`flowpy.py`), exploratory simplified flow-path model (`sflm.py`, `run_sflm.sh`) |
| `paper/R/` | Figures of the manuscript (`todas_las_figuras.R` runs all of them) |

Comments and identifiers are in Spanish.

## Requirements

- R ≥ 4.5 with `targets`, `tarchetypes`, `terra`, `sf`, `classInt`, `yaml`, `data.table`, `dbscan`,
  `cluster`, `aricode`, `e1071`, `pROC`, `lme4`, `r5r` (Java 21) and, for figures, `ggplot2`,
  `patchwork`, `tidyterra`, `ggspatial`, `scico`, `ggnewscale`, `ragg`.
- [SurtGIS](https://github.com/franciscoparrao/surtgis) 1.5.1 (terrain analysis: filling, D8,
  flow accumulation, HAND); DOI 10.5281/zenodo.23042861.
- Python 3 with AvaFrame 2.1 (Flow-Py), `rasterio` and `numba`.
- `osmium` and GDAL command-line tools.

Copy `.Renviron.example` to `.Renviron` and set the paths to the SurtGIS binary and to the Python
environment with AvaFrame.

## Input data

Inputs are not redistributed here; place them under `datos/entrada/` with the paths given in
`config.yml`:

| Input | Source |
|---|---|
| Commune boundary; critical infrastructure (health, police, schools, power lines) | Official boundary; IDE Chile |
| 2024 Census microdata and cartography | Instituto Nacional de Estadísticas (INE) |
| OpenStreetMap extract (Chile, 2026-09-28) | Geofabrik |
| Building footprints | Google Open Buildings |
| Sentinel-2 Level-2A composite, 2024 (vulnerability typology only) | Copernicus, via SurtGIS STAC |
| DEM (Copernicus GLO-30) and FABDEM v1.2 | Copernicus; University of Bristol |
| Initiation susceptibility maps (XGBoost, RF, LightGBM ensemble) | Companion study; Zenodo DOI in the manuscript |
| Areas affected in March 2015 (six communes) | Sernageomin; the digitised Taltal footprint is available from the corresponding author on request |
| Official risk zoning of the Communal Regulatory Plan of Taltal (2012) | MINVU Geoportal Open Data |

The derived data used in the analysis (indices per census unit, block-level predictors for the six
communes, reach index raster and all result tables) are archived at Zenodo (DOI in the manuscript).

## Running

```bash
Rscript run.R corregido                 # full pipeline; results in salidas/corregido/
Rscript paper/R/todas_las_figuras.R     # figures in paper/figs/
```

`targets` only rebuilds what changed. The first run takes about one hour; heavy steps (isochrones,
SurtGIS) run under memory limits set in `.Renviron`.

## License

Code: MIT (see `LICENSE`). Please cite the article and the archived software release
(`CITATION.cff`).
