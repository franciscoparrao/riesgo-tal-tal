#!/usr/bin/env python3
"""HAND con WhiteboxTools (control independiente de SurtGIS) sobre los mismos DEM por comuna.

Mismo protocolo que R/validacion_externa.R: DEM Copernicus GLO-30 a 30 m (dem.tif), relleno de
depresiones (Planchon-Darboux), D8, cauces iniciados a 2000 celdas y HAND; la regla costera (HAND = elevación donde el
drenaje llega al mar sin tocar un cauce) se aplica después, en R.

Uso: hand_whitebox.py <dir_validacion_externa> <comuna> [<comuna> ...]
Salida: <dir>/<comuna>/wbt/hand_wbt_2000.tif
"""
import os
import sys

import numpy as np
import rasterio
import whitebox

UMBRAL = 2000


def main(base, comunas):
    wbt = whitebox.WhiteboxTools()
    wbt.set_verbose_mode(False)
    for c in comunas:
        d = os.path.abspath(os.path.join(base, c, "wbt"))
        os.makedirs(d, exist_ok=True)
        f = lambda n: os.path.join(d, n)
        # WhiteboxTools no reconoce NaN como NoData (el DEM de SurtGIS lo usa): sin esto no ve el
        # mar ni el borde como salida y rellena cuencas costeras enteras.
        dem = f("dem_nd.tif")
        with rasterio.open(os.path.join(base, c, "dem.tif")) as r:
            z = r.read(1); prof = r.profile
        z = np.where(np.isfinite(z), z, -9999.0).astype("float32")
        prof.update(nodata=-9999.0, compress="deflate")
        with rasterio.open(dem, "w", **prof) as w:
            w.write(z, 1)
        # fill_depressions (priority-flood) de WBT 2.4 eleva incluso el mar (~71 m) en Chañaral con
        # este DEM; Planchon-Darboux da el relleno esperado y es el análogo de fill-sinks de SurtGIS.
        wbt.fill_depressions_planchon_and_darboux(dem, f("dem_fill.tif"))
        wbt.d8_pointer(f("dem_fill.tif"), f("d8.tif"))
        wbt.d8_flow_accumulation(f("d8.tif"), f("facc.tif"), out_type="cells", pntr=True)
        wbt.extract_streams(f("facc.tif"), f("cauces.tif"), threshold=UMBRAL, zero_background=False)
        wbt.elevation_above_stream(f("dem_fill.tif"), f("cauces.tif"), f(f"hand_wbt_{UMBRAL}.tif"))
        print(c, "ok", flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
