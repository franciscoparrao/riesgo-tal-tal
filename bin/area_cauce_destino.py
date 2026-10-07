#!/usr/bin/env python3
"""Área aportante del cauce al que drena cada celda (para el Geomorphic Flood Index, GFI;
Samela et al. 2017): A_r = acumulación de la primera celda de cauce (>= umbral celdas) aguas abajo
por D8. Celdas que llegan al mar o al borde sin tocar un cauce: A_r = umbral (nivel base marino,
coherente con HAND = elevación en esas celdas).

Uso: area_cauce_destino.py <dem_fill.tif> <fdir.tif> <facc.tif> <salida.tif> [umbral=2000]
fdir: códigos D8 de SurtGIS (1-8 antihorario desde el este; 0 = sin salida).
Salida: A_r en km2 (float32).
"""
import sys
import numpy as np
import rasterio
from numba import njit

DR = np.array([0, 0, -1, -1, -1, 0, 1, 1, 1], dtype=np.int64)
DC = np.array([0, 1, 1, 0, -1, -1, -1, 0, 1], dtype=np.int64)


@njit(cache=True)
def destino(fd, acc, umbral):
    nr, nc = fd.shape
    n = nr * nc
    dn = -np.ones(n, dtype=np.int64)
    for i in range(nr):
        for j in range(nc):
            k = fd[i, j]
            if 1 <= k <= 8:
                ii = i + DR[k]; jj = j + DC[k]
                if 0 <= ii < nr and 0 <= jj < nc:
                    dn[i * nc + j] = ii * nc + jj
    a = acc.ravel()
    orden = np.argsort(-a)                 # aguas abajo primero (acumulación decreciente)
    res = np.full(n, -1.0)
    for c in orden:
        if a[c] >= umbral:
            res[c] = a[c]
        else:
            d = dn[c]
            res[c] = res[d] if d >= 0 and res[d] > 0 else umbral
    return res


if __name__ == "__main__":
    dem, fdir, facc, out = sys.argv[1:5]
    umbral = float(sys.argv[5]) if len(sys.argv) > 5 else 2000.0
    with rasterio.open(fdir) as r:
        fd = r.read(1).astype(np.int64)
    with rasterio.open(facc) as r:
        acc = np.nan_to_num(r.read(1).astype(np.float64), nan=0.0); prof = r.profile; tr = r.transform
    with rasterio.open(dem) as r:
        z = r.read(1)
    ar = destino(fd, acc, umbral).reshape(fd.shape) * abs(tr.a * tr.e) / 1e6
    ar = np.where(np.isfinite(z), ar, np.nan).astype(np.float32)
    prof.update(dtype="float32", nodata=np.nan, count=1, compress="deflate")
    with rasterio.open(out, "w", **prof) as w:
        w.write(ar, 1)
    print(out, float(np.nanmin(ar)), float(np.nanmax(ar)))
