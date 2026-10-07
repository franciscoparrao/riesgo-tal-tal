#!/usr/bin/env python3
"""Modelo de fricción simple con tope de velocidad (SFLM de Flow-R; Horton et al. 2013, Eq. 10)
propagado a lo largo de D8, y utilidades de cuencas. Entorno: ~/.local/opt/avaframe-env (numba).

  K_i = min(K_0 + Δz - Δx tan(α), h_max),  h_max = v_max^2 / (2 g)   (altura cinética, m)

Un flujo sigue mientras K >= 0 a lo largo de D8; lateralmente alcanza las celdas cuya altura sobre
el cauce al que drenan (HAND, cauces >= 2000 celdas) no supera la altura de energía K en ese cauce. Para cada celda se guarda el MAYOR ángulo α (de una grilla) con el
que la alcanza el flujo de alguna fuente: "ángulo crítico". Sin tope (h_max = inf) es el ángulo de
alcance clásico (H/L) de la línea de energía, con la condición de que se mantenga sobre el terreno.

Subcomandos:
  cuencas  <dem_fill> <fdir> <salida_etiquetas.tif> <afectada1.gpkg> [...]
           etiqueta (1..n) las celdas aguas arriba de cada área afectada; imprime bbox (UTM y 4326)
  angulo   <dem_fill> <fdir> <facc> <mascara_cuenca.tif|-> <salida.tif> <vmax|inf> [fuentes.tif]
           fuentes: criterio Flow-R de eventos extremos (Horton et al. 2013, Eq. 2) si no se da raster
"""
import json
import math
import sys

import numpy as np
import rasterio
from numba import njit
from rasterio import features
import geopandas as gpd
from pyproj import Transformer

# SurtGIS D8: 0 = sin salida; 1-8 en sentido antihorario desde el este (fila crece hacia el sur)
DR = np.array([0, 0, -1, -1, -1, 0, 1, 1, 1], dtype=np.int64)
DC = np.array([0, 1, 1, 0, -1, -1, -1, 0, 1], dtype=np.int64)
DIAG = np.array([0, 0, 1, 0, 1, 0, 1, 0, 1], dtype=np.int64)


@njit(cache=True)
def aguas_abajo(fd, nr, nc):
    n = nr * nc
    dn = -np.ones(n, dtype=np.int64)
    for i in range(nr):
        for j in range(nc):
            k = fd[i, j]
            if k < 1 or k > 8:
                continue
            ii = i + DR[k]; jj = j + DC[k]
            if 0 <= ii < nr and 0 <= jj < nc:
                dn[i * nc + j] = ii * nc + jj
    return dn


@njit(cache=True)
def aguas_arriba(dn, semillas, n):
    # BFS inverso: lista de padres por celda (CSR)
    cnt = np.zeros(n + 1, dtype=np.int64)
    for c in range(n):
        d = dn[c]
        if d >= 0:
            cnt[d + 1] += 1
    for c in range(n):
        cnt[c + 1] += cnt[c]
    padres = np.empty(cnt[n], dtype=np.int64)
    pos = cnt.copy()
    for c in range(n):
        d = dn[c]
        if d >= 0:
            padres[pos[d]] = c
            pos[d] += 1
    marca = np.zeros(n, dtype=np.bool_)
    cola = np.empty(n, dtype=np.int64)
    h = 0; t = 0
    for s in semillas:
        if not marca[s]:
            marca[s] = True; cola[t] = s; t += 1
    while h < t:
        c = cola[h]; h += 1
        for q in range(cnt[c], cnt[c + 1]):
            p = padres[q]
            if not marca[p]:
                marca[p] = True; cola[t] = p; t += 1
    return marca


@njit(cache=True)
def cauce_destino(dn, orden_desc, es_cauce):
    # celda de cauce a la que drena cada celda (la misma si es cauce); -1 si no llega a uno
    n = dn.shape[0]
    t = -np.ones(n, dtype=np.int64)
    for c in orden_desc:
        if es_cauce[c]:
            t[c] = c
        else:
            d = dn[c]
            if d >= 0:
                t[c] = t[d]
    return t


@njit(cache=True)
def angulo_critico(z, dn, largo, orden, fuente, activo, angulos, hmax, destino):
    n = z.shape[0]
    crit = np.full(n, np.nan, dtype=np.float32)
    K = np.empty(n, dtype=np.float64)
    for a in angulos:
        t = math.tan(math.radians(a))
        for c in range(n):
            K[c] = -1.0
        for c in orden:
            if fuente[c] and K[c] < 0.0:
                K[c] = 0.0
            if K[c] < 0.0:
                continue
            d = dn[c]
            if d < 0 or not activo[d]:
                continue
            k = K[c] + (z[c] - z[d]) - largo[c] * t
            if k > hmax:
                k = hmax
            if k >= 0.0 and k > K[d]:
                K[d] = k
        # alcance lateral: la celda c se alcanza si la altura de energía del flujo en el cauce al
        # que drena supera su altura sobre ese cauce (HAND): K[t] >= z[c] - z[t]
        for c in range(n):
            if not activo[c]:
                continue
            tt = destino[c]
            if tt < 0 or K[tt] < 0.0:
                continue
            if K[tt] >= z[c] - z[tt]:
                if math.isnan(crit[c]) or a > crit[c]:
                    crit[c] = a
    return crit


def leer(ruta):
    with rasterio.open(ruta) as r:
        return r.read(1), r.profile, r.transform


def cmd_cuencas(dem, fdir, salida, afectadas):
    z, prof, tr = leer(dem)
    fd, _, _ = leer(fdir)
    nr, nc = z.shape
    dn = aguas_abajo(fd.astype(np.int64), nr, nc)
    etiq = np.zeros(nr * nc, dtype=np.uint8)
    t4326 = Transformer.from_crs(prof["crs"], 4326, always_xy=True)
    info = {}
    for i, g in enumerate(afectadas, start=1):
        a = gpd.read_file(g).to_crs(prof["crs"])
        m = features.rasterize(((geom, 1) for geom in a.geometry), out_shape=z.shape, transform=tr,
                               all_touched=True, dtype="uint8").ravel()
        cu = aguas_arriba(dn, np.flatnonzero(m), nr * nc)
        etiq[cu & (etiq == 0)] = i
        filas, cols = np.divmod(np.flatnonzero(cu), nc)
        x0, y1 = tr * (cols.min(), filas.min()); x1, y0 = tr * (cols.max() + 1, filas.max() + 1)
        lon0, lat0 = t4326.transform(x0, y0); lon1, lat1 = t4326.transform(x1, y1)
        info[g] = {"celdas": int(cu.sum()), "km2": round(cu.sum() * abs(tr.a * tr.e) / 1e6, 1),
                   "bbox_utm": [x0, y0, x1, y1], "bbox_4326": [lon0, lat0, lon1, lat1],
                   "toca_borde": bool(filas.min() == 0 or cols.min() == 0 or filas.max() == nr - 1 or cols.max() == nc - 1)}
    prof.update(dtype="uint8", nodata=0, count=1, compress="deflate")
    with rasterio.open(salida, "w", **prof) as w:
        w.write(etiq.reshape(nr, nc), 1)
    print(json.dumps(info, indent=1))


def fuentes_flowr(z, fd, acc, celda):
    """Criterio de fuentes de Flow-R para eventos extremos (Horton et al. 2013, Eq. 2):
    tan(beta) >= 0.31 A^-0.15 si A < 2.5 km2; >= 0.26 si A >= 2.5 km2; A >= 0.01 km2.
    beta: pendiente local hacia la celda aguas abajo (D8)."""
    nr, nc = z.shape
    A = acc.astype(np.float64) * celda * celda / 1e6
    tanb = np.zeros_like(z, dtype=np.float64)
    for k in range(1, 9):
        m = fd == k
        ii, jj = np.nonzero(m)
        i2 = np.clip(ii + DR[k], 0, nr - 1); j2 = np.clip(jj + DC[k], 0, nc - 1)
        L = celda * (math.sqrt(2) if DIAG[k] else 1.0)
        tanb[ii, jj] = (z[ii, jj] - z[i2, j2]) / L
    umbral = np.where(A < 2.5, 0.31 * np.power(np.maximum(A, 1e-9), -0.15), 0.26)
    return (A >= 0.01) & (tanb >= umbral)


def cmd_angulo(dem, fdir, facc, mascara, salida, vmax, fuentes=None, umbral_cauce=2000):
    z, prof, tr = leer(dem)
    fd, _, _ = leer(fdir)
    acc, _, _ = leer(facc)
    nr, nc = z.shape
    celda = abs(tr.a)
    zf = np.where(np.isfinite(z), z, np.nan).astype(np.float64)
    activo = np.isfinite(zf)
    if mascara != "-":
        m, _, _ = leer(mascara)
        activo &= m > 0
    if fuentes:
        f, _, _ = leer(fuentes)
        fu = (f > 0) & activo
    else:
        fu = fuentes_flowr(zf, fd.astype(np.int64), acc, celda) & activo
    dn = aguas_abajo(fd.astype(np.int64), nr, nc)
    largo = np.full(nr * nc, celda, dtype=np.float64)
    diag = np.isin(fd.ravel(), [2, 4, 6, 8])
    largo[diag] = celda * math.sqrt(2)
    accf = np.where(activo, acc, np.inf).ravel()
    orden = np.argsort(accf, kind="stable")
    orden = orden[np.isfinite(accf[orden])]
    angulos = np.round(np.arange(0.25, 15.01, 0.25), 2)
    es_cauce = (acc.ravel() > umbral_cauce) & activo.ravel()
    destino = cauce_destino(dn, orden[::-1].copy(), es_cauce)
    hmax = math.inf if vmax == "inf" else float(vmax) ** 2 / (2 * 9.81)
    crit = angulo_critico(np.nan_to_num(zf.ravel(), nan=-1e9), dn, largo, orden, fu.ravel(),
                          activo.ravel(), angulos, hmax, destino)
    prof.update(dtype="float32", nodata=np.nan, count=1, compress="deflate")
    with rasterio.open(salida, "w", **prof) as w:
        w.write(crit.reshape(nr, nc), 1)
    print(json.dumps({"fuentes": int(fu.sum()), "celdas_activas": int(activo.sum()),
                      "alcanzadas": int(np.isfinite(crit).sum()), "hmax_m": hmax}))


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "cuencas":
        cmd_cuencas(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5:])
    elif cmd == "angulo":
        cmd_angulo(*sys.argv[2:8], *(sys.argv[8:9] or [None]))
    else:
        sys.exit(__doc__)
