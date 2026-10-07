#!/usr/bin/env python3
"""Flow-Py (AvaFrame com4FlowPy) con rutas propias: propagación tipo Flow-R desde celdas fuente.

Uso: flowpy.py <dem.tif> <release.tif> <workdir> <alpha> <exp> <max_z> <cpus>
Salidas en <workdir>/res_<uid>/: zDelta, cellCounts, travelLengthMax, fpTravelAngleMax.
Se ejecuta con el entorno ~/.local/opt/avaframe-env (AvaFrame 2.1).
"""
import configparser
import pathlib
import sys
from datetime import datetime

import avaframe.com4FlowPy as pkg
from avaframe.com4FlowPy import com4FlowPy

dem, rel, work, alpha, expo, maxz, cpus = sys.argv[1:8]
cfg = configparser.ConfigParser()
cfg.read(pathlib.Path(pkg.__file__).parent / "com4FlowPyCfg.ini")
g = cfg["GENERAL"]
g["alpha"], g["exp"], g["max_z"], g["cpuCount"] = alpha, expo, maxz, cpus
g["tileSize"], g["tileOverlap"] = "15000", "5000"

uid = f"a{alpha}_e{expo}_z{maxz}_{pathlib.Path(rel).stem}"
work = pathlib.Path(work)
res = work / f"res_{uid}"
(res / "temp").mkdir(parents=True, exist_ok=True)
vacio = pathlib.Path("")
p = {
    "workDir": work, "outDir": res, "resDir": res, "tempDir": res / "temp",
    "demPath": pathlib.Path(dem), "releasePath": pathlib.Path(rel), "relIdPath": vacio,
    "infraPath": vacio, "forestPath": vacio, "varUmaxPath": vacio, "varAlphaPath": vacio,
    "varExponentPath": vacio, "deleteTemp": "True", "outputFileFormat": ".tif",
    "outputFiles": "zDelta|cellCounts|travelLengthMax|fpTravelAngleMax",
    "outputNoDataValue": -9999.0, "useCompression": True, "customDirs": "True",
    "uid": uid, "timeString": datetime.now().strftime("%Y%m%d_%H%M%S"),
}
com4FlowPy.com4FlowPyMain(p, g)
print(res)
