#!/bin/bash
# Descarga GLO-30 por región, hidrología y SFLM (vmax 15 m/s e inf). Log en run_sflm.log
set -euo pipefail
cd "$(dirname "$0")"
S="prlimit --as=7000000000 $HOME/.local/opt/surtgis-1.5.1/surtgis --compress --output-dtype f32"
PY=$HOME/.local/opt/avaframe-env/bin/python
B=/home/franciscoparrao/proyectos/riesgo_tal_tal/bin/sflm.py
declare -A BBOX=( [taltal]="-70.55,-25.70,-68.94,-24.91"
                  [salado]="-70.74,-27.86,-68.16,-25.52"
                  [copiapo]="-70.59,-28.69,-69.03,-26.61"
                  [huasco]="-70.65,-29.80,-69.64,-28.52" )
declare -A AFECT=( [taltal]="afectada_taltal.gpkg"
                   [salado]="afectada_chanaral.gpkg afectada_diego_de_almagro.gpkg"
                   [copiapo]="afectada_copiapo.gpkg afectada_tierra_amarilla.gpkg"
                   [huasco]="afectada_alto_del_carmen.gpkg" )
for r in taltal huasco copiapo salado; do
  echo "== $r $(date)"
  mkdir -p $r
  [ -f $r/dem30_4326.tif ] || $S stac fetch-mosaic --catalog pc --collection cop-dem-glo-30 --asset data --bbox=${BBOX[$r]} --max-items 60 --retries 8 $r/dem30_4326.tif
  [ -f $r/dem_fill.tif ] || { $S reproject --to 32719 --pixel-size 30 --method bilinear $r/dem30_4326.tif $r/dem.tif
    $S hydrology fill-sinks $r/dem.tif $r/dem_fill.tif
    $S hydrology flow-direction $r/dem_fill.tif $r/fdir.tif
    $S hydrology flow-accumulation $r/fdir.tif $r/facc.tif; }
  $PY $B cuencas $r/dem_fill.tif $r/fdir.tif $r/cuencas.tif ${AFECT[$r]} > $r/cuencas.json
  for v in 15 inf; do
    [ -f $r/angulo_v$v.tif ] || $PY $B angulo $r/dem_fill.tif $r/fdir.tif $r/facc.tif $r/cuencas.tif $r/angulo_v$v.tif $v
  done
done
echo "LISTO $(date)"
