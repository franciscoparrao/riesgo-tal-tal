# Modelos de propagación como alternativa al índice de alcance (HAND), evaluados en Taltal con el
# mismo protocolo que el resto de las capas (AUC por unidad poblada, R/validacion_2015.R).
#   susc_aguas_arriba   susceptibilidad de inicio promedio del área que drena a la celda
#   fuente_susceptible  área fuente susceptible total aguas arriba, log(1 + suma de S)
#   laharz_graduado     nº de volúmenes LAHARZ (1e4-1e7 m3, flujo de detritos) que inundan la celda
#   cono_graduado       nº de conos de energía (H/L 1,5-4°) que alcanzan la celda
#   hand, dist_cauce    referencias (con signo invertido: mayor = más amenaza)
# LAHARZ y el cono parten de UNA fuente: el cauce principal de la Quebrada de Taltal donde entra
# al marco del mapa de 2015. Conocer ese punto es información del evento, así que ambos son una
# cota superior de lo que lograría un modelo de propagación sin conocer el evento.

capas_propagacion <- function(amenaza, marco, p, umbral = 2000,
                              volumenes = c(1e4, 3e4, 1e5, 3e5, 1e6, 3e6, 1e7),
                              angulos = c(1.5, 2, 2.5, 3, 4), entrada = c(354000, 7185000)) {
  stopifnot(length(amenaza) > 0)                      # dependencia: los intermedios de amenaza_flujo
  dir_a <- dirname(salida(p, "amenaza", "x"))
  dir <- dirname(salida(p, "propagacion", "x"))
  f <- function(n) file.path(dir, n)
  wo <- list(datatype = "FLT4S", gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
  terra::terraOptions(memfrac = 0.1, progress = 0, tempdir = dirname(salida(p, "tmp", "terra", "x")))

  # Referencias y propagación de la susceptibilidad (grilla de la comuna)
  terra::app(terra::rast(file.path(dir_a, sprintf("hand_costa_%d.tif", umbral))), function(x) -x,
             filename = f("hand_neg.tif"), overwrite = TRUE, wopt = wo)
  terra::app(terra::rast(file.path(dir_a, sprintf("dist_cauce_%d.tif", umbral))), function(x) -x,
             filename = f("dist_cauce_neg.tif"), overwrite = TRUE, wopt = wo)
  terra::app(terra::rast(file.path(dir_a, "acc_s.tif")), log1p,
             filename = f("fuente_susceptible.tif"), overwrite = TRUE, wopt = wo)

  # LAHARZ y cono de energía en un recorte de la ciudad (mismo DEM Copernicus rellenado)
  m <- sf::st_read(marco, quiet = TRUE)
  caja <- terra::ext(sf::st_bbox(sf::st_buffer(sf::st_union(m), 3000))[c("xmin", "xmax", "ymin", "ymax")])
  terra::crop(terra::rast(file.path(dir_a, "dem_fill.tif")), caja, filename = f("dem_recorte.tif"),
              overwrite = TRUE, wopt = wo)
  surtgis("hydrology", "fill-sinks", f("dem_recorte.tif"), f("dem_fill.tif"))
  surtgis("hydrology", "flow-direction", f("dem_fill.tif"), f("fdir.tif"))
  surtgis("hydrology", "flow-accumulation", f("fdir.tif"), f("facc.tif"))
  acc <- terra::rast(f("facc.tif"))
  zona <- terra::vect(sf::st_buffer(sf::st_sfc(sf::st_point(entrada), crs = p$crs), 150))
  cand <- terra::mask(terra::crop(acc, zona), zona)
  cel <- terra::cellFromXY(acc, terra::xyFromCell(cand, which.max(terra::values(cand))))
  rc <- terra::rowColFromCell(acc, cel)
  fuente <- sprintf("%d,%d", rc[1, 1] - 1, rc[1, 2] - 1)          # SurtGIS: fila,columna base 0

  suma <- function(archivos, destino) {
    s <- Reduce(`+`, lapply(archivos, function(a) terra::classify(terra::rast(a) > 0, cbind(NA, 0))))
    terra::writeRaster(s, destino, overwrite = TRUE, wopt = wo)
  }
  lz <- vapply(volumenes, function(v) {
    o <- f(sprintf("laharz_%g.tif", v))
    surtgis("hydrology", "laharz", "--flow-dir", f("fdir.tif"), "--source", fuente,
            "--volume", format(v, scientific = FALSE), "--flow-type", "debris-flow", f("dem_fill.tif"), o)
    o
  }, "")
  co <- vapply(angulos, function(a) {
    o <- f(sprintf("cono_%g.tif", a))
    surtgis("hydrology", "energy-cone", "--source", fuente, "--cone-angle", a, f("dem_fill.tif"), o)
    o
  }, "")
  suma(lz, f("laharz_graduado.tif")); suma(co, f("cono_graduado.tif"))
  existe(c(file.path(dir_a, "susc_aguas_arriba_30m.tif"),
           f(c("fuente_susceptible.tif", "laharz_graduado.tif", "cono_graduado.tif", "hand_neg.tif",
               "dist_cauce_neg.tif"))))
}
