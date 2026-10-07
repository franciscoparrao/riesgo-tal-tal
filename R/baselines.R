# Líneas base topográficas triviales (revisión ciega, Landslides): ¿HAND aporta algo sobre la
# elevación, la pendiente o la distancia a la costa? En Taltal HAND = elevación en la plataforma
# costera (el mar es el nivel base), así que la elevación es la línea base obligada.
#   elev_neg        -elevación (m s.n.m., DEM Copernicus rellenado)
#   pendiente_neg   -pendiente (grados)
#   dist_costa_neg  -distancia a la costa (celdas con elevación <= 0 en el DEM sin rellenar)
# Mayor valor = más amenaza en todas.

rasters_triviales <- function(dir_dem, dir_out, wo = list(datatype = "FLT4S", gdal = c("COMPRESS=DEFLATE", "TILED=YES"))) {
  dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
  f <- function(n) file.path(dir_out, n)
  dem <- terra::rast(file.path(dir_dem, "dem_fill.tif"))
  terra::app(dem, function(x) -x, filename = f("elev_neg.tif"), overwrite = TRUE, wopt = wo)
  terra::app(terra::terrain(dem, "slope", unit = "degrees"), function(x) -x,
             filename = f("pendiente_neg.tif"), overwrite = TRUE, wopt = wo)
  mar <- terra::classify(terra::rast(file.path(dir_dem, "dem.tif")) <= 0, cbind(0, NA))
  if (isTRUE(terra::global(mar, "notNA")[[1]] > 0)) {
    terra::app(terra::distance(mar), function(x) -x, filename = f("dist_costa_neg.tif"), overwrite = TRUE, wopt = wo)
  }
  f(c("elev_neg.tif", "pendiente_neg.tif", "dist_costa_neg.tif"))
}

# Taltal (grilla de la comuna, la misma de la amenaza), para validar por unidad poblada.
capas_triviales_taltal <- function(amenaza, p) {
  stopifnot(length(amenaza) > 0)
  terra::terraOptions(memfrac = 0.1, progress = 0, tempdir = dirname(salida(p, "tmp", "terra", "x")))
  existe(rasters_triviales(dirname(salida(p, "amenaza", "x")), dirname(salida(p, "baselines", "x"))))
}

# Por manzana en las 6 comunas del evento (mismas manzanas y DEM que la validación externa).
auc_triviales_comunas <- function(catastro, censo_gpkg, area_2015, manzanas_validacion, p, umbral = 2000) {
  stopifnot(nrow(manzanas_validacion) > 0)      # dependencia: DEM por comuna ya descargado
  terra::terraOptions(memfrac = 0.1, progress = 0, tempdir = dirname(salida(p, "tmp", "terra", "x")))
  do.call(rbind, lapply(seq_len(nrow(comunas_validacion)), function(i) {
    comuna <- comunas_validacion$comuna[i]
    a <- afectada_comuna(catastro, comuna, p$crs, area_2015)
    m <- manzanas_comuna(censo_gpkg, comunas_validacion$cut[i], p$crs)
    if (comuna != "TALTAL") m <- m[as.numeric(sf::st_distance(m, a)) <= 5000, ]
    dir <- dirname(salida(p, "validacion_externa", carpeta_ascii(comuna), "x"))
    r <- rasters_triviales(dir, file.path(dir, "triviales"))
    capas <- c(HAND = file.path(dir, sprintf("hand_costa_%d.tif", umbral)),
               "distance to channel" = file.path(dir, sprintf("dist_cauce_%d.tif", umbral)),
               elevation = r[1], slope = r[2], "distance to coast" = r[3])
    signo <- c(-1, -1, 1, 1, 1)                  # hand y dist crudos; triviales ya negados
    pts <- terra::vect(sf::st_point_on_surface(sf::st_geometry(m)))
    clase <- lengths(sf::st_intersects(m, a)) > 0
    do.call(rbind, lapply(seq_along(capas), function(k) {
      if (!file.exists(capas[[k]])) return(NULL)
      v <- signo[k] * terra::extract(terra::rast(capas[[k]]), pts, ID = FALSE)[[1]]
      if (sum(is.finite(v)) < 10) return(NULL)
      ic <- auc_ic(v, clase, semilla = p$semilla)
      data.frame(comuna = comuna, predictor = names(capas)[k], n = sum(is.finite(v)), afectadas = sum(clase),
                 auc = round(ic[["auc"]], 3), ic95_inf = round(ic[["ic_inf"]], 3), ic95_sup = round(ic[["ic_sup"]], 3))
    }))
  }))
}
