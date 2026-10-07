# Validación externa (mismo evento, 25-26/03/2015, región de Atacama): ¿la cercanía al cauce
# discrimina manzanas afectadas en otras cuencas, sin recalibrar? Área afectada: catastro
# Sernageomin (clase "Área afectada"); DEM: Copernicus GLO-30 (Planetary Computer, vía SurtGIS);
# unidades: manzanas pobladas del Censo 2024. Mismos parámetros que en Taltal.

comunas_validacion <- data.frame(
  comuna = c("TALTAL", "CHAÑARAL", "DIEGO DE ALMAGRO", "TIERRA AMARILLA", "COPIAPO", "ALTO DEL CARMEN"),
  cut = c(2104, 3201, 3202, 3103, 3101, 3302)
)

afectada_comuna <- function(catastro, comuna, epsg, area_2015 = NULL) {
  if (comuna == "TALTAL") return(sf::st_union(sf::st_read(area_2015, quiet = TRUE)))
  x <- sf::st_read(catastro, quiet = TRUE)
  x <- x[x$NOMBRE == "Área afectada" & x$COMUNA == comuna, ]
  sf::st_union(sf::st_transform(sf::st_make_valid(x), epsg))
}

manzanas_comuna <- function(censo_gpkg, cut, epsg) {
  m <- sf::st_read(censo_gpkg, quiet = TRUE,
                   query = sprintf("SELECT n_per, SHAPE FROM Manzanas_CPV24 WHERE CUT=%d AND n_per > 0", cut))
  sf::st_transform(sf::st_make_valid(m), epsg)
}

# DEM Copernicus -> UTM 30 m -> cauces (acumulación > umbral) -> distancia al cauce y HAND.
terreno_comuna <- function(bbox_utm, dir, p, umbral = 2000, margen_m = 15000, dem_local = NULL) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  bb <- sf::st_bbox(bbox_utm)
  bb <- sf::st_bbox(c(xmin = bb[["xmin"]] - margen_m, ymin = bb[["ymin"]] - margen_m,
                      xmax = bb[["xmax"]] + margen_m, ymax = bb[["ymax"]] + margen_m), crs = sf::st_crs(bbox_utm))
  b <- sf::st_bbox(sf::st_transform(sf::st_as_sfc(bb), 4326))
  caja <- sprintf("%.4f,%.4f,%.4f,%.4f", b[1], b[2], b[3], b[4])
  f <- function(n) file.path(dir, n)
  # Caché: DEM e hidrología base se reutilizan si la fuente y el recuadro no cambiaron.
  clave <- paste(if (is.null(dem_local)) "cop-dem-glo-30" else normalizePath(dem_local), caja)
  if (!file.exists(f("clave.txt")) || !identical(readLines(f("clave.txt")), clave) ||
      !file.exists(f("facc.tif"))) {
    if (is.null(dem_local)) {
      surtgis("stac", "fetch-mosaic", "--catalog", "pc", "--collection", "cop-dem-glo-30", "--asset", "data",
              paste0("--bbox=", caja), f("dem_4326.tif"))
    } else {
      surtgis("clip", paste0("--bbox=", caja), dem_local, f("dem_4326.tif"))
    }
    surtgis("reproject", "--to", p$crs, "--pixel-size", 30, "--method", "bilinear", f("dem_4326.tif"), f("dem.tif"))
    surtgis("hydrology", "fill-sinks", f("dem.tif"), f("dem_fill.tif"))
    surtgis("hydrology", "flow-direction", f("dem_fill.tif"), f("fdir.tif"))
    surtgis("hydrology", "flow-accumulation", f("fdir.tif"), f("facc.tif"))
    writeLines(clave, f("clave.txt"))
  }
  # Productos por umbral de cauce (celdas de área contribuyente)
  hand_u <- f(sprintf("hand_costa_%d.tif", umbral)); dist_u <- f(sprintf("dist_cauce_%d.tif", umbral))
  if (!file.exists(hand_u) || !file.exists(dist_u) || file.mtime(hand_u) < file.mtime(f("facc.tif"))) {
    surtgis("hydrology", "hand", "--threshold", umbral, f("dem_fill.tif"), f(sprintf("hand_%d.tif", umbral)))
    cauce <- terra::classify(terra::rast(f("facc.tif")) > umbral, cbind(0, NA))
    terra::writeRaster(terra::distance(cauce), dist_u, overwrite = TRUE)
    # HAND queda NA donde el drenaje llega al mar/borde sin tocar un cauce (celdas costeras):
    # ahí el nivel base es el mar -> HAND = elevación (m s.n.m.).
    dem <- terra::rast(f("dem_fill.tif"))
    terra::writeRaster(terra::cover(terra::rast(f(sprintf("hand_%d.tif", umbral))), terra::clamp(dem, lower = 0)),
                       hand_u, overwrite = TRUE)
  }
  c(dist = dist_u, hand = hand_u)
}

# Rutas solo ASCII: SurtGIS 1.5.1 no escribe la provenance si la ruta tiene ñ/tildes (issue #10).
carpeta_ascii <- function(x) gsub("[^a-z0-9]+", "_", tolower(iconv(x, to = "ASCII//TRANSLIT")))

datos_comuna <- function(comuna, cut, catastro, censo_gpkg, area_2015, p, umbral = 2000, dem_local = NULL,
                         carpeta = carpeta_ascii(comuna)) {
  a <- afectada_comuna(catastro, comuna, p$crs, area_2015)
  m <- manzanas_comuna(censo_gpkg, cut, p$crs)
  if (comuna != "TALTAL") {
    # evaluar solo localidades tocadas por el área afectada (manzanas a <= 5 km de ella)
    m <- m[as.numeric(sf::st_distance(m, a)) <= 5000, ]
  }
  t <- terreno_comuna(sf::st_union(c(sf::st_geometry(m), a)),
                      salida(p, "validacion_externa", carpeta, "x") |> dirname(), p,
                      umbral = umbral, dem_local = dem_local)
  pts <- sf::st_point_on_surface(sf::st_geometry(m))
  v <- terra::extract(terra::rast(unname(t)), terra::vect(pts), ID = FALSE)
  names(v) <- names(t)
  data.frame(comuna = comuna, umbral = umbral, poblacion = m$n_per, hand = v$hand, dist = v$dist,
             afectada = lengths(sf::st_intersects(m, a)) > 0)
}

auc_tabla <- function(d, scores, semilla) {
  do.call(rbind, lapply(split(d, d$comuna), function(x) do.call(rbind, lapply(names(scores), function(n) {
    r <- auc_ic(scores[[n]](x), x$afectada, semilla = semilla)
    data.frame(comuna = x$comuna[1], metrica = n, n_manzanas = nrow(x), n_afectadas = sum(x$afectada),
               pob_afectada = sum(x$poblacion[x$afectada]), auc = round(r[["auc"]], 3),
               ic95_inf = round(r[["ic_inf"]], 3), ic95_sup = round(r[["ic_sup"]], 3),
               n_sin_dato = r[["n_excluidos"]])
  }))))
}

# Probabilidad de alcance: logística afectada ~ log1p(HAND) + log1p(dist), ajustada SOLO con
# comunas de calibración (Atacama); se reporta validación dejando una comuna fuera (LOCO).
formula_alcance <- afectada ~ log1p(hand) + log1p(dist)

calibrar_alcance <- function(d, calibracion) {
  dc <- d[d$comuna %in% calibracion & is.finite(d$hand) & is.finite(d$dist), ]
  loco <- do.call(rbind, lapply(calibracion, function(cm) {
    fit <- stats::glm(formula_alcance, stats::binomial(), dc[dc$comuna != cm, ])
    x <- dc[dc$comuna == cm, ]
    data.frame(comuna = cm, auc_loco = round(auc(stats::predict(fit, x), x$afectada), 3))
  }))
  list(modelo = stats::glm(formula_alcance, stats::binomial(), dc), loco = loco)
}

# Datos por manzana de todas las comunas del evento (descarga DEM por comuna; ~2 min).
validacion_externa_datos <- function(catastro, censo_gpkg, area_2015, p, umbral = 2000) {
  do.call(rbind, lapply(seq_len(nrow(comunas_validacion)), function(i)
    datos_comuna(comunas_validacion$comuna[i], comunas_validacion$cut[i], catastro, censo_gpkg, area_2015, p,
                 umbral = umbral)))
}

guardar_csv <- function(d, p, ...) {
  destino <- salida(p, ...)
  utils::write.csv(d, destino, row.names = FALSE)
  existe(destino)
}
