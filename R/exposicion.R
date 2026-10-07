# Exposición, modo corregido. Tres componentes normalizados 0–1:
#  1. Densidad poblacional por zona censal (Censo 2017, taltal_geo): hab / km².
#  2. Densidad de construcciones por zona: edificios de Google Open Buildings / km².
#  3. Uso de suelo: industrial 3, comercial 2, residencial 1, resto 0 (÷ 3).
#     Dynamic World no tiene esas clases (solo "built"); se usa landuse de OSM, que sí
#     las distingue.
# Exposición = promedio de los tres, reescalado a 0–1.
# Todo es constante por pieza (zona censal × clase de uso), así que se calcula en vectorial
# y se rasteriza una sola vez con SurtGIS (exacto y sin cargar la grilla de 10 m en memoria).

minmax <- function(x) (x - min(x)) / (max(x) - min(x))

# Normaliza una densidad: "lineal" (min-max, índice convencional) o "log" (min-max de log1p).
normalizar_densidad <- function(x, transformacion = "lineal") {
  if (identical(transformacion, "log")) minmax(log1p(x)) else minmax(x)
}

# Unidades censales con población, en EPSG `epsg`. Particionan la comuna: donde no hay
# unidad (desierto) se agrega un polígono "resto" con población 0.
zonas_censales <- function(censo, epsg, comuna_utm) {
  if (identical(censo$version, 2017L) || identical(censo$version, 2017)) {
    z <- sf::st_transform(sf::st_make_valid(sf::st_read(censo$ruta, quiet = TRUE)), epsg)
    edades <- grep("^X?(0_15|16_25|26_50|51_64|64_mas)$", names(z), value = TRUE)
    if (length(edades) != 5) stop("No se encontraron las 5 columnas de edad en ", censo$ruta)
    z$poblacion <- rowSums(sf::st_drop_geometry(z)[, edades])
  } else {
    leer <- function(capa) sf::st_read(censo$ruta, quiet = TRUE,
                                       query = sprintf("SELECT * FROM %s WHERE CUT=%d", capa, censo$cut))
    z <- do.call(rbind, lapply(censo$capas, function(capa) {
      x <- leer(capa)
      sf::st_sf(capa = capa, poblacion = x$n_per, geometry = sf::st_geometry(x))
    }))
    z <- sf::st_transform(sf::st_make_valid(z), epsg)
  }
  com <- sf::st_union(sf::st_read(comuna_utm, quiet = TRUE))
  z <- suppressWarnings(sf::st_intersection(z[, "poblacion"], com))
  resto <- sf::st_difference(com, sf::st_union(z))
  if (length(resto) > 0 && !all(sf::st_is_empty(resto))) {
    z <- rbind(z, sf::st_sf(poblacion = 0, geometry = resto))
  }
  z$km2 <- as.numeric(sf::st_area(z)) / 1e6
  z$zona <- seq_len(nrow(z))
  z[, c("zona", "poblacion", "km2")]
}

contar_edificios <- function(zonas, buildings_csv_gz, confianza_min) {
  b <- data.table::fread(cmd = paste("zcat", shQuote(buildings_csv_gz)),
                         select = c("latitude", "longitude", "confidence"))
  b <- b[b$confidence >= confianza_min, ]
  p <- sf::st_transform(sf::st_as_sf(b, coords = c("longitude", "latitude"), crs = 4326),
                        sf::st_crs(zonas))
  lengths(sf::st_intersects(zonas, p))
}

# Clases de uso de suelo disjuntas (la de mayor puntaje gana donde se superponen).
uso_suelo_osm <- function(pbf, comuna_utm, puntajes) {
  x <- sf::st_read(pbf, query = "SELECT landuse FROM multipolygons WHERE landuse IS NOT NULL",
                   quiet = TRUE)
  com <- sf::st_union(sf::st_read(comuna_utm, quiet = TRUE))
  x <- sf::st_transform(sf::st_make_valid(x), sf::st_crs(com))
  x$puntaje <- unname(unlist(puntajes)[x$landuse])
  x$puntaje[is.na(x$puntaje)] <- 0
  capas <- list(); ocupado <- NULL
  for (s in sort(unique(c(x$puntaje[x$puntaje > 0])), decreasing = TRUE)) {
    g <- sf::st_intersection(sf::st_union(x[x$puntaje == s, ]), com)
    if (!is.null(ocupado) && length(g) > 0) g <- sf::st_difference(g, ocupado)
    if (length(g) == 0 || all(sf::st_is_empty(g))) next
    capas[[length(capas) + 1]] <- sf::st_sf(puntaje = s, geometry = g)
    ocupado <- if (is.null(ocupado)) g else sf::st_union(ocupado, g)
  }
  resto <- if (is.null(ocupado)) com else sf::st_difference(com, ocupado)
  do.call(rbind, c(capas, list(sf::st_sf(puntaje = 0, geometry = resto))))
}

exposicion_vector <- function(censo_ruta, buildings, pbf, comuna_utm, cfg) {
  destino <- salida(cfg, "componentes", "exposicion", "exposicion_piezas.gpkg")
  censo <- utils::modifyList(cfg$censo, list(ruta = censo_ruta))
  z <- zonas_censales(censo, cfg$crs, comuna_utm)
  z$edificios <- contar_edificios(z, buildings, cfg$exposicion$confianza_edificios)
  tr <- cfg$exposicion$transformacion_densidad
  z$dens_pob_n <- normalizar_densidad(z$poblacion / z$km2, tr)
  z$dens_edif_n <- normalizar_densidad(z$edificios / z$km2, tr)
  uso <- uso_suelo_osm(pbf, comuna_utm, cfg$exposicion$puntajes_uso)
  piezas <- suppressWarnings(sf::st_intersection(z, uso))
  piezas <- piezas[!sf::st_is_empty(piezas), ]
  piezas$uso_n <- piezas$puntaje / max(unlist(cfg$exposicion$puntajes_uso))
  bruto <- (piezas$dens_pob_n + piezas$dens_edif_n + piezas$uso_n) / 3
  piezas$exposicion <- bruto / max(bruto)
  sf::st_write(piezas, destino, delete_dsn = TRUE, quiet = TRUE)
  existe(destino)
}

exposicion_raster <- function(piezas, plantilla_fina, cfg) {
  destino <- salida(cfg, "componentes", "exposicion.tif")
  surtgis("rasterize", "--reference", plantilla_fina, "--attribute", "exposicion", piezas, destino)
  existe(destino)
}
