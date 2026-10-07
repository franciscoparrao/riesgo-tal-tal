# LOCC — falta de capacidad de afrontamiento (lack of coping capacity), modo corregido.
# Un indicador binario por tipo de infraestructura crítica: 0 si la celda está dentro de la
# isócrona (tiempo de viaje en auto <= umbral) de alguna instalación de ese tipo, 1 si no.
# LOCC = promedio de los indicadores (el promedio mantiene la escala 0–1 de los demás
# componentes).

osmium_bin <- function() {
  bin <- Sys.getenv("OSMIUM_BIN", path.expand("~/.local/opt/osmium-env/bin/osmium"))
  if (!file.exists(bin)) stop("No se encontró osmium (define OSMIUM_BIN).")
  bin
}

# Recorta el extracto OSM al bbox de la comuna + margen (las rutas pueden salir de ella).
recortar_osm <- function(pbf, comuna, margen_km, cfg) {
  destino <- salida(cfg, "red", "taltal.osm.pbf")
  bb <- sf::st_bbox(sf::st_transform(sf::st_read(comuna, quiet = TRUE), 4326))
  m <- margen_km / 111
  caja <- paste(round(c(bb["xmin"] - m, bb["ymin"] - m, bb["xmax"] + m, bb["ymax"] + m), 5),
                collapse = ",")
  res <- system2(osmium_bin(), c("extract", "-b", caja, "--overwrite", "-o", destino, pbf),
                 stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(res, "status"))) stop("osmium extract falló:\n", paste(res, collapse = "\n"))
  existe(destino)
}

# Construye la red r5r en un proceso aislado (JVM fuera del proceso de targets) y con tope
# de memoria. Reconstruye siempre: si cambió el .pbf, la network.dat en caché quedaría vieja.
red_r5r <- function(pbf_recortado, memoria_max = Sys.getenv("LOCC_MEM_MAX", "4G"), xmx = "2G") {
  dir <- dirname(pbf_recortado)
  unlink(file.path(dir, c("network.dat", "network_settings.json")))
  expr <- sprintf("options(java.parameters='-Xmx%s'); n <- r5r::build_network('%s', verbose=FALSE); r5r::stop_r5(n)",
                  xmx, dir)
  res <- suppressWarnings(system2("systemd-run", c("--user", "--scope", "-q", "-p",
                                                   paste0("MemoryMax=", memoria_max), "-p", "MemorySwapMax=0",
                                                   "Rscript", "-e", shQuote(expr)),
                                  stdout = TRUE, stderr = TRUE))
  if (!file.exists(file.path(dir, "network.dat"))) {
    stop("build_network falló:\n", paste(utils::tail(res, 15), collapse = "\n"), call. = FALSE)
  }
  existe(file.path(dir, "network.dat"))
}

# Lee una capa de infraestructura como puntos WGS84 (las líneas se muestrean cada `paso_m`).
leer_instalaciones <- function(ruta, paso_m = 500) {
  v <- sf::st_zm(sf::st_read(ruta, quiet = TRUE))
  v <- sf::st_transform(v, 32719)
  geom <- sf::st_geometry(v)
  if (any(sf::st_geometry_type(geom) %in% c("LINESTRING", "MULTILINESTRING"))) {
    geom <- sf::st_cast(sf::st_line_sample(sf::st_cast(geom, "LINESTRING"), density = 1 / paso_m),
                        "POINT")
  } else if (any(sf::st_geometry_type(geom) %in% c("POLYGON", "MULTIPOLYGON"))) {
    geom <- sf::st_point_on_surface(geom)
  }
  p <- sf::st_transform(sf::st_sf(geometry = geom), 4326)
  p <- p[!sf::st_is_empty(p), ]
  p$id <- as.character(seq_len(nrow(p)))
  p
}

# Instalaciones recreativas desde OSM (no hay capa oficial en el proyecto).
recreativos_osm <- function(pbf_recortado, comuna, cfg) {
  destino <- salida(cfg, "red", "recreativos_osm.gpkg")
  tags <- c("park", "pitch", "sports_centre", "playground", "stadium", "recreation_ground")
  q <- function(capa) sprintf("SELECT * FROM %s WHERE other_tags LIKE '%%\"leisure\"=>\"%s\"%%'",
                              capa, tags)
  leer <- function(capa, geom_ok) {
    do.call(rbind, lapply(q(capa), function(sql) {
      x <- sf::st_read(pbf_recortado, query = sql, quiet = TRUE)
      if (nrow(x) == 0) return(NULL)
      sf::st_sf(osm_id = x[[1]], geometry = sf::st_geometry(x))
    }))
  }
  pts <- leer("points")
  pol <- leer("multipolygons")
  if (!is.null(pol)) {
    pol <- sf::st_sf(osm_id = pol$osm_id,
                     geometry = sf::st_point_on_surface(sf::st_make_valid(sf::st_geometry(pol))))
  }
  r <- rbind(pts, pol)
  com <- sf::st_transform(sf::st_read(comuna, quiet = TRUE), sf::st_crs(r))
  r <- r[lengths(sf::st_intersects(r, sf::st_buffer(sf::st_union(com), 0))) > 0, ]
  sf::st_write(r, destino, delete_dsn = TRUE, quiet = TRUE)
  existe(destino)
}

# Isócrona (unión sobre todas las instalaciones de un tipo) para un umbral en minutos.
# Corre bin/isocrona.R en un proceso aparte bajo `systemd-run --user --scope -p MemoryMax=`:
# si r5r/JVM se excede, el kernel mata solo ese proceso en vez de congelar el equipo.
isocrona_tipo <- function(red, instalaciones, umbral, nombre, cfg,
                          memoria_max = Sys.getenv("LOCC_MEM_MAX", "4G"), xmx = "2G", lote = 25) {
  destino <- salida(cfg, "red", paste0("isocrona_", nombre, ".gpkg"))
  pts <- salida(cfg, "red", paste0("origenes_", nombre, ".gpkg"))
  sf::st_write(instalaciones, pts, delete_dsn = TRUE, quiet = TRUE)
  args <- c("--user", "--scope", "-q", "-p", paste0("MemoryMax=", memoria_max),
            "-p", "MemorySwapMax=0", "Rscript", "bin/isocrona.R",
            dirname(red), pts, umbral, destino, xmx, lote)
  res <- suppressWarnings(system2("systemd-run", args, stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(res, "status")) || !file.exists(destino)) {
    stop("isócrona ", nombre, " falló (¿límite de memoria ", memoria_max, "?):\n",
         paste(utils::tail(res, 15), collapse = "\n"), call. = FALSE)
  }
  existe(destino)
}

# Indicador binario en la grilla fina: 0 dentro de la isócrona, 1 fuera, NaN fuera de la comuna.
indicador_locc <- function(isocrona, comuna, plantilla_fina, nombre, cfg) {
  pol <- locc_poligonos(comuna, isocrona, cfg, paste0("pol_", nombre, ".gpkg"))
  destino <- salida(cfg, "componentes", "locc", paste0("locc_", nombre, ".tif"))
  surtgis("rasterize", "--reference", plantilla_fina, "--attribute", "locc", pol, destino)
  existe(destino)
}

# Promedio de indicadores (todos en la misma grilla).
# WORKAROUND issues_surtgis.md #8: `imagery calc` de SurtGIS 1.5.1 carga todas las bandas
# (~55 B/celda con 5 bandas -> ~16 GB a 10 m). terra procesa por bloques con memfrac acotado.
# Volver a surtgis("imagery", "calc", ...) cuando calc haga streaming.
locc_promedio <- function(indicadores, cfg) {
  destino <- salida(cfg, "componentes", "locc.tif")
  terra::terraOptions(memfrac = 0.2, progress = 0)
  r <- terra::rast(indicadores)
  terra::app(r, mean, filename = destino, overwrite = TRUE,
             wopt = list(datatype = "FLT4S", gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES")))
  existe(destino)
}
