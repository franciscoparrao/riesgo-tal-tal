# Plantilla raster (valor 0) sobre la extensión de un vector, para rasterize.
# Se crea con gdal_create (sin cargar celdas en memoria: a 10 m son 293 M celdas).
# Resolución exacta como `terra::rast(ext, res = 10)` de la implementación de referencia: el origen es la
# esquina superior izquierda y la extensión se amplía hacia la derecha y abajo.
plantilla <- function(vector, resolucion, destino) {
  v <- terra::vect(vector)
  e <- as.vector(terra::ext(v))
  ncol <- ceiling((e[["xmax"]] - e[["xmin"]]) / resolucion)
  nrow <- ceiling((e[["ymax"]] - e[["ymin"]]) / resolucion)
  xmax <- e[["xmin"]] + ncol * resolucion
  ymin <- e[["ymax"]] - nrow * resolucion
  epsg <- terra::crs(v, describe = TRUE)$code
  args <- c("-of", "GTiff", "-ot", "Byte", "-burn", "0", "-a_srs", paste0("EPSG:", epsg),
            "-outsize", ncol, nrow,
            "-a_ullr", sprintf("%.10f", c(e[["xmin"]], e[["ymax"]], xmax, ymin)),
            "-co", "COMPRESS=DEFLATE", "-co", "TILED=YES", "-q", destino)
  res <- suppressWarnings(system2("gdal_create", args, stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(res, "status"))) stop("gdal_create falló:\n", paste(res, collapse = "\n"))
  existe(destino)
}
