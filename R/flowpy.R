# Propagación tipo Flow-R con Flow-Py (AvaFrame com4FlowPy; D'Amboise et al. 2022), sembrada desde
# la susceptibilidad de inicio y sin información del evento: fuentes = 1 % más susceptible de la
# cuenca de la Quebrada de Taltal; parámetros publicados para flujos de detritos (Horton et al.
# 2013): exponente de Holmgren 4, velocidad máxima 15 m/s (altura de energía 11,5 m), ángulo de
# alcance 11° (gruesos) y 7° (finos). Corre bin/flowpy.py en el entorno FLOWPY_PYTHON (~15-50 min
# por ángulo); se reutiliza la corrida si ya existe con los mismos parámetros.

entradas_flowpy <- function(amenaza, marco, p, q = 0.99, entrada = c(353917.3, 7185196)) {
  stopifnot(length(amenaza) > 0)
  dir_a <- dirname(salida(p, "amenaza", "x")); dir <- dirname(salida(p, "flowpy", "x"))
  destino <- file.path(dir, c("dem.tif", sprintf("release_q%02d.tif", round(100 * q))))
  if (all(file.exists(destino))) return(destino)    # caché: borrar para regenerar (y re-correr Flow-Py)
  terra::terraOptions(memfrac = 0.1, progress = 0, tempdir = dirname(salida(p, "tmp", "terra", "x")))
  ws <- terra::classify(terra::watershed(terra::rast(file.path(dir_a, "flowdir_terra.tif")), rbind(entrada)),
                        cbind(0, NA))
  e <- terra::union(terra::ext(terra::trim(ws)), terra::ext(sf::st_read(marco, quiet = TRUE))) + 2000
  dem <- terra::crop(terra::rast(file.path(dir_a, "dem_fill.tif")), e)
  s <- terra::crop(terra::mask(terra::rast(file.path(dir_a, "susc_30m.tif")), ws), e)
  u <- terra::global(s, fun = function(x) stats::quantile(x, q, na.rm = TRUE))[[1]]
  terra::writeRaster(dem, destino[1], overwrite = TRUE, NAflag = -9999, datatype = "FLT4S", gdal = "COMPRESS=DEFLATE")
  terra::writeRaster(terra::classify(s >= u, cbind(NA, 0)), destino[2], overwrite = TRUE, datatype = "INT1U",
                     gdal = "COMPRESS=DEFLATE")
  existe(destino)
}

correr_flowpy <- function(entradas, alpha, expo = 4, max_z = 11.5, cpus = 3) {
  dir <- dirname(entradas[1]); rel <- entradas[2]
  res <- file.path(dir, sprintf("res_a%s_e%s_z%s_%s", alpha, expo, max_z, tools::file_path_sans_ext(basename(rel))))
  cc <- list.files(res, "cellCounts\\.tif$", full.names = TRUE)
  if (length(cc) != 1 || file.mtime(cc) < max(file.mtime(entradas))) {
    unlink(res, recursive = TRUE)
    py <- Sys.getenv("FLOWPY_PYTHON", path.expand("~/.local/opt/avaframe-env/bin/python"))
    out <- system2("systemd-run", c("--user", "--scope", "-q", "-p", "MemoryMax=7G", shQuote(py),
                                    shQuote(normalizePath("bin/flowpy.py")), shQuote(entradas[1]), shQuote(rel),
                                    shQuote(dir), alpha, expo, max_z, cpus), stdout = TRUE, stderr = TRUE)
    cc <- list.files(res, "cellCounts\\.tif$", full.names = TRUE)
    if (length(cc) != 1) stop("Flow-Py falló (alpha = ", alpha, "): ", paste(utils::tail(out, 5), collapse = "\n"))
  }
  # Celdas no alcanzadas = 0; puntaje = log(1 + nº de celdas fuente cuyo flujo alcanza la celda)
  destino <- file.path(dir, sprintf("flowpy_a%s.tif", alpha))
  terra::writeRaster(log1p(terra::classify(terra::rast(cc), cbind(NA, 0))), destino, overwrite = TRUE,
                     datatype = "FLT4S", gdal = "COMPRESS=DEFLATE")
  existe(destino)
}
