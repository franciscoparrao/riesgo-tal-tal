# Integración del riesgo, modo corregido: 10 m, solo donde la vulnerabilidad está
# definida (unidades censales pobladas: "el riesgo se materializa donde estos elementos están
# presentes"). Resumen por área y por población.

# Ensamble de susceptibilidad: promedio de los mapas del modelo final del paper CFS-TOPSIS
# (revisión R1: XGBoost, Random Forest, LightGBM; 8 factores CFS; FABDEM ~30 m, EPSG:4326).
# Mismas grillas y máscara (cuencas costeras: 72 % de la comuna, 100 % de la población).
susceptibilidad_ensamble <- function(mapas, p) {
  destino <- salida(p, "tmp", "susceptibilidad_ensamble_4326.tif")
  tmp <- dirname(salida(p, "tmp", "terra", "x"))
  terra::terraOptions(memfrac = 0.2, progress = 0, tempdir = tmp)
  terra::app(terra::rast(mapas), mean, filename = destino, overwrite = TRUE,
             wopt = list(names = "susceptibilidad", datatype = "FLT4S",
                         gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES")))
  existe(destino)
}

# Susceptibilidad (EPSG:4326) -> EPSG:32719 -> grilla de 10 m (bilineal).
susceptibilidad_fina <- function(susc, plantilla_fina, p) {
  rep <- salida(p, "tmp", "susceptibilidad_32719.tif")
  destino <- salida(p, "componentes", "susceptibilidad.tif")
  surtgis("reproject", "--to", p$crs, "--method", "bilinear", susc, rep)
  surtgis("resample", "--method", "bilinear", "--reference", plantilla_fina, rep, destino)
  existe(destino)
}

# Suma ponderada convencional con terra (por bloques; `imagery calc` de SurtGIS 1.5.1 no hace streaming, issue #8).
riesgo_corregido <- function(susc, vuln, expo, locc, p, nombre) {
  destino <- salida(p, "riesgo", paste0("riesgo_", nombre, ".tif"))
  # Una sola pasada por bloques (lapp) y temporales en disco: /tmp es tmpfs (RAM) en este
  # equipo y la expresión encadenada creaba intermedios de 293 M celdas que lo llenaban.
  tmp <- dirname(salida(p, "tmp", "terra", "x"))
  terra::terraOptions(memfrac = 0.2, progress = 0, tempdir = tmp)
  w <- p$pesos
  s <- terra::rast(c(susc, vuln, expo, locc))
  f <- function(a, b, c, d) w$susceptibilidad * a + w$vulnerabilidad * b + w$exposicion * c + w$locc * d
  terra::lapp(s, f, filename = destino, overwrite = TRUE,
              wopt = list(names = "riesgo", datatype = "FLT4S",
                          gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES")))
  existe(destino)
}

# Puntos regulares por unidad poblada (1 cada `paso_m2`, entre 1 y `max_px`).
puntos_unidades <- function(unidades, p) {
  pv <- p$vulnerabilidad
  u <- sf::st_read(unidades, quiet = TRUE)
  u <- u[u$poblacion > 0, ]
  n <- pmin(pmax(ceiling(as.numeric(sf::st_area(u)) / pv$espectral_paso_m2), 1), pv$espectral_max_px)
  set.seed(p$semilla)
  do.call(rbind, lapply(seq_len(nrow(u)), function(i) {
    g <- sf::st_sample(sf::st_geometry(u)[i], n[i], type = "regular")
    if (length(g) == 0) g <- sf::st_point_on_surface(sf::st_geometry(u)[i])
    sf::st_sf(unidad = u$unidad[i], poblacion = u$poblacion[i], geometry = g)
  }))
}

# Jenks sobre los valores de las unidades pobladas; clases por área (puntos) y por población
# (cada unidad aporta su población a la clase de su riesgo medio).
niveles_riesgo <- c("Muy bajo", "Bajo", "Medio", "Alto", "Muy alto")

resumen_riesgo <- function(riesgo, puntos, p, nombre) {
  v <- terra::extract(terra::rast(riesgo), terra::vect(puntos), ID = FALSE)[[1]]
  d <- data.frame(unidad = puntos$unidad, poblacion = puntos$poblacion, riesgo = v)
  cobertura <- mean(!is.na(d$riesgo))   # fuera de las cuencas de la susceptibilidad -> NA
  d <- d[!is.na(d$riesgo), ]
  q <- classInt::classIntervals(d$riesgo, n = p$jenks_clases, style = "fisher")$brks
  clase <- function(x) cut(x, q, labels = niveles_riesgo, include.lowest = TRUE)
  por_unidad <- stats::aggregate(riesgo ~ unidad, d, mean)
  por_unidad$poblacion <- d$poblacion[match(por_unidad$unidad, d$unidad)]
  area <- prop.table(table(clase(d$riesgo)))
  pob <- tapply(por_unidad$poblacion, clase(por_unidad$riesgo), sum)
  pob[is.na(pob)] <- 0
  list(
    clases = data.frame(modelo = nombre, clase = niveles_riesgo,
                        desde = round(q[-length(q)], 4), hasta = round(q[-1], 4),
                        pct_area_poblada = round(100 * as.vector(area), 2),
                        pct_poblacion = round(100 * as.vector(pob) / sum(pob), 2),
                        pct_area_poblada_con_riesgo = round(100 * cobertura, 2)),
    unidades = stats::setNames(por_unidad[, c("unidad", "riesgo")], c("unidad", nombre))
  )
}

comparar_riesgos <- function(resumenes, p) {
  destino <- salida(p, "riesgo", "comparacion_vulnerabilidad.csv")
  clases <- do.call(rbind, lapply(resumenes, `[[`, "clases"))
  utils::write.csv(clases, salida(p, "riesgo", "clases_jenks.csv"), row.names = FALSE)
  u <- Reduce(function(a, b) merge(a, b, by = "unidad"), lapply(resumenes, `[[`, "unidades"))
  m <- as.matrix(u[, -1])
  cor_s <- stats::cor(m, method = "spearman")
  utils::write.csv(round(cor_s, 3), destino)
  existe(destino)
}

# Riesgo = amenaza x consecuencias, consecuencias = promedio(vulnerabilidad, exposición, LOCC).
# Sin amenaza no hay riesgo (la suma ponderada convencional daba riesgo alto sin amenaza y no
# predijo el daño de 2015: AUC 0,44 vs 0,77 multiplicativo; R/validacion_2015.R).
riesgo_multiplicativo <- function(amenaza, vuln, expo, locc, p, nombre) {
  destino <- salida(p, "riesgo", paste0("riesgo_", nombre, ".tif"))
  tmp <- dirname(salida(p, "tmp", "terra", "x"))
  terra::terraOptions(memfrac = 0.1, progress = 0, tempdir = tmp)
  terra::lapp(terra::rast(c(amenaza, vuln, expo, locc)), function(h, v, e, l) h * (v + e + l) / 3,
              filename = destino, overwrite = TRUE,
              wopt = list(names = "riesgo", datatype = "FLT4S",
                          gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES")))
  existe(destino)
}
