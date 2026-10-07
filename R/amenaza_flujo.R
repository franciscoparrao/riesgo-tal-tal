# Amenaza por flujo (fuente -> propagación). La validación con el aluvión de 2015 mostró que la
# susceptibilidad de inicio no discrimina las manzanas dañadas (AUC 0,54), mientras que la altura
# sobre el cauce (HAND) sí lo hace en las 6 comunas del evento (AUC 0,73-0,89).
#   P(alcance)      = logit^-1(b0 + b1 * log1p(HAND)), calibrada SOLO en comunas de Atacama
#                     (comunas ponderadas por igual; R/validacion_externa.R)
#   S aguas arriba  = promedio de la susceptibilidad de inicio del área que drena a la celda
#   amenaza_flujo   = S aguas arriba * P(alcance)
#   amenaza_comb    = max(S local, amenaza_flujo)   (inicio en laderas o llegada por el cauce)

calibrar_alcance_final <- function(d, calibracion) {
  d <- d[d$comuna %in% calibracion & is.finite(d$hand), ]
  d$w <- as.numeric(1 / table(d$comuna)[d$comuna]); d$w <- d$w / mean(d$w)
  suppressWarnings(stats::glm(afectada ~ log1p(hand), stats::binomial(), d, weights = w))
}

amenaza_flujo <- function(comuna_utm, susceptibilidad, modelo_alcance, plantilla_fina, p,
                          umbral = 2000, margen_m = 10000) {
  dir <- dirname(salida(p, "amenaza", "x"))
  f <- function(n) file.path(dir, n)
  com <- sf::st_read(comuna_utm, quiet = TRUE)
  t <- terreno_comuna(sf::st_union(com), dir, p, umbral = umbral, margen_m = margen_m)
  tmp <- dirname(salida(p, "tmp", "terra", "x"))
  # memfrac bajo + cada intermedio a disco: con 0.3 los temporales de ~41 M celdas y
  # flowAccumulation superaron 6 GB. (SurtGIS no tiene acumulación ponderada: issue #12.)
  terra::terraOptions(memfrac = 0.1, progress = 0, tempdir = tmp)
  wo <- list(datatype = "FLT4S", gdal = c("COMPRESS=DEFLATE", "TILED=YES"))

  # P(alcance) desde HAND (con corrección costera)
  b <- stats::coef(modelo_alcance)
  terra::app(terra::rast(t[["hand"]]), function(h) stats::plogis(b[[1]] + b[[2]] * log1p(pmax(h, 0))),
             filename = f("p_alcance_30m.tif"), overwrite = TRUE, wopt = wo)

  # Susceptibilidad aguas arriba: acumulación ponderada / acumulación de celdas con dato
  dem <- terra::rast(f("dem_fill.tif"))
  terra::project(terra::rast(susceptibilidad), dem, method = "bilinear",
                 filename = f("susc_30m.tif"), overwrite = TRUE, wopt = wo)
  terra::terrain(dem, v = "flowdir", filename = f("flowdir_terra.tif"), overwrite = TRUE)
  s30 <- terra::rast(f("susc_30m.tif")); fd <- terra::rast(f("flowdir_terra.tif"))
  terra::classify(s30, cbind(NA, 0), filename = f("peso_s.tif"), overwrite = TRUE, wopt = wo)
  terra::ifel(is.na(s30), 0, 1, filename = f("peso_n.tif"), overwrite = TRUE, wopt = wo)
  terra::flowAccumulation(fd, weight = terra::rast(f("peso_s.tif")), filename = f("acc_s.tif"), overwrite = TRUE)
  invisible(gc())
  terra::flowAccumulation(fd, weight = terra::rast(f("peso_n.tif")), filename = f("acc_n.tif"), overwrite = TRUE)
  invisible(gc())
  acc_n <- terra::rast(f("acc_n.tif"))
  terra::ifel(acc_n > 0, terra::rast(f("acc_s.tif")) / acc_n, NA,
              filename = f("susc_aguas_arriba_30m.tif"), overwrite = TRUE, wopt = wo)
  s_up <- terra::rast(f("susc_aguas_arriba_30m.tif")); palc <- terra::rast(f("p_alcance_30m.tif"))

  terra::lapp(c(s_up, palc), function(a, b) a * b, filename = f("amenaza_flujo_30m.tif"),
              overwrite = TRUE, wopt = wo)
  terra::lapp(c(s30, terra::rast(f("amenaza_flujo_30m.tif"))), function(a, b) pmax(a, b, na.rm = TRUE),
              filename = f("amenaza_combinada_30m.tif"), overwrite = TRUE, wopt = wo)
  # Variante de sensibilidad: área fuente susceptible TOTAL aguas arriba (log, 0-1) x P(alcance)
  acc_s <- terra::rast(f("acc_s.tif"))
  mx <- terra::global(log1p(acc_s), "max", na.rm = TRUE)[[1]]
  terra::lapp(c(acc_s, palc), function(a, b) log1p(a) / mx * b, filename = f("amenaza_fuente_total_30m.tif"),
              overwrite = TRUE, wopt = wo)
  capas <- c(p_alcance = "p_alcance", flujo = "amenaza_flujo", combinada = "amenaza_combinada",
             fuente_total = "amenaza_fuente_total")
  out <- vapply(capas, function(n) salida(p, "componentes", paste0(n, ".tif")), "")
  for (n in names(capas)) {
    surtgis("resample", "--method", "bilinear", "--reference", plantilla_fina,
            f(paste0(capas[[n]], "_30m.tif")), out[[n]])
  }
  existe(out)
}

# Sensibilidad del modelo de alcance: variantes probadas, LOCO en Atacama (criterio de
# elección) y AUC en Taltal (validación independiente, NO usado para elegir).
sensibilidad_alcance <- function(d, calibracion, semilla) {
  dc <- d[d$comuna %in% calibracion & is.finite(d$hand), ]
  dc$w_iguales <- as.numeric(1 / table(dc$comuna)[dc$comuna]); dc$w_iguales <- dc$w_iguales / mean(dc$w_iguales)
  dc$w_uno <- 1
  x <- d[d$comuna == "TALTAL" & is.finite(d$hand), ]
  variantes <- list(
    "HAND + distancia (sin ponderar)" = list(f = afectada ~ log1p(hand) + log1p(dist), w = "w_uno"),
    "HAND + distancia (comunas iguales)" = list(f = afectada ~ log1p(hand) + log1p(dist), w = "w_iguales"),
    "solo distancia (comunas iguales)" = list(f = afectada ~ log1p(dist), w = "w_iguales"),
    "solo HAND (comunas iguales) [final]" = list(f = afectada ~ log1p(hand), w = "w_iguales"))
  ajustar <- function(v, datos) { datos$.w <- datos[[v$w]]; suppressWarnings(stats::glm(v$f, stats::binomial(), datos, weights = .w)) }
  do.call(rbind, lapply(names(variantes), function(n) {
    v <- variantes[[n]]
    loco <- vapply(calibracion, function(cm) {
      y <- dc[dc$comuna == cm, ]; auc(stats::predict(ajustar(v, dc[dc$comuna != cm, ]), y), y$afectada)
    }, 0)
    r <- auc_ic(stats::predict(ajustar(v, dc), x), x$afectada, semilla = semilla)
    data.frame(modelo = n, auc_loco_media = round(mean(loco), 3), auc_loco_min = round(min(loco), 3),
               auc_taltal = round(r[["auc"]], 3), ic95_inf = round(r[["ic_inf"]], 3), ic95_sup = round(r[["ic_sup"]], 3))
  }))
}
