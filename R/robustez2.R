# Robustez II (revisión ciega, Landslides):
#   1. Incertidumbre de la huella 2015: desplazar el borde ±10/20/30 m (buffer negativo/positivo)
#   2. Nivel base marino: HAND sin corrección costera (celdas que drenan al mar sin cauce = NA)
#   3. Buffer de selección de manzanas en las comunas de calibración: 2, 5 y 10 km
#   4. Diferencias de AUC pareadas (bootstrap por unidad y por bloques espaciales de 1 km)

# 1. Huella desplazada: AUC en Taltal por unidad con la clase recalculada.
sens_huella <- function(d, area_afectada, unidades, desplazamientos = c(-30, -20, -10, 0, 10, 20, 30), semilla = 42) {
  a0 <- sf::st_union(sf::st_read(area_afectada, quiet = TRUE))
  u <- sf::st_read(unidades, quiet = TRUE)
  u <- u[u$unidad %in% d$unidad, ]
  d$riesgo <- riesgo_de(d)
  do.call(rbind, lapply(desplazamientos, function(b) {
    a <- if (b == 0) a0 else sf::st_buffer(a0, b)
    afect <- u$unidad[lengths(sf::st_intersects(u, a)) > 0]
    clase <- d$unidad %in% afect
    do.call(rbind, lapply(c("Reach index" = "p_alcance", "Reach-based risk" = "riesgo"), function(k) {
      r <- auc_ic(d[[k]], clase, semilla = semilla)
      data.frame(desplazamiento_m = b, area_ha = round(as.numeric(sf::st_area(a)) / 1e4, 1),
                 n_afectadas = sum(clase), capa = k, auc = round(r[["auc"]], 3),
                 ic95_inf = round(r[["ic_inf"]], 3), ic95_sup = round(r[["ic_sup"]], 3))
    }))
  }))
}

# Manzanas de una comuna con su distancia al área afectada (para 2 y 3).
manzanas_con_distancia <- function(comuna, cut, catastro, censo_gpkg, area_2015, p) {
  a <- afectada_comuna(catastro, comuna, p$crs, area_2015)
  m <- manzanas_comuna(censo_gpkg, cut, p$crs)
  m$dist_area <- as.numeric(sf::st_distance(m, a))
  m$afectada <- lengths(sf::st_intersects(m, a)) > 0
  m
}

extraer <- function(m, archivo) {
  terra::extract(terra::rast(archivo), terra::vect(sf::st_point_on_surface(sf::st_geometry(m))), ID = FALSE)[[1]]
}

# 2 y 3. HAND con y sin nivel base marino; buffers de selección. AUC por comuna (= AUC LOCO de un
# modelo de un solo predictor monótono, que no depende de los coeficientes).
sens_hand_buffer <- function(catastro, censo_gpkg, area_2015, manzanas_validacion, p,
                             buffers = c(2000, 5000, 10000), umbral = 2000) {
  stopifnot(nrow(manzanas_validacion) > 0)
  do.call(rbind, lapply(seq_len(nrow(comunas_validacion)), function(i) {
    comuna <- comunas_validacion$comuna[i]
    m <- manzanas_con_distancia(comuna, comunas_validacion$cut[i], catastro, censo_gpkg, area_2015, p)
    dir <- dirname(salida(p, "validacion_externa", carpeta_ascii(comuna), "x"))
    m$hand_costa <- extraer(m, file.path(dir, sprintf("hand_costa_%d.tif", umbral)))
    m$hand_cruda <- extraer(m, file.path(dir, sprintf("hand_%d.tif", umbral)))
    bs <- if (comuna == "TALTAL") Inf else buffers
    do.call(rbind, lapply(bs, function(b) {
      x <- m[m$dist_area <= b, ]
      f <- function(v) { ok <- is.finite(v); r <- auc_ic(-v[ok], x$afectada[ok], semilla = p$semilla)
                         c(round(r[["auc"]], 3), round(r[["ic_inf"]], 3), round(r[["ic_sup"]], 3), sum(!ok)) }
      hc <- f(x$hand_costa); hr <- f(x$hand_cruda)
      data.frame(comuna = comuna, buffer_m = b, n = nrow(x), afectadas = sum(x$afectada),
                 auc_hand = hc[1], ic95_inf = hc[2], ic95_sup = hc[3], n_sin_dato = hc[4],
                 auc_hand_sin_mar = hr[1], ic95_inf_sin_mar = hr[2], ic95_sup_sin_mar = hr[3],
                 n_sin_dato_sin_mar = hr[4])
    }))
  }))
}

# 4. Diferencias de AUC pareadas en Taltal (mismas unidades en cada réplica).
auc_pareado <- function(d, unidades, pares, B = 2000, semilla = 42, bloque_m = 1000) {
  u <- sf::st_read(unidades, quiet = TRUE)
  xy <- sf::st_coordinates(sf::st_point_on_surface(sf::st_geometry(u)))
  bl <- stats::setNames(paste(floor(xy[, 1] / bloque_m), floor(xy[, 2] / bloque_m)), u$unidad)
  d$bloque <- bl[as.character(d$unidad)]
  clase <- d$frac_afectada > 0
  set.seed(semilla)
  idx_u <- replicate(B, sample(nrow(d), replace = TRUE), simplify = FALSE)
  por <- split(seq_len(nrow(d)), d$bloque)
  idx_b <- replicate(B, unlist(por[sample(names(por), replace = TRUE)], use.names = FALSE), simplify = FALSE)
  dif <- function(a, b, i) if (length(unique(clase[i])) < 2) NA_real_ else auc(d[[a]][i], clase[i]) - auc(d[[b]][i], clase[i])
  do.call(rbind, lapply(pares, function(pr) {
    obs <- dif(pr[1], pr[2], seq_len(nrow(d)))
    du <- vapply(idx_u, function(i) dif(pr[1], pr[2], i), 0)
    db <- vapply(idx_b, function(i) dif(pr[1], pr[2], i), 0)
    q <- function(x) round(stats::quantile(x, c(0.025, 0.975), na.rm = TRUE), 3)
    data.frame(a = pr[1], b = pr[2], dif_auc = round(obs, 3),
               ic95_unidad = paste(q(du), collapse = " – "), ic95_bloque_1km = paste(q(db), collapse = " – "),
               p_unilateral_bloque = round(mean(db <= 0, na.rm = TRUE), 3))
  }))
}

# Tabla por unidad con todas las capas (validación principal + propagación + triviales + Flow-Py).
tabla_pareada <- function(csvs) {
  ts <- lapply(csvs, function(f) utils::read.csv(sub("/validacion_", "/unidades_", f)))
  d <- Reduce(function(x, y) merge(x, y[, setdiff(names(y), c("poblacion", "frac_afectada"))], by = "unidad",
                                   suffixes = c("", ".dup")), ts)
  d[, !grepl("\\.dup$", names(d))]
}
