# Validación con el aluvión de marzo de 2015 (área afectada Sernageomin, R/evento_2015.R).
# Limitaciones a declarar: un solo evento; polígono digitalizado de capturas (±10-20 m); el área
# afectada sigue aguas arriba fuera del mapa (se evalúa solo dentro del marco); vulnerabilidad
# del Censo 2024 (9 años después del evento, con reconstrucción).

# Marco del mapa de 2015: unión de las extensiones georreferenciadas de ambas capturas.
marco_2015 <- function(png_norte, png_sur, p) {
  destino <- salida(p, "validacion", "marco_2015.gpkg")
  e <- lapply(c("norte", "sur"), function(n) {
    f <- if (n == "norte") png_norte else png_sur
    sf::st_as_sfc(sf::st_bbox(terra::ext(georref_captura(f, grilla_2015[[n]], p$crs)) |> as.vector() |>
                                stats::setNames(c("xmin", "xmax", "ymin", "ymax")), crs = p$crs))
  })
  sf::st_write(sf::st_sf(geometry = sf::st_union(e[[1]], e[[2]])), destino, delete_dsn = TRUE, quiet = TRUE)
  existe(destino)
}

auc <- function(score, clase) {           # Mann-Whitney: P(score afectada > score no afectada)
  r <- rank(score); n1 <- sum(clase); n0 <- sum(!clase)
  (sum(r[clase]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

auc_ic <- function(score, clase, B = 2000, semilla = 42) {
  ok <- is.finite(score) & !is.na(clase)          # rank() manda los NA al final y sesga el AUC
  score <- score[ok]; clase <- clase[ok]
  set.seed(semilla)
  i1 <- which(clase); i0 <- which(!clase)
  b <- replicate(B, { s <- c(sample(i1, replace = TRUE), sample(i0, replace = TRUE))
                      auc(score[s], clase[s]) })
  c(auc = auc(score, clase), ic_inf = unname(stats::quantile(b, 0.025)),
    ic_sup = unname(stats::quantile(b, 0.975)), n_excluidos = sum(!ok))
}

# Valor medio de cada capa por unidad poblada (puntos regulares de puntos_unidades()).
# Capa por capa: las grillas pueden diferir (p. ej. recortes de propagación a 30 m).
medias_unidades <- function(capas, puntos) {
  pv <- terra::vect(puntos)
  v <- as.data.frame(lapply(unname(capas), function(r) terra::extract(terra::rast(r), pv, ID = FALSE)[[1]]))
  names(v) <- names(capas)
  stats::aggregate(v, list(unidad = puntos$unidad), mean, na.rm = TRUE)
}

validar_2015 <- function(area_afectada, marco, unidades, puntos, capas, p, nombre = "validacion_2015") {
  destino <- salida(p, "validacion", paste0(nombre, ".csv"))
  a <- sf::st_union(sf::st_read(area_afectada, quiet = TRUE))
  m <- sf::st_union(sf::st_read(marco, quiet = TRUE))
  u <- sf::st_read(unidades, quiet = TRUE)
  u <- u[u$poblacion > 0, ]
  u <- u[lengths(sf::st_intersects(sf::st_point_on_surface(sf::st_geometry(u)), m)) > 0, ]
  inter <- suppressWarnings(sf::st_intersection(u[, "unidad"], a))
  frac <- stats::setNames(as.numeric(sf::st_area(inter)), inter$unidad)
  u$frac_afectada <- ifelse(as.character(u$unidad) %in% names(frac),
                            frac[as.character(u$unidad)] / as.numeric(sf::st_area(u)), 0)
  d <- merge(sf::st_drop_geometry(u)[, c("unidad", "poblacion", "frac_afectada")],
             medias_unidades(capas, puntos), by = "unidad")
  utils::write.csv(d, salida(p, "validacion", sub("validacion", "unidades", nombre) |> paste0(".csv")),
                  row.names = FALSE)
  filas <- list()
  for (umbral in c(0, 0.25)) {
    clase <- d$frac_afectada > umbral
    for (n in names(capas)) {
      ok <- !is.na(d[[n]])
      r <- auc_ic(d[[n]][ok], clase[ok], semilla = p$semilla)
      rho <- suppressWarnings(stats::cor.test(d[[n]][ok], (d$poblacion * d$frac_afectada)[ok],
                                              method = "spearman", exact = FALSE))
      filas[[length(filas) + 1]] <- data.frame(
        definicion = if (umbral == 0) "cualquier superposición" else "≥ 25 % de la unidad",
        capa = n, n_unidades = sum(ok), n_afectadas = sum(clase[ok]),
        pob_afectada = sum(d$poblacion[ok & clase]), auc = r[["auc"]],
        ic95_inf = r[["ic_inf"]], ic95_sup = r[["ic_sup"]],
        rho_pob_afectada = rho$estimate[[1]], p_rho = signif(rho$p.value, 2))
    }
  }
  t <- do.call(rbind, filas)
  t[, c("auc", "ic95_inf", "ic95_sup", "rho_pob_afectada")] <- round(t[, c("auc", "ic95_inf", "ic95_sup", "rho_pob_afectada")], 3)
  utils::write.csv(t, destino, row.names = FALSE)
  existe(destino)
}

# Validación del peligro a nivel de píxel dentro del marco: AUC de la susceptibilidad y curva
# de éxito (fracción del área afectada capturada vs fracción del área ordenada por susceptibilidad).
validar_susceptibilidad_2015 <- function(area_afectada, marco, susceptibilidad, p, paso = 20) {
  destino <- salida(p, "validacion", "susceptibilidad_2015.csv")
  a <- sf::st_union(sf::st_read(area_afectada, quiet = TRUE))
  m <- sf::st_union(sf::st_read(marco, quiet = TRUE))
  pts <- sf::st_sf(geometry = sf::st_make_grid(m, cellsize = paso, what = "centers"))
  pts <- pts[lengths(sf::st_intersects(pts, m)) > 0, ]
  s <- terra::extract(terra::rast(susceptibilidad), terra::vect(pts), ID = FALSE)[[1]]
  dentro <- lengths(sf::st_intersects(pts, a)) > 0
  ok <- !is.na(s)
  r <- auc_ic(s[ok], dentro[ok], B = 500, semilla = p$semilla)
  o <- order(s[ok], decreasing = TRUE)
  curva <- data.frame(frac_area = seq_along(o) / length(o), frac_afectada = cumsum(dentro[ok][o]) / sum(dentro[ok]))
  top <- function(q) curva$frac_afectada[which.min(abs(curva$frac_area - q))]
  utils::write.csv(data.frame(n_pixeles = sum(ok), pct_afectado = round(100 * mean(dentro[ok]), 2),
                              auc = round(r[["auc"]], 3), ic95_inf = round(r[["ic_inf"]], 3),
                              ic95_sup = round(r[["ic_sup"]], 3),
                              captura_top10 = round(top(0.10), 3), captura_top20 = round(top(0.20), 3)),
                   destino, row.names = FALSE)
  existe(destino)
}
