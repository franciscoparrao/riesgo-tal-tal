# Análisis de sensibilidad del resultado principal (riesgo = P(alcance) x consecuencias).

# 1. Umbral de cauce para HAND: recalcula HAND en las 6 comunas, recalibra SOLO en Atacama
#    (modelo final: solo HAND, comunas ponderadas por igual) y valida en Taltal.
sens_umbral <- function(catastro, censo_gpkg, area_2015, calibracion, p, umbrales = c(1000, 2000, 4000)) {
  do.call(rbind, lapply(umbrales, function(u) {
    d <- validacion_externa_datos(catastro, censo_gpkg, area_2015, p, umbral = u)
    s <- sensibilidad_alcance(d, calibracion, p$semilla)
    s <- s[grepl("final", s$modelo), ]
    m <- calibrar_alcance_final(d, calibracion)
    cbind(umbral_celdas = u, area_umbral_km2 = u * 900 / 1e6, b0 = round(stats::coef(m)[[1]], 3),
          b1 = round(stats::coef(m)[[2]], 3), s[, -1])
  }))
}

# 2. DEM: HAND de Taltal con FABDEM (sin edificios/vegetación) vs Copernicus GLO-30, aplicando
#    el modelo calibrado con Copernicus (sin reajuste).
sens_dem <- function(censo_gpkg, area_2015, fabdem, modelo_alcance, manzanas_validacion, p) {
  cop <- manzanas_validacion[manzanas_validacion$comuna == "TALTAL", ]
  fab <- datos_comuna("TALTAL", 2104, NULL, censo_gpkg, area_2015, p, umbral = 2000,
                      dem_local = fabdem, carpeta = "taltal_fabdem")
  fila <- function(d, nombre) {
    r <- auc_ic(stats::predict(modelo_alcance, d), d$afectada, semilla = p$semilla)
    data.frame(dem = nombre, n_manzanas = sum(is.finite(d$hand)), auc = round(r[["auc"]], 3),
               ic95_inf = round(r[["ic_inf"]], 3), ic95_sup = round(r[["ic_sup"]], 3),
               hand_mediana_afectadas = round(stats::median(d$hand[d$afectada], na.rm = TRUE), 1),
               hand_mediana_no_afectadas = round(stats::median(d$hand[!d$afectada], na.rm = TRUE), 1))
  }
  rbind(fila(cop, "Copernicus GLO-30 (calibración)"), fila(fab, "FABDEM v1.2"),
        data.frame(dem = "correlación Spearman HAND (Cop vs FABDEM)", n_manzanas = nrow(cop),
                   auc = round(stats::cor(cop$hand, fab$hand, method = "spearman", use = "complete.obs"), 3),
                   ic95_inf = NA, ic95_sup = NA, hand_mediana_afectadas = NA, hand_mediana_no_afectadas = NA))
}

# Tabla por unidad (dentro del marco del mapa 2015) con medias de las capas y afectación.
tabla_unidades_2015 <- function(area_afectada, marco, unidades, puntos, capas) {
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
  d$afectada <- d$frac_afectada > 0
  d$pob_afectada <- d$poblacion * d$frac_afectada
  d
}

metricas <- function(score, d, semilla) {
  r <- auc_ic(score, d$afectada, B = 500, semilla = semilla)
  rho <- suppressWarnings(stats::cor(score, d$pob_afectada, method = "spearman", use = "complete.obs"))
  c(auc = r[["auc"]], ic95_inf = r[["ic_inf"]], ic95_sup = r[["ic_sup"]], rho_pob_afectada = rho)
}

# 3 y 4. Formas del riesgo, variantes de vulnerabilidad y Monte Carlo de pesos de consecuencias.
sens_riesgo <- function(d, p, n_mc = 2000) {
  H <- d$p_alcance; E <- d$exposicion; L <- d$locc
  formas <- list(
    "solo amenaza P(alcance)" = H,
    "H x mean(V,E,L) [final, V = índice extendida]" = H * (d$v_extendida_indice + E + L) / 3,
    "H x mean(V,E,L), V = índice base" = H * (d$v_base_indice + E + L) / 3,
    "H x mean(V,E,L), V = fuzzy extendida" = H * (d$v_extendida_fuzzy + E + L) / 3,
    "H x mean(V,E,L), V = fuzzy base" = H * (d$v_base_fuzzy + E + L) / 3,
    "H x V (extendida)" = H * d$v_extendida_indice,
    "H x mean(V,E) (sin LOCC)" = H * (d$v_extendida_indice + E) / 2,
    "H x (V x E x L)^(1/3)" = H * (d$v_extendida_indice * E * L)^(1 / 3),
    "aditivo 0.25 (H,V,E,L)" = 0.25 * (H + d$v_extendida_indice + E + L),
    "aditivo convencional (S,V,E,L)" = 0.25 * (d$susceptibilidad + d$v_extendida_indice + E + L))
  fijas <- do.call(rbind, lapply(names(formas), function(n)
    data.frame(formulacion = n, t(round(metricas(formas[[n]], d, p$semilla), 3)))))
  set.seed(p$semilla)
  w <- matrix(stats::rexp(3 * n_mc), ncol = 3); w <- w / rowSums(w)        # Dirichlet(1,1,1)
  mc <- t(apply(w, 1, function(wi) {
    s <- H * (wi[1] * d$v_extendida_indice + wi[2] * E + wi[3] * L)
    c(auc = auc(s, d$afectada),
      rho = suppressWarnings(stats::cor(s, d$pob_afectada, method = "spearman", use = "complete.obs")))
  }))
  q <- function(x) round(stats::quantile(x, c(0.05, 0.5, 0.95)), 3)
  resumen_mc <- data.frame(formulacion = sprintf("Monte Carlo pesos (V,E,L) ~ Dirichlet(1,1,1), n=%d", n_mc),
                           metrica = c("AUC", "rho_pob_afectada"),
                           p05 = c(q(mc[, "auc"])[1], q(mc[, "rho"])[1]),
                           mediana = c(q(mc[, "auc"])[2], q(mc[, "rho"])[2]),
                           p95 = c(q(mc[, "auc"])[3], q(mc[, "rho"])[3]),
                           pct_auc_sobre_0.7 = c(round(100 * mean(mc[, "auc"] > 0.7), 1), NA))
  list(formulaciones = fijas, monte_carlo = resumen_mc)
}
