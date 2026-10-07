# Robustez de la validación en Taltal (respuesta a revisión simulada, Landslides):
#   1. IC por bootstrap espacial en bloques (las unidades afectadas son contiguas: una quebrada)
#   2. Correlaciones sin circularidad: ρ con la fracción afectada y ρ parcial con la población
#      afectada controlando la población de la unidad
#   3. Dominancia de la amenaza en el producto R = R_h · C
#   4. Sensibilidad del ranking a los coeficientes del modelo de alcance (β0, β1)
#   5. Población 2017 vs 2024 en las zonas censales afectadas (desfase censo-evento)
# Todas las métricas usan la tabla de unidades de tabla_unidades_2015() (R/sensibilidad.R).

riesgo_de <- function(d, rh = d$p_alcance) rh * (d$v_extendida_indice + d$exposicion + d$locc) / 3

auc_bloques <- function(score, clase, bloque, B = 2000, semilla = 42) {
  ok <- is.finite(score); score <- score[ok]; clase <- clase[ok]; bloque <- bloque[ok]
  ids <- unique(bloque); por <- split(seq_along(score), bloque)
  set.seed(semilla)
  b <- replicate(B, {
    i <- unlist(por[as.character(sample(ids, replace = TRUE))], use.names = FALSE)
    if (length(unique(clase[i])) < 2) NA_real_ else auc(score[i], clase[i])
  })
  c(auc = auc(score, clase), ic_inf = unname(stats::quantile(b, 0.025, na.rm = TRUE)),
    ic_sup = unname(stats::quantile(b, 0.975, na.rm = TRUE)), n_bloques = length(ids))
}

# ρ de Spearman parcial: correlación de los residuos de los rangos.
rho_parcial <- function(x, y, z) {
  r <- function(a) rank(a); ok <- is.finite(x) & is.finite(y) & is.finite(z)
  stats::cor(stats::resid(stats::lm(r(x[ok]) ~ r(z[ok]))), stats::resid(stats::lm(r(y[ok]) ~ r(z[ok]))))
}

robustez_2015 <- function(d, unidades, puntos, hand, manzanas, calibracion, modelo, p,
                          tamanos = c(250, 500, 1000)) {
  u <- sf::st_read(unidades, quiet = TRUE)
  xy <- sf::st_coordinates(sf::st_point_on_surface(sf::st_geometry(u)))
  xy <- data.frame(unidad = u$unidad, x = xy[, 1], y = xy[, 2])
  d <- merge(d, xy, by = "unidad")
  d$riesgo <- riesgo_de(d)
  d$aditivo <- 0.25 * (d$susceptibilidad + d$v_extendida_indice + d$exposicion + d$locc)
  capas <- c("Reach index" = "p_alcance", "Reach-based risk" = "riesgo",
             "Initiation susceptibility" = "susceptibilidad", "Conventional additive" = "aditivo")

  # 1. Bootstrap espacial en bloques cuadrados
  bloques <- do.call(rbind, lapply(tamanos, function(t) {
    b <- paste(floor(d$x / t), floor(d$y / t))
    do.call(rbind, lapply(names(capas), function(n) {
      r <- auc_bloques(d[[capas[[n]]]], d$afectada, b, semilla = p$semilla)
      data.frame(capa = n, bloque_m = t, n_bloques = r[["n_bloques"]], auc = r[["auc"]],
                 ic95_inf = r[["ic_inf"]], ic95_sup = r[["ic_sup"]])
    }))
  }))

  # 2. Correlaciones
  sp <- function(a, b) suppressWarnings(stats::cor(a, b, method = "spearman", use = "complete.obs"))
  correl <- do.call(rbind, lapply(names(capas), function(n) {
    s <- d[[capas[[n]]]]
    data.frame(capa = n, rho_pob_afectada = sp(s, d$pob_afectada), rho_frac_afectada = sp(s, d$frac_afectada),
               rho_parcial_pob_afectada_dado_pob = rho_parcial(s, d$pob_afectada, d$poblacion),
               rho_pob_afectada_solo_afectadas = sp(s[d$afectada], d$pob_afectada[d$afectada]),
               rho_frac_solo_afectadas = sp(s[d$afectada], d$frac_afectada[d$afectada]))
  }))

  # 3. Dominancia de la amenaza
  alc <- d$p_alcance > stats::median(d$p_alcance)
  C <- (d$v_extendida_indice + d$exposicion + d$locc) / 3
  lr <- log(d$riesgo); lh <- log(d$p_alcance); lc <- log(C)
  dominancia <- data.frame(
    metrica = c("rho(R, R_h) todas", "rho(R, R_h) R_h > mediana", "rho(R, C) R_h > mediana",
                "varianza log R_h / (var log R_h + var log C)", "CV R_h", "CV C"),
    valor = c(sp(d$riesgo, d$p_alcance), sp(d$riesgo[alc], d$p_alcance[alc]), sp(d$riesgo[alc], C[alc]),
              stats::var(lh) / (stats::var(lh) + stats::var(lc)),
              stats::sd(d$p_alcance) / mean(d$p_alcance), stats::sd(C) / mean(C)))

  # 4. Coeficientes del alcance: R_h por punto y promedio por unidad (como p_alcance)
  h <- terra::extract(terra::rast(hand), terra::vect(puntos), ID = FALSE)[[1]]
  b_final <- stats::coef(modelo)
  b_salado <- stats::coef(calibrar_alcance_final(manzanas, c("CHAÑARAL", "DIEGO DE ALMAGRO")))
  casos <- list("final (Atacama)" = b_final, "Salado valley only" = b_salado,
                "beta1 x 0.5" = c(b_final[[1]], 0.5 * b_final[[2]]), "beta1 x 2" = c(b_final[[1]], 2 * b_final[[2]]),
                "beta0 - 2" = c(b_final[[1]] - 2, b_final[[2]]), "beta0 + 2" = c(b_final[[1]] + 2, b_final[[2]]))
  r_ref <- NULL
  beta <- do.call(rbind, lapply(names(casos), function(n) {
    b <- casos[[n]]
    rh_p <- stats::plogis(b[[1]] + b[[2]] * log1p(pmax(h, 0)))
    rh <- stats::aggregate(rh_p, list(unidad = puntos$unidad), mean, na.rm = TRUE)
    rh <- rh$x[match(d$unidad, rh$unidad)]
    r <- riesgo_de(d, rh); if (is.null(r_ref)) r_ref <<- r
    a <- auc_ic(r, d$afectada, semilla = p$semilla)
    data.frame(caso = n, beta0 = round(b[[1]], 3), beta1 = round(b[[2]], 3), auc_riesgo = a[["auc"]],
               ic95_inf = a[["ic_inf"]], ic95_sup = a[["ic_sup"]], rho_pob_afectada = sp(r, d$pob_afectada),
               rho_ranking_vs_final = sp(r, r_ref), top20_comun_final = length(intersect(
                 order(-r)[1:20], order(-r_ref)[1:20])))
  }))

  num <- function(t) { i <- vapply(t, is.numeric, TRUE); t[i] <- round(t[i], 3); t }
  list(bloques = num(bloques), correlaciones = num(correl), dominancia = num(dominancia), beta = num(beta))
}

# Población 2017 (zonas censales) vs 2024 (manzanas agregadas a esas zonas), por zona afectada.
censo_2017_vs_2024 <- function(zonas_2017, unidades, area_afectada, marco, p) {
  z <- sf::st_transform(sf::st_read(zonas_2017, quiet = TRUE), p$crs)
  z$pob_2017 <- rowSums(sf::st_drop_geometry(z)[, c("X0_15", "X16_25", "X26_50", "X51_64", "X64_mas")])
  u <- sf::st_read(unidades, quiet = TRUE)
  pts <- sf::st_point_on_surface(u[, "poblacion"])
  j <- sf::st_join(pts, z[, "GEOCODIGO"])
  p24 <- stats::aggregate(poblacion ~ GEOCODIGO, sf::st_drop_geometry(j), sum)
  z <- merge(z, stats::setNames(p24, c("GEOCODIGO", "pob_2024")), by = "GEOCODIGO", all.x = TRUE)
  z$pob_2024[is.na(z$pob_2024)] <- 0
  a <- sf::st_union(sf::st_read(area_afectada, quiet = TRUE))
  m <- sf::st_union(sf::st_read(marco, quiet = TRUE))
  z <- z[lengths(sf::st_intersects(z, m)) > 0, ]
  inter <- sf::st_intersects(z, a)
  z$frac_afectada <- 0
  z$frac_afectada[lengths(inter) > 0] <- vapply(which(lengths(inter) > 0), function(i)
    as.numeric(sf::st_area(sf::st_intersection(sf::st_geometry(z)[i], a))) / as.numeric(sf::st_area(z[i, ])), 0)
  out <- sf::st_drop_geometry(z)[, c("GEOCODIGO", "pob_2017", "pob_2024", "frac_afectada")]
  out$cambio_pct <- round(100 * (out$pob_2024 / out$pob_2017 - 1), 1)
  out$frac_afectada <- round(out$frac_afectada, 3)
  out[order(-out$frac_afectada), ]
}
