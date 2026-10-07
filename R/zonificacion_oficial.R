# Zonificación oficial de riesgo del Plan Regulador Comunal de Taltal (2012), anterior al evento de
# 2015: zonas ARA (área de riesgo por aluvión), ARR (remoción en masa) y ARM (tsunami), restringidas
# al desarrollo urbano. Fuente: MINVU, Geoportal Open Data, servicio 02__PRC_Antofagasta, capa 27
# "PRC Riesgo Taltal" (descargada 2026-10-06 a datos/entrada/minvu/prc_riesgo_taltal.geojson).
# Se evalúa como término de amenaza de un planificador: fracción de cada unidad dentro de la zona.

zonificacion_oficial <- function(prc, d, unidades, area_afectada, limite_urbano, semilla = 42, bloque_m = 1000) {
  lu <- sf::st_union(sf::st_transform(sf::st_make_valid(sf::st_read(limite_urbano, quiet = TRUE)), 32719))
  z <- sf::st_transform(sf::st_make_valid(sf::st_read(prc, quiet = TRUE)), 32719)
  a <- sf::st_union(sf::st_read(area_afectada, quiet = TRUE))
  u <- sf::st_read(unidades, quiet = TRUE); u <- u[u$unidad %in% d$unidad, ]
  frac <- function(g) {
    i <- suppressWarnings(sf::st_intersection(u[, "unidad"], g))
    f <- stats::setNames(as.numeric(sf::st_area(i)), i$unidad)
    x <- ifelse(as.character(u$unidad) %in% names(f), f[as.character(u$unidad)] / as.numeric(sf::st_area(u)), 0)
    x[match(d$unidad, u$unidad)]
  }
  zonas <- list("ARA (debris-flow/flood zone)" = "ARA", "ARA + ARR" = c("ARA", "ARR"))
  sp <- function(x, y) suppressWarnings(stats::cor(x, y, method = "spearman"))
  xy <- sf::st_coordinates(sf::st_point_on_surface(sf::st_geometry(u)))[match(d$unidad, u$unidad), ]
  bl <- paste(floor(xy[, 1] / bloque_m), floor(xy[, 2] / bloque_m))
  val <- do.call(rbind, lapply(names(zonas), function(n) {
    g <- sf::st_union(z[z$ZONA %in% zonas[[n]], ])
    s <- frac(g)
    r <- auc_ic(s, d$afectada, semilla = semilla)
    rb <- auc_bloques(s, d$afectada, bl, semilla = semilla)
    area_z <- as.numeric(sf::st_area(g)); dentro <- as.numeric(sf::st_area(sf::st_intersection(g, a)))
    a_u <- sf::st_intersection(a, lu); dentro_u <- as.numeric(sf::st_area(sf::st_intersection(a_u, g)))
    data.frame(zona = n, area_ha = round(area_z / 1e4, 1),
               pct_huella_dentro = round(100 * dentro / as.numeric(sf::st_area(a)), 1),
               pct_zona_afectada = round(100 * dentro / area_z, 1),
               pct_huella_urbana_dentro = round(100 * dentro_u / as.numeric(sf::st_area(a_u)), 1),
               unidades_en_zona = sum(s > 0), afectadas_en_zona = sum(d$afectada & s > 0),
               auc = round(r[["auc"]], 3), ic95_inf = round(r[["ic_inf"]], 3), ic95_sup = round(r[["ic_sup"]], 3),
               ic95_bloques = paste(round(rb[["ic_inf"]], 2), round(rb[["ic_sup"]], 2), sep = "-"),
               rho_pob_afectada = round(sp(s, d$pob_afectada), 3))
  }))
  # diferencia pareada alcance (HAND) - ARA, bootstrap por unidades y por bloques de 1 km
  s_ara <- frac(sf::st_union(z[z$ZONA == "ARA", ]))
  set.seed(semilla)
  por <- split(seq_len(nrow(d)), bl)
  f <- function(i) if (length(unique(d$afectada[i])) < 2) NA_real_ else auc(d$p_alcance[i], d$afectada[i]) - auc(s_ara[i], d$afectada[i])
  du <- replicate(2000, f(sample(nrow(d), replace = TRUE)))
  db <- replicate(2000, f(unlist(por[sample(names(por), replace = TRUE)], use.names = FALSE)))
  q <- function(v) paste(round(stats::quantile(v, c(.025, .975), na.rm = TRUE), 3), collapse = " to ")
  # unidades afectadas fuera de ARA: ¿las capta el índice de alcance?
  fuera <- d$afectada & s_ara == 0
  en_lu <- lengths(sf::st_intersects(suppressWarnings(sf::st_point_on_surface(u)), lu))[match(d$unidad, u$unidad)] > 0
  list(validacion = val,
       pareado = data.frame(contraste = "Reach index - ARA", dif_auc = round(f(seq_len(nrow(d))), 3),
                            ic95_unidades = q(du), ic95_bloques = q(db)),
       fuera_ara = data.frame(afectadas_fuera_ara = sum(fuera), afectadas_en_limite_urbano = sum(d$afectada & en_lu),
                              afectadas_en_limite_urbano_fuera_ara = sum(fuera & en_lu),
                              alcance_mediana_fuera = round(stats::median(d$p_alcance[fuera]), 3),
                              alcance_mediana_no_afectadas = round(stats::median(d$p_alcance[!d$afectada]), 3),
                              auc_alcance_fuera_vs_no_afectadas = round(auc(c(d$p_alcance[fuera], d$p_alcance[!d$afectada]),
                                                                             c(rep(TRUE, sum(fuera)), rep(FALSE, sum(!d$afectada)))), 3)))
}
